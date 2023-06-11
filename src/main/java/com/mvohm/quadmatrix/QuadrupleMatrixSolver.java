/*

 Copyright 2021-2023 M.Vokhmentsev

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

*/

package com.mvohm.quadmatrix;

import com.mvohm.quadmatrix.Matrix.ErrorCodes;
import com.mvohm.quadruple.Quadruple;

import static com.mvohm.quadmatrix.Util.*;



class QuadrupleMatrixSolver extends MatrixSolver {

  // TODO_for_next_release 2023-05-27 11:10:12 Reorder methods accordingly with their usage

  enum SolutionMethod {NONE, CHOLESKY, LU};

  private final int size;
  private final Quadruple[][] matrix;  // Original matrix
  private Quadruple[] vector;          // Original vector, used for refinement

  private SolutionMethod lastSolvedWith = SolutionMethod.NONE;
  private ErrorCodes errorCode = ErrorCodes.OK;

  private Quadruple[][] luDecomposition;
  private boolean luDecompositionError;

  private Quadruple[][] choleskyDecomposition;
  private boolean choleskyDecompositionError;

  private int[] pivot;              // filled by decomposeLU(), used by solveLU()
  private Quadruple[] rowScales;
  private Quadruple[] solution;

  private Quadruple[][] inversion = null;

  private int detSignCorrection = 1;            // To correct the determinant sign in case of swapping rows while pivoting
  private Quadruple determinant = Quadruple.nan();      // Unless already computed
  private Quadruple norm;




  QuadrupleMatrixSolver(Quadruple[][] matrix, boolean needToScale) {
    this.size = matrix.length;
    this.matrix = matrix;
    this.needToScale = needToScale;
  }



  /* *******************************************************************************
  /***** Solutions with respect of vectors *****************************************
  /*********************************************************************************/

  /**
   * Solves a system formed by the matrix and the given vector
   * using LU-decomposition
   * @param vector
   * @return the solution
   */
  Quadruple[] solveLU(Quadruple[] vector) {
    this.vector = deepCopyOf(vector);
    solution = solveLUInternally(vector);
    lastSolvedWith = SolutionMethod.LU;
    return solution.clone();                  // To protect the internally stored solutions for possible subsequent retrievals
  }

  Quadruple[] solveLUAccurately(Quadruple[] vector) {
    this.vector = deepCopyOf(vector);
    solution = solveLUInternally(vector);
    lastSolvedWith = SolutionMethod.LU;
    refine();
    return solution.clone();
  }

  Quadruple[] solveCholesky(Quadruple[] vector) {
    this.vector = deepCopyOf(vector);
    solution = solveCholeskyInternally(vector);
    lastSolvedWith = SolutionMethod.CHOLESKY;
    return solution.clone();                  // To protect the internally stored solutions for possible subsequent retrievals
  }

  Quadruple[] solveCholeskyAccurately(Quadruple[] vector) {
    this.vector = deepCopyOf(vector);
    solution = solveCholeskyInternally(vector);
    lastSolvedWith = SolutionMethod.CHOLESKY;
    refine();
    return solution.clone();
  }

  Quadruple[] getSolution() {
    return solution;
  }

  String errorCode() {
    return errorCode.toString();
  }

  /* ********************************************************************************
  /*** Solutions with respect of matrices *******************************************
  /* ********************************************************************************/

  /** Attention -- spoils the matrixB, but returns a newly-created array */
  Quadruple[][] solve(Quadruple[][] matrixB) {
    if (luDecompositionError)
      throwNonInvertibleError();

    if (luDecomposition == null) { // The very first invocation.
      scaleAndDecompose();
    }

    // Copy right hand side with pivoting and scaling
    final Quadruple[][] matrixX = scaleAndPermute(matrixB, pivot); // Spoils matrixB, and returns a new array

    // Solve L*Y = B(piv,:)
    for (int k = 0; k < size; k++) {
      for (int i = k + 1; i < size; i++) {
        for (int j = 0; j < size; j++) {
          matrixX[i][j].subtract(Quadruple.multiply(matrixX[k][j], luDecomposition[i][k]));
        }
      }
    }

    // Solve U * X = Y;
    for (int k = size - 1; k >= 0; k--) {
      for (int j = 0; j < size; j++) {
        matrixX[k][j].divide(luDecomposition[k][k]);
      }
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < size; j++) {
          matrixX[i][j].subtract(Quadruple.multiply(matrixX[k][j], luDecomposition[i][k]));
        }
      }
    }

    return matrixX;
  }

  /* ********************************************************************************
  /*** Inverse **********************************************************************
  /* ********************************************************************************/

  Quadruple[][] inverse() {
    if (inversion == null)
      inversion = solve(unityMatrix());
    return inversion;
  }

  Quadruple[][] inverseAccurately() {
    return solveAccurately(unityMatrix());
  }

  /**
   * Creates a new two-dimentional double array and fills it with the transposition of the source array
   * @param matrix
   * @return
   */
  Quadruple[][] transpose(Quadruple[][] matrix) {
    final int size = matrix.length;
    final Quadruple[][] result = new Quadruple[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        result[i][j] = new Quadruple(matrix[j][i]);
      }
    }
    return result;
  }

  Quadruple[][] unityMatrix() {
    final Quadruple[][] result = new Quadruple[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        result[i][j] = new Quadruple();
      }
      result[i][i] = Quadruple.one();
    }
    return result;
  }

  /* *******************************************************************************
  /* **** Multiplications **********************************************************
  /*********************************************************************************/


  Quadruple[][] multiply(Quadruple[][] factor) {
    final Quadruple[][] product = new Quadruple[size][size];
    final Quadruple[] prodVector = new Quadruple[size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        for (int k = 0; k < size; k++) {
          prodVector[k] = Quadruple.multiply(matrix[i][k], factor[k][j]);
        }
        product[i][j] = sumOfVector(prodVector);
      }
    }
    return product;
  }

  Quadruple[] multiply(Quadruple[] factor) {
    final int size = factor.length;
    final Quadruple[] product = new Quadruple[size];
    for (int i = 0; i < size; i++) {
      final Quadruple[] rowProduct = multiplyElements(matrix[i], factor);
      product[i] = sumOfVector(rowProduct);
    }
    return product;
  }

  Quadruple[][] multiply(Quadruple scalar) {
    final Quadruple[][] product = new Quadruple[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        product[i][j] = Quadruple.multiply(scalar, matrix[i][j]);
      }
    }
    return product;
  }

  /* *******************************************************************************
  /* **** Additions and subtractions ***********************************************
  /*********************************************************************************/

  static Quadruple[][] addMatrices(Quadruple[][] summandA, Quadruple[][] summandB) {
    final int size = summandA.length;
    final Quadruple[][] sum = new Quadruple[size][size];
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++) {
        sum[i][j] = Quadruple.add(summandA[i][j], summandB[i][j]);
      }
    return sum;
  }

  static Quadruple[][] subtractMatrices(Quadruple[][] minuend, Quadruple[][] subtrahend) {
    final int size = minuend.length;
    final Quadruple[][] difference = new Quadruple[size][size];
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++) {
        difference[i][j] = Quadruple.subtract(minuend[i][j], subtrahend[i][j]);
      }
    return difference;
  }

  Quadruple determinant() {
    if (determinant.isNaN()) {
      if (luDecompositionError) {
        determinant = new Quadruple();
      } else {
        try {
          determinant = computeDeterminant();
        } catch (final IllegalArgumentException x) { // Non-invertible
          determinant = new Quadruple();
        }
      }
    }
    return determinant;
  }

  Quadruple norm() {
    if (norm == null)
      norm = computeNorm(matrix);
    return norm;
  }

  double cond() {
    if (errorCode == ErrorCodes.NON_INVERTIBLE)
      return Double.POSITIVE_INFINITY;
    final double matrixNorm = computeNorm(matrix).doubleValue();
    try {
      final double inversionNorm = computeNorm(inverse()).doubleValue();
      return matrixNorm * inversionNorm;
    } catch (final IllegalArgumentException x) { // Matrix is non-invertible
      return Double.POSITIVE_INFINITY;
    }
  }

  /* ************************************************************************************
  /*** Private methods ******************************************************************
  /**************************************************************************************/

  private static final double MIN_CORRECTION_FACTOR = 1.0 / 8;
  private static final int MAX_REFINEMENT_ITERATIONS = 20;

  /**
   * @param vector
   * @return a newly-created vector containing the solution
   */
  private Quadruple[] solveLUInternally(Quadruple[] vector) {
    if (matrix == null || vector == null || size != vector.length)
      throw new IllegalArgumentException("Matrix and vector must be of equal size");

    if (luDecompositionError)
      throwNonInvertibleError();

    errorCode = ErrorCodes.OK;                                    // May remain after failed SPD decomposition
    vector = deepCopyOf(vector);                                  // preserve the data passed-in

    if (luDecomposition == null) {                                // The very first invocation.
      vector = scaleMatrixAndVector(vector);          // creates luDecomposition and a scaled copy of matrix
      decomposeLU();                                              // modifies luDecomposition and creates pivot.
    } else {
      vector = scaleVector(vector);                               // scales and returns vector
    }

    final Quadruple[] solution = getPermutted(vector, pivot);     // returns a new array

    // Solve L * Y = B(piv,:)
    for (int k = 0; k < size; k++)
     for (int i = k + 1; i < size; i++) {
       solution[i].subtract(Quadruple.multiply(solution[k], luDecomposition[i][k]));    // (N^2 - N) / 2
     }

    // Solve U * X = Y;
    for (int k = size - 1; k >= 0; k--) {
      solution[k].divide(luDecomposition[k][k]);                  // N
      for (int i = 0; i < k; i++) {
        solution[i].subtract(Quadruple.multiply(solution[k], luDecomposition[i][k]));   // (N^2 - N) / 2
      }
    }

    // (N^2 - N) multiplications and subtractions, N divisions
    return solution;
  }

  private void throwNonInvertibleError() {
    errorCode = ErrorCodes.NON_INVERTIBLE;
    this.luDecompositionError = true;
    throw new IllegalArgumentException(String.format("Matrix is unsolvable (%s)",
                                                     errorCode.name()));
  }

  /**
   * Creates luDecomposition and fills it with the scaled matrix
   * Scales and returns vector
   * @param vector
   * @return
   */
  private Quadruple[] scaleMatrixAndVector(Quadruple[] vector) {
    scaleMatrix();
    if (needToScale)
      vector = scaleVector(vector);
    return vector;
  }

  /**
   * Solves a system formed by the matrix and the given vector
   * using Cholesky decomposition
   * @param vector
   * @return the solution
   */
  private Quadruple[] solveCholeskyInternally(Quadruple[] vector) {
    if (matrix == null || vector == null || size != vector.length)
      throw new IllegalArgumentException("Matrix and vector must be of equal size");

    if (choleskyDecompositionError)
      throwCholeskyError(errorCode);

    if (choleskyDecomposition == null)
      decomposeCholesky();

    final Quadruple[] solution = vector.clone();

    // Solve L*Y = B;
    for (int i = 0; i < size; i++) {
      for (int k = 0; k < i ; k++) {
        solution[i].subtract(Quadruple.multiply(solution[k], choleskyDecomposition[i][k]));
      }
      solution[i].divide(choleskyDecomposition[i][i]);
    }

    // Solve L'*X = Y;
    for (int k = size - 1; k >= 0; k--) {
      for (int i = k + 1; i < size ; i++) {
        solution[k].subtract(Quadruple.multiply(solution[i], choleskyDecomposition[i][k]));
      }
      solution[k].divide(choleskyDecomposition[k][k]);
    }

    return solution;
  }

  /**
   * Iteratively refines the last found solution while the iteration count doesn't exceed
   * {@code iterationsToDo} and iterations continue to improve the accuracy.
   * The refined solution can be obtained via {@linkplain #getSolution()},
   * @return the number of iterations actually done
   */
  private void refine() {
    Quadruple[] bestSolution = solution;
    Quadruple[] deltaX = null;               // It will be filled by the result of the solution at the next iteration
    Quadruple bestError = Quadruple.positiveInfinity();
    double correctionFactor = 1.0;
    final int iterationsToDo = MAX_REFINEMENT_ITERATIONS;

    for (int i = 0; i < iterationsToDo; i++) {
      final Quadruple[] computedVector = multiply(solution);
      final Quadruple[] deltaB = subtractVectors(computedVector, vector);  // B1 - B0
      final Quadruple currentError = findError(deltaB);

      if (currentError.compareTo(bestError) < 0) { // The best solution so far
        bestError = currentError;
        bestSolution = solution;
        if (i >= iterationsToDo || currentError.isZero()) {
          return;
        }
      } else {
        correctionFactor /= 2;
        if (correctionFactor < MIN_CORRECTION_FACTOR || i >=  iterationsToDo) {
          solution = bestSolution;
          return;
        }
      }

      deltaX = solve(deltaB);
      solution = correctSolutionBy(solution, deltaX, correctionFactor);
    }
  }


  /**
   * Subtracts elements of the subtrahend from corresponding elements of the minuend
   * @return minuend with the new values of the elements
   */
  private static Quadruple[] subtractVectors(Quadruple[] minuend, Quadruple[] subtrahrend) {
    for (int i = 0; i < minuend.length; i++)
      minuend[i].subtract(subtrahrend[i]);
    return minuend;
  }

  private static Quadruple findError(Quadruple[] deviations) {
    final Quadruple mse = new Quadruple();
    for (final Quadruple d: deviations) {
      mse.add(Quadruple.multiply(d, d));
    }
    return mse;
  }

  private Quadruple[] solve(Quadruple[] deltaB) {
    if (lastSolvedWith == SolutionMethod.LU)
      return solveLUInternally(deltaB);
    else
      return solveCholeskyInternally(deltaB);
  }

  /**
   * Subtracts elements of the {@code errors} each multiplied by {@code correctionFactor}
   * from corresponding elements of the the previous solution
   * @return the previous solution with the corrected values
   */
  private static Quadruple[] correctSolutionBy(Quadruple[] prevSolution, Quadruple[] errors, double correctionFactor) {
    final Quadruple[] result = new Quadruple[errors.length];
    for (int i = 0; i < errors.length; i++) {
      result[i] = Quadruple.subtract(prevSolution[i], Quadruple.multiply(errors[i], correctionFactor));
    }
    return result;
  }

  private static Quadruple[] multiplyElements(Quadruple[] a, Quadruple[] b) {
    final int length = a.length;
    final Quadruple[] result = new Quadruple[length];
    for (int i = 0; i < length; i++)
      result[i] = Quadruple.multiply(a[i], b[i]);
    return result;
  }

  // Bare Kahan summation
  private static Quadruple sumOfVector(Quadruple[] vector) {
    Quadruple sum = new Quadruple();
    Quadruple c = new Quadruple();
    for (int i = 0; i < vector.length; i++) {
      final Quadruple y = Quadruple.subtract(vector[i], c);
      final Quadruple t = Quadruple.add(sum, y);
      c = Quadruple.subtract(t, sum).subtract(y);
      sum = t;
    }
    return sum;
  }

  private void scaleAndDecompose() {
    scaleMatrix(); // creates luDecomposition and a scaled copy of matrix
    decomposeLU(); // modifies luDecomposition and creates pivot.
  }

  Quadruple[][] solveAccurately(Quadruple[][] matrixB) {
    Quadruple[][] matrixX = solve(deepCopyOf(matrixB));

    Quadruple[][] deltaX = null;
    Quadruple[][] bestSolution = matrixX;
    Quadruple bestError = Quadruple.positiveInfinity();
    double correctionFactor = 1.0;
    final double min_correction_factor =  // 1.0/16;
                                    1.0/8;

    final int iterationsToDo = 20;
    for (int i = 0; i < iterationsToDo; i++) {
      final Quadruple[][] matrixB_1 = multiply(matrixX);
      final Quadruple currentError = subtractBfromA(matrixB_1, matrixB);

      if (currentError.compareTo(bestError) < 0) {
        bestError = currentError;
        bestSolution = matrixX;
        if (i >= iterationsToDo || currentError.isZero()) {
          return bestSolution;
        }
      } else {
        correctionFactor /= 2;
        if (correctionFactor < min_correction_factor || i >=  iterationsToDo) {
          return bestSolution;
        }
      }

      deltaX = solve(matrixB_1);
      matrixX = subtractMatrices(matrixX, deltaX, correctionFactor);
    }

    return matrixX;
  }

  /**
   * Subtracts subtrahend from minuend and computes the MSE of the difference D as sqrt(sum(D[i][j]^2)/N)
   * @param minuend
   * @param subtrahend
   * @return
   */
  private static Quadruple subtractBfromA(Quadruple[][] minuend, Quadruple[][] subtrahend) {
    final int size = minuend.length;
    final Quadruple mse = new Quadruple();
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        final Quadruple d = minuend[i][j].subtract(subtrahend[i][j]);
        mse.add(Quadruple.multiply(d, d));
      }
    }
    return mse.divide(size * size).sqrt();
  }

  private static Quadruple[][] subtractMatrices(Quadruple[][] minuend, Quadruple[][] subtrahend, double correctionFactor) {
    final int size = minuend.length;
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        minuend[i][j].subtract(Quadruple.multiply(subtrahend[i][j], correctionFactor));
      }
    }
    return minuend;
  }

  /* *****************************************************************************
   * Private methods used by solveLU *********************************************
   *******************************************************************************/

  private void decomposeLU() {
    pivot = new int[size];
    for (int i = 0; i < size; i++) {
      pivot [i] = i;
    }

    for (int i = 0; i < size; i++) {
      int p = i;
      for (int j = i + 1; j < size; j++) {
        if (luDecomposition[j][i].abs().compareTo(luDecomposition[p][i].abs()) > 0)
          p = j;
      }

      if (p != i) {                     // swap rows p and i
        final Quadruple[] t = luDecomposition[p]; luDecomposition[p] = luDecomposition[i]; luDecomposition[i] = t;
        final int tt = pivot[p]; pivot[p] = pivot[i]; pivot[i] = tt;
        detSignCorrection *= -1;        // Swapped rows, the determinant changes sign
      }

      final Quadruple[] row_i = luDecomposition[i];
      if (!row_i[i].isZero()) { // To compute some threshold here?
        final Quadruple row_i_i = Quadruple.one().divide(row_i[i]);

        for (int j = i + 1; j < size; j++) {
          final Quadruple[] row_j = luDecomposition[j];
          if (!row_j[i].isZero()) {
            final Quadruple row_j_i = row_j[i].multiply(row_i_i);
            for (int k = i + 1; k < size; k++) {
              row_j[k].subtract(Quadruple.multiply(row_i[k], row_j_i));
            }
          }
        }
      } else {
        throwNonInvertibleError();
      }
    }

  }

  private void decomposeCholesky() {
    choleskyDecomposition = new Quadruple[size][size];   // If we are here, decomosition == null (called from within if statement)

    for (int i = 0; i < size; i++) {
      final Quadruple sum2 = new Quadruple();
      for (int j = 0; j < i; j++) {
        if (matrix[i][j].compareTo(matrix[j][i]) != 0)
          throwCholeskyError(ErrorCodes.ASYMMETRIC);

        final Quadruple sum = new Quadruple();
        for (int k = 0; k < j; k++)
          sum.add(Quadruple.multiply(choleskyDecomposition[i][k], choleskyDecomposition[j][k]));

        final Quadruple s = choleskyDecomposition[i][j] = Quadruple.subtract(matrix[i][j], sum).divide(choleskyDecomposition[j][j]);
        sum2.add(Quadruple.multiply(s, s));
      }

      final Quadruple dd = Quadruple.subtract(matrix[i][i], sum2);
      if (dd.isNegative() || dd.isInfinite() || dd.isNaN()) {
        throwCholeskyError(ErrorCodes.NON_SPD);
      }

      choleskyDecomposition[i][i] = dd.sqrt();
    }
  }

  /**
   */
  private void throwCholeskyError(ErrorCodes errorCode) {
    this.errorCode = errorCode;
    choleskyDecompositionError = true;
    throw new IllegalArgumentException(String.format("holeskyDecompose: matrix must be %s",
            errorCode == ErrorCodes.NON_SPD ? "positively defined" : "symmetric"));
  }


  /**
   * Puts the scaled matrix to the luDecomposition array.
   * Direct or indirect call to it is always followed by call to decomposeLU(),
   * so that luDecomposition is always either null or contains the LU-decomposition
   */
  private void scaleMatrix() {
    rowScales = new Quadruple[size];
    luDecomposition = new Quadruple[size][size];
    if (needToScale) {
      for (int i = 0; i < size; i++) {
        final Quadruple rowSum = sumOfAbsValues(matrix[i]);
        final Quadruple scale = (rowSum.compareTo(new Quadruple()) > 0) ?
            Quadruple.one().divide(rowSum) :
            Quadruple.one();
        rowScales[i] = scale;
        for (int j = 0; j < size; j++)
          luDecomposition[i][j] = Quadruple.multiply(matrix[i][j], scale);
      }
    } else {
      for (int i = 0; i < size; i++) {
        rowScales[i] = Quadruple.one();
        luDecomposition[i] = matrix[i].clone() ;
      }
    }
  }

  /** Attention -- spoils the matrix, returns a new array */
  private Quadruple[][] scaleAndPermute(Quadruple[][] matrix, int[] indices) {
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++)
        matrix[i][j].multiply(rowScales[i]);
    final Quadruple[][] result = new Quadruple[size][size];
    for (int i = 0; i < size; i++)
      result[i] = matrix[indices[i]];
    return result;
  }

  private Quadruple[] scaleVector(Quadruple[] vector) {
    for (int i = 0; i < size; i++)
      vector[i].multiply(rowScales[i]);
    return vector;
  }

  private Quadruple sumOfAbsValues(Quadruple[] vector) {
    final Quadruple[] absValues = new Quadruple[size];
    for (int i = 0; i < size; i++)
      absValues[i] = vector[i].abs();
    return sumOfVector(absValues);
  }

  private static Quadruple[] getPermutted(Quadruple[] vector, int[] indices) {
    final Quadruple[] result = new Quadruple[vector.length];
    for (int i = 0; i < vector.length; i++)
      result[i] = vector[indices[i]];
    return result;
  }

  private Quadruple computeDeterminant() {
    if (luDecomposition == null) { // The very first invocation.
      scaleAndDecompose();
    }

    final Quadruple result = new Quadruple(1.0);
    final Quadruple scale = new Quadruple(1.0);;

    if (needToScale) {
      for (int i = 0; i < size; i++) {
        result.multiply(luDecomposition[i][i]);
        scale.multiply(rowScales[i]);
      }
      result.divide(scale);
    } else {
      for (int i = 0; i < size; i++) {
        result.multiply(luDecomposition[i][i]);
      }
    }
    return result.multiply(detSignCorrection);
  }

  private Quadruple computeNorm(Quadruple[][] matrix) { // May be used to find the norm of the inversion
    Quadruple result = new Quadruple();
    for (int i = 0; i < size; i++) {
      final Quadruple rowNorm = new Quadruple();
      for (int j = 0; j < size; j++) {
        rowNorm.add(matrix[i][j].abs());
      }
      result = Quadruple.max(result, rowNorm);
    }
    return result;
  }

}
