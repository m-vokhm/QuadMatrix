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
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

import static com.mvohm.quadmatrix.Util.*;

class BigDecimalMatrixSolver extends MatrixSolver {

  // TODO_for_next_release 2023-05-27 11:10:12 Reorder methods accordingly with their usage

  enum SolutionMethod {NONE, CHOLESKY, LU};

  private final int size;
  private final BigDecimal[][] matrix;  // Original matrix
  private BigDecimal[] vector;          // Original vector, used for refinement

  private SolutionMethod lastSolvedWith = SolutionMethod.NONE;
  private ErrorCodes errorCode = ErrorCodes.OK;

  private BigDecimal[][] luDecomposition;
  private boolean luDecompositionError;

  private BigDecimal[][] choleskyDecomposition;
  private boolean choleskyDecompositionError;

  private int[] pivot;              // filled by decomposeLU(), used by solveLU()
  private BigDecimal[] rowScales;
  private BigDecimal[] solution;

  private BigDecimal[][] inversion = null;

  private int detSignCorrection = 1;            // To correct the determinant sign in case of swapping rows while pivoting
  private BigDecimal determinant = null;      // Unless already computed
  private BigDecimal norm;

  private final MathContext mc;
  private final MathContext mc2;

  BigDecimalMatrixSolver(BigDecimal[][] matrix, boolean needToScale, MathContext mc) {
    this.size = matrix.length;
    this.matrix = matrix;
    this.needToScale = needToScale;
    this.mc = mc;
    this.mc2 = new MathContext(mc.getPrecision() + 10, RoundingMode.HALF_EVEN); // For cases where the error can accumulate
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
  BigDecimal[] solveLU(BigDecimal[] vector) {
    this.vector = deepCopyOf(vector);
    solution = solveLUInternally(vector);
    lastSolvedWith = SolutionMethod.LU;
    return solution.clone();                  // To protect the internally stored solutions for possible subsequent retrievals
  }

  BigDecimal[] solveLUAccurately(BigDecimal[] vector) {
    this.vector = deepCopyOf(vector);
    solution = solveLUInternally(vector);
    lastSolvedWith = SolutionMethod.LU;
    refine();
    return solution.clone();
  }

  BigDecimal[] solveCholesky(BigDecimal[] vector) {
    this.vector = deepCopyOf(vector);
    solution = solveCholeskyInternally(vector);
    lastSolvedWith = SolutionMethod.CHOLESKY;
    return solution.clone();                  // To protect the internally stored solutions for possible subsequent retrievals
  }

  BigDecimal[] solveCholeskyAccurately(BigDecimal[] vector) {
    this.vector = deepCopyOf(vector);
    solution = solveCholeskyInternally(vector);
    lastSolvedWith = SolutionMethod.CHOLESKY;
    refine();
    return solution.clone();
  }

  BigDecimal[] getSolution() {
    return solution;
  }

  String errorCode() {
    return errorCode.toString();
  }

  /* ********************************************************************************
  /*** Solutions with respect of matrices *******************************************
  /* ********************************************************************************/

  /** Attention -- spoils the matrixB, but returns a newly-created array */
  BigDecimal[][] solve(BigDecimal[][] matrixB) {
    if (luDecompositionError)
      throwNonInvertibleError();

    if (luDecomposition == null) { // The very first invocation.
      scaleAndDecompose();
    }

    // Copy right hand side with pivoting and scaling
    final BigDecimal[][] matrixX = scaleAndPermute(matrixB, pivot); // Spoils matrixB, and returns a new array

    // Solve L*Y = B(piv,:)
    for (int k = 0; k < size; k++) {
      for (int i = k + 1; i < size; i++) {
        for (int j = 0; j < size; j++) {
          matrixX[i][j] = matrixX[i][j].subtract(matrixX[k][j].multiply(luDecomposition[i][k], mc), mc);
        }
      }
    }

    // Solve U * X = Y;
    for (int k = size - 1; k >= 0; k--) {
      for (int j = 0; j < size; j++) {
        matrixX[k][j] = matrixX[k][j].divide(luDecomposition[k][k], mc);
      }
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < size; j++) {
          matrixX[i][j] = matrixX[i][j].subtract(matrixX[k][j].multiply(luDecomposition[i][k], mc), mc);
        }
      }
    }

    return matrixX;
  }

  /* ********************************************************************************
  /*** Inverse **********************************************************************
  /* ********************************************************************************/

  BigDecimal[][] inverse() {
    if (inversion == null)
      inversion = solve(unityMatrix());
    return inversion;
  }

  BigDecimal[][] inverseAccurately() {
    return solveAccurately(unityMatrix());
  }

  /**
   * Creates a new two-dimentional double array and fills it with the transposition of the source array
   * @param matrix
   * @return
   */
  BigDecimal[][] transpose(BigDecimal[][] matrix) {
    final int size = matrix.length;
    final BigDecimal[][] result = new BigDecimal[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        result[i][j] = matrix[j][i];
      }
    }
    return result;
  }

  BigDecimal[][] unityMatrix() {
    final BigDecimal[][] result = new BigDecimal[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        result[i][j] = BigDecimal.ZERO;
      }
      result[i][i] = BigDecimal.ONE;
    }
    return result;
  }

  /* *******************************************************************************
  /* **** Multiplications **********************************************************
  /*********************************************************************************/


  BigDecimal[][] multiply(BigDecimal[][] factor) {
    final BigDecimal[][] product = new BigDecimal[size][size];
    final BigDecimal[] prodVector = new BigDecimal[size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        for (int k = 0; k < size; k++) {
          prodVector[k] = matrix[i][k].multiply(factor[k][j], mc);
        }
        product[i][j] = sumOfVector(prodVector);
      }
    }
    return product;
  }

  BigDecimal[] multiply(BigDecimal[] factor) {
    final int size = factor.length;
    final BigDecimal[] product = new BigDecimal[size];
    for (int i = 0; i < size; i++) {
      final BigDecimal[] rowProduct = multiplyElements(matrix[i], factor);
      product[i] = sumOfVector(rowProduct);
    }
    return product;
  }

  BigDecimal[][] multiply(BigDecimal scalar) {
    final BigDecimal[][] product = new BigDecimal[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        product[i][j] = scalar.multiply(matrix[i][j], mc);
      }
    }
    return product;
  }

  /* *******************************************************************************
  /* **** Additions and subtractions ***********************************************
  /*********************************************************************************/

  BigDecimal[][] addMatrices(BigDecimal[][] summandA, BigDecimal[][] summandB) {
    final int size = summandA.length;
    final BigDecimal[][] sum = new BigDecimal[size][size];
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++) {
        sum[i][j] = summandA[i][j].add(summandB[i][j], mc);
      }
    return sum;
  }

  BigDecimal[][] subtractMatrices(BigDecimal[][] minuend, BigDecimal[][] subtrahend) {
    final int size = minuend.length;
    final BigDecimal[][] difference = new BigDecimal[size][size];
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++) {
        difference[i][j] = minuend[i][j].subtract(subtrahend[i][j], mc);
      }
    return difference;
  }

  BigDecimal determinant() {
    if (determinant == null) {
      if (luDecompositionError) {
        determinant = BigDecimal.ZERO;
      } else {
        try {
          determinant = computeDeterminant();
        } catch (final IllegalArgumentException x) { // Non-invertible
          determinant = BigDecimal.ZERO;
        }
      }
    }
    return determinant;
  }

  BigDecimal norm() {
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
  private BigDecimal[] solveLUInternally(BigDecimal[] vector) {
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

    final BigDecimal[] solution = getPermutted(vector, pivot);     // returns a new array

    // Solve L * Y = B(piv,:)
    for (int k = 0; k < size; k++)
     for (int i = k + 1; i < size; i++) {
       solution[i] = solution[i].subtract(solution[k].multiply(luDecomposition[i][k], mc), mc);    // (N^2 - N) / 2
     }

    // Solve U * X = Y;
    for (int k = size - 1; k >= 0; k--) {
      solution[k] = solution[k].divide(luDecomposition[k][k], mc);                  // N
      for (int i = 0; i < k; i++) {
        solution[i] = solution[i].subtract(solution[k].multiply(luDecomposition[i][k], mc), mc);    // (N^2 - N) / 2
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
  private BigDecimal[] scaleMatrixAndVector(BigDecimal[] vector) {
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
  private BigDecimal[] solveCholeskyInternally(BigDecimal[] vector) {
    if (matrix == null || vector == null || size != vector.length)
      throw new IllegalArgumentException("Matrix and vector must be of equal size");

    if (choleskyDecompositionError)
      throwCholeskyError(errorCode);

    if (choleskyDecomposition == null)
      decomposeCholesky();

    final BigDecimal[] solution = vector.clone();

    // Solve L*Y = B;
    for (int i = 0; i < size; i++) {
      for (int k = 0; k < i ; k++) {
        solution[i] = solution[i].subtract(solution[k].multiply(choleskyDecomposition[i][k], mc), mc);
      }
      solution[i] = solution[i].divide(choleskyDecomposition[i][i], mc);
    }

    // Solve L'*X = Y;
    for (int k = size - 1; k >= 0; k--) {
      for (int i = k + 1; i < size ; i++) {
        solution[k] = solution[k].subtract(solution[i].multiply(choleskyDecomposition[i][k], mc), mc);
      }
      solution[k] = solution[k].divide(choleskyDecomposition[k][k], mc);
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
    BigDecimal[] bestSolution = solution;
    BigDecimal[] deltaX = null;               // It will be filled by the result of the solution at the next iteration
    BigDecimal bestError = BigDecimal.valueOf(Double.MAX_VALUE);
    double correctionFactor = 1.0;
    final int iterationsToDo = MAX_REFINEMENT_ITERATIONS;

    for (int i = 0; i < iterationsToDo; i++) {
      final BigDecimal[] computedVector = multiply(solution);
      final BigDecimal[] deltaB = subtractVectors(computedVector, vector);  // B1 - B0
      final BigDecimal currentError = findError(deltaB);

      if (currentError.compareTo(bestError) < 0) { // The best solution so far
        bestError = currentError;
        bestSolution = solution;
        if (i >= iterationsToDo || (currentError.compareTo(BigDecimal.ZERO) == 0)) {
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
  private BigDecimal[] subtractVectors(BigDecimal[] minuend, BigDecimal[] subtrahrend) {
    for (int i = 0; i < minuend.length; i++)
      minuend[i] = minuend[i].subtract(subtrahrend[i], mc2);
    return minuend;
  }

  private BigDecimal findError(BigDecimal[] deviations) {
    BigDecimal mse = BigDecimal.ZERO;
    for (final BigDecimal d: deviations) {
      mse = mse.add(d.multiply(d, mc), mc);
    }
    return mse;
  }

  private BigDecimal[] solve(BigDecimal[] deltaB) {
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
  private BigDecimal[] correctSolutionBy(BigDecimal[] prevSolution, BigDecimal[] errors, double correctionFactor) {
    final BigDecimal[] result = new BigDecimal[errors.length];
    for (int i = 0; i < errors.length; i++) {
      result[i] = prevSolution[i].subtract(errors[i].multiply(BigDecimal.valueOf(correctionFactor), mc2), mc2);
    }
    return result;
  }

  private BigDecimal[] multiplyElements(BigDecimal[] a, BigDecimal[] b) {
    final int length = a.length;
    final BigDecimal[] result = new BigDecimal[length];
    for (int i = 0; i < length; i++)
      result[i] = a[i].multiply(b[i], mc);
    return result;
  }

  // Bare Kahan summation
  private BigDecimal sumOfVector(BigDecimal[] vector) {
    BigDecimal sum = BigDecimal.ZERO;
    BigDecimal c = BigDecimal.ZERO;
    for (int i = 0; i < vector.length; i++) {
      final BigDecimal y = vector[i].subtract(c, mc2);
      final BigDecimal t = sum.add(y, mc2);
      c = t.subtract(sum, mc2).subtract(y, mc2);
      sum = t;
    }
    return sum;
  }

  private void scaleAndDecompose() {
    scaleMatrix(); // creates luDecomposition and a scaled copy of matrix
    decomposeLU(); // modifies luDecomposition and creates pivot.
  }

  BigDecimal[][] solveAccurately(BigDecimal[][] matrixB) {
    BigDecimal[][] matrixX = solve(deepCopyOf(matrixB));

    BigDecimal[][] deltaX = null;
    BigDecimal[][] bestSolution = matrixX;
    BigDecimal bestError = BigDecimal.valueOf(Double.MAX_VALUE);
    double correctionFactor = 1.0;
    final double min_correction_factor =  // 1.0/16;
                                    1.0/8;

    final int iterationsToDo = 20;
    for (int i = 0; i < iterationsToDo; i++) {
      final BigDecimal[][] matrixB_1 = multiply(matrixX);
      final BigDecimal currentError = subtractBfromA(matrixB_1, matrixB);

      if (currentError.compareTo(bestError) < 0) {
        bestError = currentError;
        bestSolution = matrixX;
        if (i >= iterationsToDo || currentError.doubleValue() == 0) {
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
  private BigDecimal subtractBfromA(BigDecimal[][] minuend, BigDecimal[][] subtrahend) {
    final int size = minuend.length;
    BigDecimal mse = BigDecimal.ZERO;
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        final BigDecimal d = minuend[i][j] = minuend[i][j].subtract(subtrahend[i][j], mc2);
        mse = mse.add(d.multiply(d, mc2), mc2);
      }
    }
    return BigDecimal.valueOf(Math.sqrt(mse.doubleValue() / (size * size)));
  }

  private BigDecimal[][] subtractMatrices(BigDecimal[][] minuend, BigDecimal[][] subtrahend, double correctionFactor) {
    final int size = minuend.length;
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        minuend[i][j] = minuend[i][j].subtract(subtrahend[i][j].multiply(new BigDecimal(correctionFactor), mc), mc);
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
        final BigDecimal[] t = luDecomposition[p]; luDecomposition[p] = luDecomposition[i]; luDecomposition[i] = t;
        final int tt = pivot[p]; pivot[p] = pivot[i]; pivot[i] = tt;
        detSignCorrection *= -1;        // Swapped rows, the determinant changes its sign
      }

      final BigDecimal[] row_i = luDecomposition[i];
      if (row_i[i].doubleValue() != 0) {                // Do we better compute some threshold here?
                                                        // Some practically unsolvable matrices may have non-zero here
        final BigDecimal row_i_i = BigDecimal.ONE.divide(row_i[i], mc);

        for (int j = i + 1; j < size; j++) {
          final BigDecimal[] row_j = luDecomposition[j];
          if (row_j[i].doubleValue() != 0) {
            final BigDecimal row_j_i = row_j[i] = row_j[i].multiply(row_i_i, mc); // !!!
            for (int k = i + 1; k < size; k++) {
              row_j[k] = row_j[k].subtract(row_i[k].multiply(row_j_i, mc), mc);
            }
          }
        }
      } else {
        throwNonInvertibleError();
      }
    }

  }

  private void decomposeCholesky() {
    choleskyDecomposition = new BigDecimal[size][size];   // If we are here, decomosition == null (called from within if statement)

    for (int i = 0; i < size; i++) {
      BigDecimal sum2 = BigDecimal.ZERO;
      for (int j = 0; j < i; j++) {
        if (matrix[i][j].compareTo(matrix[j][i]) != 0)
          throwCholeskyError(ErrorCodes.ASYMMETRIC);

        BigDecimal sum = BigDecimal.ZERO;
        for (int k = 0; k < j; k++)
          sum = sum.add(choleskyDecomposition[i][k].multiply(choleskyDecomposition[j][k], mc), mc);

        final BigDecimal s = choleskyDecomposition[i][j] = matrix[i][j].subtract(sum, mc).divide(choleskyDecomposition[j][j], mc);
        sum2 = sum2.add(s.multiply(s, mc), mc);
      }

      final BigDecimal dd = matrix[i][i].subtract(sum2, mc);
      if (dd.signum() < 0 || Double.isInfinite(dd.doubleValue())) {
        throwCholeskyError(ErrorCodes.NON_SPD);
      }

      choleskyDecomposition[i][i] = dd.sqrt(mc);
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
    rowScales = new BigDecimal[size];
    luDecomposition = new BigDecimal[size][size];
    if (needToScale) {
      for (int i = 0; i < size; i++) {
        final BigDecimal rowSum = sumOfAbsValues(matrix[i]);
        final BigDecimal scale = (rowSum.signum() > 0) ?
            BigDecimal.ONE.divide(rowSum, mc) :
            BigDecimal.ONE;
        rowScales[i] = scale;
        for (int j = 0; j < size; j++)
          luDecomposition[i][j] = matrix[i][j].multiply(scale, mc);
      }
    } else {
      for (int i = 0; i < size; i++) {
        rowScales[i] = BigDecimal.ONE;
        luDecomposition[i] = matrix[i].clone() ;
      }
    }
  }

  /** Attention -- spoils the matrix, returns a new array */
  private BigDecimal[][] scaleAndPermute(BigDecimal[][] matrix, int[] indices) {
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++)
        matrix[i][j] = matrix[i][j].multiply(rowScales[i], mc2);
    final BigDecimal[][] result = new BigDecimal[size][size];
    for (int i = 0; i < size; i++)
      result[i] = matrix[indices[i]];
    return result;
  }

  private BigDecimal[] scaleVector(BigDecimal[] vector) {
    for (int i = 0; i < size; i++)
      vector[i] = vector[i].multiply(rowScales[i], mc2);
    return vector;
  }

  private BigDecimal sumOfAbsValues(BigDecimal[] vector) {
    final BigDecimal[] absValues = new BigDecimal[size];
    for (int i = 0; i < size; i++)
      absValues[i] = vector[i].abs();
    return sumOfVector(absValues);
  }

  private static BigDecimal[] getPermutted(BigDecimal[] vector, int[] indices) {
    final BigDecimal[] result = new BigDecimal[vector.length];
    for (int i = 0; i < vector.length; i++)
      result[i] = vector[indices[i]];
    return result;
  }

  private BigDecimal computeDeterminant() {
    if (luDecomposition == null) { // The very first invocation.
      scaleAndDecompose();
    }

    BigDecimal result = new BigDecimal(1.0);
    BigDecimal scale = new BigDecimal(1.0);;

    if (needToScale) {
      for (int i = 0; i < size; i++) {
        result = result.multiply(luDecomposition[i][i], mc);
        scale = scale.multiply(rowScales[i], mc);
      }
      result = result.divide(scale, mc);
    } else {
      for (int i = 0; i < size; i++) {
        result = result.multiply(luDecomposition[i][i], mc);
      }
    }
    return result.multiply(new BigDecimal(detSignCorrection));
  }

  private BigDecimal computeNorm(BigDecimal[][] matrix) { // May be used to find the norm of the inversion
    BigDecimal result = BigDecimal.ZERO;
    for (int i = 0; i < size; i++) {
      BigDecimal rowNorm = BigDecimal.ZERO;;
      for (int j = 0; j < size; j++) {
        rowNorm = rowNorm.add(matrix[i][j].abs(), mc);
      }
      result = result.max(rowNorm);
    }
    return result;
  }

}
