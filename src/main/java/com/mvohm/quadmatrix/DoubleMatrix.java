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

import java.math.BigDecimal;
import java.util.Arrays;

import com.mvohm.quadruple.Quadruple;

import static com.mvohm.quadmatrix.DoubleMatrixSolver.*;
import static com.mvohm.quadmatrix.Util.*;

/**
 * An implementation of the abstract class {@link Matrix} that uses an array of primitive
 * {@code double} values to store the matrix elements and standard {@code double} arithmetic
 * to perform calculations.
 * <br>
 * @author M.Vokhmentsev
 */
public class DoubleMatrix extends Matrix {

  /* *******************************************************************************
  /***** Internal data ****************************************************************
  /*********************************************************************************/

  private double[][] matrix;

  private double[] solution;
  private double[][] matrixSolution;

  private DoubleMatrixSolver solver;

  /**
   * A default value of the {@code needToScale} flag that is set in new instances by
   * constructors that don't have the corresponding parameter.
   * <br>This flag determines the need to scale the rows of the matrix before LU-decomposing it
   * so that the values of their norms are close enough to each other. This may increase the accuracy
   * of the solution and the inversion, especially for matrices with a significantly non-uniform
   * distribution of element values, at the cost of a insignificant increase in computation time.
   * The initial default value is {@code true} but it can be changed using {@link #setDefaultScaling(boolean)}
   * static method.
   * @see #setDefaultScaling(boolean)
   * @see #getDefaultScaling()
   * @see #DoubleMatrix(Matrix, boolean)
   * @see #DoubleMatrix(double[][], boolean)
   * @see #DoubleMatrix(Number[][], boolean)
   */
  private static boolean scaleByDefault = true;

  /* *******************************************************************************
  /***** Constructors **************************************************************
  /*********************************************************************************/

  /**
   * Creates a new {@code DoubleMatrix} with a copy of the data of the given {@code source} matrix and
   * the default value of the {@code needToScale} flag.
   * <br>Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag designating
   * necessity to scale matrix rows before solving equation systems and inverting.<br>
   * @param source a matrix whose data get copied into the new Matrix's internal storage.
   * @see #setDefaultScaling(boolean)
   */
  public DoubleMatrix(Matrix source) {
    this(source, scaleByDefault);
  }

  /**
   * Creates a new {@code DoubleMatrix} with a copy of the data of the given {@code source} matrix and
   * the value of the {@code needToScale} flag that is passed in as the {@code needToScale} argument.
   * <br>Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag designating
   * necessity to scale matrix rows while solving equation systems and inverting.
   * The value of this flag can be obtained via {@link #getScaling()} instance method.
   * @param source a matrix whose data get copied into the new Matrix's internal storage.
   * @param needToScale the value defining necessity to scale rows while solving equation systems and inverting.
   * @see #getScaling()
   */
  public DoubleMatrix(Matrix source, boolean needToScale) {
    checkMatrixData(source);
    this.size = source.size;
    this.matrix = source.getDoubleData();
    solver = makeSolver(needToScale);
  }

  /**
   * Creates a new {@code DoubleMatrix} with a copy of the data of the given {@code source} array and
   * the default value of the {@code needToScale} flag.
   * <br>Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag designating
   * necessity to scale matrix rows before solving equation systems and inverting.
   *
   * <br>Throws {@linkplain IllegalArgumentException} if the {@code source} argument is a non-square
   * or empty array, or contains non-numeric values, {@link Double#NaN}, {@link Double#NEGATIVE_INFINITY},
   * or {@link Double#POSITIVE_INFINITY}.
   *
   * @param source a two-dimentional square array of {@code double}s
   *    to be copied into the internal storage of this instance.
   * @see #setDefaultScaling(boolean)
   * @see #scaleByDefault
   *
   */
  public DoubleMatrix(double[][] source) {
    this(source, scaleByDefault);
  }

  /**
   * Creates a new {@code DoubleMatrix} with a copy of the data of the given {@code source} array and
   * the value of the {@code needToScale} flag that is passed in as the {@code needToScale} argument.
   * <br>Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag designating
   * necessity to scale matrix rows while solving equation systems and inverting.
   * The value of this flag can be obtained via {@link #getScaling()} instance method.
   *
   * <br>Throws {@linkplain IllegalArgumentException} if the {@code source} argument is a non-square or empty array,
   * or contains non-numeric values, {@link Double#NaN}, {@link Double#NEGATIVE_INFINITY},
   * or {@link Double#POSITIVE_INFINITY}.
   *
   * @param source a two-dimentional square array of {@code double}s
   *    to be copied into the internal storage of this instance.
   * @param needToScale
   *   the value defining necessity to scale rows while solving equation
   *   systems and inverting.
   * @see #getScaling()
   */
  public DoubleMatrix(double[][] source, boolean needToScale) {
    checkMatrixData(source);
    copyDataFrom(source);
    solver = makeSolver(needToScale);
  }

  /**
   * Creates a new {@code DoubleMatrix} with a copy of the data of the given {@code source} array and
   * the default value of the {@code needToScale} flag.
   * <br>Translates the values of the given {@code source} array into {@code double} values
   * using their {@linkplain Number#doubleValue()} methods.
   * <br>Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag designating
   * necessity to scale matrix rows before solving equation systems and inverting.<br>
   *
   * <br>Throws {@linkplain IllegalArgumentException} if the {@code source} argument is a non-square or empty array,
   * or contains non-numeric values, {@link Double#NaN}, {@link Double#NEGATIVE_INFINITY}
   * or {@link Double#POSITIVE_INFINITY}), or {@code null}.
   *
   * @param source a two-dimentional square array of {@linkplain Number} values
   *    to be copied into the internal storage of this instance.
   * @see #setDefaultScaling(boolean)
   */
  public DoubleMatrix(Number[][] source) {
    this(source, scaleByDefault);
  }

  /**
   * Creates a new {@code DoubleMatrix} with a copy of the data of the given {@code source} array and
   * the value of the {@code needToScale} flag that is passed in as the {@code needToScale} argument.
   * <br>Translates the values of the given {@code source} array into {@code double} values
   * using their {@linkplain Number#doubleValue()} methods.
   * <br>Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag designating
   * necessity to scale matrix rows while solving equation systems and inverting.
   * The value of this flag can be obtained via {@link #getScaling()} instance method.
   *
   * <br>Throws {@linkplain IllegalArgumentException} if the {@code source} argument is a non-square or empty array,
   * or contains non-numeric values, {@link Double#NaN}, {@link Double#NEGATIVE_INFINITY}
   * or {@link Double#POSITIVE_INFINITY}), or {@code null,} .
   *
   * @param source a two-dimentional square array of {@linkplain Number} values
   *    to be copied into the internal storage of this instance.
   * @param needToScale the flag defining necessity to scale rows while solving equation
   *   systems and inverting.
   * @see #setDefaultScaling(boolean)
   */
  public DoubleMatrix(Number[][] source, boolean needToScale) {
    checkMatrixData(source);
    copyDataFrom(source);
    solver = makeSolver(needToScale);
  }

  //***** Private constructor and fabric method creating a new instance with the same data array *************

  /** Private constructor that leaves the object fields blank. Used by and factory method below
   * to make a new Matrix with the same instance of the given data array without allocating a new array
   * to avoid unnecessary memory fragmentation and GC load */
  private DoubleMatrix() {
    super();
  }

  /** A factory method creating a new Matrix
  * with the same instance of the given data array, without allocating a new array.
  * Used by methods returning new DoubleMatrix instances to avoid unnecessary memory fragmentation and GC load */
  private static DoubleMatrix newWithThisArray(double[][] array) {
    final DoubleMatrix matrix = new DoubleMatrix();
    matrix.size = array.length;
    matrix.matrix = array;
    matrix.solver = new DoubleMatrixSolver(array, scaleByDefault);
    return matrix;
  }

  //***** Static methods *************

  /**
   * Sets the value of the static {@link #scaleByDefault} flag that is used to set corresponding
   * instance flags when creating a new instances by constructors without {@code needToScale} parameter.
   * <br>This flag determines the need to scale the rows of the matrix before LU-decomposing it
   * so that the values of their norms are close enough to each other. This may increase the accuracy
   * of the solution and the inversion, especially for matrices with a significantly non-uniform
   * distribution of element values, at the cost of a insignificant increase in computation time.
   * @param scaleByDefault the value of {@code needToScale} flag to set for new matrices
   * created with constructors without {@code needToScale} parameter.
   */
  public static void setDefaultScaling(boolean scaleByDefault) {
    DoubleMatrix.scaleByDefault = scaleByDefault;
  };

  /**
   * Returns the value of the static {@code #scaleByDefault} flag that is used to set corresponding
   * instance flags when creating a new instances by constructors without {@code needToScale} parameter.
   * @see #setDefaultScaling(boolean)
   * @return the value of static {@code scaleByDefault} flag
   */
  public static boolean getDefaultScaling() {
    return scaleByDefault;
  };

  /* *******************************************************************************
  /***** Getting data **************************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean getScaling() {
    return solver.getScaling();
  }

  /**
   * {@inheritDoc}
   * <br>DoubleMatrix returns a copy of the internal matrix data as a two-dimentional
   * array of {@link Double}.
   */
  @Override
  public Number[][] getData() {
    return convertToDoublesAsNumbers(matrix);
  };

  /**
   * {@inheritDoc}
   * <br>DoubleMatrix returns a copy of the internal matrix data containing exact values
   * of the corresponding elements of the internal storage.
   * @return a two-dimentional array of primitive <b>{@code double}</b> values containing exact
   * copy of the internal data array */
  @Override
  public double[][] getDoubleData() {
    return deepCopyOf(matrix);
  };

  /**
   * {@inheritDoc}
   * <br>DoubleMatrix returns an array containing the exact values of the corresponding elements of the internal storage,
   * expanded to {@linkplain Quadruple} values, without precision loss.
   */
  @Override
  public Quadruple[][] getQuadrupleData() {
    final Quadruple[][] result = new Quadruple[size][size];
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++)
        result[i][j] = new Quadruple(matrix[i][j]);
    return result;
  };

  /**
   * {@inheritDoc}
   * <br>DoubleMatrix returns an array containing the values of the corresponding elements of the internal storage,
   * translated to {@linkplain BigDecimal} values using {@linkplain BigDecimal#valueOf(double)}.
   */
  @Override
  public BigDecimal[][] getBigDecimalData() {
    final BigDecimal[][] result = new BigDecimal[size][size];
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++)
        result[i][j] = BigDecimal.valueOf(matrix[i][j]);
    return result;
  };

  /** {@inheritDoc} */
  @Override
  public boolean equals(Object anotherOne) {
    if (this == anotherOne)
      return true;
    if (!(anotherOne instanceof DoubleMatrix))
      return false;
    return  (Arrays.deepEquals(matrix, ((DoubleMatrix)anotherOne).matrix))
            && (solver.needToScale == ((DoubleMatrix)anotherOne).solver.needToScale);
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    int result = Arrays.deepHashCode(matrix);
    return result = result * 31 + Boolean.hashCode(solver.needToScale);
  };

  /* *******************************************************************************
  /***** Solutions with respect of vectors *****************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code DoubleMatrix} is {@code Double}.
   */
  @Override
  public Number[] solve(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solve(double[])");
    solution = solver.solveLU(vector);
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * Before solving the equation, the values of the {@code vector} argument get rounded to the nearest
   * {@code double} values.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code DoubleMatrix} is {@code Double}.
   */
  @Override
  public Number[] solve(Number[] vector) throws IllegalArgumentException, NullPointerException{
    checkVector(vector, "solve(Number[])");
    return solve(convertToDoubles(vector));
  }

  /**
   * {@inheritDoc}
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code DoubleMatrix} is {@code Double}.
   */
  @Override
  public Number[] solveSPD(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPD(double[])");
    solution = solver.solveCholesky(vector);
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * Before solving the equation, the values of the {@code vector} argument get rounded to the nearest
   * {@code double} values.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code DoubleMatrix} is {@code Double}.
   */
  @Override
  public Number[] solveSPD(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPD(Number[])");
    return solveSPD(convertToDoubles(vector));
  }

  /**
   * {@inheritDoc}
   * <br>
   * The execution time is about twice longer than that of the simple
   * {@linkplain DoubleMatrix#solve(double[])}.
   * In a typical case, the result error is reduced, compared to the simple {@code solve()},
   * by a factor from several units to several tens or a few hundreds, depending on the properties
   * of the specific matrix data.<br>
   * For matrices 500 x 500 with uniformly-distributed random values, average square root of MSE
   * is reduced by approximately 50 times, from about 2.0e-13 down to 4.0e-15.<br>
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code DoubleMatrix} is {@code Double}.
   */
  @Override
  public Number[] solveAccurately(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveAccurately(double[])");
    solution = solver.solveLUAccurately(vector);
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get rounded to the nearest {@code double} values.
   * <br>
   * The execution time is about twice longer than that of the simple
   * {@linkplain DoubleMatrix#solve(double[])}.
   * In a typical case, the result error is reduced, compared to the simple {@code solve()},
   * by a factor from several units to several tens or a few hundreds, depending on the properties
   * of the specific matrix data.<br>
   * For matrices 500 x 500 with uniformly-distributed random values, average square root of MSE
   * is reduced by approximately 50 times, from about 2.0e-13 down to 4.0e-15.<br>
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code DoubleMatrix} is {@code Double}.
   */
  @Override
  public Number[] solveAccurately(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveAccurately(Number[])");
    return solveAccurately(convertToDoubles(vector));
  }

  /**
   * {@inheritDoc}
   * <br>
   * The execution time is about 70% longer than that of the simple
   * {@linkplain DoubleMatrix#solveSPD(double[])}.
   * In a typical case, the result error is reduced, compared to the simple solveSPD(),
   * by a factor of 6-8 times, depending on specific values of the matrix elements.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code DoubleMatrix} is {@code Double}.
   */
  @Override
  public Number[] solveSPDAccurately(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPDAccurately(double[])");
    solution = solver.solveCholeskyAccurately(vector);
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get rounded to the nearest {@code double} values.
   * <br>
   * The execution time is about 70% longer than that of the simple
   * {@linkplain DoubleMatrix#solveSPD(double[])}.
   * In a typical case, the result error is reduced, compared to the simple solveSPD(),
   * by a factor of 6-8 times, depending on specific values of the matrix elements.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code DoubleMatrix} is {@code Double}.
   */
  @Override
  public Number[] solveSPDAccurately(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPDAccurately(Number[])");
    return solveSPDAccurately(convertToDoubles(vector));
  }

  /**
   * {@inheritDoc}
   * @return the last of the previously found vector solutions, <b>x</b>, to a system of of linear equations
   *  of form <b><i>A</i>x = b</b>, as an array of {@link Double}s, or {@code null},
   *  if no system was solved with the matrix.
   */
  @Override
  public Number[] getSolution() {
    return convertToDoublesAsNumbers(solution);
  }

  /**
   * {@inheritDoc}
   * <br>{@code DoubleMatrix} returns the exact values that was found by the corresponding solution method.
   */
  @Override
  public double[] getDoubleSolution() {
    if (solution == null) return null;
    return solution.clone();
  }

  /**
   * {@inheritDoc}
   * <br>{@code DoubleMatrix} translates double values of the internally stored solution to exactly equal {@link Quadruple} values.
   */
  @Override
  public Quadruple[] getQuadrupleSolution() {
    return convertToQuadruples(solution);
  }

  /**
   * {@inheritDoc}
   * <br>{@code DoubleMatrix} translates double values of the internally stored solution to
   * {@link BigDecimal} values, using {@link BigDecimal#valueOf(double)} method.
   */
  @Override
  public BigDecimal[] getBigDecimalSolution() {
    return convertToBigDecimals(solution);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String getErrorCode() {
    return solver.errorCode();
  }

  /* *******************************************************************************
  /***** Solutions with respect of matrices ****************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   * as a new instance of {@code DoubleMatrix}.
   */
  @Override
  public Matrix solve(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solve(Matrix)");
    matrixSolution = solver.solve(matrixB.getDoubleData());
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   * as a new instance of {@code DoubleMatrix}.
   */
  @Override
  public Matrix solve(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solve(double[][])");
    matrixSolution = solver.solve(deepCopyOf(matrixB));
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   * as a new instance of {@code DoubleMatrix}.
   */
  @Override
  public Matrix solve(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solve(Number[][])");
    matrixSolution = solver.solve(convertToDoubles(matrixB));
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   * <br>The execution time is about 70 times longer than that of the simple
   * {@linkplain DoubleMatrix#solve(Matrix)}. In a typical case, the result error is reduced,
   * compared to the simple {@code solve()}, by a factor from several units to several tens,
   * depending on specific values of the matrix elements. For matrices 200 x 200
   * with uniformly-distributed random values, average square root of MSE is reduced
   * by approximately 30 times, from about 7.5e-14 down to about 2.5e-15<br>
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   *  as a new instance of {@code DoubleMatrix}.
   */
  @Override
  public Matrix solveAccurately(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solveAccurately(Matrix)");
    return solveAccurately(matrixB.getDoubleData());
  }

  /**
   * {@inheritDoc}
   * <br>The execution time is about 70 times longer than that of the simple {@linkplain DoubleMatrix#solve(Matrix)}.
   * In a typical case, the result error is reduced, compared to the simple {@code solve()}
   * by a factor from several units to several tens, depending on specific values of the matrix elements.
   * For matrices 200 x 200 with uniformly-distributed random values, average square root of MSE
   * is reduced by approximately 30 times, from about 7.5e-14 down to about 2.5e-15<br>
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   *  as a new instance of {@code DoubleMatrix}.
   */
  @Override
  public Matrix solveAccurately(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solveAccurately(double[][])");
    matrixSolution = solver.solveAccurately(deepCopyOf(matrixB));
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   * <br>The execution time is about 70 times longer than that of the simple {@linkplain DoubleMatrix#solve(Matrix)}.
   * In a typical case, the result error is reduced, compared to the simple {@code solve()}
   * by a factor from several units to several tens, depending on specific values of the matrix elements.
   * For matrices 200 x 200 with uniformly-distributed random values, average square root of MSE
   * is reduced by approximately 30 times, from about 7.5e-14 down to about 2.5e-15<br>
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   *  as a new instance of {@code DoubleMatrix}.
   */
  @Override
  public Matrix solveAccurately(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solveAccurately(Number[][])");
    return solveAccurately(convertToDoubles(matrixB));
  }

  /**
   * {@inheritDoc}
   *
   * @return the last of the previously found matrix solution <b><i>X</i></b>, to a system of
   * linear equations of form <b><i>AX = B</i></b>, as a {@code DoubleMatrix}, or {@code null},
   * if no matrix equation was solved with the matrix.
   */
  @Override
  public Matrix getMatrixSolution() {
    if (matrixSolution == null) return null;
    return new DoubleMatrix(matrixSolution); // Sic! The constructor makes a copy to preserve the solution
  }

  /**
   * {@inheritDoc}
   * <br>{@code DoubleMatrix#getNumberMatrixSolution()} returns the exact values of the elements
   * of the matrix solution found by the call to one of the {@code DoubleMatrix#solve()} or
   * {@code DoubleMatrix#solveAccurately()} methods without rounding or any other conversion.
   *
   * @return the last of the previously found matrix solution <b><i>X</i></b>, to a system of
   *  linear equations of form <b><i>AX = B</i></b>, as a two-dimentional array of {@linkplain Number},
   *  or {@code null}, if no matrix equation was solved with the matrix.
   */
  @Override
  public Number[][] getNumberMatrixSolution() {
    return convertToDoublesAsNumbers(matrixSolution);
  }

  /**
   * {@inheritDoc}
   * <br>The returned array contains the exact values of the elements of the matrix solution found by a call
   * to one of the {@code DoubleMatrix#solve()} or {@code DoubleMatrix#solveAccurately()} methods,
   * without rounding or any other conversion.
   */
  @Override
  public double[][] getDoubleMatrixSolution() {
    return deepCopyOf(matrixSolution);
  }

  /**
   * {@inheritDoc}
   * <br>The returned array contains the exact values of the elements of the matrix solution found by a call
   * to one of the {@code DoubleMatrix#solve()} or {@code DoubleMatrix#solveAccurately()} methods,
   * translated into corresponding {@link Quadruple} values, without any precision loss.
   */
  @Override
  public Quadruple[][] getQuadrupleMatrixSolution() {
    return convertToQuadruples(matrixSolution);
  }

  /**
   * {@inheritDoc}
   * <br>The returned array contains the exact values of the elements of the matrix solution found by a call
   * to one of the {@code DoubleMatrix#solve()} or {@code DoubleMatrix#solveAccurately()} methods,
   * translated to {@linkplain BigDecimal} values using {@linkplain BigDecimal#valueOf(double)}.
   * If the internal storage contains element that can't be translated  to {@linkplain BigDecimal} values (NaN or Infinity),
   * throws {@linkplain NumberFormatException}.
   */
  @Override
  public BigDecimal[][] getBigDecimalMatrixSolution() {
    return convertToBigDecimals(matrixSolution);
  }

  /* *******************************************************************************
  /***** Inversions ****************************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   * <br>The particular subtype of {@code Matrix} returned by {@code DoubleMatrix#inverse()}
   * is {@code DoubleMatrix}.
   */
  @Override
  public Matrix inverse() throws IllegalArgumentException {
    return newWithThisArray(solver.inverse());
  }

  /**
   * {@inheritDoc}
   * <br>The execution time is about 70 times longer than that of the simple
   * {@linkplain DoubleMatrix#inverse()}.
   * In a typical case, the result error is reduced, compared to the simple {@code inverse()}
   * by a factor from several units to several tens, depending on specific values of the matrix elements.
   * For matrices 200 x 200 with uniformly-distributed random values, average square root of MSE
   * is reduced by approximately 30 times, from about 1.2e-14 down to about 4.0e-16<br>
   *
   * @return a new {@code DoubleMatrix} containing the inversion of the given matrix.
   */
  @Override
  public Matrix inverseAccurately() throws IllegalArgumentException {
    return newWithThisArray(solver.inverseAccurately());
  };

  /**
   * {@inheritDoc}
   *
   * @return a new {@code DoubleMatrix} containing the transposition of the given matrix.
   */
  @Override
  public Matrix transpose() {
    return newWithThisArray(solver.transpose(matrix));
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code DoubleMatrix} containing the unity matrix of the same size as this {@code Matrix}
   */
  @Override
  public Matrix unity() {
    return newWithThisArray(solver.unityMatrix());
  }

  /* *******************************************************************************
  /* **** Multiplications ***********************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   *
   * @return a new {@code DoubleMatrix} representing the product
   */
  @Override
  public Matrix multiply(Matrix factor) throws IllegalArgumentException, NullPointerException {
    checkMatrix(factor, "multiply(Matrix)");
    return multiply(factor.getDoubleData());
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code DoubleMatrix} representing the product
   */
  @Override
  public Matrix multiply(double[][] factor) throws IllegalArgumentException, NullPointerException {
    checkMatrix(factor, "multiply(double[][])");
    final double[][] product = solver.multiply(factor);
    return newWithThisArray(product);
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code DoubleMatrix} representing the product
   */
  @Override
  public Matrix multiply(Number[][] factor) throws IllegalArgumentException, NullPointerException {
    checkMatrix(factor, "multiply(Number[][])");
    return multiply(convertToDoubles(factor));
  }

  /**
   * {@inheritDoc}
   *
   * @return an array of {@linkplain Double} values containing the product.
   */
  @Override
  public Number[] multiply(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "multiply(double[])");
    final double[] product = solver.multiply(vector);
    return convertToDoublesAsNumbers(product);
  }

  /**
   * {@inheritDoc}
   *
   * @return an array of {@linkplain Double} values containing the product.
   */
  @Override
  public Number[] multiply(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "multiply(Number[])");
    return multiply(convertToDoubles(vector));
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code DoubleMatrix} containing the product of the source matrix and the given scalar.
   */
  @Override
  public Matrix multiply(double scalar) throws IllegalArgumentException {
    checkScalar(scalar, "multiply(double)");
    final double[][] product = deepCopyOf(matrix);
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        product[i][j] *= scalar;
      }
    }
    return newWithThisArray(product);
  }

  /**
   * {@inheritDoc}
   * Before the multiplication, the value of the argument is translated
   * into the corresponding {@code double} value by {@linkplain Number#doubleValue()} method.
   *
   * @return a new {@code DoubleMatrix} containing the product of the source matrix and the given scalar.
   */
  @Override
  public Matrix multiply(Number scalar) throws IllegalArgumentException, NullPointerException {
    checkScalar(scalar, "multiply(Number)");
    return multiply(scalar.doubleValue());
  }

  /* *******************************************************************************
  /* **** Additions and subtractions ***********************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   * @return a new {@code DoubleMatrix} containing the sum of the two matrices.
   */
  @Override
  public Matrix add(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "add(Matrix)");
    return newWithThisArray(addMatrices(this.matrix, matrixB.getDoubleData()));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code DoubleMatrix} containing the sum of the two matrices.
   */
  @Override
  public Matrix add(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "add(double[][])");
    return newWithThisArray(addMatrices(this.matrix, matrixB));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code DoubleMatrix} containing the sum of the two matrices.
   */
  @Override
  public Matrix add(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "add(Number[][])");
    return newWithThisArray(addMatrices(this.matrix, convertToDoubles(matrixB)));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code DoubleMatrix} containing the difference of the two matrices.
   */
  @Override
  public Matrix subtract(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "subtract(Matrix)");
    return newWithThisArray(subtractMatrices(this.matrix, matrixB.getDoubleData()));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code DoubleMatrix} containing the difference of the two matrices.
   */
  @Override
  public Matrix subtract(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "subtract(double[][])");
    return newWithThisArray(subtractMatrices(this.matrix, matrixB));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code DoubleMatrix} containing the difference of the two matrices.
   */
  @Override
  public Matrix subtract(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "subtract(Number[][])");
    return newWithThisArray(subtractMatrices(this.matrix, convertToDoubles(matrixB)));
  }

  /* *******************************************************************************
  /* **** Determinant **************************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   * @return the value of the determinant of the given {@code Matrix} as an instance
   * of class {@code Double}
   */
  @Override
  public Number determinant() {
    return Double.valueOf(solver.determinant());
  }

  /**
   * {@inheritDoc}
   * {@code DoubleMatrix#determinantAsDouble()} returns the exact result of the corresponding
   * computations without any type conversions.
   */
  @Override
  public double determinantAsDouble() {
    return solver.determinant();
  }

  /**
   * {@inheritDoc}
   * {@code DoubleMatrix#determinantAsQuadruple()} returns a {@code Quadruple} value exactly
   * equal to the internally calculated double value, with no loss or increase in precision.
   */
  @Override
  public Quadruple determinantAsQuadruple() {
    return new Quadruple(solver.determinant());
  }

  /**
   * {@inheritDoc}
   * For {@code DoubleMatrix}, it is limited to the precision provided by calculations with
   * {@code double}s.
   */
  @Override
  public BigDecimal determinantAsBigDecimal() {
    return BigDecimal.valueOf(solver.determinant());
  }

  /* *******************************************************************************
  /* **** Norm *********************************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   * @return the value of the norm of the given {@code Matrix} as an instance
   * of class {@code Double}
   */
  @Override
  public Number norm() {
    return Double.valueOf(solver.norm());
  }

  /**
   * {@inheritDoc}
   * {@code DoubleMatrix#normAsDouble()} returns the exact result of the corresponding
   * computations without any type conversions.
   */
  @Override
  public double normAsDouble() {
    return solver.norm();
  }

  /**
   * {@inheritDoc}
   * {@code DoubleMatrix#normAsQuadruple()} returns a {@code Quadruple} value exactly
   * equal to the internally calculated double value, with no loss or increase in precision.
   */
  @Override
  public Quadruple normAsQuadruple() {
    return new Quadruple(solver.norm());
  }

  /**
   * {@inheritDoc}
   * For {@code DoubleMatrix}, it is limited to the precision provided by calculations with
   * {@code double}s.
   */
  @Override
  public BigDecimal normAsBigDecimal() {
    return BigDecimal.valueOf(solver.norm());
  }

  /* *******************************************************************************
  /* **** Condition number *********************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   */
  @Override
  public double cond() {
    return solver.cond();
  }

  /* *****************************************************************************************
   *** Support constructors of the superclass ************************************************
  /* *****************************************************************************************/

  private final void copyDataFrom(double[][] source) {
    this.size = source.length;
    matrix = deepCopyOf(source);
  }

  private final void copyDataFrom(Number[][] source) {
    this.size = source.length;
    matrix = checkAndCopy(source);
  }

  private final DoubleMatrixSolver makeSolver(boolean needToScale) {
    solver = new DoubleMatrixSolver(matrix, needToScale);
    return solver;
  }

  /* *******************************************************************************
  /***** Private methods ***********************************************************
  /*********************************************************************************/

  private double[][] checkAndCopy(Number[][] source) {
    final double[][] data = new double[size][];
    for (int i = 0; i < size; i++) {
      final Number[] row = source[i];
      if (row.length != size) {
        throw new IllegalArgumentException("Can't create matrix: data array must be square");
      }
      data[i] = convertToDoubles(row);
    }
    return data;
  }

}
