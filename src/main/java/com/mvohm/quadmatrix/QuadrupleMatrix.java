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

import static com.mvohm.quadmatrix.QuadrupleMatrixSolver.*;
import static com.mvohm.quadmatrix.Util.*;

/**
 * An implementation of the abstract class {@link Matrix} that uses an array of
 * {@link Quadruple} values to store the internal data and arithmetic provided by {@link Quadruple}
 * to perform calculations.
 * <br>
 * @author M.Vokhmentsev
 */
public class QuadrupleMatrix extends Matrix {

  /* *******************************************************************************
  /***** Internal data *************************************************************
  /*********************************************************************************/

  private Quadruple[][] matrix;

  private Quadruple[] solution;
  private Quadruple[][] matrixSolution;

  private QuadrupleMatrixSolver solver;

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
   * @see #QuadrupleMatrix(Matrix, boolean)
   * @see #QuadrupleMatrix(double[][], boolean)
   * @see #QuadrupleMatrix(Number[][], boolean)
   */
  private static boolean scaleByDefault = true;

  /* *******************************************************************************
  /***** Constructors **************************************************************
  /*********************************************************************************/

  /**
   * Creates a new {@code QuadrupleMatrix} with a copy of the data of the given {@code source} matrix and
   * the default value of the {@code needToScale} flag.
   * <br>Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows before solving equation systems and inverting.<br>
   * @param source a matrix whose data get copied into the new Matrix's internal storage.
   * @see #setDefaultScaling(boolean)
   */
  public QuadrupleMatrix(Matrix source) {
    this(source, scaleByDefault);
  }

  /**
   * Creates a new {@code QuadrupleMatrix} with a copy of the data of the given {@code source} matrix and
   * the value of the {@code needToScale} flag that is passed in as the {@code needToScale} argument.
   * <br>Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows while solving equation systems and inverting.
   * The value of this flag can be obtained via {@link #getScaling()} instance method.
   * @param source a matrix whose data get copied into the new Matrix's internal storage.
   * @param needToScale the value defining necessity to scale rows while solving equation systems and inverting.
   * @see #getScaling()
   */
  public QuadrupleMatrix(Matrix source, boolean needToScale) {
    checkMatrixData(source);
    this.size = source.size;
    this.matrix = source.getQuadrupleData();
    solver = makeSolver(needToScale);
  }

  /**
   * Creates a new {@code QuadrupleMatrix} whose inner data are obtained by converting
   * the values of the given {@code source} array into corresponding {@link Quadruple}
   * values without precision loss, and the default value of the {@code needToScale} flag.
   * <br>Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows before solving equation systems and inverting.
   *
   * <br>Throws {@linkplain IllegalArgumentException} if the {@code source} argument is a non-square
   * or empty array, or contains non-numeric values, {@link Double#NaN}, {@link Double#NEGATIVE_INFINITY},
   * or {@link Double#POSITIVE_INFINITY}.
   *
   * @param source a two-dimentional square array of {@code double}s
   *    to be converted into the {@link Quadruple} values of the internal storage of this instance.
   * @see #setDefaultScaling(boolean)
   * @see #scaleByDefault
   *
   */
  public QuadrupleMatrix(double[][] source) {
    this(source, scaleByDefault);
  }

  /**
   * Creates a new {@code QuadrupleMatrix} whose inner data are obtained by converting
   * the values of the given {@code source} array into corresponding {@link Quadruple}
   * values without precision loss, and the specified value of the {@code needToScale} flag
   * that is passed in as the {@code needToScale} argument.
   * <br>
   * Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows before solving equation systems and inverting.
   * The value of this flag can be obtained via {@link #getScaling()} instance method.
   * <br>
   * Throws {@linkplain IllegalArgumentException} if the {@code source} argument is a non-square
   * or empty array, or contains non-numeric values, {@link Double#NaN}, {@link Double#NEGATIVE_INFINITY},
   * or {@link Double#POSITIVE_INFINITY}.
   *
   * @param source a two-dimentional square array of {@code double}s
   *    to be converted into the {@link Quadruple} values of the internal storage of this instance.
   * @param needToScale the value defining necessity to scale rows while solving equation systems and inverting.
   * @see #setDefaultScaling(boolean)
   * @see #scaleByDefault
   *
   */
  public QuadrupleMatrix(double[][] source, boolean needToScale) {
    checkMatrixData(source);
    copyDataFrom(source);
    solver = makeSolver(needToScale);
  }

  /**
   * Creates a new {@code QuadrupleMatrix} with a copy of the data of the given {@code source} array and
   * the default value of the {@code needToScale} flag.
   *
   * <br>Translates the values of the specified {@code source} array to {@code Quadruple} values
   * using {@code Quadruple} constructors depending on the actual element type of the source array.
   *
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
  public QuadrupleMatrix(Number[][] source) {
    this(source, scaleByDefault);
  }

  /**
   * Creates a new {@code QuadrupleMatrix} with a copy of the data of the given {@code source} array and
   * the specified value of the {@code needToScale} flag.
   * <br>
   * Translates the values of the specified {@code source} array to {@code Quadruple} values
   * using {@code Quadruple} constructors depending on the actual element type of the source array.
   * <br>
   * Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows while solving equation systems and inverting.
   * The value of this flag can be obtained via {@link #getScaling()} instance method.
   * <br>
   * Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag designating
   * necessity to scale matrix rows before solving equation systems and inverting.
   * <br>
   * Throws {@linkplain IllegalArgumentException} if the {@code source} argument is a non-square or empty array,
   * or contains non-numeric values, {@link Double#NaN}, {@link Double#NEGATIVE_INFINITY}
   * or {@link Double#POSITIVE_INFINITY}), or {@code null}.
   *
   * @param source a two-dimentional square array of {@linkplain Number} values
   *    to be copied into the internal storage of this instance.
   * @param needToScale the value defining necessity to scale rows while solving equation systems and inverting.
   */
  public QuadrupleMatrix(Number[][] source, boolean needToScale) {
    checkMatrixData(source);
    copyDataFrom(source);
    solver = makeSolver(needToScale);
  }

  //***** Private constructor and fabric method creating a new instance with the same data array *************

  /** Private constructor that leaves the fields of the object blank. Used by the factory method
   * below */
  private QuadrupleMatrix() {
    super();
  }

  /** A factory method creating a new Matrix
  * with the same instance of the given data array, without allocating a new array.
  * Used by methods returning new QuadrupleMatrix instances to avoid unnecessary memory fragmentation and GC load */
  private static QuadrupleMatrix newWithThisArray(Quadruple[][] array) {
    final QuadrupleMatrix matrix = new QuadrupleMatrix();
    matrix.size = array.length;
    matrix.matrix = array;
    matrix.solver = new QuadrupleMatrixSolver(array, scaleByDefault);
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
    QuadrupleMatrix.scaleByDefault = scaleByDefault;
  };

  /**
   * Returns the value of the static {@code #scaleByDefault} flag that is used to set corresponding
   * instance flags when creating a new instances by constructors without {@code needToScale} parameter.
   * @see #setDefaultScaling(boolean)
   * @return the value of the static {@code #scaleByDefault} flag
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
   * <br>QuadrupleMatrix returns a copy of the internal matrix data as a two-dimentional
   * array of {@link Quadruple}.
   */
  @Override
  public Number[][] getData() {
    return deepCopyOf(matrix);
  };

  /**
   * {@inheritDoc}
   * <br>QuadrupleMatrix returns an array containing values of the internal matrix data
   * rounded to nearest <b>{@code double}</b> values.
   * @return a two-dimentional array of primitive <b>{@code double}</b> values containing
   * rounded values of the corresponding elements of the internal data array */
  @Override
  public double[][] getDoubleData() {
    return convertToDoubles(matrix);
  };

  /**
   * {@inheritDoc}
   * <br>QuadrupleMatrix returns an array containing the exact values of the corresponding elements
   * of the internal storage.
   */
  @Override
  public Quadruple[][] getQuadrupleData() {
    return deepCopyOf(matrix);
  };

  /**
   * {@inheritDoc}
   * <br>QuadrupleMatrix returns an array containing the values of the corresponding elements
   * of the internal storage, translated to {@linkplain BigDecimal} values
   * using {@linkplain Quadruple#bigDecimalValue()}.
   */
  @Override
  public BigDecimal[][] getBigDecimalData() {
    return convertToBigDecimals(matrix);
  };

  /** {@inheritDoc} */
  @Override
  public boolean equals(Object anotherOne) {
    if (this == anotherOne)
      return true;
    if (!(anotherOne instanceof QuadrupleMatrix))
      return false;
    return  (Arrays.deepEquals(matrix, ((QuadrupleMatrix)anotherOne).matrix))
            && (solver.needToScale == ((QuadrupleMatrix)anotherOne).solver.needToScale);
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

  /** {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * {@code Quadruple} values, using {@link Quadruple#Quadruple(double)} constructor.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array
   *  returned by {@code QuadrupleMatrix} is {@code Quadruple}.
   */
  @Override
  public  Number[] solve(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solve(double[])");
    solution = solver.solveLU(convertToQuadruples(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * exactly equal or nearest possible {@code Quadruple} values, depending on the particular type
   * of the {@code vector} elements.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code QuadrupleMatrix} is {@code Quadruple}.
   */
  @Override
  public  Number[] solve(Number[] vector) throws IllegalArgumentException, NullPointerException{
    checkVector(vector, "solve(Number[])");
    solution = solver.solveLU(convertToQuadruples(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * {@code Quadruple} values, using {@link Quadruple#Quadruple(double)} constructor.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code QuadrupleMatrix} is {@code Quadruple}.
   */
  @Override
  public  Number[] solveSPD(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPD(double[])");
    solution = solver.solveCholesky(convertToQuadruples(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * exactly equal or nearest possible {@code Quadruple} values, depending on the particular type
   * of the {@code vector} elements.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code QuadrupleMatrix} is {@code Quadruple}.
   */
  @Override
  public  Number[] solveSPD(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPD(Number[])");
    solution = solver.solveCholesky(convertToQuadruples(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * The execution time is about 50% longer than that of the simple
   * {@linkplain QuadrupleMatrix#solve(double[])}.
   * In a typical case, the result error is reduced, compared to the simple {@code solve()},
   * by a factor from several units to several tens or a few hundreds, depending on the properties
   * of the specific matrix data.<br>
   * For matrices 150 x 150 with uniformly-distributed random values, average square root of MSE
   * is reduced by approximately 40 times, from about 2.5e-37 down to 6.3e-39.<br>
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   * of {@linkplain Number} instances. The particular type of the elements of the array returned by
   * {@code QuadrupleMatrix} is {@code Quadruple}.
   */
  @Override
  public  Number[] solveAccurately(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveAccurately(double[])");
    solution = solver.solveLUAccurately(convertToQuadruples(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * exactly equal or nearest possible {@code Quadruple} values, depending on the particular type
   * of the {@code vector} elements.
   * <br>
   * The execution time is about 50% longer than that of the simple
   * {@linkplain QuadrupleMatrix#solve(Number[])}.
   * In a typical case, the result error is reduced, compared to the simple {@code solve()},
   * by a factor from several units to several tens or a few hundreds, depending on the properties
   * of the specific matrix data.<br>
   * For matrices 150 x 150 with uniformly-distributed random values, average square root of MSE
   * is reduced by approximately 40 times, from about 2.5e-37 down to 6.3e-39.<br>
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code QuadrupleMatrix} is {@code Quadruple}.
   */
  @Override
  public  Number[] solveAccurately(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveAccurately(Number[])");
    solution = solver.solveLUAccurately(convertToQuadruples(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * The execution time is about 75% longer than that of the simple
   * {@linkplain QuadrupleMatrix#solveSPD(double[])}.
   * In a typical case, the result error is reduced, compared to the simple {@code solveSPD()},
   * by a factor of several units, depending on the properties of the specific matrix data.<br>
   * For matrices 150 x 150, average square root of MSE
   * is reduced by approximately 5 times, from about 3.5e-38 down to 7.0e-39.<br>
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code QuadrupleMatrix} is {@code Quadruple}.
   */
  @Override
  public  Number[] solveSPDAccurately(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPDAccurately(double[])");
    solution = solver.solveCholeskyAccurately(convertToQuadruples(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * exactly equal or nearest possible {@code Quadruple} values, depending on the particular type
   * of the {@code vector} elements.
   * <br>
   * The execution time is about 75% longer than that of the simple
   * {@linkplain QuadrupleMatrix#solveSPD(double[])}.
   * In a typical case, the result error is reduced, compared to the simple {@code solveSPD()},
   * by a factor of several units, depending on the properties of the specific matrix data.<br>
   * For matrices 150 x 150, average square root of MSE
   * is reduced by approximately 5 times, from about 3.5e-38 down to 7.0e-39.<br>
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code QuadrupleMatrix} is {@code Quadruple}.
   */
  @Override
  public  Number[] solveSPDAccurately(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPDAccurately(Number[])");
    solution = solver.solveCholeskyAccurately(convertToQuadruples(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * @return the last of the previously found vector solutions, <b>x</b>, to a system of of linear equations
   *  of form <b><i>A</i>x = b</b>, as an array of {@link Quadruple}s, or {@code null},
   *  if no system was solved with the matrix.
   */
  @Override
  public Number[] getSolution() {
    if (solution == null) return null;

    final Quadruple[] result = new Quadruple[size];
    for (int i = 0; i < size; i++) {
      result[i] = new Quadruple(solution[i]);
    }
    return result;
  }

  /**
   * {@inheritDoc}
   * <br>The elements of the array returned by {@code QudrupleMatrix} are corresponded values of the found
   * solution, rounded to the closest {@code double} values
   */
  @Override
  public double[] getDoubleSolution() {
    return convertToDoubles(solution);
  }

  /**
   * {@inheritDoc}
   * <br>{@code QudrupleMatrix} returns the exact values that was found by the corresponding solution method.
   */
  @Override
  public Quadruple[] getQuadrupleSolution() {
    return deepCopyOf(solution);
  }

  /**
   * {@inheritDoc}
   * <br>{@code QuadrupleMatrix} translates {@link Quadruple} values of the internally stored solution to
   * {@link BigDecimal} values, using {@link Quadruple#bigDecimalValue()} method.
   */
  @Override
  public BigDecimal[] getBigDecimalSolution() {
    return convertToBigDecimals(solution);
  }

  /** {@inheritDoc} */
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
   * as a new instance of {@code QuadrupleMatrix}.
   */
  @Override
  public Matrix solve(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solve(Matrix)");
    matrixSolution = solver.solve(matrixB.getQuadrupleData());
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   * as a new instance of {@code QuadrupleMatrix}.
   */
  @Override
  public Matrix solve(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solve(double[][])");
    matrixSolution = solver.solve(convertToQuadruples(matrixB));
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   * as a new instance of {@code QuadrupleMatrix}.
   */
  @Override
  public Matrix solve(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solve(Number[][])");
    matrixSolution = solver.solve(convertToQuadruples(matrixB));
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   * <br>The execution time is about 30 times longer than that of the simple
   * {@linkplain QuadrupleMatrix#solve(Matrix)}. In a typical case, the result error is reduced,
   * compared to the simple {@code solve()}, by a factor from several units to several tens,
   * depending on specific values of the matrix elements. For matrices 200 x 200
   * with uniformly-distributed random values, average square root of MSE is reduced
   * by approximately 40 times, from about 1.2e-36 down to about 3.0e-38.<br>
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   *  as a new instance of {@code QuadrupleMatrix}.
   */
  @Override
  public Matrix solveAccurately(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solveAccurately(Matrix)");
    return solveAccurately(matrixB.getQuadrupleData());
  }

  /**
   * {@inheritDoc}
   * <br>The execution time is about 30 times longer than that of the simple
   * {@linkplain QuadrupleMatrix#solve(Matrix)}. In a typical case, the result error is reduced,
   * compared to the simple {@code solve()}, by a factor from several units to several tens,
   * depending on specific values of the matrix elements. For matrices 200 x 200
   * with uniformly-distributed random values, average square root of MSE is reduced
   * by approximately 40 times, from about 1.2e-36 down to about 3.0e-38.<br>
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   *  as a new instance of {@code QuadrupleMatrix}.
   */
  @Override
  public Matrix solveAccurately(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solveAccurately(double[][])");
    matrixSolution = solver.solveAccurately(convertToQuadruples(matrixB));
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   * <br>The execution time is about 30 times longer than that of the simple
   * {@linkplain QuadrupleMatrix#solve(Matrix)}. In a typical case, the result error is reduced,
   * compared to the simple {@code solve()}, by a factor from several units to several tens,
   * depending on specific values of the matrix elements. For matrices 200 x 200
   * with uniformly-distributed random values, average square root of MSE is reduced
   * by approximately 40 times, from about 1.2e-36 down to about 3.0e-38.<br>
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   *  as a new instance of {@code QuadrupleMatrix}.
   */
  @Override
  public Matrix solveAccurately(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solveAccurately(Number[][])");
    matrixSolution = solver.solveAccurately(convertToQuadruples(matrixB));
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   *
   * @return the last of the previously found matrix solution <b><i>X</i></b> to a system of
   * linear equations of form <b><i>AX = B</i></b>, as a {@code QuadrupleMatrix}, or {@code null},
   * if no matrix equation was solved with the matrix.
   */
  @Override
  public Matrix getMatrixSolution() {
    if (matrixSolution == null) return null;
    return new QuadrupleMatrix(matrixSolution); // Sic! The constructor makes a copy to preserve the solution
  }

  /**
   * {@inheritDoc}
   * <br>{@code QuadrupleMatrix#getNumberMatrixSolution()} returns the exact values of the elements
   * of the matrix solution found by the call to one of the {@code QuadrupleMatrix#solve()} or
   * {@code QuadrupleMatrix#solveAccurately()} methods.
   *
   * @return the last of the previously found matrix solution <b><i>X</i></b>, to a system of
   *  linear equations of form <b><i>AX = B</i></b>, as a two-dimentional array of {@linkplain Quadruple},
   *  or {@code null}, if no matrix equation was solved with the matrix.
   */
  @Override
  public Number[][] getNumberMatrixSolution() {
    if (matrixSolution == null) return null;
    return deepCopyOf(matrixSolution);
  }

  /**
   * {@inheritDoc}
   * <br>The returned array contains the values of the elements of the matrix solution found by a call
   * to one of the {@code QuadrupleMatrix#solve()} or {@code QuadrupleMatrix#solveAccurately()} methods,
   * rounded to nearest {@code double} values.
   */
  @Override
  public double[][] getDoubleMatrixSolution() {
    if (matrixSolution == null) return null;
    return convertToDoubles(matrixSolution);
  }

  /**
   * {@inheritDoc}
   * <br>The returned array contains the exact values of the elements of the matrix solution found by a call
   * to one of the {@code QuadrupleMatrix#solve()} or {@code QuadrupleMatrix#solveAccurately()} methods.
   */
  @Override
  public Quadruple[][] getQuadrupleMatrixSolution() {
    if (matrixSolution == null) return null;
    return deepCopyOf(matrixSolution);
  }

  /**
   * {@inheritDoc}
   * <br>The returned array contains the values of the elements of the matrix solution found by a call
   * to one of the {@code QuadrupleMatrix#solve()} or {@code QuadrupleMatrix#solveAccurately()} methods,
   * translated to {@linkplain BigDecimal} values using {@linkplain Quadruple#bigDecimalValue()} method.
   * If the internal storage contains element that can't be translated to {@linkplain BigDecimal} values (NaN or Infinity),
   * throws {@linkplain NumberFormatException}.
   */
  @Override
  public BigDecimal[][] getBigDecimalMatrixSolution() {
    if (matrixSolution == null) return null;
    return convertToBigDecimals(matrixSolution);
  }

/* *******************************************************************************
  /***** Inversions ****************************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   * <br>The particular subtype of {@code Matrix} returned by {@code QuadrupleMatrix#inverse()}
   * is {@code QuadrupleMatrix}.
   */
  @Override
  public Matrix inverse() throws IllegalArgumentException {
    return newWithThisArray(solver.inverse());
  }

  /**
   * {@inheritDoc}
   * <br>The execution time is about 30 times longer than that of the simple
   * {@linkplain QuadrupleMatrix#inverse()}.
   * In a typical case, the result error is reduced, compared to the simple {@code inverse()},
   * by a factor of several units, depending on specific values of the matrix elements.
   * For matrices 200 x 200 with uniformly-distributed random values, average square root of MSE
   * is reduced by approximately 5 times, from about 2.0e-37 down to about 4.0e-38<br>
   *
   * @return a new {@code QuadrupleMatrix} containing the inversion of the given matrix.
   */
  @Override
  public Matrix inverseAccurately() throws IllegalArgumentException {
    return newWithThisArray(solver.inverseAccurately());
  };

  /**
   * {@inheritDoc}
   *
   * @return a new {@code QuadrupleMatrix} containing the transposition of the given matrix.
   */
  @Override
  public Matrix transpose() {
    return newWithThisArray(solver.transpose(matrix));
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code QuadrupleMatrix} containing the unity matrix of the same size as this {@code Matrix}
   */
  @Override
  public Matrix unity() {
    return newWithThisArray(solver.unityMatrix());
  }

  /* *******************************************************************************
  /***** Multiplications ***********************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   *
   * @return a new {@code QuadrupleMatrix} representing the product
   */
  @Override
  public Matrix multiply(Matrix factor) throws IllegalArgumentException, NullPointerException {
    checkMatrix(factor, "multiply(Matrix)");
    final Quadruple[][] product = solver.multiply(factor.getQuadrupleData());
    return newWithThisArray(product);
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code QuadrupleMatrix} representing the product
   */
  @Override
  public Matrix multiply(double[][] factor) throws IllegalArgumentException, NullPointerException {
    checkMatrix(factor, "multiply(double[][])");
    final Quadruple[][] product = solver.multiply(convertToQuadruples(factor));
    return newWithThisArray(product);
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code QuadrupleMatrix} representing the product
   */
  @Override
  public Matrix multiply(Number[][] factor) throws IllegalArgumentException, NullPointerException {
    checkMatrix(factor, "multiply(Number[][])");
    final Quadruple[][] product = solver.multiply(convertToQuadruples(factor));
    return newWithThisArray(product);
  }

  /**
   * {@inheritDoc}
   *
   * @return an array of {@linkplain Quadruple} values containing the product.
   */
  @Override
  public Number[] multiply(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "multiply(double[])");
    return solver.multiply(convertToQuadruples(vector));
  }

  /**
   * {@inheritDoc}
   *
   * @return an array of {@linkplain Quadruple} values containing the product.
   */
  @Override
  public Number[] multiply(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "multiply(Number[])");
    return solver.multiply(convertToQuadruples(vector));
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code QuadrupleMatrix} containing the product of the source matrix and the given scalar.
   */
  @Override
  public Matrix multiply(double scalar) throws IllegalArgumentException {
    checkScalar(scalar, "multiply(double)");
    return newWithThisArray(solver.multiply(new Quadruple(scalar)));
  }

  /**
   * {@inheritDoc}
   * Before the multiplication, the value of the argument is translated
   * into the corresponding {@code Quadruple} by one of the {@code Quadruple} constructors,
   * depending on the actual type of the argument.
   *
   * @return a new {@code QuadrupleMatrix} containing the product of the source matrix and the given scalar.
   */
  @Override
  public Matrix multiply(Number scalar) throws IllegalArgumentException, NullPointerException {
    checkScalar(scalar, "multiply(Number)");
    Quadruple factor;
    if (scalar instanceof Double)
      factor = new Quadruple((double)scalar);
    else if (scalar instanceof Quadruple)
      factor = (Quadruple)scalar;
    else if (scalar instanceof BigDecimal)
      factor = new Quadruple((BigDecimal)scalar);
    else
      throw new IllegalArgumentException("only Double, Quadruple and BigDecimal values are allowed for Number arguments");
    return newWithThisArray(solver.multiply(factor));
  }

  /* *******************************************************************************
  /* **** Additions and subtractions ***********************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   * @return a new {@code QuadrupleMatrix} containing the sum of the two matrices.
   */
  @Override
  public Matrix add(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "add(Matrix)");
    return newWithThisArray(addMatrices(this.matrix, matrixB.getQuadrupleData()));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code QuadrupleMatrix} containing the sum of the two matrices.
   */
  @Override
  public Matrix add(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "add(double[][])");
    return newWithThisArray(addMatrices(this.matrix, convertToQuadruples(matrixB)));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code QuadrupleMatrix} containing the sum of the two matrices.
   */
  @Override
  public Matrix add(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "add(Number[][])");
    return newWithThisArray(addMatrices(this.matrix, convertToQuadruples(matrixB)));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code QuadrupleMatrix} containing the difference of the two matrices.
   */
  @Override
  public Matrix subtract(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "subtract(Matrix)");
    return newWithThisArray(subtractMatrices(this.matrix, matrixB.getQuadrupleData()));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code QuadrupleMatrix} containing the difference of the two matrices.
   */
  @Override
  public Matrix subtract(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "subtract(double[][])");
    return newWithThisArray(subtractMatrices(this.matrix, convertToQuadruples(matrixB)));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code QuadrupleMatrix} containing the difference of the two matrices.
   */
  @Override
  public Matrix subtract(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "subtract(Number[][])");
    return newWithThisArray(subtractMatrices(this.matrix, convertToQuadruples(matrixB)));
  }

  /* *******************************************************************************
  /* **** Determinant **************************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   * @return the value of the determinant of the given {@code Matrix} as an instance
   * of class {@code Quadruple}
   */
  @Override
  public Number determinant() {
    return new Quadruple(solver.determinant());
  }

  /**
   * {@inheritDoc}
   * {@code QuadrupleMatrix#determinantAsDouble()} returns the result of the corresponding
   * computations rounded to the nearest {@code double} value.
   */
  @Override
  public double determinantAsDouble() {
    return solver.determinant().doubleValue();
  }

  /**
   * {@inheritDoc}
   * {@code QuadrupleMatrix#determinantAsQuadruple()} returns a {@code Quadruple} value exactly
   * equal to the internally calculated {@code Quadruple} value.
   */
  @Override
  public Quadruple determinantAsQuadruple() {
    return new Quadruple(solver.determinant());
  }

  /**
   * {@inheritDoc}
   * For {@code QuadrupleMatrix}, it is limited to the precision provided by calculations with
   * {@code Quadruple}s.
   */
  @Override
  public BigDecimal determinantAsBigDecimal() {
    return solver.determinant().bigDecimalValue();
  }

  /* *******************************************************************************
  /* **** Norm *********************************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   * @return the value of the norm of the given {@code Matrix} as an instance
   * of class {@code Quadruple}
   */
  @Override
  public Number norm() {
    return new Quadruple(solver.norm());
  }

  /**
   * {@inheritDoc}
   * {@code QuadrupleMatrix#normAsDouble()} returns the result of the corresponding
   * computations rounded to the nearest {@code double} value.
   */
  @Override
  public double normAsDouble() {
    return solver.norm().doubleValue();
  }

  /**
   * {@inheritDoc}
   * {@code QuadrupleMatrix#normAsQuadruple()} returns a {@code Quadruple} value exactly
   * equal to the internally calculated value.
   */
  @Override
  public Quadruple normAsQuadruple() {
    return new Quadruple(solver.norm());
  }

  /**
   * {@inheritDoc}
   * For {@code QuadrupleMatrix}, it is limited to the precision provided by calculations with
   * {@code Quadruple}s.
   */
  @Override
  public BigDecimal normAsBigDecimal() {
    return solver.norm().bigDecimalValue();
  }

  /* *******************************************************************************
  /* **** Condition number *********************************************************
  /*********************************************************************************/

  /** {@inheritDoc} */
  @Override
  public double cond() {
    return solver.cond();
  }

  /* *******************************************************************************
  /* **** Private methods **********************************************************
  /*********************************************************************************/

  private final QuadrupleMatrixSolver makeSolver(boolean needToScale) {
    solver = new QuadrupleMatrixSolver(matrix, needToScale);
    return solver;
  }

  private final void copyDataFrom(double[][] source) {
    this.size = source.length;
    matrix = convertToQuadruples(source);
  }

  private final void copyDataFrom(Number[][] source) {
    this.size = source.length;
    matrix = convertToQuadruples(source);
  }

}
