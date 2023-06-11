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
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;

import org.apache.commons.math3.util.Precision;

import com.mvohm.quadruple.Quadruple;

import static com.mvohm.quadmatrix.BigDecimalMatrixSolver.*;
import static com.mvohm.quadmatrix.Util.*;

/**
 * An implementation of the abstract class {@link Matrix} that uses an array of
 * {@link BigDecimal} values to store the internal data and arithmetic provided by {@link BigDecimal}
 * to perform calculations. <br>
 * The precision of the calculations to be performed by a new instance of {@code BigDecimalMatrix}
 * can be controlled via the {@code precision} parameter of a constructor,
 * or by setting {@code defaultPrecision} wti {@link #setDefaultPrecision(int)} method.
 * The default precision value is used to set the precision for all new instances which
 * are created by constructors that don't have {@code precision} parameter.
 * If the default value of {@code defaultPrecision} is 40, this value will be used
 * if no special actions to control the precision was done.<br>
 * All arithmetic operations performed internally are done with a {@link MathContext} instance
 * which is created with the precision set for the given instance and {@link RoundingMode#HALF_EVEN} mode.
 *
 * <br>
 * @author M.Vokhmentsev
 */
public class BigDecimalMatrix extends Matrix {

  /* *******************************************************************************
  /***** Internal data *************************************************************
  /*********************************************************************************/

  private BigDecimal[][] matrix;

  private BigDecimal[] solution;
  private BigDecimal[][] matrixSolution;

  private BigDecimalMatrixSolver solver;

  private int precision;
  private MathContext mc;

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
   * @see #BigDecimalMatrix(Matrix, boolean)
   * @see #BigDecimalMatrix(double[][], boolean)
   * @see #BigDecimalMatrix(Number[][], boolean)
   */
  private static boolean scaleByDefault = true;

  /**
   * A default value of the precision that is set in new instances by
   * constructors that don't have the corresponding parameter.
   */
  private static int defaultPrecision = 40;

  /* *******************************************************************************
  /***** Constructors **************************************************************
  /*********************************************************************************/

  /**
   * Creates a new {@code BigDecimalMatrix} with a copy of the data of the given {@code source}
   * matrix and the default values of the {@code needToScale} flag and {@code precision} field.
   * <br>
   * Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows before solving equation systems and inverting.
   * <br>
   * Uses static {@linkplain defaultPrecision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * @param source a matrix whose data get copied into the new Matrix's internal storage.
   * @see #setDefaultScaling(boolean)
   * @see #setDefaultPrecision(int)
   */
  public BigDecimalMatrix(Matrix source) {
    this(source, scaleByDefault, defaultPrecision);
  }

  /**
   * Creates a new {@code BigDecimalMatrix} with a copy of the data of the given {@code source} matrix,
   * the specified value of the {@code needToScale} flag, and the default precision set for the class.
   * <br>
   * Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows while solving equation systems and inverting.
   * The value of this flag can be obtained via {@link #getScaling()} instance method.
   * <br>
   * Uses static {@linkplain defaultPrecision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * @param source a matrix whose data get copied into the new Matrix's internal storage.
   * @param needToScale the value defining necessity to scale rows while solving equation systems and inverting.
   * @see #getScaling()
   * @see #setDefaultPrecision(int)
   */
  public BigDecimalMatrix(Matrix source, boolean needToScale) {
    this(source, needToScale, defaultPrecision);
  }

  /**
   * Creates a new {@code BigDecimalMatrix} with a copy of the data of the given {@code source} matrix,
   * the default value of the {@code needToScale} flag, and the specified precision.
   * <br>
   * Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows before solving equation systems and inverting.
   * <br>
   * Uses the specified {@code precision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * @see #setDefaultScaling(boolean)
   * @see #getPrecision()
   * @param source a matrix whose data get copied into the new Matrix's internal storage.
   * @param precision the value defining the precision of the calculations that the
   *  newly-created instance will do.
   */
  public BigDecimalMatrix(Matrix source, int precision) {
    this(source, scaleByDefault, precision);
  }

  /**
   * Creates a new {@code BigDecimalMatrix} with a copy of the data of the given {@code source}
   * matrix and the specified values of the {@code needToScale} flag and {@code precision} field.
   * <br>
   * Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows while solving equation systems and inverting.
   * <br>
   * Uses the specified {@code precision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * @see #getPrecision()
   * @see #getScaling()
   * @param source a matrix whose data get copied into the new Matrix's internal storage.
   * @param needToScale the value defining necessity to scale rows while solving equation systems and inverting.
   * @param precision the value defining the precision of the calculations that the
   *  newly-created instance will do.
   */
  public BigDecimalMatrix(Matrix source, boolean needToScale, int precision) {
    checkMatrixData(source);
    this.size = source.size;
    this.matrix = source.getBigDecimalData();
    setPrecision(precision);
    solver = makeSolver(needToScale);
  }

  /**
   * Creates a new {@code BigDecimalMatrix} whose inner data are obtained by converting
   * the values of the given {@code source} array into corresponding {@link BigDecimal}
   * values by {@link BigDecimal#valueOf(double)} method, and the default values
   * of the {@code needToScale} flag and {@code precision} field.
   * <br>
   * Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows before solving equation systems and inverting.
   * <br>
   * Uses static {@linkplain defaultPrecision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * Throws {@linkplain IllegalArgumentException} if the {@code source} argument is a non-square
   * or empty array, or contains non-numeric values, {@link Double#NaN}, {@link Double#NEGATIVE_INFINITY},
   * or {@link Double#POSITIVE_INFINITY}.
   *
   * @see #setDefaultScaling(boolean)
   * @see #setDefaultPrecision(int)
   * @param source a two-dimentional square array of {@code double}s
   *    to be converted into the {@link BigDecimal} values of the internal storage of this instance.
   */
  public BigDecimalMatrix(double[][] source) {
    this(source, scaleByDefault, defaultPrecision);
  }

  /**
   * Creates a new {@code BigDecimalMatrix} whose inner data are obtained by converting
   * the values of the given {@code source} array into corresponding {@link BigDecimal}
   * values by {@link BigDecimal#valueOf(double)} method, the specified value of the
   * {@code needToScale} flag, and the default precision set for the class.
   * <br>
   * Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows before solving equation systems and inverting.
   * The value of this flag can be obtained via {@link #getScaling()} instance method.
   * <br>
   * Uses static {@linkplain defaultPrecision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * Throws {@linkplain IllegalArgumentException} if the {@code source} argument is a non-square
   * or empty array, or contains non-numeric values, {@link Double#NaN}, {@link Double#NEGATIVE_INFINITY},
   * or {@link Double#POSITIVE_INFINITY}.
   * <br>
   * @see #getScaling()
   * @see #setDefaultPrecision(int)
   * @param source a two-dimentional square array of {@code double}s
   *  to be converted into the {@link BigDecimal} values of the internal storage of this instance.
   * @param needToScale the value defining necessity to scale rows while solving equation systems and inverting.
   */
  public BigDecimalMatrix(double[][] source, boolean needToScale) {
    this(source, needToScale, defaultPrecision);
  }

  /**
   * Creates a new {@code BigDecimalMatrix} whose inner data are obtained by converting
   * the values of the given {@code source} array into corresponding {@link BigDecimal}
   * values by {@link BigDecimal#valueOf(double)} method,
   * the default value of the {@code needToScale} flag and the specified precision.
   * <br>
   * Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows before solving equation systems and inverting.
   * <br>
   * Uses the specified {@code precision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * @see #setDefaultScaling(boolean)
   * @see #getPrecision()
   * @param source a two-dimentional square array of {@code double}s
   *  to be converted into the {@link BigDecimal} values of the internal storage of this instance.
   * @param precision the value defining the precision of the calculations that the
   *  newly-created instance will do.
   */
  public BigDecimalMatrix(double[][] source, int precision) {
    this(source, scaleByDefault, precision);
  }

  /**
   * Creates a new {@code BigDecimalMatrix} whose inner data are obtained by converting
   * the values of the given {@code source} array into corresponding {@link BigDecimal}
   * values by {@link BigDecimal#valueOf(double)} method,
   * and the specified values of the {@code needToScale} flag and {@code precision} field.
   * <br>
   * Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows while solving equation systems and inverting.
   * <br>
   * Uses the specified {@code precision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * @see #getPrecision()
   * @see #getScaling()
   * @param source a two-dimentional square array of {@code double}s
   *  to be converted into the {@link BigDecimal} values of the internal storage of this instance.
   * @param needToScale the value defining necessity to scale rows while solving equation systems and inverting.
   * @param precision the value defining the precision of the calculations that the
   *  newly-created instance will do.
   */
  public BigDecimalMatrix(double[][] source, boolean needToScale, int precision) {
    checkMatrixData(source);
    copyDataFrom(source);
    setPrecision(precision);
    solver = makeSolver(needToScale);
  }

  /**
   * Creates a new {@code BigDecimalMatrix} with a copy of the data of the given
   * {@code source} array and the specified value of the {@code needToScale} flag.
   * <br>
   * Translates the values of the specified {@code source} array to {@code BigDecimal} values
   * with methods depending on the actual element type of the source array --
   * it's {@link BigDecimal#valueOf(double)} for Doubles,
   * {@link Quadruple#bigDecimalValue()} for Quadruples,
   * and just copying for {@link BigDecimal} values.
   * <br>
   * Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows before solving equation systems and inverting.
   * <br>
   * Uses static {@linkplain defaultPrecision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * Throws {@linkplain IllegalArgumentException} if the {@code source} argument is a non-square
   * or empty array, or contains non-numeric values, {@link Double#NaN}, {@link Double#NEGATIVE_INFINITY},
   * or {@link Double#POSITIVE_INFINITY}.
   *
   * @see #setDefaultScaling(boolean)
   * @see #setDefaultPrecision(int)
   * @param source a two-dimentional square array of {@code double}s
   *    to be converted into the {@link BigDecimal} values of the internal storage of this instance.
   */
  public BigDecimalMatrix(Number[][] source) {
    this(source, scaleByDefault, defaultPrecision);
  }

  /**
   * Creates a new {@code BigDecimalMatrix} with a copy of the data of the given
   * {@code source} array, the specified value of the {@code needToScale} flag,
   * and the default precision set for the class.
   * <br>
   * Translates the values of the specified {@code source} array to {@code BigDecimal} values
   * with methods depending on the actual element type of the source array --
   * it's {@link BigDecimal#valueOf(double)} for Doubles,
   * {@link Quadruple#bigDecimalValue()} for Quadruples,
   * and just copying for {@link BigDecimal} values.
   * <br>
   * Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows before solving equation systems and inverting.
   * The value of this flag can be obtained via {@link #getScaling()} instance method.
   * <br>
   * Uses static {@linkplain defaultPrecision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * Throws {@linkplain IllegalArgumentException} if the {@code source} argument is a non-square
   * or empty array, or contains non-numeric values, {@link Double#NaN}, {@link Double#NEGATIVE_INFINITY},
   * or {@link Double#POSITIVE_INFINITY}.
   * <br>
   * @see #getScaling()
   * @see #setDefaultPrecision(int)
   * @param source a two-dimentional square array of {@code double}s
   *  to be converted into the {@link BigDecimal} values of the internal storage of this instance.
   * @param needToScale the value defining necessity to scale rows while solving equation systems and inverting.
   */
  public BigDecimalMatrix(Number[][] source, boolean needToScale) {
    this(source, needToScale, defaultPrecision);
  }

  /**
   * Creates a new {@code BigDecimalMatrix} with a copy of the data of the given
   * {@code source} array, the default value of the {@code needToScale}
   * flag and the specified precision.
   * <br>
   * Translates the values of the specified {@code source} array to {@code BigDecimal} values
   * with methods depending on the actual element type of the source array --
   * it's {@link BigDecimal#valueOf(double)} for Doubles,
   * {@link Quadruple#bigDecimalValue()} for Quadruples,
   * and just copying for {@link BigDecimal} values.
   * <br>
   * Uses static {@linkplain scaleByDefault} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows before solving equation systems and inverting.
   * <br>
   * Uses the specified {@code precision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * @see #setDefaultScaling(boolean)
   * @see #getPrecision()
   * @param source a two-dimentional square array of {@code double}s
   *  to be converted into the {@link BigDecimal} values of the internal storage of this instance.
   * @param precision the value defining the precision of the calculations that the
   *  newly-created instance will do.
   */
  public BigDecimalMatrix(Number[][] source, int precision) {
    this(source, scaleByDefault, precision);
  }

  /**
   * Creates a new {@code BigDecimalMatrix} with a copy of the data of the given
   * {@code source} array,
   * and the specified values of the {@code needToScale} flag and {@code precision} field.
   * <br>
   * <br>
   * Translates the values of the specified {@code source} array to {@code BigDecimal} values
   * with methods depending on the actual element type of the source array --
   * it's {@link BigDecimal#valueOf(double)} for Doubles,
   * {@link Quadruple#bigDecimalValue()} for Quadruples,
   * and just copying for {@link BigDecimal} values.
   * <br>
   * <br>
   * Uses the given {@code needToScale} value to set the new matrix's {@code needToScale} flag
   * designating necessity to scale matrix rows while solving equation systems and inverting.
   * <br>
   * Uses the specified {@code precision} value to set the new matrix's {@code precision} field
   * that defines the precision of the computations that the instance methods perform.
   * <br>
   * @see #getPrecision()
   * @see #getScaling()
   * @param source a two-dimentional square array of {@code double}s
   *  to be converted into the {@link BigDecimal} values of the internal storage of this instance.
   * @param precision the value defining the precision of the calculations that the
   * @param needToScale the value defining necessity to scale rows while solving equation systems and inverting.
   *  newly-created instance will do.
   */
  public BigDecimalMatrix(Number[][] source, boolean needToScale, int precision) {
    checkMatrixData(source);
    copyDataFrom(source);
    setPrecision(precision);
    solver = makeSolver(needToScale);
  }

  //***** Private constructor and fabric method creating a new instance with the same data array *************

  /** Private constructor that leaves the fields of the object blank. Used by the factory method
   * below  */
  private BigDecimalMatrix() {
    super();
  }

  /** A factory method creating a new Matrix
  * with the same instance of the given data array, without allocating a new array.
  * Used by methods returning new BigDecimalMatrix instances to avoid unnecessary memory fragmentation and GC load */
  private static BigDecimalMatrix newWithThisArray(BigDecimal[][] array) {
    final BigDecimalMatrix matrix = new BigDecimalMatrix();
    matrix.size = array.length;
    matrix.matrix = array;
    matrix.setPrecision(defaultPrecision);
    matrix.solver = new BigDecimalMatrixSolver(array, scaleByDefault, matrix.mc);
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
   *   created with constructors without {@code needToScale} parameter.
   */
  public static void setDefaultScaling(boolean scaleByDefault) {
    BigDecimalMatrix.scaleByDefault = scaleByDefault;
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

  /**
   * Sets the value of the static {@link #defaultPrecision} variable that is used to set corresponding
   * instance fields when creating a new instances by constructors without {@code precision} parameter.
   * <br>This value is used to create a {@link MathContext} object that is used by all arithmetic operations
   * performed by the instance and thus defines the precision of the calculations.
   * @param precision the value of the precision to be used when
   *  creating new instances by constructors without {@code precision}parameter
   */
  public static void setDefaultPrecision(int precision) {
    defaultPrecision = precision;
  }

  /**
   * Returns the current value of the precision used by default when creating new instances with
   * constructors that don't have {@code precision} parameter.
   * @see #setDefaultPrecision(int)
   * @return the default precision value currently set for the {@code BigDecimalMatrix} class
   */
  public static int getDefaultPrecision() {
    return defaultPrecision;
  }

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
   * Returns the value of {@code precision} used by this instance of {@code BigDecimalMatrix}.
   * @see #setDefaultPrecision(int)
   * @see #BigDecimalMatrix(Matrix, int)
   * @see #BigDecimalMatrix(double[][], int)
   * @see #BigDecimalMatrix(Number[][], int)
   * @return the value of {@code precision} used by this instance.
   */
  public int getPrecision() {
    return precision;
  }

  /**
   * {@inheritDoc}
   * <br>BigDecimalMatrix returns an array containing the exact values of the corresponding elements
   * of the internal storage
   */
  @Override
  public Number[][] getData() {
    return deepCopyOf(matrix);
  };

  /**
   * {@inheritDoc}
   * <br>BigDecimalMatrix returns an array containing values of the internal matrix data
   * rounded to nearest <b>{@code double}</b> values by {@link BigDecimal#doubleValue()} method.
   * @return a two-dimentional array of primitive <b>{@code double}</b> values containing
   * rounded values of the corresponding elements of the internal data array.
   */
  @Override
  public double[][] getDoubleData() {
    return convertToDoubles(matrix);
  };

  /**
   * {@inheritDoc}
   * <br>BigDecimalMatrix returns an array containing {@link Quadruple} values nearest to the values
   * of the corresponding elements of the internal array of {@link BigDecimal}
   * of the internal storage.
   */
  @Override
  public Quadruple[][] getQuadrupleData() {
    return convertToQuadruples(matrix);
  };

  /**
   * {@inheritDoc}
   * <br>BigDecimalMatrix returns an array containing the exact values of the corresponding elements
   * of the internal storage.
   */
  @Override
  public BigDecimal[][] getBigDecimalData() {
    return deepCopyOf(matrix);
  };

  /** {@inheritDoc} */
  @Override
  public boolean equals(Object anotherOne) {
    if (this == anotherOne)
      return true;
    if (!(anotherOne instanceof BigDecimalMatrix))
      return false;
    return  (Arrays.deepEquals(matrix, ((BigDecimalMatrix)anotherOne).matrix))
            && (solver.needToScale == ((BigDecimalMatrix)anotherOne).solver.needToScale
            && (precision == ((BigDecimalMatrix)anotherOne).precision));
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    final int result = Arrays.deepHashCode(matrix);
    return (result * 31 + Boolean.hashCode(solver.needToScale)) * 31 + Integer.hashCode(precision);
  };

  /* *******************************************************************************
  /***** Solutions with respect of vectors *****************************************
  /*********************************************************************************/

  /** {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * {@code BigDecimal} values, using {@link BigDecimal#valueOf(double)} method.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array
   *  returned by {@code BigDecimalMatrix} is {@code BigDecimal}.
   */
  @Override
  public  Number[] solve(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solve(double[])");
    solution = solver.solveLU(convertToBigDecimals(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * exactly equal or nearest possible {@code BigDecimal} values, depending on the particular type
   * of the {@code vector} elements.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code BigDecimalMatrix} is {@code BigDecimal}.
   */
  @Override
  public  Number[] solve(Number[] vector) throws IllegalArgumentException, NullPointerException{
    checkVector(vector, "solve(Number[])");
    solution = solver.solveLU(convertToBigDecimals(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * {@code BigDecimal} values, using {@link BigDecimal#valueOf(double)} method.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code BigDecimalMatrix} is {@code BigDecimal}.
   */
  @Override
  public  Number[] solveSPD(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPD(double[])");
    solution = solver.solveCholesky(convertToBigDecimals(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * exactly equal or nearest possible {@code BigDecimal} values, depending on the particular type
   * of the {@code vector} elements.
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code BigDecimalMatrix} is {@code BigDecimal}.
   */
  @Override
  public  Number[] solveSPD(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPD(Number[])");
    solution = solver.solveCholesky(convertToBigDecimals(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * {@code BigDecimal} values, using {@link BigDecimal#valueOf(double)} method.
   * <br>
   * The execution time is about 35% longer than that of the simple
   * {@linkplain BigDecimalMatrix#solve(double[])}.
   * In a typical case, the result error is reduced, compared to the simple {@code solve()},
   * by a factor from several units to several tens or a few hundreds, depending on the properties
   * of the specific matrix data.<br>
   * For matrices 150 x 150 with uniformly-distributed random values and {@code precision} value of 40,
   * average square root of MSE is reduced by approximately 70 times, from about 7.0e-38 down to 1.0e-39.<br>
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   * of {@linkplain Number} instances. The particular type of the elements of the array returned by
   * {@code BigDecimalMatrix} is {@code BigDecimal}.
   */
  @Override
  public  Number[] solveAccurately(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveAccurately(double[])");
    solution = solver.solveLUAccurately(convertToBigDecimals(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * exactly equal or nearest possible {@code BigDecimal} values, depending on the particular type
   * of the {@code vector} elements.
   * <br>
   * The execution time is about 35% longer than that of the simple
   * {@linkplain BigDecimalMatrix#solve(double[])}.
   * In a typical case, the result error is reduced, compared to the simple {@code solve()},
   * by a factor from several units to several tens or a few hundreds, depending on the properties
   * of the specific matrix data.<br>
   * For matrices 150 x 150 with uniformly-distributed random values and {@code precision} value of 40,
   * average square root of MSE is reduced by approximately 70 times, from about 7.0e-38 down to 1.0e-39.<br>
   *
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code BigDecimalMatrix} is {@code BigDecimal}.
   */
  @Override
  public  Number[] solveAccurately(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveAccurately(Number[])");
    solution = solver.solveLUAccurately(convertToBigDecimals(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * {@code BigDecimal} values, using {@link BigDecimal#valueOf(double)} method.
   * <br>
   * The execution time is about 65% longer than that of the simple
   * {@linkplain BigDecimalMatrix#solveSPD(double[])}.
   * In a typical case, the result error is reduced, compared to the simple {@code solveSPD()},
   * by a factor from several units to several tens, depending on the properties of the specific matrix data.<br>
   * For matrices 150 x 150 with {@code precision} value of 40, average square root of MSE
   * is reduced by approximately 11 times, from about 9.0e-38 down to 7.7e-39.<br>
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code BigDecimalMatrix} is {@code BigDecimal}.
   */
  @Override
  public  Number[] solveSPDAccurately(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPDAccurately(double[])");
    solution = solver.solveCholeskyAccurately(convertToBigDecimals(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * Before solving the equation, the values of the {@code vector} argument get converted to
   * exactly equal or nearest possible {@code BigDecimal} values, depending on the particular type
   * of the {@code vector} elements.
   * <br>
   * The execution time is about 65% longer than that of the simple
   * {@linkplain BigDecimalMatrix#solveSPD(double[])}.
   * In a typical case, the result error is reduced, compared to the simple {@code solveSPD()},
   * by a factor from several units to several tens, depending on the properties of the specific matrix data.<br>
   * For matrices 150 x 150 with {@code precision} value of 40, average square root of MSE
   * is reduced by approximately 11 times, from about 9.0e-38 down to 7.7e-39.<br>
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances. The particular type of the elements of the array returned by
   *  {@code BigDecimalMatrix} is {@code BigDecimal}.
   */
  @Override
  public  Number[] solveSPDAccurately(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "solveSPDAccurately(Number[])");
    solution = solver.solveCholeskyAccurately(convertToBigDecimals(vector));
    return getSolution();
  }

  /**
   * {@inheritDoc}
   * @return the last of the previously found vector solutions, <b>x</b>, to a system of of linear equations
   *  of form <b><i>A</i>x = b</b>, as an array of {@link BigDecimal}s, or {@code null},
   *  if no system was solved with the matrix.
   */
  @Override
  public Number[] getSolution() {
    if (solution == null) return null;

    final BigDecimal[] result = new BigDecimal[size];
    for (int i = 0; i < size; i++) {
      result[i] = solution[i];
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
   * <br>{@code BigDecimalMatrix} translates {@link BigDecimal} values of the internally stored solution to
   * {@link Quadruple} values, using {@link Quadruple#Quadruple(BigDecimal)} constructor.
   */
  @Override
  public Quadruple[] getQuadrupleSolution() {
    return convertToQuadruples(solution);
  }

  /**
   * {@inheritDoc}
   * {@code BigDecimalMatrix} returns the exact values that was found by the corresponding solution method.
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
   * as a new instance of {@code BigDecimalMatrix}.
   */
  @Override
  public Matrix solve(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solve(Matrix)");
    matrixSolution = solver.solve(matrixB.getBigDecimalData());
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   * as a new instance of {@code BigDecimalMatrix}.
   */
  @Override
  public Matrix solve(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solve(double[][])");
    matrixSolution = solver.solve(convertToBigDecimals(matrixB));
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   * as a new instance of {@code BigDecimalMatrix}.
   */
  @Override
  public Matrix solve(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solve(Number[][])");
    matrixSolution = solver.solve(convertToBigDecimals(matrixB));
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * The execution time is about 60 times longer than that of the simple
   * {@linkplain BigDecimalMatrix#solve(Matrix)}. In a typical case, the result error is reduced,
   * compared to the simple {@code solve()}, by a factor from units to several tens or more,
   * depending on specific values of the matrix elements. For matrices 150 x 150
   * with uniformly-distributed random values and {@code precision} value of 40,
   * average square root of MSE is reduced by approximately 100 times,
   * from about 2.0e-37 down to 2.0e-39.<br>
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   *  as a new instance of {@code BigDecimalMatrix}.
   */
  @Override
  public Matrix solveAccurately(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solveAccurately(Matrix)");
    return solveAccurately(matrixB.getBigDecimalData());
  }

  /**
   * {@inheritDoc}
   * <br>
   * The execution time is about 60 times longer than that of the simple
   * {@linkplain BigDecimalMatrix#solve(Matrix)}. In a typical case, the result error is reduced,
   * compared to the simple {@code solve()}, by a factor from units to several tens or more,
   * depending on specific values of the matrix elements. For matrices 150 x 150
   * with uniformly-distributed random values and {@code precision} value of 40,
   * average square root of MSE is reduced by approximately 100 times,
   * from about 2.0e-37 down to 2.0e-39.<br>
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   *  as a new instance of {@code BigDecimalMatrix}.
   */
  @Override
  public Matrix solveAccurately(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solveAccurately(double[][])");
    matrixSolution = solver.solveAccurately(convertToBigDecimals(matrixB));
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   * <br>
   * The execution time is about 60 times longer than that of the simple
   * {@linkplain BigDecimalMatrix#solve(Matrix)}. In a typical case, the result error is reduced,
   * compared to the simple {@code solve()}, by a factor from units to several tens or more,
   * depending on specific values of the matrix elements. For matrices 150 x 150
   * with uniformly-distributed random values and {@code precision} value of 40,
   * average square root of MSE is reduced by approximately 100 times,
   * from about 2.0e-37 down to 2.0e-39.<br>
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   *  as a new instance of {@code BigDecimalMatrix}.
   */
  @Override
  public Matrix solveAccurately(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "solveAccurately(Number[][])");
    matrixSolution = solver.solveAccurately(convertToBigDecimals(matrixB));
    return getMatrixSolution();
  }

  /**
   * {@inheritDoc}
   *
   * @return the last of the previously found matrix solution <b><i>X</i></b> to a system of
   * linear equations of form <b><i>AX = B</i></b>, as a {@code BigDecimalMatrix}, or {@code null},
   * if no matrix equation was solved with the matrix.
   */
  @Override
  public Matrix getMatrixSolution() {
    if (matrixSolution == null) return null;
    return new BigDecimalMatrix(matrixSolution); // Sic! The constructor makes a copy to preserve the solution
  }

  /**
   * {@inheritDoc}
   * <br>{@code BigDecimalMatrix#getNumberMatrixSolution()} returns the exact values of the elements
   * of the matrix solution found by a call to one of the {@code BigDecimalMatrix#solve()} or
   * {@code BigDecimalMatrix#solveAccurately()} methods.
   *
   * @return the last of the previously found matrix solution <b><i>X</i></b>, to a system of
   *  linear equations of form <b><i>AX = B</i></b>, as a two-dimentional array of {@linkplain BigDecimal},
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
   * to one of the {@code BigDecimalMatrix#solve()} or {@code BigDecimalMatrix#solveAccurately()} methods,
   * rounded to nearest {@code double} values.
   */
  @Override
  public double[][] getDoubleMatrixSolution() {
    if (matrixSolution == null) return null;
    return convertToDoubles(matrixSolution);
  }

  /**
   * {@inheritDoc}
   * <br>The returned array contains {@link Quadruple} values nearest to the corresponding
   * {@link BigDecimal} elements of the matrix solution found by a call to one of the
   * {@code BigDecimalMatrix#solve()} or {@code BigDecimalMatrix#solveAccurately()} methods.
   */
  @Override
  public Quadruple[][] getQuadrupleMatrixSolution() {
    if (matrixSolution == null) return null;
    return convertToQuadruples(matrixSolution);
  }

  /**
   * {@inheritDoc}
   * <br>The returned array contains the exact values of the elements of the matrix solution
   * found by a call to one of the {@code BigDecimalMatrix#solve()} or
   * {@code BigDecimalMatrix#solveAccurately()} methods.
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
   * <br>The particular subtype of {@code Matrix} returned by {@code BigDecimalMatrix#inverse()}
   * is {@code BigDecimalMatrix}.
   */
  @Override
  public Matrix inverse() throws IllegalArgumentException {
    return newWithThisArray(solver.inverse());
  }

  /**
   * {@inheritDoc}
   * <br>The execution time is about 30 times longer than that of the simple
   * {@linkplain BigDecimalMatrix#inverse()}.
   * In a typical case, the result error is reduced, compared to the simple {@code inverse()},
   * by a factor from several units to several tens, depending on specific values of the matrix elements.
   * For matrices 150 x 150 with uniformly-distributed random values, average square root of MSE
   * is reduced by approximately 25 times, from about 4.5e-38 down to about 1.8e-39<br>
   *
   * @return a new {@code BigDecimalMatrix} containing the inversion of the given matrix.
   */
  @Override
  public Matrix inverseAccurately() throws IllegalArgumentException {
    return newWithThisArray(solver.inverseAccurately());
  };

  /**
   * {@inheritDoc}
   *
   * @return a new {@code BigDecimalMatrix} containing the transposition of the given matrix.
   */
  @Override
  public Matrix transpose() {
    return newWithThisArray(solver.transpose(matrix));
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code BigDecimalMatrix} containing the unity matrix of the same size as this {@code Matrix}
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
   * @return a new {@code BigDecimalMatrix} representing the product
   */
  @Override
  public Matrix multiply(Matrix factor) throws IllegalArgumentException, NullPointerException {
    checkMatrix(factor, "multiply(Matrix)");
    final BigDecimal[][] product = solver.multiply(factor.getBigDecimalData());
    return newWithThisArray(product);
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code BigDecimalMatrix} representing the product
   */
  @Override
  public Matrix multiply(double[][] factor) throws IllegalArgumentException, NullPointerException {
    checkMatrix(factor, "multiply(double[][])");
    final BigDecimal[][] product = solver.multiply(convertToBigDecimals(factor));
    return newWithThisArray(product);
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code BigDecimalMatrix} representing the product
   */
  @Override
  public Matrix multiply(Number[][] factor) throws IllegalArgumentException, NullPointerException {
    checkMatrix(factor, "multiply(Number[][])");
    final BigDecimal[][] product = solver.multiply(convertToBigDecimals(factor));
    return newWithThisArray(product);
  }

  /**
   * {@inheritDoc}
   *
   * @return an array of {@linkplain BigDecimal} values containing the product.
   */
  @Override
  public Number[] multiply(double[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "multiply(double[])");
    return solver.multiply(convertToBigDecimals(vector));
  }

  /**
   * {@inheritDoc}
   *
   * @return an array of {@linkplain BigDecimal} values containing the product.
   */
  @Override
  public Number[] multiply(Number[] vector) throws IllegalArgumentException, NullPointerException {
    checkVector(vector, "multiply(Number[])");
    return solver.multiply(convertToBigDecimals(vector));
  }

  /**
   * {@inheritDoc}
   *
   * @return a new {@code BigDecimalMatrix} containing the product of the source matrix and the given scalar.
   */
  @Override
  public Matrix multiply(double scalar) throws IllegalArgumentException {
    checkScalar(scalar, "multiply(double)");
    return newWithThisArray(solver.multiply(BigDecimal.valueOf(scalar)));
  }

  /**
   * {@inheritDoc}
   * Before the multiplication, the value of the argument is translated
   * into the corresponding {@code BigDecimal} value by a method depending on the actual
   * type of the argument.
   *
   * @return a new {@code BigDecimalMatrix} containing the product of the source matrix and the given scalar.
   */
  @Override
  public Matrix multiply(Number scalar) throws IllegalArgumentException, NullPointerException {
    checkScalar(scalar, "multiply(Number)");
    BigDecimal factor;
    if (scalar instanceof Double)
      factor = new BigDecimal((double)scalar);
    else if (scalar instanceof Quadruple)
      factor = ((Quadruple)scalar).bigDecimalValue();
    else if (scalar instanceof BigDecimal)
      factor = (BigDecimal)scalar;
    else
      throw new IllegalArgumentException("only Double, Quadruple and BigDecimal values are allowed for Number arguments");
    return newWithThisArray(solver.multiply(factor));
  }

  /* *******************************************************************************
  /* **** Additions and subtractions ***********************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   * @return a new {@code BigDecimalMatrix} containing the sum of the two matrices.
   */
  @Override
  public Matrix add(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "add(Matrix)");
    return newWithThisArray(solver.addMatrices(this.matrix, matrixB.getBigDecimalData()));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code BigDecimalMatrix} containing the sum of the two matrices.
   */
  @Override
  public Matrix add(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "add(double[][])");
    return newWithThisArray(solver.addMatrices(this.matrix, convertToBigDecimals(matrixB)));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code BigDecimalMatrix} containing the sum of the two matrices.
   */
  @Override
  public Matrix add(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "add(Number[][])");
    return newWithThisArray(solver.addMatrices(this.matrix, convertToBigDecimals(matrixB)));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code BigDecimalMatrix} containing the difference of the two matrices.
   */
  @Override
  public Matrix subtract(Matrix matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "subtract(Matrix)");
    return newWithThisArray(solver.subtractMatrices(this.matrix, matrixB.getBigDecimalData()));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code BigDecimalMatrix} containing the difference of the two matrices.
   */
  @Override
  public Matrix subtract(double[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "subtract(double[][])");
    return newWithThisArray(solver.subtractMatrices(this.matrix, convertToBigDecimals(matrixB)));
  }

  /**
   * {@inheritDoc}
   * @return a new {@code BigDecimalMatrix} containing the difference of the two matrices.
   */
  @Override
  public Matrix subtract(Number[][] matrixB) throws IllegalArgumentException, NullPointerException {
    checkMatrix(matrixB, "subtract(Number[][])");
    return newWithThisArray(solver.subtractMatrices(this.matrix, convertToBigDecimals(matrixB)));
  }

  /* *******************************************************************************
  /* **** Determinant **************************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   * @return the value of the determinant of the given {@code Matrix} as an instance
   * of class {@code BigDecimal}
   */
  @Override
  public Number determinant() {
    return solver.determinant();
  }

  /**
   * {@inheritDoc}
   * {@code BigDecimalMatrix#determinantAsDouble()} returns the result of the corresponding
   * computations rounded to the nearest {@code double} value.
   */
  @Override
  public double determinantAsDouble() {
    return solver.determinant().doubleValue();
  }

  /**
   * {@inheritDoc}
   * {@code BigDecimalMatrix#determinantAsQuadruple()} returns a {@code Quadruple} value
   * closest to the internally calculated {@code BigDecimal} value.
   */
  @Override
  public Quadruple determinantAsQuadruple() {
    return new Quadruple(solver.determinant());
  }

  /**
   * {@inheritDoc}
   * For {@code BigDecimalMatrix}, it is defined by the {@code precision} value of this
   * {@code BigDecimalMatrix} instance.
   */
  @Override
  public BigDecimal determinantAsBigDecimal() {
    return solver.determinant();
  }

  /* *******************************************************************************
  /* **** Norm *********************************************************************
  /*********************************************************************************/

  /**
   * {@inheritDoc}
   * @return the value of the norm of the given {@code Matrix} as an instance
   * of class {@code BigDecimal}
   */
  @Override
  public Number norm() {
    return solver.norm();
  }

  /**
   * {@inheritDoc}
   * {@code BigDecimalMatrix#normAsDouble()} returns the result of the corresponding
   * computations rounded to the nearest {@code double} value.
   */
  @Override
  public double normAsDouble() {
    return solver.norm().doubleValue();
  }

  /**
   * {@inheritDoc}
   * {@code BigDecimalMatrix#normAsQuadruple()} returns a {@code Quadruple} value closest to
   * the internally calculated {@link BigDecimal} value.
   */
  @Override
  public Quadruple normAsQuadruple() {
    return new Quadruple(solver.norm());
  }

  /**
   * {@inheritDoc}
   * {@code BigDecimalMatrix} returns a {@code BigDecimal} value exactly
   * equal to the internally calculated value.
   */
  @Override
  public BigDecimal normAsBigDecimal() {
    return solver.norm();
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

  /**
   * Sets the precision field and creates an appropriate MathContext to be used by calculations
   * @param precision
   */
  private void setPrecision(int precision) {
    this.precision = precision;
    mc = new MathContext(precision, RoundingMode.HALF_EVEN);
  }

  private final BigDecimalMatrixSolver makeSolver(boolean needToScale) {
    solver = new BigDecimalMatrixSolver(matrix, needToScale, mc);
    return solver;
  }

  private final void copyDataFrom(double[][] source) {
    this.size = source.length;
    matrix = convertToBigDecimals(source);
  }

  private final void copyDataFrom(Number[][] source) {
    this.size = source.length;
    matrix = convertToBigDecimals(source);
  }

}
