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
import java.util.Objects;

import com.mvohm.quadruple.Quadruple;

// Started 21.06.20 10:05:16


/**
 * A generic abstract class which defines the set of operations that all its subclasses must implement.<br>
 * Instances of subclasses are capable of:
 * <ul>
 *  <li>solving systems of linear equations of the forms <b>A * X = B</b> and <b>A * x = b</b>,
 *    including versions with enhanced accuracy that use an iterative refinement of the solution
 *    (see {@linkplain #solve(Number[])}, {@linkplain #solveAccurately(Number[])}, {@linkplain #solve(Matrix)} and alike);
 *  <li>inversion (including a version with enhanced accuracy) and transposition of the matrix;
 *  <li>multiplying this matrix by another matrix, by a vector, and by a scalar;
 *  <li>addition and subtraction of a matrix;
 *  <li>computation of the determinant.
 * </ul>
 * Instances of the matrix are immutable; all operations whose results are matrices
 * create and return new instances representing the resulting matrices.
 * <br>
 * @author M.Vokhmentev
 */
public abstract class Matrix {

  /**
   * Enumeration including error codes that may be set in case of unsuccessful attempt to solve a system of linear equations
   * or to inverse the matrix.
   * <br>A string designation of the Matrix's error code after an attempt to solve a system or to inverse a matrix
   * can be obtained using {@linkplain Matrix#getErrorCode()} method of the instance. This method can return the following values:
   * <ul>
   *  <li><b>{@code "OK"}</b> -- The last solution or inversion was successful;
   *  <li><b>{@code "ASYMMETRIC"}</b> -- There was an attempt to solve or inverse the matrix using Cholesky decomposition, but the matrix is asymmetric;
   *  <li><b>{@code "NON_SPD"}</b> -- There was an attempt to solve or inverse the matrix using Cholesky decomposition, but the matrix is not positively-defined;
   *  <li><b>{@code "NON_INVERTIBLE"}</b> -- There was an attempt to solve or inverse the matrix using LU decomposition, but the matrix is inconsistent or underdetermined
   * </ul>   *
   * @see Matrix#getErrorCode()
   */
  protected enum ErrorCodes {
    /** The last solution or inversion was successful */
    OK,
    /** There was an attempt to solve or inverse the matrix using Cholesky decomposition, but the matrix is asymmetric */
    ASYMMETRIC,
    /** There was an attempt to solve or inverse the matrix using Cholesky decomposition, but the matrix is not positively-defined */
    NON_SPD,
    /** There was an attempt to solve or inverse the matrix using LU decomposition, but the matrix is inconsistent or underdetermined */
    NON_INVERTIBLE
    };

  protected int size;

  /** Used by private constructors of the subclasses. Creates an instance with empty fields. */
  Matrix() {
  }

//==== Getting data =============================

  /**
   * Returns the size of the matrix, <code>m</code> for a matrix <code>m x m</code>
   * @return the size of the matrix
   */
  public int getSize() {
    return size;
  }

  /**
   * Returns the value of an internal flag that defines whether row scaling will be applied while solving a system by this instance of Matrix.
   * <br>Scaling of the rows of the matrix may be used while solving systems of linear equations.
   * When an internal flag controlling the scaling is set to <code>true</code>,
   * the rows of the matrix along with the corresponding elements of the vector <b>b</b> or the matrix <b>B</b>
   * are scaled so that the norms of the rows are all be 1.0. In most cases this improves the accuracy of the solution,
   * especially for matrices including both very large and very small elements. The default value of the flag signifying
   * the necessity of the scaling for a certain subclass may be set via a static method <code>setScaling()</code> of the subclass,
   * and the necessity of the scaling for a newly-created instance of Matrix may be controlled via
   * the corresponding parameter of a constructor.
   * @return the value of the flag that defines whether the scaling will be applied while solving a system by this instance of Matrix
   */
  public abstract boolean         getScaling();

  /**
   * Returns a two-dimentional array of {@link Number} containing the values of the corresponding matrix elements.
   * <br>The exact type of the array elements depends on the particular subclass of this instance
   * (it is Double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * The returned array is a copy of the corresponding internal array and can be safely modified.
   * There is no loss of precision when copying, so the returned values are exactly equal
   * to the values of the corresponding elements of the internal storage.
   * @return a two-dimentional array of {@link Number} containing the values of the corresponding matrix elements.
   */
  public abstract Number[][]      getData();

  /**
   * Returns a two-dimentional array of primitive <b>{@code double}</b> values containing the values
   * of the corresponding matrix elements, perhaps rounded.
   * <br>If the type of the internal data of the particular class enables for higher precision
   * than that provided by <b>{@code double}</b>, the values get rounded to the nearest possible
   * <b>{@code double}</b> value.<br>
   * If the value of an internal data element exceeds the range of <b>{@code double}</b> values,
   * it gets converted to {@linkplain Double#POSITIVE_INFINITY} or {@linkplain Double#NEGATIVE_INFINITY}, depending on its sign.
   * @return a two-dimentional array of primitive <b>{@code double}</b> values containing the values
   * of the corresponding matrix elements, perhaps rounded.
   */
  public abstract double[][]      getDoubleData();

  /**
   * Returns a two-dimentional array of {@link Quadruple} instances containing the values
   * of the corresponding matrix elements, perhaps rounded.
   * <br>If the type of the internal data of the particular class enables for higher precision
   * than that provided by <b>{@code Quadruple}</b>, the values get rounded to the
   * nearest possible <b>{@code Quadruple}</b> value.<br>
   * If the value of an internal data element exceeds the range of <b>{@code Quadruple}</b> values,
   * it gets converted to {@linkplain Quadruple#POSITIVE_INFINITY} or {@linkplain Quadruple#NEGATIVE_INFINITY}, depending on its sign.
   * @return a two-dimentional array of {@link Quadruple} values containing the values
   * of the corresponding matrix elements, perhaps rounded.
   */
  public abstract Quadruple[][]   getQuadrupleData();

  /**
   * Returns a two-dimentional array of {@link BigDecimal} instances containing the values
   * of the corresponding matrix elements.
   * <br>If the type of the internal data of the particular class is {@code double}
   * then the corresponding {@code BigDecimal} values are obtained using {@linkplain BigDecimal#valueOf(double)} method.
   * If the type of the internal data of the particular class is {@code Quadruple}
   * then the corresponding {@code BigDecimal} values are obtained using {@linkplain Quadruple#bigDecimalValue()} method.
   * If an internal data element is not convertible to {@link BigDecimal} (i.e. it is NaN or Infinity), throws
   * {@link NumberFormatException}.
   * @return a two-dimentional array of {@link BigDecimal} values containing the values
   * of the corresponding matrix elements translated to {@link BigDecimal} instances.
   */
  public abstract BigDecimal[][]  getBigDecimalData();

  /**
   * Indicates whether the other {@code Matrix} is equal to this one.
   * <br>Matrices are considered to be equal if they belong to the same subtype,
   * and their internal arrays containing the elements of the matrices are equal,
   * and their {@code needToScale} flags are equal.<br> Under those condition,
   * the matrices yield equal results for all operations performed on them.
   * @return true, if the matrices are equal.
   */
  @Override
  public abstract boolean equals(Object anotherOne);

  /**
   * Returns a hash code value for the {@code Matrix}.
   * <br>It is guaranteed, that for two matrices returning different hashcodes
   * their {@link #equals(Object)} methods return {@code false},
   * and two matrices considered to be equal return equal hashcodes.
   * For two different matrices, the probability of the equality of their hashcodes is reasonably low.
   * @return a hash code value for this object.
   * @see #equals(Object)
   */
  @Override
  public abstract int hashCode();

  /* *******************************************************************************
  /***** Solutions with respect of vectors *****************************************
  /*********************************************************************************/

  /**
   * Solves a system of linear equations of form <b><i>A</i>x = b</b> and returns the found solution.
   * <br>
   * Solves the system <b><i>A</i>x = b</b>, formed by the inner matrix data <b><i>A</i></b> and the
   * vector <b>b</b> passed in as the {@code vector} argument, using LU decomposition.<br>
   * The exact type of the elements of the returned array depends on the particular subclass of this instance
   * (it is Double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * In case of an inconsistent or underdefined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the length of the {@code vector} does not match the size of the matrix,
   * or the argument contains non-numeric values (NaN or Infinity),
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the length of the {@code vector} does not match the size
   *  of the matrix,, or the system has no solution or has infinitely many solutions,
   *  or the argument contains non-numeric values (NaN or Infinity).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param vector the column vector <b>b</b> of the equation to be solved, <b><i>A</i>x = b</b>, as an array of {@code double}s
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of {@linkplain Number} instances
   *  whose particular type depends on the particular {@code Matrix} subtype.
   *  @see #getErrorCode()
   */
  public abstract Number[]        solve(double[] vector);

  /**
   * Solves a system of linear equations of form <b><i>A</i>x = b</b> and returns the found solution.
   * <br>
   * Solves the system <b><i>A</i>x = b</b>,
   * formed by the inner matrix data <b><i>A</i></b> and the
   * vector <b>b</b> passed in as the {@code vector} argument, using LU decomposition.
   * Before solving, the given vector gets translated to an array of values of the type that is used
   * to store internal matrix data (double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * The exact type of the elements of the returned array depends on the particular subclass of this instance
   * (it is Double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * In case of an inconsistent or underdefined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the length of the {@code vector} does not match the size of the matrix,
   * or {@code vector} contains non-numeric values ({@code NaN} or {@code Infinity}) or {@code null},
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the length of the {@code vector} does not match the size
   *  of the matrix, or the system has no solution or has infinitely many solutions,
   *  or the argument contains invalid values (NaN, Infinity, or null).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param vector the column vector <b>b</b> of the equation to be solved, <b><i>A</i>x = b</b>,
   *  as an array of {@code Number}s
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances, whose particular type depends on the particular
   *  {@code Matrix} subtype.
   *
   * @see #getErrorCode()
   */
  public abstract Number[]        solve(Number[] vector);

  /**
   * Solves a system of linear equations of form <b><i>A</i>x = b</b> for a symmetric positively-defined matrix of coefficients
   * and returns the found solution.
   * <br>
   * Solves the system <b><i>A</i>x = b</b>,
   * formed by the inner matrix data <b><i>A</i></b> and the
   * vector <b>b</b> passed in as the {@code vector} argument, using Cholesky decomposition.
   * The exact type of the elements of the returned array depends on the particular subclass of this instance
   * (it is Double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * In case of an asymmetric or non-positively-defined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the length of the {@code vector} does not match the size of the matrix,
   * or the argument contains non-numeric values (NaN or Infinity),
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the length of the {@code vector} does not match the size
   *  of the matrix, or the matrix is asymmetric or not positively-defined,
   *  or the argument contains non-numeric values (NaN or Infinity).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param vector the column vector <b>b</b> of the equation to be solved, <b><i>A</i>x = b</b>,
   *  as an array of {@code Number}s
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances, whose particular type depends on the particular
   *  {@code Matrix} subtype.
   *
   * @see #getErrorCode()
   */
  public abstract Number[]        solveSPD(double[] vector);

  /**
   * Solves a system of linear equations of form <b><i>A</i>x = b</b> for a symmetric positively-defined matrix of coefficients
   * and returns the found solution.
   * <br>
   * Solves the system <b><i>A</i>x = b</b>,
   * formed by the inner matrix data <b><i>A</i></b> and the
   * vector <b>b</b> passed in as the {@code vector} argument, using Cholesky decomposition.
   * Before solving, the given vector gets translated to an array of values of the type that is used
   * to store internal matrix data (double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).
   * The exact type of the elements of the returned array depends on the particular subclass of this instance
   * (it is Double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * In case of an asymmetric or non-positively-defined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the length of the {@code vector} does not match the size of the matrix,
   * or {@code vector} contains non-numeric values ({@code NaN} or {@code Infinity}) or {@code null},
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the length of the {@code vector} does not match
   *  the size of the matrix, or the matrix is asymmetric or not positively-defined,
   *  or the argument contains invalid values (NaN, Infinity, or null)
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param vector the column vector <b>b</b> of the equation to be solved, <b><i>A</i>x = b</b>,
   *  as an array of {@code Number}s
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances, whose particular type depends on the particular
   *  {@code Matrix} subtype.
   *
   * @see #getErrorCode()
   */
  public abstract Number[]        solveSPD(Number[] vector);

  /**
   * Solves a system of linear equations of form <b><i>A</i>x = b</b> with increased accuracy and returns the found solution.
   * <br>
   * Solves the system <b><i>A</i>x = b</b>,
   * formed by the inner matrix data <b><i>A</i></b> and the
   * vector <b>b</b> passed in as the {@code vector} argument, using LU decomposition.
   * Uses an iterative refinement to find a more accurate solution.
   * The exact type of the elements of the returned array depends on the particular subclass of this instance
   * (it is Double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * In case of an inconsistent or underdefined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the length of the {@code vector} does not match the size of the matrix,
   * or {@code vector} contains non-numeric values (NaN or Infinity),
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the length of the {@code vector} does not match
   *  the size of the matrix, or the system has no solution or has infinitely many solutions,
   *  or the argument contains non-numeric values (NaN or Infinity).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param vector the column vector <b>b</b> of the equation to be solved, <b><i>A</i>x = b</b>,
   *  as an array of {@code double}s
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances, whose particular type depends on the particular
   *  {@code Matrix} subtype.
   *
   * @see #getErrorCode()
   */
  public abstract Number[]        solveAccurately(double[] vector);

  /**
   * Solves a system of linear equations of form <b><i>A</i>x = b</b> with increased accuracy and returns the found solution.
   * <br>
   * Solves the system <b><i>A</i>x = b</b>,
   * formed by the inner matrix data <b><i>A</i></b> and the
   * vector <b>b</b> passed in as the {@code vector} argument, using LU decomposition.
   * Before solving, the given vector gets translated to an array of values of the type that is used
   * to store internal matrix data (double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * Uses an iterative refinement to find a more accurate solution.
   * The exact type of the elements of the returned array depends on the particular subclass of this instance
   * (it is Double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * In case of an inconsistent or underdefined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the length of the {@code vector} does not match the size of the matrix,
   * or {@code vector} contains non-numeric values ({@code NaN} or {@code Infinity}) or {@code null},
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the length of the {@code vector} does not match
   *  the size of the matrix, or the system has no solution or has infinitely many solutions,
   *  or the argument contains invalid values (NaN, Infinity, or null)
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param vector the column vector <b>b</b> of the equation to be solved, <b><i>A</i>x = b</b>, as an array of {@code Number}s
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array of
   *  {@linkplain Number} instances, whose particular type depends on the particular
   *  {@code Matrix} subtype.
   *
   * @see #getErrorCode()
   */
  public abstract Number[]        solveAccurately(Number[] vector);

  /**
   * Solves a system of linear equations of form <b><i>A</i>x = b</b> for a symmetric positively-defined matrix
   * of coefficients, using an iterative refinement algorithm to achieve higher solution accuracy, and returns the found solution.
   * <br>
   * Solves the system <b><i>A</i>x = b</b>,
   * formed by the inner matrix data <b><i>A</i></b> and the
   * vector <b>b</b> passed in as the {@code vector} argument, using Cholesky decomposition.
   * Uses an iterative refinement to find a more accurate solution.
   * The exact type of the elements of the returned array depends on the particular subclass of this instance
   * (it is Double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * In case of an asymmetric or non-SPD matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the length of the {@code vector} does not match the size of the matrix,
   * or {@code vector} contains non-numeric values (NaN or Infinity),
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the length of the {@code vector} does not match
   *  the size of the matrix, or the matrix is asymmetric or not positively-defined,
   *  or the argument contains non-numeric values (NaN or Infinity).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param vector the column vector <b>b</b> of the equation to be solved, <b><i>A</i>x = b</b>, as an array of {@code double}s
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances, whose particular type depends on the particular
   *  {@code Matrix} subtype.
   *
   *  @see #getErrorCode()
   */
  public abstract Number[]        solveSPDAccurately(double[] vector);

  /**
   * Solves a system of linear equations of form <b><i>A</i>x = b</b> for a symmetric positively-defined matrix
   * of coefficients, using an iterative refinement algorithm to achieve higher solution accuracy, and returns the found solution.
   * <br>
   * Solves the system <b><i>A</i>x = b</b>,
   * formed by the inner matrix data <b><i>A</i></b> and the
   * vector <b>b</b> passed in as the {@code vector} argument, using Cholesky decomposition.
   * Uses an iterative refinement to find a more accurate solution.
   * The exact type of the elements of the returned array depends on the particular subclass of this instance
   * (it is Double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * In case of an asymmetric or non-SPD matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.
   *
   * If the length of the {@code vector} does not match the size of the matrix, or {@code vector} contains
   * non-numeric values ({@code NaN} or {@code Infinity}) or {@code null} throws {@code IllegalArgumentException}
   * with a relevant message.
   *
   * @throws IllegalArgumentException if the length of the {@code vector} does not match
   *  the size of the matrix, or the matrix is asymmetric or not positively-defined,
   *  or the argument contains invalid values (NaN, Infinity, or null)
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param vector the column vector <b>b</b> of the equation to be solved, <b><i>A</i>x = b</b>,
   * as an array of {@code double}s
   *
   * @return the found solution <b>x</b> to the equation <b><i>A</i>x = b</b> as an array
   *  of {@linkplain Number} instances, whose particular type depends on the particular
   *  {@code Matrix} subtype.
   *
   *  @see #getErrorCode()
   */
  public abstract Number[]        solveSPDAccurately(Number[] vector);

  /**
   * Returns a copy of the last previously found vector solutions, <b>x</b>, to a system of linear equations
   * of form <b><i>A</i>x = b</b>, as an array of {@linkplain Number}s.
   * <br>
   * The exact type of the elements of the returned array depends on the particular subclass of this instance
   * (it is Double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * If no system was solved with this instance of the {@code Matrix}, returns {@code null}.
   *
   * @return the last of the previously found vector solutions, <b>x</b>, to a system of linear equations
   *  of form <b><i>A</i>x = b</b>, as an array of {@linkplain Number}s,
   *  or {@code null}, if no system was solved with the matrix.
   */
  public abstract Number[]        getSolution();

  /**
   * Returns a copy of the last previously found vector solutions, <b>x</b>, to a system of linear equations
   * of form <b><i>A</i>x = b</b>, as an array of {@code double}s.
   * <br>
   * If the particular subtype of the instance
   * stores internal data with higher precision than that provided by {@code double}, the values
   * of the elements of the returned vector get rounded to the nearest {@code double} values.
   * <br>
   * If no system was solved with this instance of the {@code Matrix}, returns {@code null}.
   *
   * @return the last of the previously found vector solutions, <b>x</b>, to a system of linear equations
   *  of form <b><i>A</i>x = b</b>, as an array of primitive {@code double}s.
   *  or {@code null}, if no system was solved with the matrix.
   */
  public abstract double[]        getDoubleSolution();

  /**
   * Returns the last of the previously found vector solutions, <b>x</b>, to a system of linear equations
   * of form <b><i>A</i>x = b</b>, as an array of {@linkplain Quadruple}s.
   * <br>
   * If the particular subtype of the instance
   * stores internal data with higher precision than that provided by {@link Quadruple}, the values
   * of the elements of the returned vector get rounded to the nearest {@code double} values.
   * <br>
   * If no system was solved with this instance of the {@code Matrix}, returns {@code null}.
   *
   * @return the last of the previously found vector solutions, <b>x</b>, to a system of
   *  linear equations of form <b><i>A</i>x = b</b>, as an array of {@linkplain Quadruple}s,
   *  or {@code null}, if no system was solved with the matrix.
   */
  public abstract Quadruple[]     getQuadrupleSolution();

  /**
   * Returns the last of the previously found vector solutions, <b>x</b>, to a system of linear equations
   * of form <b><i>A</i>x = b</b>, as an array of {@linkplain BigDecimal}s.
   * <br>
   * If no system was solved with this instance of the {@code Matrix}, returns {@code null}.<br>
   * If the solution contains values that can't be translated to {@code BigDeciaml} ({@code NaN} or {@code Infinity}),
   * throws {@link NumberFormatException}.
   *
   * @return the last of the previously found vector solutions, <b>x</b>, to a system of linear equations
   * of form <b><i>A</i>x = b</b>, as an array of {@linkplain BigDecimal}s.
   * or {@code null}, if no system was solved with the matrix.
   */
  public abstract BigDecimal[]    getBigDecimalSolution();

  /**
   * Returns a string designation of the error code if an error has occurred
   * during the solving or inversion of the matrix.
   * <br>
   * This method can return the following values:
   * <ul>
   *  <li><b>{@code "OK"}</b> -- The last solution or inversion was successful;
   *  <li><b>{@code "ASYMMETRIC"}</b> -- There was an attempt to solve or inverse the matrix using
   *          Cholesky decomposition, but the matrix is asymmetric;
   *  <li><b>{@code "NON_SPD"}</b> -- There was an attempt to solve or inverse the matrix using Cholesky
   *          decomposition, but the matrix is not positively-defined;
   *  <li><b>{@code "NON_INVERTIBLE"}</b> -- There was an attempt to solve or inverse the matrix using
   *          LU decomposition, but the matrix is inconsistent or underdetermined
   * </ul>
   * If no operations requiring decompositions were performed with this instance yet, returns OK.
   *
   * @return a string designation of the error code if an error has occurred during the solving
   *  or inversion of the matrix, or "OK" in case of success.
   *
   * @see ErrorCodes */
  public abstract String          getErrorCode();

  /* *******************************************************************************
  /***** Solutions with respect of matrices ****************************************
  /*********************************************************************************/

  /**
   * Solves a matrix equation of form <b><i>AX = B</i></b> and returns the found solution.
   * <br>
   * Solves a matrix equation of form <b><i>AX = B</i></b>
   * formed by the inner matrix data <b><i>A</i></b> and the matrix <b><i>B</i></b> passed in as the {@code Matrix matrixB} argument,
   * using LU decomposition.<br>
   * In case of an inconsistent or underdefined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.
   * If the size of the {@code matrixB} does not match the size of this instance of the matrix,
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the size of the {@code matrixB} does not match
   *  the size of this instance or the system has no solution or has infinitely many solutions.
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix <b><i>B</i></b> of the equation to be solved, <b><i>AX = B</i></b>
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b> as a new instance
   *  of Matrix of the same particular subtype as this instance.
   *
   *  @see #getErrorCode()
   */
  public abstract Matrix          solve(Matrix matrixB);

  /**
   * Solves a matrix equation of form <b><i>AX = B</i></b> and returns the found solution.
   * <br>
   * Solves a matrix equation of form <b><i>AX = B</i></b>
   * formed by the inner matrix data <b><i>A</i></b> and the matrix <b><i>B</i></b> passed in as the {@code double[][] matrixB} argument,
   * using LU decomposition.<br>
   * In case of an inconsistent or underdefined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or the argument contains non-numeric values (NaN or Infinity),
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the passed-in array is non-square,
   *  or its size does not match the size of this instance,
   *  or the system has no solution or has infinitely many solutions,
   *  or the argument contains non-numeric values ({@code NaN} or {@code Infinity}).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix <b><i>B</i></b> of the equation to be solved, <b><i>AX = B</i></b>,
   *  as a two-dimentional array of primitive {@code double}s
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b> as a new instance
   *  of Matrix of the same particular subtype as this instance.
   *
   *  @see #getErrorCode()
   */
  public abstract Matrix          solve(double[][] matrixB);

  /**
   * Solves a matrix equation of form <b><i>AX = B</i></b> and returns the found solution.
   * <br>
   * Solves a matrix equation of form <b><i>AX = B</i></b> formed by the inner matrix data <b><i>A</i></b>
   * and the matrix <b><i>B</i></b> passed in as the {@code Number[][] matrixB} argument,
   * using LU decomposition.<br>
   *
   * Before solving, the given {@code matrixB} array gets translated to an array of values
   * of the type that is used to store internal matrix data (double for {@linkplain DoubleMatrix},
   * {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   *
   * In case of an inconsistent or underdefined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or {@code matrixB} contains non-numeric values ({@code NaN} or {@code Infinity}) or {@code null},
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the passed-in array is non-square,
   *  or its size does not match the size of this instance,
   *  or the system has no solution or has infinitely many solutions,
   *  or the argument contains invalid values (NaN, Infinity, or null).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix <b><i>B</i></b> of the equation to be solved, <b><i>AX = B</i></b>,
   *  as a two-dimentional array of {@linkplain Number}
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b> as a new instance
   *  of Matrix of the same particular subtype as this instance.
   *
   *  @see #getErrorCode()
   */
  public abstract Matrix          solve(Number[][] matrixB);

  /**
   * Solves a matrix equation of form <b><i>AX = B</i></b> and returns the found solution.
   * <br>
   * Solves a matrix equation of form <b><i>AX = B</i></b> formed by the inner matrix data
   * <b><i>A</i></b> and the matrix <b><i>B</i></b> passed in as the {@code matrixB} argument,
   * using LU decomposition.<br>
   * Uses an iterative refinement to find a more accurate solution.<br>
   *
   * In case of an inconsistent or underdefined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the size of the {@code matrixB} does not match the size of this instance of the matrix,
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the size of the {@code matrixB} does not match
   * the size of this instance or the system has no solution or has infinitely many solutions.
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix <b><i>B</i></b> of the equation to be solved, <b><i>AX = B</i></b>
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   *  as a new instance of Matrix of the same particular subtype as this instance.
   *
   *  @see #getErrorCode()
   */
  public abstract Matrix          solveAccurately(Matrix matrixB);

  /**
   * Solves a matrix equation of form <b><i>AX = B</i></b> and returns the found solution.
   * <br>
   * Solves a matrix equation of form <b><i>AX = B</i></b>
   * formed by the inner matrix data <b><i>A</i></b> and the matrix <b><i>B</i></b> passed in
   * as the {@code double[][] matrixB} argument, using LU decomposition.<br>
   * Uses an iterative refinement to find a more accurate solution.
   * In case of an inconsistent or underdefined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or the argument contains non-numeric values (NaN or Infinity),
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the size of the {@code matrixB} does not match the size of this instance
   *  or the system has no solution or has infinitely many solutions,
   *  or the argument contains non-numeric values (NaN or Infinity).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix <b><i>B</i></b> of the equation to be solved, <b><i>AX = B</i></b>,
   *  as a two-dimentional array of primitive {@code double}s
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b> as a new instance
   *  of Matrix of the same particular subtype as this instance.
   *
   *  @see #getErrorCode()
   */
  public abstract Matrix          solveAccurately(double[][] matrixB);

  /**
   * Solves a matrix equation of form <b><i>AX = B</i></b> and returns the found solution.
   * <br>
   * Solves a matrix equation of form <b><i>AX = B</i></b>
   * formed by the inner matrix data <b><i>A</i></b> and the matrix <b><i>B</i></b> passed in as the {@code Number[][] matrixB} argument,
   * using LU decomposition.<br>
   * Before solving, the given {@code matrixB} array gets translated to an array of values of the type that is used
   * to store internal matrix data (double for {@linkplain DoubleMatrix}, {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, etc).<br>
   * Uses an iterative refinement to find a more accurate solution.
   * In case of an inconsistent or underdefined matrix throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.<br>
   *
   * If the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or {@code matrixB} contains non-numeric values ({@code NaN} or {@code Infinity}) or {@code null},
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the size of the {@code matrixB} does not match the size of this instance
   *  or the system has no solution or has infinitely many solutions, or the argument contains invalid
   *  values (NaN, Infinity, or null).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix <b><i>B</i></b> of the equation to be solved, <b><i>AX = B</i></b>,
   *  as a two-dimentional array of {@linkplain Number}
   *
   * @return the found solution <b><i>X</i></b> to the equation <b><i>AX = B</i></b>
   *  as a new instance of Matrix of the same particular subtype as this instance.
   *
   *  @see #getErrorCode()
   */
  public abstract Matrix          solveAccurately(Number[][] matrixB);

  /**
   * Returns the last of the previously found matrix solutions, <b><i>X</i></b>, to a matrix equation
   * of form <b><i>AX = B</i></b>, as a {@code Matrix}.
   * If no matrix equation was solved with this instance of the {@code Matrix}, returns {@code null}.
   *
   * @return the last of the previously found matrix solution <b><i>X</i></b>, to a system of
   * linear equations of form <b><i>AX = B</i></b>, as a {@code Matrix},
   * or {@code null}, if no matrix equation was solved with the matrix.
   */
  public abstract Matrix          getMatrixSolution();

  /**
   * Returns the last of the previously found matrix solutions, <b><i>X</i></b>, to a matrix equation
   * of form <b><i>AX = B</i></b>, as a two-dimentional array of {@linkplain Number}.
   * <br>
   * The particular type of the elements of the returned vector depends on the particular
   * subtype of {@code Matrix} (it is Double for {@linkplain DoubleMatrix}, {@linkplain Quadruple}
   * for {@linkplain QuadrupleMatrix}, etc).<br>
   * If no matrix equation was solved with this instance of the {@code Matrix}, returns {@code null}.
   *
   * @return the last of the previously found matrix solution <b><i>X</i></b>, to a system of
   *  linear equations of form <b><i>AX = B</i></b>, as two-dimentional array of {@linkplain Number},
   *  or {@code null}, if no matrix equation was solved with the matrix.
   */
  public abstract Number[][]      getNumberMatrixSolution();

  /**
   * Returns the last of the previously found matrix solutions, <b><i>X</i></b>, to a matrix equation
   * of form <b><i>AX = B</i></b>, as a two-dimentional array of primitive {@code double}s.
   * <br>
   * If no matrix equation was solved with this instance of the {@code Matrix}, returns {@code null}.
   *
   * @return the last of the previously found matrix solution <b><i>X</i></b>, to a system of linear equations
   *  of form <b><i>AX = B</i></b>, as two-dimentional array of {@code double}s, or {@code null},
   *  if no matrix equation was solved with the matrix.
   */
  public abstract double[][]      getDoubleMatrixSolution();

  /**
   * Returns the last of the previously found matrix solutions, <b><i>X</i></b>, to a matrix equation
   * of form <b><i>AX = B</i></b>, as a two-dimentional array of {@linkplain Quadruple}s.
   * <br>
   * If no matrix equation was solved with this instance of the {@code Matrix}, returns {@code null}.
   *
   * @return the last of the previously found matrix solution <b><i>X</i></b>, to a system of linear equations
   *  of form <b><i>AX = B</i></b>, as two-dimentional array of {@linkplain Quadruple}s, or {@code null},
   *  if no matrix equation was solved with the matrix.
   */
  public abstract Quadruple[][]   getQuadrupleMatrixSolution();

  /**
   * Returns the last of the previously found matrix solutions, <b><i>X</i></b>, to a matrix equation
   * of form <b><i>AX = B</i></b>, as a two-dimentional array of {@linkplain BigDecimal}s.
   * <br>
   * If no matrix equation was solved with this instance of the {@code Matrix}, returns {@code null}.
   *
   * @return the last of the previously found matrix solution <b><i>X</i></b>, to a system of linear equations
   *  of form <b><i>AX = B</i></b>, as two-dimentional array of {@linkplain BigDecimal}, or {@code null},
   *  if no matrix equation was solved with the matrix.
   */
  public abstract BigDecimal[][]  getBigDecimalMatrixSolution();

  /* *******************************************************************************
  /***** Inversions ****************************************************************
  /*********************************************************************************/

  /**
   * Creates and returns a new Matrix instance containing the inversion of this instance.
   * <br>Computes the inversion of the matrix by solving the equation <b><i>AX = E</i></b>,
   * creates and returns a new {@code Matrix} containing the found inversion.
   * The exact subtype of the returned matrix is the same as that of this instance.
   * If the matrix is not invertible (i.e. inconsistent or underdefined),
   * throws {@code IllegalArgumentException}
   * with a relevant message and sets internal variable {@code errorCode},
   * whose value can be obtained using method {@linkplain #getErrorCode()}.
   *
   * @return a new {@code Matrix} containing the inversion of the given matrix.
   *
   * @throws IllegalArgumentException if the matrix is not invertible
   *  (i.e. inconsistent or underdefined)
   */
  public abstract Matrix          inverse();

  /**
   * Creates and returns a new Matrix instance containing the inversion of this instance.
   * <br>Computes the inversion of the matrix by solving the equation <b><i>AX = E</i></b>,
   * creates and returns a new {@code Matrix} containing the found inversion.
   * Uses an iterative refinement to achieve a more accurate solution to the equation.
   * The exact subtype of the returned matrix is the same as that of this instance.
   * If the matrix is not invertible (i.e. inconsistent or underdefined),
   * throws {@code IllegalArgumentException} with a relevant message and sets
   * internal variable {@code errorCode}, whose value can be obtained using method {@linkplain #getErrorCode()}.
   *
   * @return a new {@code Matrix} containing the inversion of the given matrix.
   *
   * @throws IllegalArgumentException if the matrix is not invertible
   *  (i.e. inconsistent or underdefined)
   */
  public abstract Matrix          inverseAccurately();


  /**
   * Creates and returns a new Matrix instance containing the transposition of this instance.
   * <br>Computes the transposition of the matrix, creates and returns a new {@code Matrix}
   * containing the found transposition.
   * <br>The exact subtype of the returned matrix is the same as that of this instance.
   *
   * @return a new {@code Matrix} containing the transposition of the given matrix.
   */
  public abstract Matrix          transpose();

  /**
   * Creates and returns a new Matrix instance containing a unity matrix of the same size as the source matrix.
   * <br>The exact subtype of the returned matrix is the same as that of this instance.
   *
   * @return a new {@code Matrix} containing the unity matrix of the same size as this {@code Matrix}
   */
  public abstract Matrix          unity();

  /* *******************************************************************************
  /* **** Multiplications ***********************************************************
  /*********************************************************************************/

  /**
   * Multiplies this instance of {@code Matrix} by the matrix passed in as the {@code factor} argument,
   * creates and returns a new matrix containing the product.
   * The exact subtype of the returned matrix is the same as that of this instance.
   * If the size of the {@code factor} does not match the size of the source matrix,
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @param factor a matrix to multiply this instance by
   *
   * @throws IllegalArgumentException if the size of the argument does not match the size of
   * this instance
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @return a new {@code Matrix} representing the product
   */
  public abstract Matrix          multiply(Matrix factor);

  /**
   * Multiplies this instance of {@code Matrix} by the {@code factor} passed in
   * as a two-dimentional double array, and creates and returns a new matrix containing the product.
   * The exact subtype of the returned matrix is the same as that of this instance.<br>
   *
   * If the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or {@code matrixB} contains non-numeric values ({@code NaN} or {@code Infinity}),
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or {@code matrixB} contains non-numeric values ({@code NaN} or {@code Infinity})
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param factor two-dimentional double array representing a matrix to multiply this instance by
   *
   * @return a new {@code Matrix} representing the product
   */
  public abstract Matrix          multiply(double[][] factor);

  /**
   * Multiplies this instance of {@code Matrix} by the {@code factor} passed in
   * as a two-dimentional array of {@link Number}s, and creates and returns a new matrix containing the product.
   * The exact subtype of the returned matrix is the same as that of this instance.<br>
   *
   * If the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or {@code matrixB} contains non-numeric values ({@code NaN} or {@code Infinity}) or {@code null},
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or {@code matrixB} contains non-numeric values ({@code NaN} or {@code Infinity}) or {@code null}
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param factor two-dimentional array of {@linkplain Number} representing the matrix to multiply this instance by
   * @return a new {@code Matrix} representing the product
   */
  public abstract Matrix          multiply(Number[][] factor);

  /**
   * Multiplies this instance of {@code Matrix} by a vector
   * passed in as an array of {@code double}s, and returns an array of {@linkplain Number} values
   * containing the product.
   * The exact type of the elements of the returned array depends
   * on the particular subtype of the {@code Matrix}.<br>
   *
   * If the length of the {@code vector} does not match the size of the source matrix,
   * or the argument contains non-numeric values (NaN or Infinity),
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the length of the array
   *  does not match the size of the matrix,
   *  or the argument contains non-numeric values (NaN or Infinity).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @return an array of {@linkplain Number} values containing the product.
   * @param vector the vector to multiply the matrix by.
   */
  public abstract Number[]        multiply(double[]   vector);

  /**
   * Multiplies this instance of {@code Matrix} by a vector
   * passed in as an array of {@linkplain Number} values and returns an array of {@linkplain Number}
   * values containing the product.
   * The exact type of the elements of the returned array depends
   * on the particular subtype of the {@code Matrix}.<br>
   *
   * If the length of the {@code vector} does not match the size of the source matrix,
   * or the argument contains non-numeric values (NaN or Infinity) or {@code null,}
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the length of the array
   *  does not match the size of the matrix,
   *  or the argument contains non-numeric values (NaN or Infinity) or {@code null}
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @return an array of {@linkplain Number} values containing the product.
   * @param vector the vector to multiply the matrix by.
   */
  public abstract Number[]        multiply(Number[]   vector);

  /**
   * Multiplies this instance of {@code Matrix} by a scalar factor {@code scalar}
   * passed in as a {@code double} parameter, and creates and returns a new matrix containing the product.
   * The exact subtype of the returned matrix is the same as that of this instance.
   *
   * If the argument is a non-numeric value ({@code NaN} or {@code Infinity}),
   * throws {@linkplain IllegalArgumentException}.
   *
   * @throws IllegalArgumentException if the argument is {@code NaN} or {@code Infinity}.
   *
   * @param scalar a {@code double} value to multiply this matrix by.
   * @return a new {@code Matrix} containing the product of the source matrix and the given scalar.
   */
  public abstract Matrix          multiply(double     scalar);

  /**
   * Multiplies this instance of {@code Matrix} by a scalar factor {@code scalar}
   * passed in as a {@linkplain Number} parameter, and creates and returns a new matrix containing the product.
   * The exact subtype of the returned matrix is the same as that of this instance.<br>
   *
   * If the argument is a non-numeric value ({@code NaN} or {@code Infinity}) or null,
   * throws {@linkplain IllegalArgumentException}.
   *
   * @throws IllegalArgumentException if the argument is {@code NaN}, {@code Infinity} or null.
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param scalar a {@linkplain Number} value to multiply this matrix by.
   * @return a new {@code Matrix} containing the product of the source matrix and the given scalar.
   */
  public abstract Matrix          multiply(Number     scalar);

  /* *******************************************************************************
  /* **** Additions and subtractions ***********************************************
  /*********************************************************************************/

  /**
   * Adds the given {@code matrixB} to this matrix and returns the sum.
   * <br>Computes the sum of this matrix and the {@code Matrix} passed in,
   * and creates and returns a new {@code Matrix} containing the found sum.
   * The exact subtype of the returned matrix is the same as that of this instance.
   *
   * If the size of the {@code matrixB} does not match the size of this instance of the matrix,
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the size of the {@code matrixB} does not match
   *  the size of this instance.
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix to add to the source matrix.
   * @return a new {@code Matrix} containing the sum of the two matrices.
   */
  public abstract Matrix          add(Matrix matrixB);

  /**
   * Adds the given {@code matrixB} passed in as a two-dimentional array of {@code double}s
   * to this matrix and returns the sum.
   * <br>Computes the sum of this matrix and the matrix passed in,
   * and creates and returns a new {@code Matrix} containing the found sum.
   * The exact subtype of the returned matrix is the same as that of this instance.<br>
   *
   * If the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or the argument contains non-numeric values (NaN or Infinity),
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the passed-in array is non-square,
   *  or its size does not match the size of this instance,
   *  or the argument contains non-numeric values ({@code NaN} or {@code Infinity}).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix to add to the source matrix.
   * @return a new {@code Matrix} containing the sum of the two matrices.
   */
  public abstract Matrix          add(double[][] matrixB);

  /**
   * Adds the given {@code matrixB} passed in as a two-dimentional array of {@linkplain Number}
   * to this matrix and returns the sum.
   * <br>Computes the sum of this matrix and the matrix passed in,
   * and creates and returns a new {@code Matrix} containing the found sum.
   * The exact subtype of the returned matrix is the same as that of this instance.<br>
   *
   * If the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or the argument contains non-numeric values (NaN or Infinity) or {@code null},
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the passed-in array is non-square,
   *  or its size does not match the size of this instance,
   *  or the argument contains non-numeric values ({@code NaN} or {@code Infinity}), or {@code null}.
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix to add to the source matrix.
   * @return a new {@code Matrix} containing the sum of the two matrices.
   */
  public abstract Matrix          add(Number[][] matrixB);

  /**
   * Subtracts the given {@code matrixB} from this matrix and returns the difference.
   * <br>Computes the difference of this matrix and the {@code Matrix} passed in,
   * and creates and returns a new {@code Matrix} containing the found difference.
   * The exact subtype of the returned matrix is the same as that of this instance.
   *
   * If the size of the {@code matrixB} does not match the size of this instance of the matrix,
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the size of the {@code matrixB} does not match
   *  the size of this instance.
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix to subtract from the source matrix.
   * @return a new {@code Matrix} containing the difference of the two matrices.
   */
  public abstract Matrix          subtract(Matrix matrixB);

  /**
   * Subtracts the given {@code matrixB} passed in as a two-dimentional array of {@code double}s
   * from this matrix and returns the difference.
   * <br>Computes the difference of this matrix and the matrix passed in,
   * creates and returns a new {@code Matrix} containing the found difference.
   * The exact subtype of the returned matrix is the same as that of this instance.<br>
   *
   * If the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or the argument contains non-numeric values (NaN or Infinity),
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the passed-in array is non-square,
   *  or its size does not match the size of this instance,
   *  or the argument contains non-numeric values ({@code NaN} or {@code Infinity}).
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix to subtract from the source matrix.
   * @return a new {@code Matrix} containing the difference of the two matrices.
   */
  public abstract Matrix          subtract(double[][] matrixB);

  /**
   * Subtracts the given {@code matrixB} passed in as a two-dimentional array of {@linkplain Number}
   * from this matrix and returns the difference.
   * <br>Computes the difference of this matrix and the matrix passed in,
   * creates and returns a new {@code Matrix} containing the found difference.
   * The exact subtype of the returned matrix is the same as that of this instance.<br>
   *
   * If the passed-in array is non-square,
   * or its size does not match the size of this instance of the matrix,
   * or the argument contains non-numeric values (NaN or Infinity) or {@code null},
   * throws {@code IllegalArgumentException} with a relevant message.
   *
   * @throws IllegalArgumentException if the passed-in array is non-square,
   *  or its size does not match the size of this instance,
   *  or the argument contains non-numeric values ({@code NaN} or {@code Infinity}) or {@code null}.
   *
   * @throws NullPointerException if the argument is {@code null}
   *
   * @param matrixB the matrix to subtract from the source matrix.
   * @return a new {@code Matrix} containing the difference of the two matrices.
   */
  public abstract Matrix          subtract(Number[][] matrixB);

  /* *******************************************************************************
  /* **** Determinant **************************************************************
  /*********************************************************************************/

  /**
   * Computes the determinant of the matrix and returns its value as a {@linkplain Number} value.
   * <br>The particular type of the returned {@linkplain Number} depends on the exact
   * subtype of this instance of {@code Matrix} (it is {@code double} for {@linkplain DoubleMatrix},
   * {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, and {@linkplain BigDecimal}
   * for {@linkplain BigDecimalMatrix}).
   * @return the value of the determinant of the given {@code Matrix}
   */
  public abstract Number          determinant();

  /**
   * Computes the determinant of the matrix and returns its value as a {@code double} value.
   * @return the value of the determinant of the given {@code Matrix}
   */
  public abstract double          determinantAsDouble();

  /**
   * Computes the determinant of the matrix and returns its value as a {@linkplain Quadruple} value.
   * The precision of the result depends on the particular subtype of the instance.
   * @return the value of the determinant of the given {@code Matrix}
   */
  public abstract Quadruple       determinantAsQuadruple();

  /**
   * Computes the determinant of the matrix and returns its value as a {@linkplain BigDecimal} value.
   * The precision of the result depends on the particular subtype of the instance.
   * @return the value of the determinant of the given {@code Matrix}
   */
  public abstract BigDecimal      determinantAsBigDecimal();

  /**
   * Computes the row-based norm of the matrix, ║<b>A</b>║<sub>∞</sub>, and returns
   * its value as a {@linkplain Number} value.
   * <br>The particular type of the returned {@linkplain Number} depends on the exact
   * subtype of this instance of {@code Matrix} (it is {@code double} for {@linkplain DoubleMatrix},
   * {@linkplain Quadruple} for {@linkplain QuadrupleMatrix}, and {@linkplain BigDecimal}
   * for {@linkplain BigDecimalMatrix}).
   * @return the value of the norm of the given {@code Matrix}
   */
  public abstract Number          norm();

  /**
   * Computes the row-based norm of the matrix, ║<b>A</b>║<sub>∞</sub>, and returns
   * its value as a {@code double} value.
   * @return the value of the norm of the given {@code Matrix}
   */
  public abstract double          normAsDouble();

  /**
   * Computes the row-based norm of the matrix, ║<b>A</b>║<sub>∞</sub>, and returns
   * its value as a {@linkplain Quadruple} value.
   * The precision of the result depends on the particular subtype of the instance.
   * @return the value of the norm of the given {@code Matrix}
   */
  public abstract Quadruple       normAsQuadruple();

  /**
   * Computes the row-based norm of the matrix, ║<b>A</b>║<sub>∞</sub>, and returns its
   * value as a {@linkplain BigDecimal} value.
   * The precision of the result depends on the particular subtype of the instance.
   * @return the value of the norm of the given {@code Matrix}
   */
  public abstract BigDecimal      normAsBigDecimal();

  /**
   * Computes the condition number of the matrix, ║<b>A</b>║•║<b>A<sup>-1</sup></b>║, and returns
   * its value as a {@code double} value.
   * <br>For non-invertible matrices returns {@code Double.POSITIVE_INFINITY}.
   * @return the value of the condition number of the given {@code Matrix},
   * or {@code Double.POSITIVE_INFINITY} for non-invertible matrices
   */
  public abstract double          cond();

  /*#######################################################################
  /*### Protected methods used by descendant classes ######################
  /*#######################################################################*/

  /**
   * Check the argument of a constructor. Throws IllegalArgumentException if it's invalid
   * @param data -- the argument to check
   */
  protected final static void checkMatrixData(Matrix data) {
    checkThatArgumentIsNotNull(data, "Can't create matrix: data");
  }

  protected final static void checkMatrixData(double[][] data) {
    checkThatArgumentIsNotNull(data, "Can't create matrix: data");
    if (data.length == 0)
      throwMatrixDataIsEmpty();
    for (int i = 0; i < data.length; i++) {
      if (data == null || data[i].length != data.length)
        throwMatrixDataIsNonSquare();
      for (int j = 0; j < data.length; j++) {
        if (!Double.isFinite(data[i][j]))
          throwNaNOrInfinityInMatrixData();
      }
    }
  }

  protected final static void checkMatrixData(Number[][] data) {
    checkThatArgumentIsNotNull(data, "Can't create matrix: data");
    if (data.length == 0)
      throwMatrixDataIsEmpty();
    for (int i = 0; i < data.length; i++) {
      if (data == null || data[i].length != data.length)
        throwMatrixDataIsNonSquare();
      for (int j = 0; j < data.length; j++) {
        if (data[i][j] == null)
          throwNullInMatrixData();
        if (data[i][j] instanceof Quadruple) {
          final Quadruple q = (Quadruple)data[i][j];
          if (q.isInfinite() || q.isNaN())
            throwNaNOrInfinityInMatrixData();
        }
      }
    }
  }

  protected final static void checkScalar(double value, String callerName) {
    checkThatArgumentIsANumber(value, callerName);
  }

  protected final static void checkScalar(Number value, String callerName) {
    checkThatArgumentIsNotNull(value, "The value passed to " + callerName);
    if (value instanceof Quadruple || value instanceof Double) {
      checkThatArgumentIsANumber(value.doubleValue(), callerName);
    }
  }

  protected final void checkVector(double[] vector, String callerName) {
    checkThatVectorIsNotNull(vector, callerName);
    if (vector.length != size) {
      trowVectorOfWrongLength(callerName);
    }
    for (int i = 0; i < size; i++) {
      if (!Double.isFinite(vector[i])) {
        throwNaNOrInfinityInVectot(callerName);
      }
    }
  }

  protected final void checkVector(Number[] vector, String callerName) {
    checkThatVectorIsNotNull(vector, callerName);
    if (vector.length != size) {
      trowVectorOfWrongLength(callerName);
    }
    for (int i = 0; i < size; i++) {
      if (vector[i] == null)
        throwNullInVector(callerName);
      if (vector[i] instanceof Quadruple) {
        final Quadruple q = (Quadruple)vector[i];
        if (q.isInfinite() || q.isNaN())
          throwNaNOrInfinityInVectot(callerName);
      }
    }
  }

  protected final void checkMatrix(Matrix matrix, String callerName) {
    checkThatMatrixIsNotNull(matrix, callerName);
    if (size != matrix.size) {
      throwMatrixOfWrongSize(callerName);
    }
  }

  protected final void checkMatrix(double[][] matrix, String callerName) {
    checkThatMatrixIsNotNull(matrix, callerName);
    if (size != matrix.length) {
      throwMatrixOfWrongSize(callerName);
    }
    for (int i = 0; i < size; i++) {
      if (matrix[i] == null || matrix[i].length != this.size)
        throwNonSquareArray(callerName);
      for (int j = 0; j < size; j++) {
        if (!Double.isFinite(matrix[i][j]))
          throwNaNOrInfinityInMatrix(callerName);
      }
    }
  }

  protected final void checkMatrix(Number[][] matrix, String callerName) {
    checkThatMatrixIsNotNull(matrix, callerName);
    if (size != matrix.length) {
      throwMatrixOfWrongSize(callerName);
    }
    for (int i = 0; i < size; i++) {
      if (matrix[i] == null || matrix[i].length != this.size)
        throwNonSquareArray(callerName);
      for (int j = 0; j < size; j++) {
        if (matrix[i][j] == null)
          throwNullInMatrix(callerName);
        if (matrix[i][j] instanceof Quadruple) {
          final Quadruple q = (Quadruple)matrix[i][j];
          if (q.isInfinite() || q.isNaN())
            throwNaNOrInfinityInMatrix(callerName);
        }
      }
    }
  }

  /*###############################################################################
   * Private members
   ###############################################################################*/

  private static void checkThatVectorIsNotNull(Object vector, String callerName) {
    checkThatArgumentIsNotNull(vector, "The vector passed to " + callerName);
  }

  private static void checkThatMatrixIsNotNull(Object matrix, String callerName) {
    checkThatArgumentIsNotNull(matrix, "The matrix passed to " + callerName);
  }

  private static void checkThatArgumentIsNotNull(Object data, String message) {
    Objects.requireNonNull(data, message + " is null");
  }

  private static void checkThatArgumentIsANumber(double value, String callerName) {
    if (!Double.isFinite(value))
      throw new IllegalArgumentException("The value passed to " + callerName + " must not be NaN or Infinity");
  }

  private static void throwMatrixDataIsNonSquare() {
    throw new IllegalArgumentException("Can't create matrix: argument must be square");
  }

  private static void throwMatrixDataIsEmpty() {
    throw new IllegalArgumentException("Can't create matrix: argument is empty");
  }

  private static void throwNullInMatrixData() {
    throw new IllegalArgumentException("Can't create matrix: data must not contain null");
  }

  private static void throwNaNOrInfinityInMatrixData() {
    throw new IllegalArgumentException("Can't create matrix: data must not contain NaN or Infinity");
  }

  private static void trowVectorOfWrongLength(String callerName) {
    throw new IllegalArgumentException("The vector passed to " + callerName + " must have the same size as the matrix");
  }

  private static void throwNaNOrInfinityInVectot(String callerName) {
    throw new IllegalArgumentException("The vector passed to " + callerName + " must not contain NaN or Infinity");
  }

  private static void throwNullInVector(String callerName) {
    throw new IllegalArgumentException("The vector passed to " + callerName + " must not contain null");
  }

  private static void throwMatrixOfWrongSize(String callerName) {
    throw new IllegalArgumentException("The matrix passed to " + callerName + " must have the same size as the callee");
  }

  private static void throwNaNOrInfinityInMatrix(String callerName) {
    throw new IllegalArgumentException("The matrix passed to " + callerName + " must not contain NaN or Infinity");
  }

  private static void throwNullInMatrix(String callerName) {
    throw new IllegalArgumentException("The matrix passed to " + callerName + " must not contain null");
  }

  private static void throwNonSquareArray(String callerName) {
    throw new IllegalArgumentException("The array passed to " + callerName + " must be square");
  }

}
