/**
 * A simple library for matrix calculations on square matrices of real numbers, with different levels of precision.
 * Includes <ul>
 * <li><b>{@linkplain com.mvohm.quadmatrix.Matrix}</b> -- a generic abstract class which defines
 *  the set of operations that all its subclasses must implement;
 * <li><b>{@linkplain com.mvohm.quadmatrix.DoubleMatrix}</b> -- a subclass of the Matrix that
 *  stores the matrix data and implements operations on them using primitive {@code double} type;
 * <li><b>{@linkplain com.mvohm.quadmatrix.BigDecimalMatrix}</b> -- a subclass of the Matrix
 *  that stores the matrix data and implements operations on them using {@link java.math.BigDecimal} class,
 *  which allows to gain (theoretically) unlimited precision;
 * <li><b>{@linkplain com.mvohm.quadmatrix.QuadrupleMatrix}</b> -- a subclass of the Matrix
 *  that stores the matrix data and implements operations on them using {@link com.mvohm.quadruple.Quadruple}
 *  class (see <a href=https://github.com/m-vokhm/Quadruple>Github repository</a>),
 *  which allows to gain an accuracy of about 36 decimal places (mat vary depending on the matrix size and the nature
 *  of performed operation) and performs calculations much faster than <code>BigDecimalMatrix</code>.
 * </ul>
 * The current version is limited to square matrices and does not provide any optimization for sparse matrices.
 */

package com.mvohm.quadmatrix;

