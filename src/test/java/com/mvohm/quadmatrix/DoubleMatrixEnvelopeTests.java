/*

 Copyright 2021 M.Vokhmentsev

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

import java.util.Random;

import java.lang.IllegalArgumentException;
import java.math.BigDecimal;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import com.mvohm.quadruple.Quadruple;

import static org.junit.jupiter.api.Assertions.*;
import static com.mvohm.quadmatrix.test.AuxMethods.*;
import static org.assertj.core.api.Assertions.*;

/**
 * Tests for {@link DoubleMatrix} constructors and matrix data getters
 * <br>
 * The methods to be tested:
 *   public static void setDefaultScaling(boolean scaleByDefault)
 *   public static boolean getDefaultScaling()
 *   public DoubleMatrix(Matrix source)
 *   public DoubleMatrix(Matrix source, boolean needToScale)
 *   public DoubleMatrix(double[][] source)
 *   public DoubleMatrix(double[][] source, boolean needToScale)
 *   public DoubleMatrix(Number[][] source)
 *   public DoubleMatrix(Number[][] source, boolean needToScale)
 *   public boolean getScaling()
 *   public int getSize()
 *   public Number[][] getData()
 *   public double[][] getDoubleData()
 *   public Quadruple[][] getQuadrupleData()
 *   public BigDecimal[][] getBigDecimalData()
 *   public abstract boolean equals(Object anotherOne);
 *   public abstract int hashCode();
 */
public class DoubleMatrixEnvelopeTests {

  /* *******************************************************************************
  /***** Testing methods that don't involve DoubleMatrixSolver ********************
  /***** Just data manipulations ***************************************************
  /*********************************************************************************/

  private static final int SOME_SIZE                            = 1234;
  private static final int ANY_NUMBER_BUT_ZERO                  = 123;
  private static final int SAMPLE_MATRIX_SIZE                   = 10;
  private static final Random RANDOM                            = new Random(12345);
  private static final double[][] SAMPLE_MATRIX_DATA            = new double[SAMPLE_MATRIX_SIZE][SAMPLE_MATRIX_SIZE];
  private static final double[][] SAMPLE_TOO_LARGE_MATRIX_DATA  = new double[SAMPLE_MATRIX_SIZE + 1][SAMPLE_MATRIX_SIZE + 1];
  private static final double[][] NONSQUARE_ARRAY               = new double[SAMPLE_MATRIX_SIZE][SAMPLE_MATRIX_SIZE + 1];
  private static final Quadruple[][] SAMPLE_QUADRUPLE_DATA      = new Quadruple[SAMPLE_MATRIX_SIZE][SAMPLE_MATRIX_SIZE];
  private static final BigDecimal[][] SAMPLE_BIGDECIMAL_DATA    = new BigDecimal[SAMPLE_MATRIX_SIZE][SAMPLE_MATRIX_SIZE];
  private static final DoubleMatrix SAMPLE_DOUBLE_MATRIX;
  private static final double[][] SAMPLE_ANOTHER_MATRIX_DATA;

  static {
    for (int i = 0; i < SAMPLE_MATRIX_SIZE; i++) {
      for (int j = 0; j < SAMPLE_MATRIX_SIZE; j++) {
        SAMPLE_MATRIX_DATA[i][j] = RANDOM.nextDouble();
        NONSQUARE_ARRAY[i][j] = SAMPLE_MATRIX_DATA[i][j];
        SAMPLE_QUADRUPLE_DATA[i][j] = new Quadruple(SAMPLE_MATRIX_DATA[i][j]);
        SAMPLE_BIGDECIMAL_DATA[i][j] = BigDecimal.valueOf(SAMPLE_MATRIX_DATA[i][j]);
        NONSQUARE_ARRAY[i][j] = SAMPLE_MATRIX_DATA[i][j];
      }
      NONSQUARE_ARRAY[i][SAMPLE_MATRIX_SIZE] = 1.0;
    }
    SAMPLE_DOUBLE_MATRIX = new DoubleMatrix(SAMPLE_MATRIX_DATA);
    SAMPLE_ANOTHER_MATRIX_DATA = deepCopyOf(SAMPLE_MATRIX_DATA);
    SAMPLE_ANOTHER_MATRIX_DATA[1][1] += 0.123;
  }

  /*################################################################################
   * Test methods for
   *   public static void setDefaultScaling(boolean scaleByDefault)
   *   public static boolean getDefaultScaling()
   * Tested behavior:
   *  -- getDefaultScaling() returns true after true was set with setDefaultScaling()
   *  -- getDefaultScaling() returns false after false was set with setDefaultScaling()
   #################################################################################*/

  @Test
  @DisplayName("getDefaultScaling() returns true after true was set with setDefaultScaling()")
  void testGetDefaultScalingReturnsTrueAfterSettingTrue() {
    final boolean expectedValue = true;

    DoubleMatrix.setDefaultScaling(!expectedValue);
    DoubleMatrix.setDefaultScaling(expectedValue);
    final boolean actualValue = DoubleMatrix.getDefaultScaling();

    assertThat(actualValue).
        withFailMessage("getDefaultScaling() returns false when true was set with setDefaultScaling()").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("getDefaultScaling() returns false after false was set with setDefaultScaling()")
  void testGetDefaultScalingReturnsFalseAfterSettingFalse() {
    final boolean expectedValue = false;

    DoubleMatrix.setDefaultScaling(!expectedValue);
    DoubleMatrix.setDefaultScaling(expectedValue);
    final boolean actualValue = DoubleMatrix.getDefaultScaling();

    assertThat(actualValue).
        withFailMessage("getDefaultScaling() returns true when false was set with setDefaultScaling()").
        isEqualTo(expectedValue);
  }

  /*################################################################################
   * Test methods for
   *   public DoubleMatrix(Matrix source)
   * Tested behavior:
   *  -- DoubleMatrix(Matrix) throws NullPointerException with a relevant message when called with null
   *  -- DoubleMatrix(Matrix) creates a matrix whose size equal the source's size
   *  -- DoubleMatrix(Matrix) sets scaling to true, when default scaling is set to true
   *  -- DoubleMatrix(Matrix) sets scaling to false, when default scaling is set to false
   *  -- DoubleMatrix(Matrix, true) creates a matrix whose neetToScale flag is set to true
   *  -- DoubleMatrix(Matrix, false) creates a matrix whose neetToScale flag is set to false
   *  -- DoubleMatrix(Matrix) creates a matrix whose data equal the source's data
   #################################################################################*/

  @Test
  @DisplayName("DoubleMatrix(Matrix) throws NullPointerException when called with null")
  void testConstructorWithMatrixThrowsExceptionGivenNull() {

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> new DoubleMatrix((Matrix)null),
        "DoubleMatrix(Matrix) should throw NullPointerException when called with null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is null'").
        contains("can't create matrix").contains("is null");
  }

  @Test
  @DisplayName("DoubleMatrix(Matrix) creates a matrix whose size equal the source's size")
  void testConstructorWithMatrixMakesMatrixOfEqualSize() {
    final int expectedSize = SAMPLE_DOUBLE_MATRIX.getSize();

    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_DOUBLE_MATRIX);
    final int actualSize = newMatrix.getSize();

    assertThat(actualSize).
        withFailMessage("Newly created matrix's size (%d) is not equal to the size of the source matrix (%d)",
            actualSize, expectedSize).
        isEqualTo(expectedSize);
  }

  @Test
  @DisplayName("DoubleMatrix(Matrix) sets scaling to true, when default scaling is set to true")
  void testConstructorWithMatrixProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    DoubleMatrix.setDefaultScaling(expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_DOUBLE_MATRIX);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the default scaling is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(Matrix) sets scaling to false, when default scaling is set to false")
  void testConstructorWithMatrixProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    DoubleMatrix.setDefaultScaling(expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_DOUBLE_MATRIX);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the default scaling is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(Matrix, true) creates a matrix whose neetToScale flag is set to true ")
  void testConstructorWithMatrixAndNeedToScaleProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    DoubleMatrix.setDefaultScaling(!expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_DOUBLE_MATRIX, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the needToScale argument is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(Matrix, false) creates a matrix whose neetToScale flag is set to false ")
  void testConstructorWithMatrixAndNeedToScaleProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    DoubleMatrix.setDefaultScaling(!expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_DOUBLE_MATRIX, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the needToScale argument is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(Matrix) creates a matrix whose data equal the source's data")
  void testConstructorWithMatrixMakesMatrixWithEqualData() {
    final double[][] expectedData = SAMPLE_DOUBLE_MATRIX.getDoubleData();

    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_DOUBLE_MATRIX);
    final double[][] actualData = newMatrix.getDoubleData();

    assertThat(actualData).
        withFailMessage("Newly created matrix contains data not equal to the data of the source matrix").
        isEqualTo(expectedData);
  }

  /*################################################################################
   * public DoubleMatrix(Matrix source, boolean needToScale)
   * Doesn't actually require special tests since it's used by the already tested constructor
   #################################################################################*/

  /*################################################################################
   * Test methods for
   *   public DoubleMatrix(double[][])
   * Tested behavior:
   *  -- DoubleMatrix(double[][]) throws NullPointerException when called with null
   *  -- DoubleMatrix(double[][]) throws IllegalArgumentException when given an empty array
   *  -- DoubleMatrix(double[][]) throws IllegalArgumentException when given a non-square array
   *  -- DoubleMatrix(double[][]) throws IllegalArgumentException if the argument contains NaN
   *  -- DoubleMatrix(double[][]) throws IllegalArgumentException if the argument contains infinity
   *  -- DoubleMatrix(double[][]) creates a matrix whose size equal the source's size
   *  -- DoubleMatrix(double[][]) sets scaling to true, when default scaling is set to true
   *  -- DoubleMatrix(double[][]) sets scaling to false, when default scaling is set to false
   *  -- DoubleMatrix(double[][], true) creates a matrix whose neetToScale flag is set to true
   *  -- DoubleMatrix(double[][], false) creates a matrix whose neetToScale flag is set to false
   *  -- DoubleMatrix(double[][]) creates a matrix with data equal to the source data
   #################################################################################*/

  @Test
  @DisplayName("DoubleMatrix(double[][]) throws NullPointerException when called with null")
  void testConstructorWithDoubleArrayThrowsExceptionGivenNull() {

    final Throwable expectedException = assertThrows( NullPointerException.class,
        () -> new DoubleMatrix((double[][])null),
        "Constructor DoubleMatrix(double[][]) should throw NullPointerException when called with null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is null'").
        contains("can't create matrix").contains("is null");
  }

  @Test
  @DisplayName("DoubleMatrix(double[][]) throws IllegalArgumentException when given an empty array")
  void testConstructorWithDoubleArrayThrowsExceptionGivenEmptyArray() {
    final double[][] sourceData = new double[0][];

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(sourceData),
        "Constructor DoubleMatrix(double[][]) should throw IllegalArgumentException when given an empty array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is empty'").
        contains("can't create matrix").contains("is empty");
  }

  @Test
  @DisplayName("DoubleMatrix(double[][]) throws IllegalArgumentException when given a non-square array")
  void testConstructorWithDoubleArrayThrowsExceptionGivenNonSquareArray() {
    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(NONSQUARE_ARRAY),
        "Constructor DoubleMatrix(double[][]) should throw IllegalArgumentException when given a non-square array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'must be square'").
        contains("can't create matrix").contains("must be square");
  }

  @Test
  @DisplayName("DoubleMatrix(double[][]) throws IllegalArgumentException if the argument contains NaN")
  void testConstructorWithDoubleArrayWithNaNThrowsException() {
    final double[][] badMatrixData = deepCopyOf(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = Double.NaN;

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(badMatrixData),
        "Constructor DoubleMatrix(double[][]) should throw IllegalArgumentException if the argument contains NaN");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain NaN or Infinity'").
        contains("can't create matrix").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("DoubleMatrix(double[][]) throws IllegalArgumentException if the argument contains infinity")
  void testConstructorWithDoubleArrayWithInfinityThrowsException() {
    final double[][] badMatrixData = deepCopyOf(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = Double.POSITIVE_INFINITY;

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(badMatrixData),
        "Constructor DoubleMatrix(double[][]) should throw IllegalArgumentException if the argument contains infinity");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain NaN or Infinity'").
        contains("can't create matrix").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("DoubleMatrix(double[][]) creates a matrix whose size equal the source's size")
  void testConstructorWithDoubleArrayMakesMatrixOfEqualSize() {
    final int expectedSize = SAMPLE_MATRIX_DATA.length;

    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_MATRIX_DATA);
    final int actualSize = newMatrix.getSize();

    assertThat(actualSize).
        withFailMessage("Newly created matrix's size (%d) is not equal to the size of the source matrix (%d)",
            actualSize, expectedSize).
        isEqualTo(expectedSize);
  }

  @Test
  @DisplayName("DoubleMatrix(double[][]) sets scaling to true, when default scaling is set to true")
  void testConstructorWithDoubleArrayProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    DoubleMatrix.setDefaultScaling(expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_MATRIX_DATA);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the default scaling is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(double[][]) sets scaling to false, when default scaling is set to false")
  void testConstructorWithDoubleArrayProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    DoubleMatrix.setDefaultScaling(expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_MATRIX_DATA);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the default scaling is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(double[][], true) creates a matrix whose neetToScale flag is set to true ")
  void testConstructorWithDoubleArrayAndNeedToScaleProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    DoubleMatrix.setDefaultScaling(!expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_MATRIX_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the needToScale argument is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(double[][], false) creates a matrix whose neetToScale flag is set to false ")
  void testConstructorWithDoubleArrayAndNeedToScaleProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    DoubleMatrix.setDefaultScaling(!expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_MATRIX_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the needToScale argument is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(double[][]) creates a matrix with data equal to the source data")
  void testConstructorWithDoubleArrayCreatesInstanceWithCorrectData() {
    final double[][] expectedData = SAMPLE_MATRIX_DATA;

    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_MATRIX_DATA);
    final double[][] actualData = newMatrix.getDoubleData();

    assertThat(actualData).
      withFailMessage("Constructor DoubleMatrix(double[][]) creates a matrix containing data not equal to the source data").
      isEqualTo(expectedData);
  }

  /*################################################################################
   * public DoubleMatrix(double[][] source, boolean needToScale)
   * Doesn't actually require special tests since it's used by the already tested constructor
   #################################################################################*/

  /*################################################################################
   * Test methods for
   *   public DoubleMatrix(Number[][]) with an array of Quadruples as the parameter
   * Tested behavior:
   *  -- DoubleMatrix(Quadruple[][]) throws NullPointerException when called with null
   *  -- DoubleMatrix(Quadruple[][]) throws IllegalArgumentException when given an empty array
   *  -- DoubleMatrix(Quadruple[][]) throws IllegalArgumentException when given a non-square array
   *  -- DoubleMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains NaN
   *  -- DoubleMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains infinity
   *  -- DoubleMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains null
   *  -- DoubleMatrix(Quadruple[][]) creates a matrix whose size equal the source's size
   *  -- DoubleMatrix(Quadruple[][]) sets scaling to true, when default scaling is set to true
   *  -- DoubleMatrix(Quadruple[][]) sets scaling to false, when default scaling is set to false
   *  -- DoubleMatrix(Quadruple[][], true) creates a matrix whose neetToScale flag is set to true
   *  -- DoubleMatrix(Quadruple[][], false) creates a matrix whose neetToScale flag is set to false
   *  -- DoubleMatrix(Quadruple[][]) creates a matrix with data equal to the source data
   #################################################################################*/

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][]) throws NullPointerException when called with null")
  void testConstructorWithQuadrupleArrayThrowsExceptionGivenNull() {
    final Throwable expectedException = assertThrows( NullPointerException.class,
        () -> new DoubleMatrix((Quadruple[][])null),
        "Constructor DoubleMatrix(Quadruple[][]]) should throws NullPointerException when called with null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is null'").
        contains("can't create matrix").contains("is null");
  }

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][]) throws IllegalArgumentException when given an empty array")
  void testConstructorWithQuadrupleArrayThrowsExceptionGivenEmptyArray() {
    final Quadruple[][] sourceData = new Quadruple[0][];

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(sourceData),
        "Constructor DoubleMatrix(double[][]) should throw IllegalArgumentException when given an empty array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is empty'").
        contains("can't create matrix").contains("is empty");
  }

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][]) throws IllegalArgumentException when given a non-square array")
  void testConstructorWithQuadrupleArrayThrowsExceptionGivenNonSquareArray() {
    final Quadruple[][] sourceData = convertToQuadruples(NONSQUARE_ARRAY);

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(sourceData),
        "Constructor DoubleMatrix(Quadruple[][]) should throw IllegalArgumentException when given a non-square array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'must be square'").
        contains("can't create matrix").contains("must be square");
  }

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains NaN")
  void testConstructorWithQuadrupleArrayWithNaNThrowsException() {
    final Quadruple[][] badMatrixData = convertToQuadruples(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = Quadruple.nan();

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(badMatrixData),
        "Constructor DoubleMatrix(Quadruple[][]) should throw IllegalArgumentException if the argument contains NaN");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain NaN or Infinity'").
        contains("can't create matrix").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains infinity")
  void testConstructorWithQuadrupleArrayWithInfinityThrowsException() {
    final Quadruple[][] badMatrixData = convertToQuadruples(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = Quadruple.negativeInfinity();

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(badMatrixData),
        "Constructor DoubleMatrix(Quadruple[][]) should throw IllegalArgumentException if the argument contains infinity");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain NaN or Infinity'").
        contains("can't create matrix").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains null")
  void testConstructorWithQuadrupleArrayWithNullThrowsException() {
    final Quadruple[][] badMatrixData = convertToQuadruples(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = null;

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(badMatrixData),
        "Constructor DoubleMatrix(Quadruple[][]) should throw IllegalArgumentException if the argument contains null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain null'").
        contains("can't create matrix").contains("must not contain null");
  }

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][]) creates a matrix whose size equal the source's size")
  void testConstructorWithQuadrupleArrayMakesMatrixOfEqualSize() {
    final int expectedSize = SAMPLE_QUADRUPLE_DATA.length;

    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_QUADRUPLE_DATA);
    final int actualSize = newMatrix.getSize();

    assertThat(actualSize).
        withFailMessage("Newly created matrix's size (%d) is not equal to the size of the source matrix (%d)",
            actualSize, expectedSize).
        isEqualTo(expectedSize);
  }

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][]) sets scaling to true, when default scaling is set to true")
  void testConstructorWithQuadrupleArrayProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    DoubleMatrix.setDefaultScaling(expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_QUADRUPLE_DATA);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the default scaling is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][]) sets scaling to false, when default scaling is set to false")
  void testConstructorWithQuadrupleArrayProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    DoubleMatrix.setDefaultScaling(expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_QUADRUPLE_DATA);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the default scaling is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][], true) creates a matrix whose neetToScale flag is set to true ")
  void testConstructorWithQuadrupleArrayAndNeedToScaleProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    DoubleMatrix.setDefaultScaling(!expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_QUADRUPLE_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the needToScale argument is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][], false) creates a matrix whose neetToScale flag is set to false ")
  void testConstructorWithQuadrupleArrayAndNeedToScaleProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    DoubleMatrix.setDefaultScaling(!expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_QUADRUPLE_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the needToScale argument is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(Quadruple[][]) creates a matrix with data equal to the source data")
  void testConstructorWithQuadrupleArrayCreatesInstanceWithCorrectData() {
    final double[][] expectedData = convertToDoubles(SAMPLE_QUADRUPLE_DATA);

    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_QUADRUPLE_DATA);
    final double[][] actualData = newMatrix.getDoubleData();

    assertThat(actualData).
        withFailMessage("Constructor DoubleMatrix(Quadruple[][]) creates a matrix containing data not equal to the source data").
        isEqualTo(expectedData);
  }

  /*################################################################################
   * public DoubleMatrix(Quadruple[][] source, boolean needToScale)
   * Doesn't actually require special tests since it's used by the already tested constructor
   #################################################################################*/

  /*################################################################################
   * Test methods for
   *   public DoubleMatrix(Number[][]) with an array of BigDecimals as the parameter
   * Tested behavior:
   *  -- DoubleMatrix(BigDecimal[][]) throws NullPointerException when called with null
   *  -- DoubleMatrix(BigDecimal[][]) throws IllegalArgumentException when given an empty array
   *  -- DoubleMatrix(BigDecimal[][]) throws IllegalArgumentException when given a non-square array
   *  -- DoubleMatrix(BigDecimal[][]) throws IllegalArgumentException if the argument contains null
   *  -- DoubleMatrix(BigDecimal[][]) creates a matrix whose size equal the source's size
   *  -- DoubleMatrix(BigDecimal[][]) sets scaling to true, when default scaling is set to true
   *  -- DoubleMatrix(BigDecimal[][]) sets scaling to false, when default scaling is set to false
   *  -- DoubleMatrix(BigDecimal[][], true) creates a matrix whose neetToScale flag is set to true
   *  -- DoubleMatrix(BigDecimal[][], false) creates a matrix whose neetToScale flag is set to false
   *  -- DoubleMatrix(BigDecimal[][]) creates a matrix with data equal to the source matrix's data
   #################################################################################*/

  @Test
  @DisplayName("DoubleMatrix(BigDecimal[][]) throws NullPointerException when called with null")
  void testConstructorWithBigDecimalArrayThrowsExceptionGivenNull() {

    final Throwable expectedException = assertThrows( NullPointerException.class,
        () -> new DoubleMatrix((BigDecimal[][])null),
        "Constructor DoubleMatrix(BigDecimal[][]]) should throw NullPointerException when called with null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is null'").
        contains("can't create matrix").contains("is null");
  }

  @Test
  @DisplayName("DoubleMatrix(BigDecimal[][]) throws IllegalArgumentException when given an empty array")
  void testConstructorWithBigDecimalArrayThrowsExceptionGivenEmptyArray() {
    final BigDecimal[][] sourceData = new BigDecimal[0][];

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(sourceData),
        "Constructor DoubleMatrix(double[][]) should throw IllegalArgumentException when given an empty array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is empty'").
        contains("can't create matrix").contains("is empty");
  }

  @Test
  @DisplayName("DoubleMatrix(BigDecimal[][]) throws IllegalArgumentException when given a non-square array")
  void testConstructorWithBigDecimalArrayThrowsExceptionGivenNonSquareArray() {
    final BigDecimal[][] sourceData = convertToBigDecimals(NONSQUARE_ARRAY);

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(sourceData),
        "Constructor DoubleMatrix(BigDecimal[][]) should throw IllegalArgumentException when given a non-square array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'must be square'").
        contains("can't create matrix").contains("must be square");
  }

  @Test
  @DisplayName("DoubleMatrix(BigDecimal[][]) throws IllegalArgumentException if the argument contains null")
  void testConstructorWithBigDecimalArrayWithNullThrowsException() {
    final BigDecimal[][] badMatrixData = convertToBigDecimals(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = null;

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new DoubleMatrix(badMatrixData),
        "Constructor DoubleMatrix(Quadruple[][]) should throw IllegalArgumentException if the argument contains null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain null'").
        contains("can't create matrix").contains("must not contain null");
  }

  @Test
  @DisplayName("DoubleMatrix(BigDecimal[][]) creates a matrix whose size equal the source's size")
  void testConstructorWithBigDecimalArrayMakesMatrixOfEqualSize() {
    final int expectedSize = SAMPLE_BIGDECIMAL_DATA.length;

    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_BIGDECIMAL_DATA);
    final int actualSize = newMatrix.getSize();

    assertThat(actualSize).
        withFailMessage("Newly created matrix's size (%d) is not equal to the size of the source matrix (%d)",
            actualSize, expectedSize).
        isEqualTo(expectedSize);
  }

  @Test
  @DisplayName("DoubleMatrix(BigDecimal[][]) sets scaling to true, when default scaling is set to true")
  void testConstructorWithBigDecimalArrayProperlySetsScalingToTrue() {
    final boolean expectedValue = true;
    DoubleMatrix.setDefaultScaling(expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_BIGDECIMAL_DATA);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the default scaling is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(BigDecimal[][]) sets scaling to false, when default scaling is set to false")
  void testConstructorWithBigDecimalArrayProperlySetsScalingToFalse() {
    final boolean expectedValue = false;
    DoubleMatrix.setDefaultScaling(expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_BIGDECIMAL_DATA);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the default scaling is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(BigDecimal[][], true) creates a matrix whose neetToScale flag is set to true ")
  void testConstructorWithBigDecimalArrayAndNeedToScaleProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    DoubleMatrix.setDefaultScaling(!expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_BIGDECIMAL_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the needToScale argument is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(BigDecimal[][], false) creates a matrix whose neetToScale flag is set to false ")
  void testConstructorWithBigDecimalArrayAndNeedToScaleProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    DoubleMatrix.setDefaultScaling(!expectedValue);
    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_BIGDECIMAL_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the needToScale argument is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("DoubleMatrix(BigDecimal[][]) creates a matrix with data equal to the source matrix's data")
  void testConstructorWithBigDecimalArrayCreatesInstanceWithCorrectData() {
    final double[][] expectedData = convertToDoubles(SAMPLE_BIGDECIMAL_DATA);

    final DoubleMatrix newMatrix = new DoubleMatrix(SAMPLE_BIGDECIMAL_DATA);
    final double[][] actualData = newMatrix.getDoubleData();

    assertThat(actualData).
        withFailMessage("Constructor DoubleMatrix(BigDecimal[][]) creates a matrix containing data not equal to the source data").
        isEqualTo(expectedData);
  }

  // public boolean getScaling() is actually already tested by the previous tests

  /*#########################################################################################
   * Tests for the methods:
   *   public boolean getSize()
   * Tested behavior:
   *  -- getSize() returns the correct value
   #################################################################################*/

  @Test
  @DisplayName("getSize() returns the correct value")
  void testGetSizeReturnsCorrectValue() {
    final int expectedSize = SOME_SIZE;

    final DoubleMatrix matrix = makeAMatrixOfSize(expectedSize);
    final int actualSize = matrix.getSize();

    assertThat(actualSize).
        withFailMessage("getSize() returns wrong value: (%d) instead of expected (%d)", actualSize, expectedSize).
        isEqualTo(expectedSize);
  }

  /*#########################################################################################
   * Tests for the methods:
   *   public Number[][] getData()
   *  -- getData() returns the correct values
   *  -- getData() returns a copy of the internal data, not a reference to it
   #################################################################################*/

  @Test
  @DisplayName("getData() returns the correct values")
  void testGetDataReturnsCorrectValues() {
    final Number[][] expectedData = convertToDoublesAsNumbers(SAMPLE_MATRIX_DATA);
    final Number[][] actualData =  SAMPLE_DOUBLE_MATRIX.getData();

    assertThat(actualData).
        withFailMessage("getData() returns wrong values").
        isEqualTo(expectedData);
  }

  @Test
  @DisplayName("getData() returns a copy of the internal data, not a reference to it")
  void testGetDataReturnsACopyOfInternalDataNotAReferenceToIt() {
    final Number[][] expectedData = convertToDoublesAsNumbers(SAMPLE_MATRIX_DATA);

    final Number[][] tmpData =  SAMPLE_DOUBLE_MATRIX.getData();
    tmpData[0][0] = (Double)tmpData[0][0] + ANY_NUMBER_BUT_ZERO;
    final Number[][] actualData =  SAMPLE_DOUBLE_MATRIX.getData();

    assertThat(actualData).
        withFailMessage("getData() returns a reference to the internal data").
        isEqualTo(expectedData);
  }

  /*#########################################################################################
   * Tests for the methods:
   *   public Number[][] getDoubleData()
   *  -- getDoubleData() returns the correct values
   *  -- getDoubleData() returns a copy of the internal data, not a reference to it
   #################################################################################*/

  @Test
  @DisplayName("getDoubleData() returns the correct values")
  void testGetDoubleDataReturnsCorrectValues() {
    final double[][] expectedData = SAMPLE_MATRIX_DATA;

    final double[][] actualData =  SAMPLE_DOUBLE_MATRIX.getDoubleData();

    assertThat(actualData).
        withFailMessage("getDoubleData() returns wrong values").
        isEqualTo(expectedData);
  }

  @Test
  @DisplayName("getDoubleData() returns a copy of the internal data, not a reference to it")
  void testGetDoubleDataReturnsACopyOfInternalDataNotAReferenceToIt() {
    final double[][] expectedData = SAMPLE_MATRIX_DATA;

    final double[][] tmpData =  SAMPLE_DOUBLE_MATRIX.getDoubleData();
    tmpData[0][0] += ANY_NUMBER_BUT_ZERO;
    final double[][] actualData =  SAMPLE_DOUBLE_MATRIX.getDoubleData();

    assertThat(actualData).
        withFailMessage("getDoubleData() returns a reference to the internal data").
        isEqualTo(expectedData);
  }


  /*#########################################################################################
   * Tests for the methods:
   *   public Quadruple[][] getQuadrupleData()
   *  -- getQuadrupleData() returns the correct data
   *  -- getQuadrupleData() returns a copy of the internal data, not a reference to it
   #################################################################################*/

  @Test
  @DisplayName("getQuadrupleData() returns the correct data")
  void testGetQuadrupleDataReturnsCorrectData() {
    final Quadruple[][] expectedData = convertToQuadruples(SAMPLE_MATRIX_DATA);

    final Quadruple[][] actualData =  SAMPLE_DOUBLE_MATRIX.getQuadrupleData();

    assertThat(actualData).
        withFailMessage("getQuadrupleData() returns wrong values").
        isEqualTo(expectedData);
  }

  @Test
  @DisplayName("getQuadrupleData() returns a copy of the internal data, not a reference to it")
  void testGetQuadrupleDataReturnsACopyOfInternalDataNotAReferenceToIt() {
    final Quadruple[][] expectedData = convertToQuadruples(SAMPLE_MATRIX_DATA);

    final Quadruple[][] tmpData =  SAMPLE_DOUBLE_MATRIX.getQuadrupleData();
    tmpData[0][0].add(ANY_NUMBER_BUT_ZERO);
    final Quadruple[][] actualData =  SAMPLE_DOUBLE_MATRIX.getQuadrupleData();

    assertThat(actualData).
        withFailMessage("getQuadrupleData() returns a reference to the internal data").
        isEqualTo(expectedData);
  }

  /*#########################################################################################
   * Tests for the methods:
   *   public BigDecimal[][] getBigDecimalData()
   *  -- getBigDecimalData() returns the correct data
   *  -- getBigDecimalData() returns a copy of the internal data, not a reference to it
   #################################################################################*/

  @Test
  @DisplayName("getBigDecimalData() returns the correct data")
  void testGetBigDecimalDataReturnsCorrectData() {
    final BigDecimal[][] expectedData = convertToBigDecimals(SAMPLE_MATRIX_DATA);

    final BigDecimal[][] actualData =  SAMPLE_DOUBLE_MATRIX.getBigDecimalData();

    assertThat(actualData).
        withFailMessage("getBigDecimaleData() returns wrong values").
        isEqualTo(expectedData);
  }

  @Test
  @DisplayName("getBigDecimalData() returns a copy of the internal data, not a reference to it")
  void testGetBigDecimalDataReturnsACopyOfInternalDataNotAReferenceToIt() {
    final BigDecimal[][] expectedData = convertToBigDecimals(SAMPLE_MATRIX_DATA);

    final BigDecimal[][] tmpData =  SAMPLE_DOUBLE_MATRIX.getBigDecimalData();
    tmpData[0][0] = tmpData[0][0].add(BigDecimal.valueOf(ANY_NUMBER_BUT_ZERO));
    final BigDecimal[][] actualData =  SAMPLE_DOUBLE_MATRIX.getBigDecimalData();

    assertThat(actualData).
        withFailMessage("getBigDecimaleData() returns a reference to the internal data").
        isEqualTo(expectedData);
  }

  /*#########################################################################################
   * Tests for the method:
   *   public boolean equals(Object anotherOne)
   * Tested behavior:
   *  -- equals(Object) returns false if the argument is null
   *  -- equals(Object) with Matrix returns false if the argument is not a DoubleMatrix
   *  -- equals(Object) with Matrix returns false if the argument has a different size
   *  -- equals(Matrix) with Matrix returns false if the argument has a different value of the needToScale flag
   *  -- equals(Matrix) with Matrix returns false if the argument has a different matrix data
   *  -- equals(Matrix) with Matrix returns true if both matrix data and needToScale are equal
   #################################################################################*/


  @Test
  @DisplayName("equals(Object) returns false if the argument is null")
  void testEqualsReturnsFalseIfArgumentIsNull() {
    final boolean expectedResult = false;

    final boolean actualResult = SAMPLE_DOUBLE_MATRIX.equals(null);

    assertThat(actualResult).
        withFailMessage("equals() returns true when called with null argument").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("equals(Object) with Matrix returns false if the argument is not a DoubleMatrix")
  void testEqualsReturnsFalseIfArgumentIsNotDoubleMatrix() {
    final boolean expectedResult = false;

    @SuppressWarnings("unlikely-arg-type")
    final boolean actualResult = SAMPLE_DOUBLE_MATRIX.equals(new BigDecimalMatrix(SAMPLE_MATRIX_DATA));

    assertThat(actualResult).
        withFailMessage("equals() returns true when called with argument of another type").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("equals(Object) with Matrix returns false if the argument has a different size")
  void testEqualsReturnsFalseIfArgumentHasAnotherSize() {
    final boolean expectedResult = false;

    final boolean actualResult = SAMPLE_DOUBLE_MATRIX.equals(new DoubleMatrix(SAMPLE_TOO_LARGE_MATRIX_DATA));

    assertThat(actualResult).
        withFailMessage("equals() returns true when the argument has a different size").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("equals(Object) with Matrix returns false if the argument has a different value of the needToScale flag")
  void testEqualsReturnsFalseIfArgumentHasAnotherNeedToScale()  {
    final boolean expectedResult = false;

    final DoubleMatrix matrix1 = new DoubleMatrix(SAMPLE_MATRIX_DATA, false);
    final DoubleMatrix matrix2 = new DoubleMatrix(SAMPLE_MATRIX_DATA, true);
    final boolean actualResult = matrix1.equals(matrix2);

    assertThat(actualResult).
        withFailMessage("equals() returns true when called with matrix with another value of needToScale flag").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("equals(Object) with Matrix returns false if the argument has a different matrix data")
  void testEqualsReturnsFalseIfArgumentHasAnotherData() {
    final boolean expectedResult = false;

    final DoubleMatrix matrix1 = new DoubleMatrix(SAMPLE_MATRIX_DATA);
    final DoubleMatrix matrix2 = new DoubleMatrix(SAMPLE_ANOTHER_MATRIX_DATA);
    final boolean actualResult = matrix1.equals(matrix2);

    assertThat(actualResult).
        withFailMessage("equals() returns true when called with matrix with another data").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("equals(Object) with Matrix returns true if both matrix data and needToScale are equal")
  void testEqualsReturnsTrueIfArgumentsAreEqual() {
    final boolean expectedResult = true;

    final DoubleMatrix matrix1 = new DoubleMatrix(SAMPLE_MATRIX_DATA);
    final DoubleMatrix matrix2 = new DoubleMatrix(SAMPLE_MATRIX_DATA);
    final boolean actualResult = matrix1.equals(matrix2);

    assertThat(actualResult).
        withFailMessage("equals() returns false when called with equal matrix").
        isEqualTo(expectedResult);
  }

  /*#########################################################################################
   * Tests for the methods:
   *   public abstract int hashCode()
   * Tested behaviors:
   *  -- hashCode() returns different values for matrices with different values of the needToScale flag
   *  -- hashCode() returns different values for matrices with different matrix data
   *  -- hashCode() returns equal values for equal matrices
   #################################################################################*/

  @Test
  @DisplayName("hashCode() returns different values for matrices with different values of the needToScale flag")
  void testHashCodeReturnsDifferentResultsForDifferentNeedToScaleValues()  {
    final int result1 = new DoubleMatrix(SAMPLE_MATRIX_DATA, false).hashCode();
    final int result2 = new DoubleMatrix(SAMPLE_MATRIX_DATA, true).hashCode();

    assertThat(result1).
        withFailMessage("hashCode() returns equal values for matrices with different values of the needToScale flag").
        isNotEqualTo(result2);
  }

  @Test
  @DisplayName("hashCode() returns different values for matrices with different matrix data")
  void testHashCodeReturnsDifferentResultsForDifferentMatrixData()  {

    final int result1 = new DoubleMatrix(SAMPLE_MATRIX_DATA).hashCode();
    final int result2 = new DoubleMatrix(SAMPLE_ANOTHER_MATRIX_DATA).hashCode();

    assertThat(result1).
        withFailMessage("hashCode() returns equal values for matrices with different data").
        isNotEqualTo(result2);
  }

  @Test
  @DisplayName("hashCode() returns equal values for equal matrices")
  void testHashCodeReturnsEqualResultsForEqualMatrices()  {

    final int result1 = new DoubleMatrix(SAMPLE_MATRIX_DATA).hashCode();
    final int result2 = new DoubleMatrix(SAMPLE_MATRIX_DATA).hashCode();

    assertThat(result1).
        withFailMessage("hashCode() returns different values for equal matrices").
        isEqualTo(result2);
  }

  /* *******************************************************************************
  /***** Private methods ***********************************************************
  /*********************************************************************************/

  private static DoubleMatrix makeAMatrixOfSize(int expectedSize) {
    final double[][] data = new double[expectedSize][expectedSize];
    return new DoubleMatrix(data);
  }

}
