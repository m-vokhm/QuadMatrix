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
 * Tests for {@link BigDecimalMatrix} constructors and matrix data getters
 * <br>
 * The methods to be tested:
 *   public static void setDefaultScaling(boolean scaleByDefault)
 *   public static boolean getDefaultScaling()
 *   public static void setDefaultPrecision(int precision)
 *   public static int getDefaultPrecision()
 *
 *   public BigDecimalMatrix(Matrix source)
 *   public BigDecimalMatrix(Matrix source, boolean needToScale)
 *   public BigDecimalMatrix(double[][] source)
 *   public BigDecimalMatrix(double[][] source, boolean needToScale)
 *   public BigDecimalMatrix(Number[][] source)
 *   public BigDecimalMatrix(Number[][] source, boolean needToScale)
 *
 *   public BigDecimalMatrix(Matrix source, int precision)
 *   public BigDecimalMatrix(Matrix source, boolean needToScale, int precision)
 *   public BigDecimalMatrix(double[][] source, int precision)
 *   public BigDecimalMatrix(double[][] source, boolean needToScale, int precision)
 *   public BigDecimalMatrix(Number[][] source, int precision)
 *   public BigDecimalMatrix(Number[][] source, boolean needToScale, int precision)
 *
 *   public boolean getScaling()
 *   public int getSize()
 *   public Number[][] getData()
 *   public double[][] getDoubleData()
 *   public Quadruple[][] getQuadrupleData()
 *   public BigDecimal[][] getBigDecimalData()
 *   public abstract boolean equals(Object anotherOne);
 *   public abstract int hashCode();
 */
public class BigDecimalMatrixEnvelopeTests {

  /* *******************************************************************************
  /***** Testing methods that don't involve BigDecimalMatrixSolver ********************
  /***** Just data manipulations ***************************************************
  /*********************************************************************************/

  private static final int SOME_SIZE                            = 1234;
  private static final int SAMPLE_MATRIX_SIZE                   = 10;
  private static final Random RANDOM                            = new Random(12345);
  private static final double[][] SAMPLE_MATRIX_DATA            = new double[SAMPLE_MATRIX_SIZE][SAMPLE_MATRIX_SIZE];
  private static final double[][] SAMPLE_TOO_LARGE_MATRIX_DATA  = new double[SAMPLE_MATRIX_SIZE + 1][SAMPLE_MATRIX_SIZE + 1];
  private static final double[][] NONSQUARE_ARRAY               = new double[SAMPLE_MATRIX_SIZE][SAMPLE_MATRIX_SIZE + 1];
  private static final Quadruple[][] SAMPLE_QUADRUPLE_DATA      = new Quadruple[SAMPLE_MATRIX_SIZE][SAMPLE_MATRIX_SIZE];
  private static final BigDecimal[][] SAMPLE_BIGDECIMAL_DATA    = new BigDecimal[SAMPLE_MATRIX_SIZE][SAMPLE_MATRIX_SIZE];
  private static final BigDecimalMatrix SAMPLE_BIGDECIMAL_MATRIX;
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
    SAMPLE_BIGDECIMAL_MATRIX = new BigDecimalMatrix(SAMPLE_MATRIX_DATA);
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

    BigDecimalMatrix.setDefaultScaling(!expectedValue);
    BigDecimalMatrix.setDefaultScaling(expectedValue);
    final boolean actualValue = BigDecimalMatrix.getDefaultScaling();

    assertThat(actualValue).
        withFailMessage("getDefaultScaling() returns false when true was set with setDefaultScaling()").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("getDefaultScaling() returns false after false was set with setDefaultScaling()")
  void testGetDefaultScalingReturnsFalseAfterSettingFalse() {
    final boolean expectedValue = false;

    BigDecimalMatrix.setDefaultScaling(!expectedValue);
    BigDecimalMatrix.setDefaultScaling(expectedValue);
    final boolean actualValue = BigDecimalMatrix.getDefaultScaling();

    assertThat(actualValue).
        withFailMessage("getDefaultScaling() returns true when false was set with setDefaultScaling()").
        isEqualTo(expectedValue);
  }

  /*################################################################################
   * Test methods for
   *   public static void setDefaultPrecision(int precision)
   *   public static int getDefaultPrecision()
   * tested behavior:
   *   getDefaultPrecision() returns the value that was set with setDefaultPrecision()
   #################################################################################*/

  @Test
  @DisplayName("getDefaultPrecision() returns the value that was set with setDefaultPrecision()")
  void testGetDefaultPrecisionReturnsValueSetBySetDefaultPrecision() {
    final int expectedPrecision = BigDecimalMatrix.getDefaultPrecision() + 10;

    BigDecimalMatrix.setDefaultPrecision(expectedPrecision);
    final int actualPrecision = BigDecimalMatrix.getDefaultPrecision();

    assertThat(actualPrecision).
        withFailMessage("getDefaultPrecision() returns value other than the value set with setDefaultPrecision()").
        isEqualTo(expectedPrecision);
  }

  /*################################################################################
   * Test methods for
   *   public BigDecimalMatrix(Matrix source)
   * Tested behavior:
   *  -- BigDecimalMatrix(Matrix) throws NullPointerException with a relevant message when called with null
   *  -- BigDecimalMatrix(Matrix) creates a matrix whose size equal the source's size
   *  -- BigDecimalMatrix(Matrix) sets scaling to true, when default scaling is set to true
   *  -- BigDecimalMatrix(Matrix) sets scaling to false, when default scaling is set to false
   *  -- BigDecimalMatrix(Matrix, true) creates a matrix whose neetToScale flag is set to true
   *  -- BigDecimalMatrix(Matrix, false) creates a matrix whose neetToScale flag is set to false
   *  -- BigDecimalMatrix(Matrix) creates a matrix whose data equal the source's data
   *  -- BigDecimalMatrix(Matrix) creates a matrix with the default precision
   *  -- BigDecimalMatrix(Matrix, precision) creates a matrix with specified precision
   #################################################################################*/

  @Test
  @DisplayName("BigDecimalMatrix(Matrix) throws NullPointerException when called with null")
  void testConstructorWithMatrixThrowsExceptionGivenNull() {

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> new BigDecimalMatrix((Matrix)null),
        "BigDecimalMatrix(Matrix) should throw NullPointerException when called with null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is null'").
        contains("can't create matrix").contains("is null");
  }

  @Test
  @DisplayName("BigDecimalMatrix(Matrix) creates a matrix whose size equal the source's size")
  void testConstructorWithMatrixMakesMatrixOfEqualSize() {
    final int expectedSize = SAMPLE_BIGDECIMAL_MATRIX.getSize();

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_MATRIX);
    final int actualSize = newMatrix.getSize();

    assertThat(actualSize).
        withFailMessage("Newly created matrix's size (%d) is not equal to the size of the source matrix (%d)",
            actualSize, expectedSize).
        isEqualTo(expectedSize);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Matrix) sets scaling to true, when default scaling is set to true")
  void testConstructorWithMatrixProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    BigDecimalMatrix.setDefaultScaling(expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_MATRIX);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the default scaling is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Matrix) sets scaling to false, when default scaling is set to false")
  void testConstructorWithMatrixProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    BigDecimalMatrix.setDefaultScaling(expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_MATRIX);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the default scaling is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Matrix, true) creates a matrix whose neetToScale flag is set to true ")
  void testConstructorWithMatrixAndNeedToScaleProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    BigDecimalMatrix.setDefaultScaling(!expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_MATRIX, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the needToScale argument is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Matrix, false) creates a matrix whose neetToScale flag is set to false ")
  void testConstructorWithMatrixAndNeedToScaleProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    BigDecimalMatrix.setDefaultScaling(!expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_MATRIX, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the needToScale argument is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Matrix) creates a matrix whose data equal the source's data")
  void testConstructorWithMatrixMakesMatrixWithEqualData() {
    final double[][] expectedData = SAMPLE_BIGDECIMAL_MATRIX.getDoubleData();

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_MATRIX);
    final double[][] actualData = newMatrix.getDoubleData();

    assertThat(actualData).
        withFailMessage("Newly created matrix contains data not equal to the data of the source matrix").
        isEqualTo(expectedData);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Matrix) creates a matrix with the default precision")
  void testConstructorWithMatrixMakesMatrixWithDefaultPrecision() {
    final int expectedPrecision = BigDecimalMatrix.getDefaultPrecision();

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_MATRIX);
    final int actualPrecision = newMatrix.getPrecision();

    assertThat(actualPrecision).
        withFailMessage("Newly created matrix has precision other than the default precision").
        isEqualTo(expectedPrecision);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Matrix, int) creates a matrix with the speciafied precision")
  void testConstructorWithMatrixAndPrecisionMakesMatrixWithSpecifiedPrecision() {
    final int expectedPrecision = BigDecimalMatrix.getDefaultPrecision() + 123;

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_MATRIX, expectedPrecision);
    final int actualPrecision = newMatrix.getPrecision();

    assertThat(actualPrecision).
        withFailMessage("Newly created matrix has precision other than the default precision").
        isEqualTo(expectedPrecision);
  }

  /*################################################################################
   * public BigDecimalMatrix(Matrix source, boolean needToScale)
   * Doesn't actually require special tests since it's used by the already tested constructor
   #################################################################################*/

  /*################################################################################
   * Test methods for
   *   public BigDecimalMatrix(double[][])
   * Tested behavior:
   *  -- BigDecimalMatrix(double[][]) throws NullPointerException when called with null
   *  -- BigDecimalMatrix(double[][]) throws IllegalArgumentException when given an empty array
   *  -- BigDecimalMatrix(double[][]) throws IllegalArgumentException when given a non-square array
   *  -- BigDecimalMatrix(double[][]) throws IllegalArgumentException if the argument contains NaN
   *  -- BigDecimalMatrix(double[][]) throws IllegalArgumentException if the argument contains infinity
   *  -- BigDecimalMatrix(double[][]) creates a matrix whose size equal the source's size
   *  -- BigDecimalMatrix(double[][]) sets scaling to true, when default scaling is set to true
   *  -- BigDecimalMatrix(double[][]) sets scaling to false, when default scaling is set to false
   *  -- BigDecimalMatrix(double[][], true) creates a matrix whose neetToScale flag is set to true
   *  -- BigDecimalMatrix(double[][], false) creates a matrix whose neetToScale flag is set to false
   *  -- BigDecimalMatrix(double[][]) creates a matrix with data equal to the source data
   *  -- BigDecimalMatrix(double[][]) creates a matrix with the default precision
   *  -- BigDecimalMatrix(double[][], precision) creates a matrix with specified precision
   #################################################################################*/

  @Test
  @DisplayName("BigDecimalMatrix(double[][]) throws NullPointerException when called with null")
  void testConstructorWithDoubleArrayThrowsExceptionGivenNull() {

    final Throwable expectedException = assertThrows( NullPointerException.class,
        () -> new BigDecimalMatrix((double[][])null),
        "Constructor BigDecimalMatrix(double[][]) should throw NullPointerException when called with null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is null'").
        contains("can't create matrix").contains("is null");
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][]) throws IllegalArgumentException when given an empty array")
  void testConstructorWithDoubleArrayThrowsExceptionGivenEmptyArray() {
    final double[][] sourceData = new double[0][];

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(sourceData),
        "Constructor BigDecimalMatrix(double[][]) should throw IllegalArgumentException when given an empty array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is empty'").
        contains("can't create matrix").contains("is empty");
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][]) throws IllegalArgumentException when given a non-square array")
  void testConstructorWithDoubleArrayThrowsExceptionGivenNonSquareArray() {
    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(NONSQUARE_ARRAY),
        "Constructor BigDecimalMatrix(double[][]) should throw IllegalArgumentException when given a non-square array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'must be square'").
        contains("can't create matrix").contains("must be square");
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][]) throws IllegalArgumentException if the argument contains NaN")
  void testConstructorWithDoubleArrayWithNaNThrowsException() {
    final double[][] badMatrixData = deepCopyOf(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = Double.NaN;

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(badMatrixData),
        "Constructor BigDecimalMatrix(double[][]) should throw IllegalArgumentException if the argument contains NaN");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain NaN or Infinity'").
        contains("can't create matrix").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][]) throws IllegalArgumentException if the argument contains infinity")
  void testConstructorWithDoubleArrayWithInfinityThrowsException() {
    final double[][] badMatrixData = deepCopyOf(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = Double.POSITIVE_INFINITY;

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(badMatrixData),
        "Constructor BigDecimalMatrix(double[][]) should throw IllegalArgumentException if the argument contains infinity");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain NaN or Infinity'").
        contains("can't create matrix").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][]) creates a matrix whose size equal the source's size")
  void testConstructorWithDoubleArrayMakesMatrixOfEqualSize() {
    final int expectedSize = SAMPLE_MATRIX_DATA.length;

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_MATRIX_DATA);
    final int actualSize = newMatrix.getSize();

    assertThat(actualSize).
        withFailMessage("Newly created matrix's size (%d) is not equal to the size of the source matrix (%d)",
            actualSize, expectedSize).
        isEqualTo(expectedSize);
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][]) sets scaling to true, when default scaling is set to true")
  void testConstructorWithDoubleArrayProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    BigDecimalMatrix.setDefaultScaling(expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_MATRIX_DATA);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the default scaling is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][], true) creates a matrix whose neetToScale flag is set to true ")
  void testConstructorWithDoubleArrayAndNeedToScaleProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    BigDecimalMatrix.setDefaultScaling(!expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the needToScale argument is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][], false) creates a matrix whose neetToScale flag is set to false ")
  void testConstructorWithDoubleArrayAndNeedToScaleProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    BigDecimalMatrix.setDefaultScaling(!expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the needToScale argument is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][], boolean) sets scaling to the specified value")
  void testConstructorWithWithDoubleArraySetsScalingToTheSpecifiedValue() {
    final boolean expectedValue = !BigDecimalMatrix.getDefaultScaling();

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the default scaling is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][]) creates a matrix with data equal to the source data")
  void testConstructorWithDoubleArrayCreatesInstanceWithCorrectData() {
    final double[][] expectedData = SAMPLE_MATRIX_DATA;

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_MATRIX_DATA);
    final double[][] actualData = newMatrix.getDoubleData();

    assertThat(actualData).
      withFailMessage("Constructor BigDecimalMatrix(double[][]) creates a matrix containing data not equal to the source data").
      isEqualTo(expectedData);
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][]) creates a matrix with the default precision")
  void testConstructorWithDoubleArrayMakesMatrixWithDefaultPrecision() {
    final int expectedPrecision = BigDecimalMatrix.getDefaultPrecision();

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_MATRIX_DATA);
    final int actualPrecision = newMatrix.getPrecision();

    assertThat(actualPrecision).
        withFailMessage("Newly created matrix has precision other than the default precision").
        isEqualTo(expectedPrecision);
  }

  @Test
  @DisplayName("BigDecimalMatrix(double[][], int) creates a matrix with the speciafied precision")
  void testConstructorWithDoubleArrayAndPrecisionMakesMatrixWithSpecifiedPrecision() {
    final int expectedPrecision = BigDecimalMatrix.getDefaultPrecision() + 123;

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, expectedPrecision);
    final int actualPrecision = newMatrix.getPrecision();

    assertThat(actualPrecision).
        withFailMessage("Newly created matrix has precision other than the default precision").
        isEqualTo(expectedPrecision);
  }

  /*################################################################################
   * public BigDecimalMatrix(double[][] source, boolean needToScale)
   * Doesn't actually require special tests since it's used by the already tested constructor
   #################################################################################*/

  /*################################################################################
   * Test methods for
   *   public BigDecimalMatrix(Number[][]) with an array of Quadruples as the parameter
   * Tested behavior:
   *  -- BigDecimalMatrix(Quadruple[][]) throws NullPointerException when called with null
   *  -- BigDecimalMatrix(Quadruple[][]) throws IllegalArgumentException when given an empty array
   *  -- BigDecimalMatrix(Quadruple[][]) throws IllegalArgumentException when given a non-square array
   *  -- BigDecimalMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains NaN
   *  -- BigDecimalMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains infinity
   *  -- BigDecimalMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains null
   *  -- BigDecimalMatrix(Quadruple[][]) creates a matrix whose size equal the source's size
   *  -- BigDecimalMatrix(Quadruple[][]) sets scaling to true, when default scaling is set to true
   *  -- BigDecimalMatrix(Quadruple[][]) sets scaling to false, when default scaling is set to false
   *  -- BigDecimalMatrix(Quadruple[][], true) creates a matrix whose neetToScale flag is set to true
   *  -- BigDecimalMatrix(Quadruple[][], false) creates a matrix whose neetToScale flag is set to false
   *  -- BigDecimalMatrix(Quadruple[][]) creates a matrix with data equal to the source data
   *  -- BigDecimalMatrix(Quadruple[][]) creates a matrix with the default precision
   *  -- BigDecimalMatrix(Quadruple[][], precision) creates a matrix with specified precision
   #################################################################################*/

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][]) throws NullPointerException when called with null")
  void testConstructorWithQuadrupleArrayThrowsExceptionGivenNull() {
    final Throwable expectedException = assertThrows( NullPointerException.class,
        () -> new BigDecimalMatrix((Quadruple[][])null),
        "Constructor BigDecimalMatrix(Quadruple[][]]) should throws NullPointerException when called with null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is null'").
        contains("can't create matrix").contains("is null");
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][]) throws IllegalArgumentException when given an empty array")
  void testConstructorWithQuadrupleArrayThrowsExceptionGivenEmptyArray() {
    final Quadruple[][] sourceData = new Quadruple[0][];

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(sourceData),
        "Constructor BigDecimalMatrix(double[][]) should throw IllegalArgumentException when given an empty array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is empty'").
        contains("can't create matrix").contains("is empty");
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][]) throws IllegalArgumentException when given a non-square array")
  void testConstructorWithQuadrupleArrayThrowsExceptionGivenNonSquareArray() {
    final Quadruple[][] sourceData = convertToQuadruples(NONSQUARE_ARRAY);

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(sourceData),
        "Constructor BigDecimalMatrix(Quadruple[][]) should throw IllegalArgumentException when given a non-square array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'must be square'").
        contains("can't create matrix").contains("must be square");
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains NaN")
  void testConstructorWithQuadrupleArrayWithNaNThrowsException() {
    final Quadruple[][] badMatrixData = convertToQuadruples(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = Quadruple.nan();

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(badMatrixData),
        "Constructor BigDecimalMatrix(Quadruple[][]) should throw IllegalArgumentException if the argument contains NaN");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain NaN or Infinity'").
        contains("can't create matrix").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains infinity")
  void testConstructorWithQuadrupleArrayWithInfinityThrowsException() {
    final Quadruple[][] badMatrixData = convertToQuadruples(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = Quadruple.negativeInfinity();

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(badMatrixData),
        "Constructor BigDecimalMatrix(Quadruple[][]) should throw IllegalArgumentException if the argument contains infinity");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain NaN or Infinity'").
        contains("can't create matrix").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][]) throws IllegalArgumentException if the argument contains null")
  void testConstructorWithQuadrupleArrayWithNullThrowsException() {
    final Quadruple[][] badMatrixData = convertToQuadruples(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = null;

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(badMatrixData),
        "Constructor BigDecimalMatrix(Quadruple[][]) should throw IllegalArgumentException if the argument contains null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain null'").
        contains("can't create matrix").contains("must not contain null");
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][]) creates a matrix whose size equal the source's size")
  void testConstructorWithQuadrupleArrayMakesMatrixOfEqualSize() {
    final int expectedSize = SAMPLE_QUADRUPLE_DATA.length;

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_QUADRUPLE_DATA);
    final int actualSize = newMatrix.getSize();

    assertThat(actualSize).
        withFailMessage("Newly created matrix's size (%d) is not equal to the size of the source matrix (%d)",
            actualSize, expectedSize).
        isEqualTo(expectedSize);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][]) sets scaling to true, when default scaling is set to true")
  void testConstructorWithQuadrupleArrayProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    BigDecimalMatrix.setDefaultScaling(expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_QUADRUPLE_DATA);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the default scaling is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][]) sets scaling to false, when default scaling is set to false")
  void testConstructorWithQuadrupleArrayProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    BigDecimalMatrix.setDefaultScaling(expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_QUADRUPLE_DATA);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the default scaling is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][], true) creates a matrix whose neetToScale flag is set to true ")
  void testConstructorWithQuadrupleArrayAndNeedToScaleProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    BigDecimalMatrix.setDefaultScaling(!expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_QUADRUPLE_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the needToScale argument is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][], false) creates a matrix whose neetToScale flag is set to false ")
  void testConstructorWithQuadrupleArrayAndNeedToScaleProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    BigDecimalMatrix.setDefaultScaling(!expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_QUADRUPLE_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the needToScale argument is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][]) creates a matrix with data equal to the source data")
  void testConstructorWithQuadrupleArrayCreatesInstanceWithCorrectData() {
    final double[][] expectedData = convertToDoubles(SAMPLE_QUADRUPLE_DATA);

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_QUADRUPLE_DATA);
    final double[][] actualData = newMatrix.getDoubleData();

    assertThat(actualData).
        withFailMessage("Constructor BigDecimalMatrix(Quadruple[][]) creates a matrix containing data not equal to the source data").
        isEqualTo(expectedData);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][]) creates a matrix with the default precision")
  void testConstructorWithQuadrupleArrayMakesMatrixWithDefaultPrecision() {
    final int expectedPrecision = BigDecimalMatrix.getDefaultPrecision();

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_QUADRUPLE_DATA);
    final int actualPrecision = newMatrix.getPrecision();

    assertThat(actualPrecision).
        withFailMessage("Newly created matrix has precision other than the default precision").
        isEqualTo(expectedPrecision);
  }

  @Test
  @DisplayName("BigDecimalMatrix(Quadruple[][], int) creates a matrix with the speciafied precision")
  void testConstructorWithQuadrupleArrayAndPrecisionMakesMatrixWithSpecifiedPrecision() {
    final int expectedPrecision = BigDecimalMatrix.getDefaultPrecision() + 123;

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_QUADRUPLE_DATA, expectedPrecision);
    final int actualPrecision = newMatrix.getPrecision();

    assertThat(actualPrecision).
        withFailMessage("Newly created matrix has precision other than the default precision").
        isEqualTo(expectedPrecision);
  }

  /*################################################################################
   * public BigDecimalMatrix(Quadruple[][] source, boolean needToScale)
   * Doesn't actually require special tests since it's used by the already tested constructor
   #################################################################################*/

  /*################################################################################
   * Test methods for
   *   public BigDecimalMatrix(Number[][]) with an array of BigDecimals as the parameter
   * Tested behavior:
   *  -- BigDecimalMatrix(BigDecimal[][]) throws NullPointerException when called with null
   *  -- BigDecimalMatrix(BigDecimal[][]) throws IllegalArgumentException when given an empty array
   *  -- BigDecimalMatrix(BigDecimal[][]) throws IllegalArgumentException when given a non-square array
   *  -- BigDecimalMatrix(BigDecimal[][]) throws IllegalArgumentException if the argument contains null
   *  -- BigDecimalMatrix(BigDecimal[][]) creates a matrix whose size equal the source's size
   *  -- BigDecimalMatrix(BigDecimal[][]) sets scaling to true, when default scaling is set to true
   *  -- BigDecimalMatrix(BigDecimal[][]) sets scaling to false, when default scaling is set to false
   *  -- BigDecimalMatrix(BigDecimal[][], true) creates a matrix whose neetToScale flag is set to true
   *  -- BigDecimalMatrix(BigDecimal[][], false) creates a matrix whose neetToScale flag is set to false
   *  -- BigDecimalMatrix(BigDecimal[][]) creates a matrix with data equal to the source matrix's data
   *  -- BigDecimalMatrix(BigDecimal[][]) creates a matrix with the default precision
   *  -- BigDecimalMatrix(BigDecimal[][], precision) creates a matrix with specified precision
   #################################################################################*/

  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][]) throws NullPointerException when called with null")
  void testConstructorWithBigDecimalArrayThrowsExceptionGivenNull() {

    final Throwable expectedException = assertThrows( NullPointerException.class,
        () -> new BigDecimalMatrix((BigDecimal[][])null),
        "Constructor BigDecimalMatrix(BigDecimal[][]]) should throw NullPointerException when called with null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is null'").
        contains("can't create matrix").contains("is null");
  }

  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][]) throws IllegalArgumentException when given an empty array")
  void testConstructorWithBigDecimalArrayThrowsExceptionGivenEmptyArray() {
    final BigDecimal[][] sourceData = new BigDecimal[0][];

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(sourceData),
        "Constructor BigDecimalMatrix(double[][]) should throw IllegalArgumentException when given an empty array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'is empty'").
        contains("can't create matrix").contains("is empty");
  }

  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][]) throws IllegalArgumentException when given a non-square array")
  void testConstructorWithBigDecimalArrayThrowsExceptionGivenNonSquareArray() {
    final BigDecimal[][] sourceData = convertToBigDecimals(NONSQUARE_ARRAY);

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(sourceData),
        "Constructor BigDecimalMatrix(BigDecimal[][]) should throw IllegalArgumentException when given a non-square array");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'Can't create matrix' and 'must be square'").
        contains("can't create matrix").contains("must be square");
  }

  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][]) throws IllegalArgumentException if the argument contains null")
  void testConstructorWithBigDecimalArrayWithNullThrowsException() {
    final BigDecimal[][] badMatrixData = convertToBigDecimals(SAMPLE_MATRIX_DATA);
    badMatrixData[1][1] = null;

    final Throwable expectedException = assertThrows( IllegalArgumentException.class,
        () -> new BigDecimalMatrix(badMatrixData),
        "Constructor BigDecimalMatrix(Quadruple[][]) should throw IllegalArgumentException if the argument contains null");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'can't create matrix' and 'must not contain null'").
        contains("can't create matrix").contains("must not contain null");
  }

  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][]) creates a matrix whose size equal the source's size")
  void testConstructorWithBigDecimalArrayMakesMatrixOfEqualSize() {
    final int expectedSize = SAMPLE_BIGDECIMAL_DATA.length;

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_DATA);
    final int actualSize = newMatrix.getSize();

    assertThat(actualSize).
        withFailMessage("Newly created matrix's size (%d) is not equal to the size of the source matrix (%d)",
            actualSize, expectedSize).
        isEqualTo(expectedSize);
  }

  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][]) sets scaling to true, when default scaling is set to true")
  void testConstructorWithBigDecimalArrayProperlySetsScalingToTrue() {
    final boolean expectedValue = true;
    BigDecimalMatrix.setDefaultScaling(expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_DATA);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the default scaling is true").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][]) sets scaling to false, when default scaling is set to false")
  void testConstructorWithBigDecimalArrayProperlySetsScalingToFalse() {
    final boolean expectedValue = false;
    BigDecimalMatrix.setDefaultScaling(expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_DATA);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the default scaling is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][], true) creates a matrix whose neetToScale flag is set to true ")
  void testConstructorWithBigDecimalArrayAndNeedToScaleProperlySetsScalingToFalse() {
    final boolean expectedValue = false;

    BigDecimalMatrix.setDefaultScaling(!expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is true when the needToScale argument is false").
        isEqualTo(expectedValue);
  }

  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][], false) creates a matrix whose neetToScale flag is set to false ")
  void testConstructorWithBigDecimalArrayAndNeedToScaleProperlySetsScalingToTrue() {
    final boolean expectedValue = true;

    BigDecimalMatrix.setDefaultScaling(!expectedValue);
    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_DATA, expectedValue);
    final boolean actualValue = newMatrix.getScaling();

    assertThat(actualValue).
        withFailMessage("Newly created matrix's scaling is false when the needToScale argument is true").
        isEqualTo(expectedValue);
  }
  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][]) creates a matrix with data equal to the source matrix's data")
  void testConstructorWithBigDecimalArrayCreatesInstanceWithCorrectData() {
    final double[][] expectedData = convertToDoubles(SAMPLE_BIGDECIMAL_DATA);

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_DATA);
    final double[][] actualData = newMatrix.getDoubleData();

    assertThat(actualData).
        withFailMessage("Constructor BigDecimalMatrix(BigDecimal[][]) creates a matrix containing data not equal to the source data").
        isEqualTo(expectedData);
  }

  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][]) creates a matrix with the default precision")
  void testConstructorWithBigDecimalArrayMakesMatrixWithDefaultPrecision() {
    final int expectedPrecision = BigDecimalMatrix.getDefaultPrecision();

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_DATA);
    final int actualPrecision = newMatrix.getPrecision();

    assertThat(actualPrecision).
        withFailMessage("Newly created matrix has precision other than the default precision").
        isEqualTo(expectedPrecision);
  }

  @Test
  @DisplayName("BigDecimalMatrix(BigDecimal[][], int) creates a matrix with the speciafied precision")
  void testConstructorWithBigDecimalArrayAndPrecisionMakesMatrixWithSpecifiedPrecision() {
    final int expectedPrecision = BigDecimalMatrix.getDefaultPrecision() + 123;

    final BigDecimalMatrix newMatrix = new BigDecimalMatrix(SAMPLE_BIGDECIMAL_DATA, expectedPrecision);
    final int actualPrecision = newMatrix.getPrecision();

    assertThat(actualPrecision).
        withFailMessage("Newly created matrix has precision other than the default precision").
        isEqualTo(expectedPrecision);
  }

//  // public boolean getScaling() is actually already tested by the previous tests

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

    final BigDecimalMatrix matrix = makeAMatrixOfSize(expectedSize);
    final int actualSize = matrix.getSize();

    assertThat(actualSize).
        withFailMessage("getSize() returns wrong value: (%d) instead of expected (%d)", actualSize, expectedSize).
        isEqualTo(expectedSize);
  }

  /*#########################################################################################
   * Tests for the methods:
   *   public Number[][] getData()
   *  -- getData() returns the correct values
   *  -- getData() returns a copy of the internal data and not a reference to it
   #################################################################################*/

  @Test
  @DisplayName("getData() returns the correct values")
  void testGetDataReturnsCorrectValues() {
    final Number[][] expectedData = convertToBigDecimals(SAMPLE_MATRIX_DATA);
    final BigDecimalMatrix matrix = new BigDecimalMatrix(expectedData);

    final Number[][] actualData =  matrix.getData();

    assertThat(actualData).
        withFailMessage("getData() returns wrong values").
        isEqualTo(expectedData);
  }

  @Test
  @DisplayName("getData() returns a copy of the internal data and not a reference to it")
  void testGetDataReturnsACopyNotAReference() {
    final Number[][] expectedData = convertToBigDecimals(SAMPLE_MATRIX_DATA);
    final BigDecimalMatrix matrix = new BigDecimalMatrix(expectedData);

    final Number[][] tmpData = matrix.getData();
    tmpData[0][0] = ((BigDecimal)tmpData[0][0]).add(BigDecimal.ONE); // corrupt the copy
    final Number[][] actualData =  matrix.getData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(expectedData);

    final double[][] actualData =  matrix.getDoubleData();

    assertThat(actualData).
        withFailMessage("getData() returns wrong values").
        isEqualTo(expectedData);
  }

  @Test
  @DisplayName("getDoubleData() returns a copy of the internal data and not a reference to it")
  void testGetDoubleDataReturnsACopyNotAReference() {
    final double[][] expectedData = SAMPLE_MATRIX_DATA;
    final BigDecimalMatrix matrix = new BigDecimalMatrix(expectedData);

    final double[][] tmpData = matrix.getDoubleData();
    tmpData[0][0] +=  1; // corrupt the copy
    final double[][] actualData =  matrix.getDoubleData();

    assertThat(actualData).
        withFailMessage("getDoubleData() returns a reference to the internal data").
        isEqualTo(expectedData);
  }

  /*#########################################################################################
   * Tests for the methods:
   *   public Quadruple[][] getQuadrupleData()
   *  -- getQuadrupleData() returns the correct data
   *  -- getQuadrupleData() returns a copy of the internal data and not a reference to it
   #################################################################################*/

  @Test
  @DisplayName("getQuadrupleData() returns the correct data")
  void testGetQuadrupleDataReturnsCorrectData() {
    final Quadruple[][] expectedData = convertToQuadruples(SAMPLE_MATRIX_DATA);

    final BigDecimalMatrix matrix = new BigDecimalMatrix(expectedData);
    final Quadruple[][] actualData =  matrix.getQuadrupleData();

    assertThat(actualData).
        withFailMessage("getQuadrupleData() returns wrong values").
        isEqualTo(expectedData);
  }

  @Test
  @DisplayName("getQuadrupleData() returns a copy of the internal data and not a reference to it")
  void testGetQuadrupleDataReturnsACopyNotAReference() {
    final Quadruple[][] expectedData = SAMPLE_QUADRUPLE_DATA;
    final BigDecimalMatrix matrix = new BigDecimalMatrix(expectedData);

    final Quadruple[][] tmpData = matrix.getQuadrupleData();
    tmpData[0][0].add(1);                         // corrupt the copy
    final Quadruple[][] actualData =  matrix.getQuadrupleData();

    assertThat(actualData).
        withFailMessage("getQuadrupleData() returns a reference to the internal data").
        isEqualTo(expectedData);
  }

  /*#########################################################################################
   * Tests for the methods:
   *   public BigDecimal[][] getBigDecimalData()
   *  -- getBigDecimalData() returns the correct data
   #################################################################################*/

  @Test
  @DisplayName("getBigDecimalData() returns the correct data")
  void testGetBigDecimalDataReturnsCorrectData() {
    final BigDecimal[][] expectedData = convertToBigDecimals(SAMPLE_MATRIX_DATA);

    final BigDecimalMatrix matrix = new BigDecimalMatrix(SAMPLE_MATRIX_DATA);
    final BigDecimal[][] actualData =  matrix.getBigDecimalData();

    assertThat(actualData).
        withFailMessage("getBigDecimaleData() returns wrong values").
        isEqualTo(expectedData);
  }

  @Test
  @DisplayName("getBigDecimalData() returns a copy of the internal data and not a reference to it")
  void testGetBigDecimalDataReturnsACopyNotAReference() {
    final BigDecimal[][] expectedData = SAMPLE_BIGDECIMAL_DATA;
    final BigDecimalMatrix matrix = new BigDecimalMatrix(expectedData);

    final BigDecimal[][] tmpData = matrix.getBigDecimalData();
    tmpData[0][0] = tmpData[0][0].add(BigDecimal.ONE);                         // corrupt the copy
    final BigDecimal[][] actualData =  matrix.getBigDecimalData();

    assertThat(actualData).
        withFailMessage("getBigDecimalData() returns a reference to the internal data").
        isEqualTo(expectedData);
  }

  /*#########################################################################################
   * Tests for the method:
   *   public boolean equals(Object anotherOne)
   * Tested behavior:
   *  -- equals(Object) returns false if the argument is null
   *  -- equals(Object) with Matrix returns false if the argument is not a BigDecimalMatrix
   *  -- equals(Object) with Matrix returns false if the argument has a different size
   *  -- equals(Matrix) with Matrix returns false if the argument has a different value of the needToScale flag
   *  -- equals(Object) with Matrix returns false if the argument has a different precision
   *  -- equals(Matrix) with Matrix returns false if the argument has a different matrix data
   *  -- equals(Matrix) with Matrix returns true if both matrix data and needToScale are equal
   #################################################################################*/


  @Test
  @DisplayName("equals(Object) returns false if the argument is null")
  void testEqualsReturnsFalseIfArgumentIsNull() {
    final boolean expectedResult = false;

    final boolean actualResult = SAMPLE_BIGDECIMAL_MATRIX.equals(null);

    assertThat(actualResult).
        withFailMessage("equals() returns true when called with null argument").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("equals(Object) with Matrix returns false if the argument is not a BigDecimalMatrix")
  void testEqualsReturnsFalseIfArgumentIsNotBigDecimalMatrix() {
    final boolean expectedResult = false;

    @SuppressWarnings("unlikely-arg-type")
    final boolean actualResult = SAMPLE_BIGDECIMAL_MATRIX.equals(new QuadrupleMatrix(SAMPLE_MATRIX_DATA));

    assertThat(actualResult).
        withFailMessage("equals() returns true when called with argument of another type").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("equals(Object) with Matrix returns false if the argument has a different size")
  void testEqualsReturnsFalseIfArgumentHasAnotherSize() {
    final boolean expectedResult = false;

    final boolean actualResult = SAMPLE_BIGDECIMAL_MATRIX.equals(new BigDecimalMatrix(SAMPLE_TOO_LARGE_MATRIX_DATA));

    assertThat(actualResult).
        withFailMessage("equals() returns true when the argument has a different size").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("equals(Object) with Matrix returns false if the argument has a different value of the needToScale flag")
  void testEqualsReturnsFalseIfArgumentHasAnotherNeedToScale()  {
    final boolean expectedResult = false;

    final BigDecimalMatrix matrix1 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, false);
    final BigDecimalMatrix matrix2 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, true);
    final boolean actualResult = matrix1.equals(matrix2);

    assertThat(actualResult).
        withFailMessage("equals() returns true when called with matrix with another value of needToScale flag").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("equals(Object) with Matrix returns false if the argument has a different precision")
  void testEqualsReturnsFalseIfIfArgumentHasAnotherPrecision() {
    final boolean expectedResult = false;

    final BigDecimalMatrix matrix1 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, 20);
    final BigDecimalMatrix matrix2 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, 30);
    final boolean actualResult = matrix1.equals(matrix2);

    assertThat(actualResult).
        withFailMessage("equals() returns true when called with matrix with another value of precison").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("equals(Object) with Matrix returns false if the argument has a different matrix data")
  void testEqualsReturnsFalseIfArgumentHasAnotherData() {
    final boolean expectedResult = false;

    final BigDecimalMatrix matrix1 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA);
    final BigDecimalMatrix matrix2 = new BigDecimalMatrix(SAMPLE_ANOTHER_MATRIX_DATA);
    final boolean actualResult = matrix1.equals(matrix2);

    assertThat(actualResult).
        withFailMessage("equals() returns true when called with matrix with another data").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("equals(Object) with Matrix returns true if both matrix data and needToScale are equal")
  void testEqualsReturnsTrueIfArgumentsAreEqual() {
    final boolean expectedResult = true;

    final BigDecimalMatrix matrix1 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA);
    final BigDecimalMatrix matrix2 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA);
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
   *  -- hashCode() returns different values for matrices with different precision
   *  -- hashCode() returns different values for matrices with different matrix data
   *  -- hashCode() returns equal values for equal matrices
   #################################################################################*/

  @Test
  @DisplayName("hashCode() returns different values for matrices with different values of the needToScale flag")
  void testHashCodeReturnsDifferentResultsForDifferentNeedToScaleValues()  {
    final int result1 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, false).hashCode();
    final int result2 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, true).hashCode();

    assertThat(result1).
        withFailMessage("hashCode() returns equal values for matrices with different values of the needToScale flag").
        isNotEqualTo(result2);
  }

  @Test
  @DisplayName("hashCode() returns different values for matrices with different precision")
  void testHashCodeReturnsDifferentResultsForDifferentPrecisions()  {
    final int result1 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, 20).hashCode();
    final int result2 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA, 30).hashCode();

    assertThat(result1).
        withFailMessage("hashCode() returns equal values for matrices with different values of the needToScale flag").
        isNotEqualTo(result2);
  }

  @Test
  @DisplayName("hashCode() returns different values for matrices with different matrix data")
  void testHashCodeReturnsDifferentResultsForDifferentMatrixData()  {

    final int result1 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA).hashCode();
    final int result2 = new BigDecimalMatrix(SAMPLE_ANOTHER_MATRIX_DATA).hashCode();

    assertThat(result1).
        withFailMessage("hashCode() returns equal values for matrices with different data").
        isNotEqualTo(result2);
  }

  @Test
  @DisplayName("hashCode() returns equal values for equal matrices")
  void testHashCodeReturnsEqualResultsForEqualMatrices()  {

    final int result1 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA).hashCode();
    final int result2 = new BigDecimalMatrix(SAMPLE_MATRIX_DATA).hashCode();

    assertThat(result1).
        withFailMessage("hashCode() returns different values for equal matrices").
        isEqualTo(result2);
  }

  /* *******************************************************************************
  /***** Private methods ***********************************************************
  /*********************************************************************************/

  private static BigDecimalMatrix makeAMatrixOfSize(int expectedSize) {
    final double[][] data = new double[expectedSize][expectedSize];
    return new BigDecimalMatrix(data);
  }

}
