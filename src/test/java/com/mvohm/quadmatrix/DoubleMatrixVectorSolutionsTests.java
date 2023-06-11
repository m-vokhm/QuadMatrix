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

import java.util.Locale;
import java.util.Random;
import java.lang.IllegalArgumentException;
import java.math.BigDecimal;

import org.junit.jupiter.api.*;

import static org.junit.jupiter.api.Assertions.*;
import static org.assertj.core.api.Assertions.*;

import com.mvohm.quadruple.Quadruple;
import com.mvohm.quadmatrix.test.AuxMethods;
import com.mvohm.quadmatrix.test.MatrixData;
import com.mvohm.quadmatrix.test.OperationComparator;
import com.mvohm.quadmatrix.test.AuxMethods.ErrorSet;

import static com.mvohm.quadmatrix.test.AuxMethods.*;
import static com.mvohm.quadmatrix.test.MatrixDataSamples.*;

/**
 * DoubleMatrixVectorSolutionsTest.java
 *
 * Tests for {@link DoubleMatrix} methods for solving systems of the form <b>Ax = b</b>
 *
 * Tested methods:
 *   public Number[] solve(double[] vector)
 *   public Number[] solve(Number[] vector)
 *   public Number[] solveSPD(double[] vector)
 *   public Number[] solveSPD(Number[] vector)
 *   public Number[] solveAccurately(double[] vector)
 *   public Number[] solveAccurately(Number[] vector)
 *   public Number[] solveSPDAccurately(double[] vector)
 *   public Number[] solveSPDAccurately(Number[] vector)
 *   public Number[] getSolution()
 *   public double[] getDoubleSolution()
 *   public Quadruple[] getQuadrupleSolution()
 *   public BigDecimal[] getBigDecimalSolution()
 *   public String getErrorCode()
 * Tested behaviors:
 *   the 'needScaling' flag actually improves the accuracy of the solution
 *   Consequent calls to solve() and solveSPD() don't interfere with each other
 */

public class DoubleMatrixVectorSolutionsTests {
  /** Enables progress indication for time-consuming tests and printing values of some computational errors */
  private static final boolean  OUTPUT_ENABLED = true;
  private static final boolean  DETAILED_OUTPUT = true;

  private static final int      RAND_SEED     =
                                                   123;
//                                                 12345;
//                                             111222333;
//                                               3141527;
//                                              76543210;
//                                                123321;

  private static Random         random;

  private static final int      REFINEMENT_TEST_ITERATIONS      =  150;
  private static final int      LARGE_MATRIX_SIZE               =  600;

  private static final int      SCALING_TEST_ITERATIONS         = 1000;
  private static final int      SCALING_MATRIX_SIZE             =  100;
  private static final double   SCALING_RANGE                   =    1.0e100;


  /**
   * To estimate the effectiveness (refinement ratio) of iterative refinements used by {@link DoubleMatrix#solveAccurately(double[])}
   * and {@link DoubleMatrix#solveSPDAccurately(double[])} methods,
   * a set of large random matrices (data samples) is used. For every matrix, a corresponding system of linear equations is solved,
   * and the mean square error of the found solution <b>x'</b> (actually, a square root of the MSE) is calculated as
   * <pre>e = sqrt(sum<sub>N</sub>( (x0[i] - x'[i])^2 ) / N)</pre>,
   * where x0 is the exact solution, x' is the found solution, N is the solution length (equal to the matrix size).
   * The solution improvement ratio is estimated as the ratio of <code>e0 / e1</code>, where <code>e0</code> is
   * sqrt(MSE) of the simple solution, and <code>e1</code> is sqrt(MSE)of the refined solution.<br>
   * The refinement ratio greatly depends on the specific data and the solution method. A solution obtained using Cholesky decomposition,
   * due to the nature of the SPD matrices and the method accuracy, usually can't be refined as effectively as the solution
   * gained with the LU-decomposition, typical refinement ratio for large SPD matrices is about 5.0
   *
   * The size of the matrix also matters, the refinement of relatively large matrices (200 x 200 and greater)
   * usually significantly improves the accuracy of the solution, while the refinement of relatively small matrices (say, 20 x 20)
   * may occur useless, in some cases the refined solution may be even worse than the solutions obtained with the simple solution.
   *
   * Saying, for example, about matrices of size 500, for some matrices the improvement ratio may be as small as 3 or a little more,
   * while for some other matrices it may exceed 300, with average value of about 70.<br>
   * The following values, that are used for testing of the refinement, are suitable for dense matrices
   * with uniform distribution of values ranged from -1.0 to 1.0 of size 500. For testing refinement on other types and/or
   * sizes of matrices, these value should be changed appropriately.
   * @see #testSolveAccuratelyIsAccurate()
   */

  // With a very little probability, the improvement ratio for a specific matrix of size 500 may be less than 1.5.
  // In this case try lessen this value or change the RAND_SEED constant to generate another set of data samples.
  private static final double EXPECTED_MIN_IMPROVEMENT         =   1.50;
  private static final double EXPECTED_AVR_IMPROVEMENT         =  60.00;

  private static final double EXPECTED_MIN_SPD_IMPROVEMENT     =   2.50;
  private static final double EXPECTED_AVR_SPD_IMPROVEMENT     =   6.00;

  private static final double EXPECTED_MIN_SCALING_IMPROVEMENT =   0.30; // In very rare cases. scaling may worsen the accuracy
  private static final double EXPECTED_AVR_SCALING_IMPROVEMENT =  15.00; // 20.00;

  private static final double ERROR_TOLERANCE                  = 1e-15;

  private static final double[][]   sampleSolvableMatrixData                      = sampleSolvableMatrixData();
  private static final double[][]   sampleInconsistentMatrixData                  = sampleInconsistentMatrixData();
  private static final double[][]   sampleUnderdeterminedMatrixData               = sampleUnderdeterminedMatrixData();

  private static final double[][]   sampleSpdSolvableMatrixData                   = sampleSpdSolvableMatrixData();
  private static final double[][]   sampleNonPositivelyDefinedMatrixData          = sampleNonPositivelyDefinedMatrixData();
  private static final double[][]   sampleNonSymmetricMatrixData                  = sampleNonSymmetricDMatrixData();

  private static final double[]     sampleSolvableMatrixVector                    = sampleSolvableMatrixVector();
  private static final double[]     sampleSolvableMatrixAnotherVector             = sampleSolvableMatrixAnotherVector();
  private static final Quadruple[]  sampleSolvableMatrixVectorAsQuadruples        = convertToQuadruples(sampleSolvableMatrixVector);
  private static final BigDecimal[] sampleSolvableMatrixVectorAsBigDecimals       = convertToBigDecimals(sampleSolvableMatrixVector);

  private static final double[]     sampleSpdSolvableMatrixVector                 = sampleSpdSolvableMatrixVector();
  private static final double[]     sampleSpdSolvableMatrixAnotherVector          = sampleSpdSolvableMatrixAnotherVector();
  private static final Quadruple[]  sampleSpdSolvableMatrixVectorAsQuadruples     = sampleSpdSolvableMatrixVectorAsQuadruples();
  private static final BigDecimal[] sampleSpdSolvableMatrixVectorAsBigDecimals    = convertToBigDecimals(sampleSpdSolvableMatrixVector);

  private static final double[]     sampleTooLongMatrixVector                     = sampleTooLongMatrixVector();
  private static final Quadruple[]  sampleTooLongMatrixVectorAsQuadruples         = convertToQuadruples(sampleTooLongMatrixVector);
  private static final BigDecimal[] sampleTooLongMatrixVectorAsBigDecimals        = convertToBigDecimals(sampleTooLongMatrixVector);

  private static final double[]     sampleSpdTooLongMatrixVector                  = sampleSpdTooLongMatrixVector();
  private static final Quadruple[]  sampleSpdTooLongMatrixVectorAsQuadruples      = convertToQuadruples(sampleSpdTooLongMatrixVector);
  private static final BigDecimal[] sampleSpdTooLongMatrixVectorAsBigDecimals     = convertToBigDecimals(sampleSpdTooLongMatrixVector);

  private static final double[]     sampleSolvableMatrixSolution                  = sampleSolvableMatrixSolution();
  private static final double[]     sampleSolvableMatrixAnotherSolution           = sampleSolvableMatrixAnotherSolution();

  private static final double[]     sampleSpdSolvableMatrixSolution               = sampleSpdSolvableMatrixSolution();
  private static final double[]     sampleSpdSolvableMatrixAnotherSolution        = sampleSpdSolvableMatrixAnotherSolution();


  @BeforeAll
  static void setup() {
    Locale.setDefault(Locale.US);
    DoubleMatrix.setDefaultScaling(true);
  }

  @BeforeEach
  void initRandom() {
    random = new Random(RAND_SEED);
  }


  /*####################################################################################
  ## Tested method:  public Number[] solve(double[] vector)
  ## Behaviors to be tested:
  ## -- solve(double[]) with null argument throws NullPointerException
  ## -- solve(double[]) with wrong vector length throws IllegalArgumentException
  ## -- solve(double[]) with vector containing NaN throws IllegalArgumentException
  ## -- solve(double[]) with vector containing Infinity throws IllegalArgumentException
  ## -- solve(double[]) does not spoil the passed vector
  ## -- solve(double[]) with inconsistent data throws IllegalArgumentException
  ## -- solve(double[]) with underdetermined data throws IllegalArgumentException
  ## -- solve(double[]) returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solve(double[]) with null argument throws NullPointerException")
  void testSolveWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solve((double[])null),
        "Attempt to call solve(double[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(double[])' and 'is null'").
        contains("solve(double[])").contains("is null");
  }

  @Test
  @DisplayName("solve(double[]) with wrong vector length throws IllegalArgumentException")
  void testSolveWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleTooLongMatrixVector),
        "Attempt to solve with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(double[])' and 'is must have the same size'").
        contains("solve(double[])").contains("must have the same size");
  }

  @Test
  @DisplayName("solve(double[]) with vector containing NaN throws IllegalArgumentException")
  void testSolveWithVectorContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[] badVector = sampleSolvableMatrixVector.clone();
    badVector[1] = Double.NaN;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badVector),
        "Attempt to solve with vector containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(double[])' and 'must not contain NaN or Infinity'").
        contains("solve(double[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solve(double[]) with vector containing Infinity throws IllegalArgumentException")
  void testSolveWithVectorContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[] badVector = sampleSolvableMatrixVector.clone();
    badVector[1] = Double.POSITIVE_INFINITY;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badVector),
        "Attempt to solve with vector containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(double[])' and 'must not contain NaN or Infinity'").
        contains("solve(double[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solve(double[]) does not spoil the passed vector")
  void testSolveDoesntSpoilTheVector() {
    final double[] expectedVector = sampleSolvableMatrixVector;
    final double[] actualVector = sampleSolvableMatrixVector.clone();

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(actualVector);

    assertThat(actualVector).
        withFailMessage("solve(double[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  @Test
  @DisplayName("solve(double[]) with inconsistent data throws IllegalArgumentException")
  void testSolveWithInconsistentDataThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixVector),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(double[]) with underdetermined data throws IllegalArgumentException")
  void testSolveWithUnderdeterminedDataThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
                            () -> matrix.solve(sampleSolvableMatrixVector),
                            "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(double[]) returns a correct solution")
  void testSolveReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSolvableMatrixSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[] actualSolution = convertToDoubles(matrix.solve(sampleSolvableMatrixVector));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:  public Number[] solve(Number[] vector), with a Quadruple[] argument
  ## Behaviors to be tested:
  ## -- solve(Quadruple[]) with null argument throws NullPointerException
  ## -- solve(Quadruple[]) with wrong vector length throws IllegalArgumentException
  ## -- solve(Quadruple[]) with vector containing NaN throws IllegalArgumentException
  ## -- solve(Quadruple[]) with vector containing Infinity throws IllegalArgumentException
  ## -- solve(Quadruple[]) with vector containing null throws IllegalArgumentException
  ## -- solve(Quadruple[]) does not spoil the passed vector
  ## -- solve(Quadruple[]) with inconsistent data throws IllegalArgumentException
  ## -- solve(Quadruple[]) with underdetermined data throws IllegalArgumentException
  ## -- solve(Quadruple[]) returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solve(Quadruple[]) with null argument throws NullPointerException")
  void testSolveWithQuadruplesWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solve((Quadruple[])null),
        "Attempt to call solve(Quadruple[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Number[])' and 'is null'").
        contains("solve(number[])").contains("is null");
  }

  @Test
  @DisplayName("solve(Quadruple[]) with wrong vector length throws IllegalArgumentException")
  void testSolveWithQuadruplesWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleTooLongMatrixVectorAsQuadruples),
        "Attempt to solve with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Number[])' and 'must have the same size'").
        contains("solve(number[])").contains("must have the same size");
  }

  @Test
  @DisplayName("solve(Quadruple[]) with vector containing NaN throws IllegalArgumentException")
  void testSolveWithQuadruplesWithVectorContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = Quadruple.nan();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badVector),
        "Attempt to solve with vector containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(number[])' and 'must not contain NaN or Infinity'").
        contains("solve(number[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solve(Quadruple[]) with vector containing Infinity throws IllegalArgumentException")
  void testSolveWithQuadruplesWithVectorContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = Quadruple.positiveInfinity();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badVector),
        "Attempt to solve with vector containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(number[])' and 'must not contain NaN or Infinity'").
        contains("solve(number[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solve(Quadruple[]) with vector containing null throws IllegalArgumentException")
  void testSolveWithQuadruplesWithVectorContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badVector),
        "Attempt to solve with vector containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(number[])' and 'must not contain null'").
        contains("solve(number[])").contains("must not contain null");
  }

  @Test
  @DisplayName("solve(Quadruple[]) does not spoil the passed vector")
  void testSolveWithQuadruplesDoesntSpoilTheVector() {
    final Quadruple[] expectedVector = sampleSolvableMatrixVectorAsQuadruples;
    final Quadruple[] actualVector = deepCopyOf(expectedVector);

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(actualVector);

    assertThat(actualVector).
        withFailMessage("solve(Quadruple[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  @Test
  @DisplayName("solve(Quadruple[]) with inconsistent data throws IllegalArgumentException")
  void testSolveWithQuadruplesWithInconsistentDataThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixVectorAsQuadruples),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(Quadruple[]) with underdetermined data throws IllegalArgumentException")
  void testSolveWithQuadruplesWithUnderdeterminedDataThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixVectorAsQuadruples),
        "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(Quadruple[]) returns a correct solution")
  void testSolveWithQuadruplesReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSolvableMatrixSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[] actualSolution = convertToDoubles(matrix.solve(sampleSolvableMatrixVectorAsQuadruples));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:  public Number[] solve(Number[] vector), with a BigDecimal[] argument
  ## Behaviors to be tested:
  ## -- solve(BigDecimal[]) with null argument throws NullPointerException
  ## -- solve(BigDecimal[]) with wrong vector length throws IllegalArgumentException
  ## -- solve(BigDecimal[]) with vector containing null throws IllegalArgumentException
  ## -- solve(BigDecimal[]) does not spoil the passed vector
  ## -- solve(BigDecimal[]) with inconsistent data throws IllegalArgumentException
  ## -- solve(BigDecimal[]) with underdetermined data throws IllegalArgumentException
  ## -- solve(BigDecimal[]) returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solve(BigDecimal[]) with null argument throws NullPointerException")
  void testSolveWithBigDecimalsWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solve((BigDecimal[])null),
        "Attempt to call solve(BigDecimal[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Number[])' and 'is null'").
        contains("solve(number[])").contains("is null");
  }

  @Test
  @DisplayName("solve(BigDecimal[]) with wrong vector length throws IllegalArgumentException")
  void testSolveWithBigDecimalsWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleTooLongMatrixVectorAsBigDecimals),
        "Attempt to solve with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Number[])' and 'must have the same size'").
        contains("solve(number[])").contains("must have the same size");
  }

  @Test
  @DisplayName("solve(BigDecimal[]) with vector containing null throws IllegalArgumentException")
  void testSolveWithBigDecimalsWithVectorContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final BigDecimal[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsBigDecimals);
    badVector[1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badVector),
        "Attempt to solve with vector containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(number[])' and 'must not contain null'").
        contains("solve(number[])").contains("must not contain null");
  }

  @Test
  @DisplayName("solve(BigDecimal[]) does not spoil the passed vector")
  void testSolveWithBigDecimalsDoesntSpoilTheVector() {
    final BigDecimal[] expectedVector = sampleSolvableMatrixVectorAsBigDecimals;
    final BigDecimal[] actualVector = expectedVector.clone();

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(actualVector);

    assertThat(actualVector).
        withFailMessage("solve(Quadruple[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  @Test
  @DisplayName("solve(BigDecimal[]) with inconsistent data throws IllegalArgumentException")
  void testSolveWithBigDecimalsWithInconsistentDataThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixVectorAsBigDecimals),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(BigDecimal[]) with underdetermined data throws IllegalArgumentException")
  void testSolveWithBigDecimalsWithUnderdeterminedDataThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixVectorAsBigDecimals),
        "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(BigDecimal[]) returns a correct solution")
  void testSolveWithBigDecimalsReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSolvableMatrixSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[] actualSolution = convertToDoubles(matrix.solve(sampleSolvableMatrixVectorAsBigDecimals));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:  public Number[] solveSPD(double[] vector)
  ## Behaviors to be tested:
  ## -- solveSPD(double[]) with null argument throws NullPointerException
  ## -- solveSPD(double[]) with wrong vector length throws IllegalArgumentException
  ## -- solveSPD(double[]) with vector containing NaN throws IllegalArgumentException
  ## -- solveSPD(double[]) with vector containing Infinity throws IllegalArgumentException
  ## -- solveSPD(double[]) does not spoil the passed vector
  ## -- solveSPD(double[]) with non-positively-defined matrix throws IllegalArgumentException
  ## -- solveSPD(double[]) with asymmetric matrix throws IllegalArgumentException
  ## -- solveSPD(double[]) returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solveSPD(double[]) with null argument throws NullPointerException")
  void testSolveSPDWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveSPD((double[])null),
        "Attempt to call solveSPD(double[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(double[]' and 'is null'").
        contains("solvespd(double[])").contains("is null");
  }

  @Test
  @DisplayName("solveSPD(double[]) with wrong vector length throws IllegalArgumentException")
  void testSolveSPDWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(sampleSpdTooLongMatrixVector),
        "Attempt to solve with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(double[]' and 'the same size'").
        contains("solvespd(double[])").contains("the same size");
  }

  @Test
  @DisplayName("solveSPD(double[]) with vector containing NaN throws IllegalArgumentException")
  void testSolveSPDWithVectorContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[] badVector = sampleSolvableMatrixVector.clone();
    badVector[1] = Double.NaN;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(badVector),
        "Attempt to solveSPD with vector containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(double[])' and 'must not contain NaN or Infinity'").
        contains("solvespd(double[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveSPD(double[]) with vector containing Infinity throws IllegalArgumentException")
  void testSolveSPDWithVectorContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[] badVector = sampleSolvableMatrixVector.clone();
    badVector[1] = Double.POSITIVE_INFINITY;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(badVector),
        "Attempt to solveSPD with vector containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(double[])' and 'must not contain NaN or Infinity'").
        contains("solvespd(double[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveSPD(double[]) does not spoil the passed vector")
  void testSolveSPDDoesntSpoilTheVector() {
    final double[] expectedVector = sampleSpdSolvableMatrixVector;
    final double[] actualVector = expectedVector.clone();

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    matrix.solveSPD(actualVector);

    assertThat(actualVector).
        withFailMessage("solveSPD(double[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  @Test
  @DisplayName("solveSPD(double[]) with non-positively-defined matrix throws IllegalArgumentException")
  void testSolveSPDWithNonPdMatrixThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleNonPositivelyDefinedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(sampleSpdSolvableMatrixVector),
        "Attempt to solve with non-positively-defined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'must be positively defined'").
        contains("must be positively defined");
  }

  @Test
  @DisplayName("solveSPD(double[]) with asymmetric matrix throws IllegalArgumentException")
  void testSolveSPDWithAsymmetricMatrixThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleNonSymmetricMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(sampleSpdSolvableMatrixVector),
        "Attempt to solve with asymmetric matrix must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'must be symmetric'").
        contains("must be symmetric");
  }

  @Test
  @DisplayName("solveSPD(double[]) returns a correct solution")
  void testSolveSPDReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSpdSolvableMatrixSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    final double[] actualSolution = convertToDoubles(matrix.solveSPD(sampleSpdSolvableMatrixVector));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:  public Number[] solveSPD(Number[] vector) with Quadruple[] argument
  ## Behaviors to be tested:
  ## -- solveSPD(Quadruple[]) with null argument throws NullPointerException
  ## -- solveSPD(Quadruple[]) with wrong vector length throws IllegalArgumentException
  ## -- solveSPD(Quadruple[]) with vector containing NaN throws IllegalArgumentException
  ## -- solveSPD(Quadruple[]) with vector containing Infinity throws IllegalArgumentException
  ## -- solveSPD(Quadruple[]) with vector containing null throws IllegalArgumentException
  ## -- solveSPD(Quadruple[]) does not spoil the passed vector
  ## -- solveSPD(Quadruple[]) with non-positively-defined matrix throws IllegalArgumentException
  ## -- solveSPD(Quadruple[]) with non-positively-defined matrix throws IllegalArgumentException
  ## -- solveSPD(Quadruple[]) with asymmetric matrix throws IllegalArgumentException
  ## -- solveSPD(Quadruple[]) returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solveSPD(Quadruple[]) with null argument throws NullPointerException")
  void testSolveSPDWithQuadruplesWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveSPD((Quadruple[])null),
        "Attempt to call solveSPD(Quadruple[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(Number[])' and 'is null'").
        contains("solvespd(number[])").contains("is null");
  }

  @Test
  @DisplayName("solveSPD(Quadruple[]) with wrong vector length throws IllegalArgumentException")
  void testSolveSPDWithQuadruplesWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(sampleSpdTooLongMatrixVectorAsQuadruples),
        "Attempt to solveSPD with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(Number[])' and 'must have the same size'").
        contains("solvespd(number[])").contains("must have the same size");
  }

  @Test
  @DisplayName("solveSPD(Quadruple[]) with vector containing NaN throws IllegalArgumentException")
  void testSolveSPDWithQuadruplesWithVectorContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = Quadruple.nan();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(badVector),
        "Attempt to solveSPD with vector containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(number[])' and 'must not contain NaN or Infinity'").
        contains("solvespd(number[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveSPD(Quadruple[]) with vector containing Infinity throws IllegalArgumentException")
  void testSolveSPDWithQuadruplesWithVectorContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = Quadruple.positiveInfinity();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(badVector),
        "Attempt to solveSPD with vector containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(number[])' and 'must not contain NaN or Infinity'").
        contains("solvespd(number[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveSPD(Quadruple[]) with vector containing null throws IllegalArgumentException")
  void testSolveSPDWithQuadruplesWithVectorContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(badVector),
        "Attempt to solve with vector containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(number[])' and 'must not contain null'").
        contains("solvespd(number[])").contains("must not contain null");
  }

  @Test
  @DisplayName("solveSPD(Quadruple[]) does not spoil the passed vector")
  void testSolveSPDWithQuadruplesDoesntSpoilTheVector() {
    final Quadruple[] expectedVector = sampleSpdSolvableMatrixVectorAsQuadruples;
    final Quadruple[] actualVector = deepCopyOf(expectedVector);

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    matrix.solveSPD(actualVector);

    assertThat(actualVector).
        withFailMessage("solveSPD(Quadruple[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  @Test
  @DisplayName("solveSPD(Quadruple[]) with non-positively-defined matrix throws IllegalArgumentException")
  void testSolveSPDWithQuadruplesWithNonPdMatrixThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleNonPositivelyDefinedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(sampleSpdSolvableMatrixVectorAsQuadruples),
        "Attempt to solve with non-positively-defined matrix must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'must be positively defined'").
        contains("must be positively defined");
  }

  @Test
  @DisplayName("solveSPD(Quadruple[]) with asymmetric matrix throws IllegalArgumentException")
  void testSolveSPDWithQuadruplesWithAssymmetricMatrixThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleNonSymmetricMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(sampleSpdSolvableMatrixVectorAsQuadruples),
        "Attempt to solve with asymmetric matrix must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'must be symmetric'").
        contains("must be symmetric");
  }

  @Test
  @DisplayName("solveSPD(Quadruple[]) returns a correct solution")
  void testSolveSPDWithQuadruplesReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSpdSolvableMatrixSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    final double[] actualSolution = convertToDoubles(matrix.solveSPD(sampleSpdSolvableMatrixVectorAsQuadruples));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:  public Number[] solveSPD(Number[] vector) with BigDecimal[] argument
  ## Behaviors to be tested:
  ## -- solveSPD(BigDecimal[]) with null argument throws NullPointerException
  ## -- solveSPD(BigDecimal[]) with wrong vector length throws IllegalArgumentException
  ## -- solveSPD(BigDecimal[]) with vector containing null throws IllegalArgumentException
  ## -- solveSPD(BigDecimal[]) does not spoil the passed vector
  ## -- solveSPD(BigDecimal[]) with non-positively-defined matrix throws IllegalArgumentException
  ## -- solveSPD(BigDecimal[]) with asymmetric matrix throws IllegalArgumentException
  ## -- solveSPD(BigDecimal[]) returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solveSPD(BigDecimal[]) with null argument throws NullPointerException")
  void testSolveSPDBigDecimalWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveSPD((BigDecimal[])null),
        "Attempt to call solveSPD(BigDecimal[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(Number[])' and 'is null'").
        contains("solvespd(number[])").contains("is null");
  }

  @Test
  @DisplayName("solveSPD(BigDecimal[]) with wrong vector length throws IllegalArgumentException")
  void testSolveSPDWithBigDecimalsWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(sampleSpdTooLongMatrixVectorAsBigDecimals),
        "Attempt to solve with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(Number[])' and 'must have the same size'").
        contains("solvespd(number[])").contains("must have the same size");
  }

  @Test
  @DisplayName("solveSPD(BigDecimal[]) with vector containing null throws IllegalArgumentException")
  void testSolveSPDWithBigDecimalsWithVectorContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final BigDecimal[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsBigDecimals);
    badVector[1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(badVector),
        "Attempt to solve with vector containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPD(number[])' and 'must not contain null'").
        contains("solvespd(number[])").contains("must not contain null");
  }

  @Test
  @DisplayName("solveSPD(BigDecimal[]) does not spoil the passed vector")
  void testSolveSPDWithBigDecimalsDoesntSpoilTheVector() {
    final BigDecimal[] expectedVector = sampleSpdSolvableMatrixVectorAsBigDecimals;
    final BigDecimal[] actualVector = expectedVector.clone();

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    matrix.solveSPD(actualVector);

    assertThat(actualVector).
        withFailMessage("solveSPD(BigDecimal[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  @Test
  @DisplayName("solveSPD(BigDecimal[]) with non-positively-defined matrix throws IllegalArgumentException")
  void testSolveSPDWithBigDecimalsWithNonPdMatrixThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleNonPositivelyDefinedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(sampleSpdSolvableMatrixVectorAsBigDecimals),
        "Attempt to solve with non-positively-defined matrix must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'must be positively defined'").
        contains("must be positively defined");
  }

  @Test
  @DisplayName("solveSPD(BigDecimal[]) with asymmetric matrix throws IllegalArgumentException")
  void testSolveSPDWithBigDecimalWithAsymmetricMatrixThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleNonSymmetricMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPD(sampleSpdSolvableMatrixVectorAsBigDecimals),
        "Attempt to solve with asymmetric matrix must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'must be symmetric'").
        contains("must be symmetric");
  }

  @Test
  @DisplayName("solveSPD(BigDecimal[]) returns a correct solution")
  void testSolveSPDWithBigDecimalsReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSpdSolvableMatrixSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    final double[] actualSolution = convertToDoubles(matrix.solveSPD(sampleSpdSolvableMatrixVectorAsBigDecimals));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:  public Number[] solveAccurately(double[] vector)
  ## Behaviors to be tested:
  ## -- solveAccurately(double[]) with null argument throws NullPointerException
  ## -- solveAccurately(double[]) with wrong vector length throws IllegalArgumentException
  ## -- solveAccurately(double[]) with vector containing NaN throws IllegalArgumentException
  ## -- solveAccurately(double[]) with vector containing Infinity throws IllegalArgumentException
  ## -- solveAccurately(double[]) does not spoil the passed vector
  ## -- solveAccurately(double[]) with inconsistent data throws IllegalArgumentException
  ## -- solveAccurately(double[]) with underdetermined data throws IllegalArgumentException
  ## -- solveAccurately(double[]) returns a correct solution
  ## -- solveAccurately(double[]) returns a much more accurate solution than solve()
  ######################################################################################*/

  @Test
  @DisplayName("solveAccurately(double[]) with null argument throws NullPointerException")
  void testSolveAccuratelyWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveAccurately((double[])null),
        "Attempt to call solveAccurately(double[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(double[])' and 'is null'").
        contains("solveaccurately(double[])").contains("is null");
  }

  @Test
  @DisplayName("solveAccurately(double[]) with wrong vector length throws IllegalArgumentException")
  void testSolveAccuratelyWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleTooLongMatrixVector),
        "Attempt to solve with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(double[])' and 'is must have the same size'").
        contains("solveaccurately(double[])").contains("must have the same size");
  }

  @Test
  @DisplayName("solveAccurately(double[]) with vector containing NaN throws IllegalArgumentException")
  void testSolveAccuratelyWithVectorContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[] badVector = sampleSolvableMatrixVector.clone();
    badVector[1] = Double.NaN;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badVector),
        "Attempt to call solveAccurately with vector containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(double[])' and 'must not contain NaN or Infinity'").
        contains("solveaccurately(double[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveAccurately(double[]) with vector containing Infinity throws IllegalArgumentException")
  void testSolveAccuratelyWithVectorContainingInfinityThrowsException() {

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[] badVector = sampleSolvableMatrixVector.clone();
    badVector[1] = Double.POSITIVE_INFINITY;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badVector),
        "Attempt to call solveAccurately with vector containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(double[])' and 'must not contain NaN or Infinity'").
        contains("solveaccurately(double[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveAccurately(double[]) does not spoil the passed vector")
  void testSolveAccuratelyDoesntSpoilTheVector() {
    final double[] expectedVector = sampleSolvableMatrixVector;
    final double[] actualVector = expectedVector.clone();

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solveAccurately(actualVector);

    assertThat(actualVector).
        withFailMessage("solveAccurately(double[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  @Test
  @DisplayName("solveAccurately(double[]) with inconsistent data throws IllegalArgumentException")
  void testSolveAccuratelyWithInconsistentDataThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixVector),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solveAccurately(double[]) with underdetermined data throws IllegalArgumentException")
  void testSolveAccuratelyWithUnderdeterminedDataThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixVector),
        "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solveAccurately(double[]) returns a correct solution")
  void testSolveAccuratelyReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSolvableMatrixSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[] actualSolution = convertToDoubles(matrix.solveAccurately(sampleSolvableMatrixVector));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /**
   * Generates a series of random matrices (data samples), solves the corresponding systems of linear equations,
   * and estimates the minimum and average ratio of
   * accuracy improvement achieved by the {@link DoubleMatrix#solveAccurately(double[])} compared with the
   * {@link DoubleMatrix#solve(double[])}. Asserts that the minimum improvement ratio is greater that {@link #EXPECTED_MIN_IMPROVEMENT}
   * and the average improvement ratio is greater than {@link #EXPECTED_AVR_IMPROVEMENT}
   * @see #EXPECTED_MIN_IMPROVEMENT
   */
  // This method may take up to a few minutes to execute. See REFINEMENT_TEST_ITERATIONS
  // Enable OUTPUT_ENABLED to see intermediate results
  // @Disabled // TO DO 2023-03-12 17:54:28 It's for speed while debugging.
  @Test
  @DisplayName("solveAccurately(double[]) returns a much more accurate solution than solve()")
  void testSolveAccuratelyIsAccurate() {
    printMethodName();
    final OperationComparator comparator = makeAccurateVsSimpleLuSolutionComparator();

    for (int i = 0; i < REFINEMENT_TEST_ITERATIONS; i++) {
      comparator.performOperations();
      say("  %3d: MaxErr ratio: %7.3f, MSE Ratio: %7.3f",
          i, comparator.getCurrentMaxErrRatio(), comparator.getCurrentMseRatio());
    }

    if (DETAILED_OUTPUT)
      comparator.showDetailedReport();

    comparator.assertImprovementAsExpected(EXPECTED_MIN_IMPROVEMENT, EXPECTED_AVR_IMPROVEMENT);
  }


  /*####################################################################################
  ## Tested method:  public Number[] solveAccurately(Number[] vector) with Quadruple[] argument
  ## Behaviors to be tested:
  ## -- solveAccurately(Quadruple[]) with null argument throws NullPointerException
  ## -- solveAccurately(Quadruple[]) with wrong vector length throws IllegalArgumentException
  ## -- solveAccurately(Quadruple[]) with vector containing NaN throws IllegalArgumentException
  ## -- solveAccurately(Quadruple[]) with vector containing Infinity throws IllegalArgumentException
  ## -- solveAccurately(Quadruple[]) with vector containing null throws IllegalArgumentException
  ## -- solveAccurately(Quadruple[]) does not spoil the passed vector
  ######################################################################################*/

  @Test
  @DisplayName("solveAccurately(Quadruple[]) with null argument throws NullPointerException")
  void testSolveAccuratelyWithQuadruplesWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveAccurately((Quadruple[])null),
        "Attempt to call solveAccurately(Quadruple[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Number[])' and 'is null'").
        contains("solveaccurately(number[])").contains("is null");
  }

  @Test
  @DisplayName("solveAccurately(Quadruple[]) with wrong vector length throws IllegalArgumentException")
  void testSolveAccuratelyWithQuadruplesWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleTooLongMatrixVectorAsQuadruples),
        "Attempt to solve with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Number[])' and 'must have the same size'").
        contains("solveaccurately(number[])").contains("must have the same size");
  }

  @Test
  @DisplayName("solveAccurately(Quadruple[]) with vector containing NaN throws IllegalArgumentException")
  void testSolveAccuratelyWithQuadruplesWithVectorContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = Quadruple.nan();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badVector),
        "Attempt to solveAccurately with vector containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(number[])' and 'must not contain NaN or Infinity'").
        contains("solveaccurately(number[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveAccurately(Quadruple[]) with vector containing Infinity throws IllegalArgumentException")
  void testSolveAccuratelyWithQuadruplesWithVectorContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = Quadruple.positiveInfinity();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badVector),
        "Attempt to solve with vector containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(number[])' and 'must not contain NaN or Infinity'").
        contains("solveaccurately(number[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveAccurately(Quadruple[]) with vector containing null throws IllegalArgumentException")
  void testSolveAccuratelyWithQuadruplesWithVectorContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badVector),
        "Attempt to solveAccurately with vector containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(number[])' and 'must not contain null'").
        contains("solveaccurately(number[])").contains("must not contain null");
  }

  @Test
  @DisplayName("solveAccurately(Quadruple[]) does not spoil the passed vector")
  void testSolveAccuratelyWithQuadruplesDoesntSpoilTheVector() {
    final Quadruple[] expectedVector = sampleSolvableMatrixVectorAsQuadruples;
    final Quadruple[] actualVector = deepCopyOf(expectedVector);

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solveAccurately(actualVector);

    assertThat(actualVector).
        withFailMessage("solveAccurately(Quadruple[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  /*####################################################################################
  ## Tested method:  public Number[] solveAccurately(Number[] vector) with BigDecimal[] argument
  ## Behaviors to be tested:
  ## -- solveAccurately(BigDecimal[]) with null argument throws NullPointerException
  ## -- solveAccurately(BigDecimal[]) with wrong vector length throws IllegalArgumentException
  ## -- solveAccurately(BigDecimal[]) with vector containing null throws IllegalArgumentException
  ## -- solveAccurately(BigDecimal[]) does not spoil the passed vector
  ######################################################################################*/

  @Test
  @DisplayName("solveAccurately(BigDecimal[]) with null argument throws NullPointerException")
  void testSolveAccuratelyWithBigDecimalsWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveAccurately((BigDecimal[])null),
        "Attempt to call solveAccurately(double[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Number[])' and 'is null'").
        contains("solveaccurately(number[])").contains("is null");
  }

  @Test
  @DisplayName("solveAccurately(BigDecimal[]) with wrong vector length throws IllegalArgumentException")
  void testSolveAccuratelyWithBigDecimalsWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleTooLongMatrixVectorAsBigDecimals),
        "Attempt to solve with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Number[])' and 'must have the same size'").
        contains("solveaccurately(number[])").contains("must have the same size");
  }

  @Test
  @DisplayName("solveAccurately(BigDecimal[]) with vector containing null throws IllegalArgumentException")
  void testSolveAccuratelyWithBigDecimalsWithVectorContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final BigDecimal[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsBigDecimals);
    badVector[1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badVector),
        "Attempt to solve with vector containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(number[])' and 'must not contain null'").
        contains("solveaccurately(number[])").contains("must not contain null");
  }

  @Test
  @DisplayName("solveAccurately(BigDecimal[]) does not spoil the passed vector")
  void testSolveAccuratelyWithBigDecimalsDoesntSpoilTheVector() {
    final BigDecimal[] expectedVector = sampleSolvableMatrixVectorAsBigDecimals;
    final BigDecimal[] actualVector = expectedVector.clone();

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solveAccurately(actualVector);

    assertThat(actualVector).
        withFailMessage("solveAccurately(double[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  /*####################################################################################
  ## Tested method:  public Number[] solveSPDAccurately(double[] vector)
  ## Behaviors to be tested:
  ## -- solveSPDAccurately(double[]) with null argument throws NullPointerException
  ## -- solveSPDAccurately(double[]) with wrong vector length throws IllegalArgumentException
  ## -- solveSPDAccurately(double[]) with vector containing NaN throws IllegalArgumentException
  ## -- solveSPDAccurately(double[]) with vector containing Infinity throws IllegalArgumentException
  ## -- solveSPDAccurately(double[]) does not spoil the passed vector
  ## -- solveSPDAccurately(double[]) with non-positively-defined matrix throws IllegalArgumentException
  ## -- solveSPDAccurately(double[]) with asymmetric matrix throws IllegalArgumentException
  ## -- solveSPDAccurately(double[]) returns a correct solution
  ## -- solveSPDAccurately(double[]) returns a much more accurate solution than solveSPD()
  ######################################################################################*/

  @Test
  @DisplayName("solveSPDAccurately(double[]) with null argument throws NullPointerException")
  void testSolveSPDAccuratelyWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveSPDAccurately((double[])null),
        "Attempt to call solveSPDAccurately(double[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(double[])' and 'is null'").
        contains("solvespdaccurately(double[])").contains("is null");
  }

  @Test
  @DisplayName("solveSPDAccurately(double[]) with wrong vector length throws IllegalArgumentException")
  void testSolveSPDAccuratelyWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPDAccurately(sampleSpdTooLongMatrixVector),
        "Attempt to solve with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(double[])' and 'is must have the same size'").
        contains("solvespdaccurately(double[])").contains("must have the same size");
  }

  @Test
  @DisplayName("solveSPDAccurately(double[]) with vector containing NaN throws IllegalArgumentException")
  void testSolveSPDAccuratelyWithVectorContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[] badVector = sampleSolvableMatrixVector.clone();
    badVector[1] = Double.NaN;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPDAccurately(badVector),
        "Attempt to solveSPDAccurately with vector containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(double[])' and 'must not contain NaN or Infinity'").
        contains("solvespdaccurately(double[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveSPDAccurately(double[]) with vector containing Infinity throws IllegalArgumentException")
  void testSolveSPDAccuratelyWithVectorContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[] badVector = sampleSolvableMatrixVector.clone();
    badVector[1] = Double.POSITIVE_INFINITY;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPDAccurately(badVector),
        "Attempt to solveSPDAccurately with vector containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(double[])' and 'must not contain NaN or Infinity'").
        contains("solvespdaccurately(double[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveSPDAccurately(double[]) does not spoil the passed vector")
  void testSolveSPDAccuratelyDoesntSpoilTheVector() {
    final double[] expectedVector = sampleSpdSolvableMatrixVector;
    final double[] actualVector = expectedVector.clone();

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    matrix.solveSPDAccurately(actualVector);

    assertThat(actualVector).
        withFailMessage("solveSPDAccurately(double[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  @Test
  @DisplayName("solveSPDAccurately(double[]) with non-positively-defined matrix throws IllegalArgumentException")
  void testSolveSPDAccuratelyWithNonPdMatrixThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleNonPositivelyDefinedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPDAccurately(sampleSpdSolvableMatrixVector),
        "Attempt to solve with non-positively-defined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'must be positively defined'").
        contains("must be positively defined");
  }

  @Test
  @DisplayName("solveSPDAccurately(double[]) with asymmetric matrix throws IllegalArgumentException")
  void testSolveSPDAccuratelyWithAsymmetricMatrixThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleNonSymmetricMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPDAccurately(sampleSpdSolvableMatrixVector),
        "Attempt to solve with asymmetric matrix must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'must be symmetric'").
        contains("must be symmetric");
  }

  @Test
  @DisplayName("solveSPDAccurately(double[]) returns a correct solution")
  void testSolveSPDAccuratelyReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSpdSolvableMatrixSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    final double[] actualSolution = convertToDoubles(matrix.solveSPDAccurately(sampleSpdSolvableMatrixVector));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /**
   * Generates a series of random matrices (data samples), solves the corresponding systems of linear equations,
   * and estimates the minimum and average ratio of
   * accuracy improvement achieved by the {@link DoubleMatrix#solveAccurately(double[])} compared with the
   * {@link DoubleMatrix#solve(double[])}. Asserts that the minimum improvement ratio is greater that {@link #EXPECTED_MIN_IMPROVEMENT}
   * and the average improvement ratio is greater than {@link #EXPECTED_AVR_IMPROVEMENT}
   * @see #EXPECTED_MIN_IMPROVEMENT
   */
  // This method may take up to a few minutes to execute. See REFINEMENT_TEST_ITERATIONS
  // Enable OUTPUT_ENABLED to see intermediate results
  // @Disabled // TO DO 2023-03-12 17:54:28 It's for speed while debugging
  @Test
  @DisplayName("solveSPDAccurately(double[]) returns a much more accurate solution than solveSPD()")
  void testSolveSPDAccuratelyIsAccurate() {
    printMethodName();
    final OperationComparator comparator = makeAccurateVsSimpleSPDSolutionComparator();

    for (int i = 0; i < REFINEMENT_TEST_ITERATIONS; i++) {
      comparator.performOperations();
      say("  %3d: MaxErr ratio: %7.3f, MSE Ratio: %7.3f",
          i, comparator.getCurrentMaxErrRatio(), comparator.getCurrentMseRatio());
    }

    if (DETAILED_OUTPUT)
      comparator.showDetailedReport();

    comparator.assertImprovementAsExpected(EXPECTED_MIN_SPD_IMPROVEMENT, EXPECTED_AVR_SPD_IMPROVEMENT);
  }

  /*####################################################################################
  ## Tested method:  public Number[] solveSPDAccurately(Number[] vector) with Quadruple[] argument
  ## Behaviors to be tested:
  ## -- solveSPDAccurately(Quadruple[]) with null argument throws NullPointerException
  ## -- solveSPDAccurately(Quadruple[]) with wrong vector length throws IllegalArgumentException
  ## -- solveSPDAccurately(Quadruple[]) with vector containing NaN throws IllegalArgumentException
  ## -- solveSPDAccurately(Quadruple[]) with vector containing Infinity throws IllegalArgumentException
  ## -- solveSPDAccurately(Quadruple[]) with vector containing null throws IllegalArgumentException
  ## -- solveSPDAccurately(Quadruple[]) does not spoil the passed vector
  ######################################################################################*/

  @Test
  @DisplayName("solveSPDAccurately(Quadruple[]) with null argument throws NullPointerException")
  void testSolveSPDAccuratelyWithQuadruplesWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveSPDAccurately((Quadruple[])null),
        "Attempt to call solveSPDAccurately(Quadruple[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(Number[])' and 'is null'").
        contains("solvespdaccurately(number[])").contains("is null");
  }

  @Test
  @DisplayName("solveSPDAccurately(Quadruple[]) with wrong vector length throws IllegalArgumentException")
  void testSolveSPDAccuratelyWithQuadruplesWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPDAccurately(sampleTooLongMatrixVectorAsQuadruples),
        "Attempt to solve with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(Number[])' and 'must have the same size'").
        contains("solvespdaccurately(number[])").contains("must have the same size");
  }


  @Test
  @DisplayName("solveSPDAccurately(Quadruple[]) with vector containing NaN throws IllegalArgumentException")
  void testSolveSPDAccuratelyWithQuadruplesWithVectorContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = Quadruple.nan();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPDAccurately(badVector),
        "Attempt to call solveSPDAccurately() with vector containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(number[])' and 'must not contain NaN or Infinity'").
        contains("solvespdaccurately(number[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveSPDAccurately(Quadruple[]) with vector containing Infinity throws IllegalArgumentException")
  void testSolveSPDAccuratelyWithQuadruplesWithVectorContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = Quadruple.positiveInfinity();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPDAccurately(badVector),
        "Attempt to call solveSPDAccurately() with vector containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(number[])' and 'must not contain NaN or Infinity'").
        contains("solvespdaccurately(number[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveSPDAccurately(Quadruple[]) with vector containing null throws IllegalArgumentException")
  void testSolveSPDAccuratelyWithQuadruplesWithVectorContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsQuadruples);
    badVector[1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPDAccurately(badVector),
        "Attempt to solveAccurately with vector containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(number[])' and 'must not contain null'").
        contains("solvespdaccurately(number[])").contains("must not contain null");
  }

  @Test
  @DisplayName("solveSPDAccurately(Quadruple[]) does not spoil the passed vector")
  void testSolveSPDAccuratelyWithQuadruplesDoesntSpoilTheVector() {
    final Quadruple[] expectedVector = sampleSolvableMatrixVectorAsQuadruples;
    final Quadruple[] actualVector = deepCopyOf(expectedVector);

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    matrix.solveSPDAccurately(actualVector);

    assertThat(actualVector).
        withFailMessage("solveSPDAccurately(Quadruple[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  /*####################################################################################
  ## Tested method:  public Number[] solveSPDAccurately(Number[] vector) with BigDecimal[] argument
  ## Behaviors to be tested:
  ## -- solveSPDAccurately(BigDecimal[]) with null argument throws NullPointerException
  ## -- solveSPDAccurately(BigDecimal[]) with wrong vector length throws IllegalArgumentException
  ## -- solveSPDAccurately(BigDecimal[]) with vector containing null throws IllegalArgumentException
  ## -- solveSPDAccurately(BigDecimal[]) does not spoil the passed vector
  ######################################################################################*/

  @Test
  @DisplayName("solveSPDAccurately(BigDecimal[]) with null argument throws NullPointerException")
  void testSolveSPDAccuratelyWithBigDecimalsWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveSPDAccurately((BigDecimal[])null),
        "Attempt to call solveSPDAccurately(double[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(Number[])' and 'is null'").
        contains("solvespdaccurately(number[])").contains("is null");
  }

  @Test
  @DisplayName("solveSPDAccurately(BigDecimal[]) with wrong vector length throws IllegalArgumentException")
  void testSolveSPDAccuratelyWithBigDecimalsWithWrongVectorLengthThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPDAccurately(sampleTooLongMatrixVectorAsBigDecimals),
        "Attempt to solve with wrong vector length must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(Number[])' and 'must have the same size'").
        contains("solvespdaccurately(number[])").contains("must have the same size");
  }

  @Test
  @DisplayName("solveSPDAccurately(BigDecimal[]) with vector containing null throws IllegalArgumentException")
  void testSolveSPDAccuratelyWithBigDecimalsWithVectorContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final BigDecimal[] badVector = deepCopyOf(sampleSolvableMatrixVectorAsBigDecimals);
    badVector[1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveSPDAccurately(badVector),
        "Attempt to solve with vector containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveSPDAccurately(number[])' and 'must not contain null'").
        contains("solvespdaccurately(number[])").contains("must not contain null");
  }

  @Test
  @DisplayName("solveSPDAccurately(BigDecimal[]) does not spoil the passed vector")
  void testSolveSPDAccuratelyWithBigDecimalsDoesntSpoilTheVector() {
    final BigDecimal[] expectedVector = sampleSolvableMatrixVectorAsBigDecimals;
    final BigDecimal[] actualVector = expectedVector.clone();

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    matrix.solveSPDAccurately(actualVector);

    assertThat(actualVector).
        withFailMessage("solveSPDAccurately(double[]) spoiled the passed vector").
        isEqualTo(expectedVector);
  }

  /*####################################################################################
  ## Tested method:  public Number[] getSolution()
  ## Behaviors to be tested:
  ## -- getSolution() returns null if no solution was found so far
  ## -- getSolution() returns null if an attempt to solve an equation was unsuccessful
  ## -- getSolution() after successful solutions returns data equal to those returned by solve()
  ## -- getSolution() returns a copy of the internally stored data, not a reference of it
  ## -- getSolution() returns the last found solution
  ######################################################################################*/

  @Test
  @DisplayName("getSolution() returns null if no solution was found so far")
  void testGetSolutionReturnsNullIfNoSolutionWasFound() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final Number[] actualSolution = matrix.getSolution();

    assertThat(actualSolution).
        withFailMessage("getSolution() returns non-null when no solutions for the matrix was found").
        isNull();
  }

  @Test
  @DisplayName("getSolution() returns null if an attempt to solve an equation was unsuccessful")
  void testGetSolutionReturnsNullIfSolutionWasUnsuccessful() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);
    try {
      matrix.solve(sampleSolvableMatrixVector);
    } catch (final Exception x) {}
    final Number[] actualSolution = matrix.getSolution();

    assertThat(actualSolution).
        withFailMessage("getSolution() returns non-null after solution has failed").
        isNull();
  }

  @Test
  @DisplayName("getSolution() after successful solutions returns data equal to those returned by solve()")
  void testGetSolutionReturnsCorrectData() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final Number[] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherVector);

    final Number[] actualSolution = matrix.getSolution();

    assertThat(actualSolution).
        withFailMessage("getSolution() returns data that differ from what is returned by solve()").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getSolution() returns a copy of the internally stored data, not a reference of it")
  void testGetSolutionReturnsACopyOfInternalData() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final Number[] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherVector);

    final Number[] tmpSolution = matrix.getSolution();
    tmpSolution[0] = Double.valueOf(12345);
    final Number[] actualSolution = matrix.getSolution();

    assertThat(actualSolution).
        withFailMessage("getSolution() returns a reference of the internally stored data").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getSolution() returns the last found solution")
  void testGetSolutionReturnsTheLastSolution() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixVector);
    final Number[] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherVector);
    final Number[] actualSolution = matrix.getSolution();

    assertThat(actualSolution).
        withFailMessage("getSolution() returns data that differ from the last found solution").
        isEqualTo(expectedSolution);
  }

  /*####################################################################################
  ## Tested method:  public double[] getDoubleSolution()
  ## Behaviors to be tested:
  ## -- getDoubleSolution() returns null if no solution was found so far
  ## -- getDoubleSolution() returns null if an attempt to solve an equation was unsuccessful
  ## -- getDoubleSolution() after successful solutions returns data equal to those returned by solve()
  ## -- getDoubleSolution() returns a copy of the internally stored data, not a reference of it
  ######################################################################################*/

  @Test
  @DisplayName("getDoubleSolution() returns null if no solution was found so far")
  void testGetDoubleSolutionReturnsNullIfNoSolutionWasFound() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[] actualSolution = matrix.getDoubleSolution();

    assertThat(actualSolution).
        withFailMessage("getDoubleSolution() returns non-null when no solutions for the matrix was found").
        isNull();
  }

  @Test
  @DisplayName("getDoubleSolution() returns null if an attempt to solve an equation was unsuccessful")
  void testGetDoubleSolutionReturnsNullIfSolutionWasUnsuccessful() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    try {
      matrix.solve(sampleSolvableMatrixVector);
    } catch (final Exception x) {}
    final double[] actualSolution = matrix.getDoubleSolution();

    assertThat(actualSolution).
        withFailMessage("getDoubleSolution() returns non-null after solution has failed").
        isNull();
  }

  @Test
  @DisplayName("getDoubleSolution() after successful solutions returns data equal to those returned by solve()")
  void testGetDoubleSolutionReturnsCorrectData() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[] expectedSolution = convertToDoubles(matrix.solve(sampleSolvableMatrixAnotherVector));

    final double[] actualSolution = matrix.getDoubleSolution();

    assertThat(actualSolution).
        withFailMessage("getSolution() returns data that differ from what is returned by solve()").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getDoubleSolution() returns a copy of the internally stored data, not a reference of it")
  void testGetDoubleSolutionReturnsACopyOfInternalData() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixAnotherVector);
    final double[] expectedSolution = convertToDoubles(matrix.solve(sampleSolvableMatrixAnotherVector));

    final double[] tmpSolution = matrix.getDoubleSolution();
    tmpSolution[0] += 12345;
    final double[] actualSolution = matrix.getDoubleSolution();

    assertThat(actualSolution).
        withFailMessage("getSolution() returns a reference of the internally stored data").
        isEqualTo(expectedSolution);
  }

  /*####################################################################################
  ## Tested method:  public Quadruple[] getQuadrupleSolution()
  ## Behaviors to be tested:
  ## -- getQuadrupleSolution() returns null if no solution was found so far
  ## -- getQuadrupleSolution() returns null if an attempt to solve an equation was unsuccessful
  ## -- getQuadrupleSolution() after successful solutions returns data equal to those returned by solve()
  ## -- getQuadrupleSolution() returns a copy of the internally stored data, not a reference of it
  ######################################################################################*/

  @Test
  @DisplayName("getQuadrupleSolution() returns null if no solution was found so far")
  void testGetQuadrupleSolutionReturnsNullIfNoSolutionWasFound() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();

    assertThat(actualSolution).
        withFailMessage("getQuadrupleSolution() returns non-null when no solutions for the matrix was found").
        isNull();
  }

  @Test
  @DisplayName("getQuadrupleSolution() returns null if an attempt to solve an equation was unsuccessful")
  void testGetQuadrupleSolutionReturnsNullIfSolutionWasUnsuccessful() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    try {
      matrix.solve(sampleSolvableMatrixVector);
    } catch (final Exception x) {}
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();

    assertThat(actualSolution).
        withFailMessage("getQuadrupleSolution() returns non-null after solution has failed").
        isNull();
  }

  @Test
  @DisplayName("getQuadrupleSolution() after successful solutions returns data equal to those returned by solve()")
  void testGetQuadrupleSolutionReturnsCorrectData() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final Quadruple[] expectedSolution = convertToQuadruples(matrix.solve(sampleSolvableMatrixAnotherVector));

    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();

    assertThat(actualSolution).
        withFailMessage("getQuadrupleSolution() returns data that differ from what is returned by solve()").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getQuadrupleSolution() returns a copy of the internally stored data, not a reference of it")
  void testGetQuadrupleSolutionReturnsACopyOfInternalData() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixAnotherVector);
    final Quadruple[] expectedSolution = convertToQuadruples(matrix.solve(sampleSolvableMatrixAnotherVector));

    final Quadruple[] tmpSolution = matrix.getQuadrupleSolution();
    tmpSolution[0] = tmpSolution[0].add(12345);
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();

    assertThat(actualSolution).
        withFailMessage("getQuadrupleSolution() returns a reference of the internally stored data").
        isEqualTo(expectedSolution);
  }

  /*####################################################################################
  ## Tested method:  public BigDecimal[] getBigDecimalSolution()
  ## Behaviors to be tested:
  ## -- getBigDecimalSolution() returns null if no solution was found so far
  ## -- getBigDecimalSolution() returns null if an attempt to solve an equation was unsuccessful
  ## -- getBigDecimalSolution() after successful solutions returns data equal to those returned by solve()
  ## -- getBigDecimalSolution() returns a copy of the internally stored data, not a reference of it
  ######################################################################################*/

  @Test
  @DisplayName("getBigDecimalSolution() returns null if no solution was found so far")
  void testGetBigDecimalSolutionReturnsNullIfNoSolutionWasFound() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();

    assertThat(actualSolution).
        withFailMessage("getBigDecimalSolution() returns non-null when no solutions for the matrix was found").
        isNull();
  }

  @Test
  @DisplayName("getBigDecimalSolution() returns null if an attempt to solve an equation was unsuccessful")
  void testGetBigDecimalSolutionReturnsNullIfSolutionWasUnsuccessful() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    try {
      matrix.solve(sampleSolvableMatrixVector);
    } catch (final Exception x) {}
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();

    assertThat(actualSolution).
        withFailMessage("getBigDecimalSolution() returns non-null after solution has failed").
        isNull();
  }

  @Test
  @DisplayName("getBigDecimalSolution() after successful solutions returns data equal to those returned by solve()")
  void testGetBigDecimalSolutionReturnsCorrectData() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final BigDecimal[] expectedSolution = convertToBigDecimals(matrix.solve(sampleSolvableMatrixAnotherVector));

    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();

    assertThat(actualSolution).
        withFailMessage("getBigDecimalSolution() returns data that differ from what is returned by solve()").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getBigDecimalSolution() returns a copy of the internally stored data, not a reference of it")
  void testGetBigDecimalSolutionReturnsACopyOfInternalData() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixAnotherVector);
    final BigDecimal[] expectedSolution = convertToBigDecimals(matrix.solve(sampleSolvableMatrixAnotherVector));

    final BigDecimal[] tmpSolution = matrix.getBigDecimalSolution();
    tmpSolution[0] = tmpSolution[0].add(new BigDecimal(12345));
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();

    assertThat(actualSolution).
        withFailMessage("getBigDecimalSolution() returns a reference of the internally stored data").
        isEqualTo(expectedSolution);
  }

  /*####################################################################################
  ## Tested method:  public String getErrorCode()
  ## Behaviors to be tested:
  ## -- If no calls to solve() was done yet, getErrorCode() returns 'OK'
  ## -- After a successful call to solve(), getErrorCode() returns 'OK'
  ## -- After an attempt to find a solution for an inconsistent matrix, getErrorCode() returns 'NON_INVERTIBLE'
  ## -- After an attempt to find a solution for an underdetermined matrix, getErrorCode() returns 'NON_INVERTIBLE'
  ## -- After an attempt to find a SPD solution for an asymmetric matrix, getErrorCode() returns 'ASYMMETRIC'
  ## -- After an attempt to find a SPD solution for a non-SPD matrix, getErrorCode() returns 'NON_SPD'
  ## -- After a successful solution following unsuccessful one, getErrorCode() returns 'OK'
  ######################################################################################*/

  @Test
  @DisplayName("If no calls to solve() was done yet, getErrorCode() returns 'OK'")
  void testGetErrorCodeReturnsOkIfSolveWasNotInvoked() {
    final String expectedResult = "OK";
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final String actualResult = matrix.getErrorCode();

    assertThat(actualResult).
        withFailMessage("getErrorCode() returns a wrong result when solve() was never invoked").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("After a successful call to solve(), getErrorCode() returns 'OK'")
  void testGetErrorCodeReturnsOkIfSolutionWasSuccessful() {
    final String expectedResult = "OK";
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    matrix.solve(sampleSolvableMatrixAnotherVector);
    final String actualResult = matrix.getErrorCode();

    assertThat(actualResult).
        withFailMessage("getErrorCode() returns a wrong result after successful call to solve()").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("After an attempt to find a solution for an inconsistent matrix, getErrorCode() returns 'NON_INVERTIBLE'")
  void testGetRightErrorCodeForInconsistentMatrix() {
    final String expectedResult = "NON_INVERTIBLE";
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    try {
      matrix.solve(sampleSolvableMatrixAnotherVector);
    } catch (final Exception x) {}
    final String actualResult = matrix.getErrorCode();

    assertThat(actualResult).
        withFailMessage("getErrorCode() returns a wrong result after solve() was called for an inconsistent matrix").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("After an attempt to find a solution for an underdetermined matrix, getErrorCode() returns 'NON_INVERTIBLE'")
  void testGetRightErrorCodeForUnderdeterminedMatrix() {
    final String expectedResult = "NON_INVERTIBLE";
    final DoubleMatrix matrix = new DoubleMatrix(sampleUnderdeterminedMatrixData);

    try {
      matrix.solve(sampleSolvableMatrixAnotherVector);
    } catch (final Exception x) {}
    final String actualResult = matrix.getErrorCode();

    assertThat(actualResult).
        withFailMessage("getErrorCode() returns a wrong result after solve() was called for an underdetermined matrix").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("After an attempt to find a SPD solution for an asymmetric matrix, getErrorCode() returns 'ASYMMETRIC'")
  void testGetRightErrorCodeForAsymmetricMatrix() {
    final String expectedResult = "ASYMMETRIC";
    final DoubleMatrix matrix = new DoubleMatrix(sampleUnderdeterminedMatrixData);

    try {
      matrix.solveSPD(sampleSolvableMatrixAnotherVector);
    } catch (final Exception x) { }
    final String actualResult = matrix.getErrorCode();

    assertThat(actualResult).
        withFailMessage("getErrorCode() returns a wrong result after solveSPD() was called for an asymmetric matrix").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("After an attempt to find a SPD solution for a non-SPD matrix, getErrorCode() returns 'NON_SPD'")
  void testGetRightErrorCodeForNonSpdMatrix() {
    final String expectedResult = "NON_SPD";
    final DoubleMatrix matrix = new DoubleMatrix(sampleNonPositivelyDefinedMatrixData);

    try {
      matrix.solveSPD(sampleSolvableMatrixAnotherVector);
    } catch (final Exception x) { }
    final String actualResult = matrix.getErrorCode();

    assertThat(actualResult).
        withFailMessage("getErrorCode() returns a wrong result after solveSPD() was called for a non-SPD matrix").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("After a successful solution following unsuccessful one, getErrorCode() returns 'OK'")
  void testGetRightErrorCodeForSuccessAfrerFailure() {
    final String expectedResult = "OK";
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    try {
      matrix.solveSPD(sampleSolvableMatrixAnotherVector);
    } catch (final Exception x) { }
    matrix.solve(sampleSolvableMatrixVector);
    final String actualResult = matrix.getErrorCode();

    assertThat(actualResult).
        withFailMessage("getErrorCode() returns a wrong result after solveSPD() was called for a non-SPD matrix").
        isEqualTo(expectedResult);
  }

  /*####################################################################################
  ## Tested behaviors: the 'needScaling' flag actually improves the accuracy of the solution
  ######################################################################################*/

  // This method may take up to a few minutes to execute. See REFINEMENT_TEST_ITERATIONS
  // Enable OUTPUT_ENABLED to see intermediate results
  // @Disabled // TO DO 2023-03-12 17:54:28 It's for speed while debugging.
  @Test
  @DisplayName("the 'needScaling' flag actually improves the accuracy of the solution")
  void testScalingImprovesAccuracy() {
    printMethodName();
    final OperationComparator comparator = makeScaledVsUnscaledVectorSolutionComparator();

    for (int i = 0; i < SCALING_TEST_ITERATIONS; i++) {
      comparator.performOperations();
      say("  %3d: MaxErr ratio: %7.3f, MSE Ratio: %7.3f",
          i, comparator.getCurrentMaxErrRatio(), comparator.getCurrentMseRatio());
    }

    if (DETAILED_OUTPUT)
      comparator.showDetailedReport();

    comparator.assertImprovementAsExpected(EXPECTED_MIN_SCALING_IMPROVEMENT, EXPECTED_AVR_SCALING_IMPROVEMENT);
  }

  /*####################################################################################
  ## Make sure that solutions don't change the state of the Matrix
  ## (i.e. consequent calls to solve() and solveSPD() don't interfere with each other)
  ## Behaviors to be tested:
  ## -- Another call to solve(double[]) with another vector also returns a correct solution
  ## -- Another call to solveSPD(double[]) with another vector also returns a correct solution
  ## -- A call to solve(double[]) after successful SPD-solution returns a correct solution
  ## -- A call to solve(double[]) after unsuccessful SPD-solution returns a correct solution
  ## -- A call to solveSPD() after solving with LU-decomposition returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("Another call to solve(double[]) with another vector also returns a correct solution")
  void testRepeatedsolveReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSolvableMatrixAnotherSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixVector);
    final double[] actualSolution = convertToDoubles(matrix.solve(sampleSolvableMatrixAnotherVector));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  @Test
  @DisplayName("Another call to solveSPD(double[]) with another vector also returns a correct solution")
  void testRepeatedSolveSPDReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSpdSolvableMatrixAnotherSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    matrix.solveSPD(sampleSpdSolvableMatrixVector);
    final double[] actualSolution = convertToDoubles(matrix.solveSPD(sampleSpdSolvableMatrixAnotherVector));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  @Test
  @DisplayName("A call to solve(double[]) after successful SPD-solution returns a correct solution")
  void testSolveAfterSolveSPDReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSpdSolvableMatrixSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    matrix.solveSPD(sampleSpdSolvableMatrixAnotherVector);
    final double[] actualSolution = convertToDoubles(matrix.solve(sampleSpdSolvableMatrixVector));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  @Test
  @DisplayName("A call to solve(double[]) after unsuccessful SPD-solution returns a correct solution")
  void testSolveAfterUnsucessfulSolveSPDReturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSolvableMatrixSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    try {
      matrix.solveSPD(sampleSpdSolvableMatrixAnotherVector); // No matter which vector we use for a non-spd-solvable system.
    } catch (final Exception x) { /* As expected. Ignore it */ }
    final double[] actualSolution = convertToDoubles(matrix.solve(sampleSolvableMatrixVector));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  @Test
  @DisplayName("A call to solveSPD() after solving with LU-decomposition returns a correct solution")
  void testSolveSpdAfterSolvingLUreturnsCorrectSolution() {
    printMethodName();
    final double[] expectedSolution = sampleSpdSolvableMatrixAnotherSolution;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSpdSolvableMatrixData);
    matrix.solve(sampleSpdSolvableMatrixVector);
    final double[] actualSolution = convertToDoubles(matrix.solveSPD(sampleSpdSolvableMatrixAnotherVector));
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /* ************************************************************************************
  /*** Private methods ******************************************************************
  /**************************************************************************************/

  private static void printMethodName() {
    String s = new Exception().getStackTrace()[1].toString();
    s = s.replaceFirst("\\(.*\\)", "():").replaceFirst(".*\\.", "");
    say(s);
  }

  @SuppressWarnings("unused")
  private static OperationComparator makeAccurateVsSimpleLuSolutionComparator() {
    return new OperationComparator(
        () -> MatrixData.makeDataSetForVectorSolutions(LARGE_MATRIX_SIZE, random),
        MatrixData::luSolutionWithScalingErrors,
        MatrixData::accurateLUSolutionWithScalingErrors
        );
  }


  @SuppressWarnings("unused")
  private static OperationComparator makeAccurateVsSimpleSPDSolutionComparator() {
    return new OperationComparator(
        () -> MatrixData.makeDataSetForSPDSolutions(LARGE_MATRIX_SIZE, random),
        MatrixData::spdSolutionErrors,
        MatrixData::accurateSPDSolutionErrors
        );
  }

  @SuppressWarnings("unused")
  private static OperationComparator makeScaledVsUnscaledVectorSolutionComparator() {
    return new OperationComparator(
        () -> MatrixData.makeLargeRangeDataSetForVectorSolutions(SCALING_MATRIX_SIZE, random, SCALING_RANGE),
        MatrixData::luSolutionWithoutScalingErrors,
        MatrixData::luSolutionWithScalingErrors
        );
  }

  @SuppressWarnings("unused")
  private static void say() {
    if (OUTPUT_ENABLED) AuxMethods.say();
  }

  @SuppressWarnings("unused")
  private static void say(Object arg) {
    if (OUTPUT_ENABLED) AuxMethods.say(arg);
  }

  @SuppressWarnings("unused")
  private static void say(String format, Object... args) {
    if (OUTPUT_ENABLED) AuxMethods.say(format, args);
  }

}
