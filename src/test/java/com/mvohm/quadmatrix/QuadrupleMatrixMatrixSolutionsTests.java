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

import java.util.Locale;
import java.util.Random;
import java.lang.IllegalArgumentException;
import java.math.BigDecimal;

import org.junit.jupiter.api.*;

import static org.junit.jupiter.api.Assertions.*;
import static org.assertj.core.api.Assertions.*;

import com.mvohm.quadmatrix.test.AuxMethods;
import com.mvohm.quadmatrix.test.MatrixData;
import com.mvohm.quadmatrix.test.OperationComparator;
import com.mvohm.quadmatrix.test.AuxMethods.ErrorSet;
import com.mvohm.quadruple.Quadruple;

import static com.mvohm.quadmatrix.test.AuxMethods.*;
import static com.mvohm.quadmatrix.test.MatrixDataSamples.*;

/**
 * Tests for {@link QuadrupleMatrix} methods for solving systems of the form <b>AX = B</b>
 * <br>
 * Tested methods:
 *   public Matrix solve(Matrix)
 *   public Matrix solve(double[][])
 *   public Matrix solve(Number[][])
 *   public Matrix solveAccurately(Matrix)
 *   public Matrix solveAccurately(double[][])
 *   public Matrix solveAccurately(Number[][])
 *   public Matrix getMatrixSolution()
 *   public Number[][] getNumberMatrixSolution()
 *   public double[][] getDoubleMatrixSolution()
 *   public Quadruple[][] getQuadrupleMatrixSolution()
 *   public BigDecimal[][] getBigDecimalMatrixSolution()
 *   public String getErrorCode()
 * Tested behaviors:
 *   the 'needScaling' flag actually improves the accuracy of the solution
 *   Consequent calls to solve() don't interfere with each other
 */
public class QuadrupleMatrixMatrixSolutionsTests {

   /** Enables progress indication for time-consuming tests and printing values of some computational errors */
  private static final boolean OUTPUT_ENABLED = true;
  private static final boolean DETAILED_OUTPUT = true;

  private static final int    RAND_SEED
                                      = 123;
//                                      = 34567890;
//                                      = 12345;
//                                      = 54321;

  private static Random       random;

  private static final double ERROR_TOLERANCE                 = 1e-36;

  /**
   * To estimate the effectiveness (refinement ratio) of iterative refinements used by
   * {@link QuadrupleMatrix#solveAccurately(double[][])} method,
   * a set of large random matrices (data samples) is used. For every matrix, a corresponding system of linear equations is solved,
   * and the mean square error of the found solution <b>X'</b> (actually, a square root of the MSE) is calculated as
   * <pre>e = sqrt(sum<sub>N^2</sub>( (X[i][j] - X'[i][j])^2 ) / N^2)</pre>,
   * where X is the exact solution, X' is the found solution, N is the solution length (equal to the matrix size).
   * The solution improvement ratio is estimated as the ratio of <code>e0 / e1</code>, where <code>e0</code> is
   * the error of the simple solution, and <code>e1</code> is the error of the refined solution.<br>
   * The refinement ratio greatly depends on the specific data and the solution method.
   * The size of the matrix also matters, the refinement of relatively large matrices (200 and more rows and columns)
   * usually significantly improves the accuracy of the solution, while the refinement of relatively small matrices (say, 20 x 20)
   * may occur useless, in some cases the refined solution may be even worse than the solutions obtained with the simple solution.
   *
   * Saying, for example, about matrices of size 500, for some matrices the improvement ratio may be as small as 3 or a little more,
   * while for some other matrices it may exceed 300, with average value of about 70.<br>
   * The following values, that are used for the refinement testing, are suitable for dense matrices
   * with uniform distribution of values ranged from -1.0 to 1.0 of size 160. For testing refinement on other types and/or
   * sizes of matrices, these value should be changed appropriately.
   * @see #testSolveAccuratelyIsAccurate()
   */
  private static final int              REFINEMENT_TEST_ITERATIONS                        =  30; // ~60 s on my machine
  private static final int              LARGE_MATRIX_SIZE                                 =  80; // On my machine, it provides approximately 2 secs per sample

  private static final double           EXPECTED_MIN_IMPROVEMENT                          =  14;
  private static final double           EXPECTED_AVR_IMPROVEMENT                          =  20;

  private static final double[]         sampleSpdSolvableMatrixAnotherVector              = sampleSpdSolvableMatrixAnotherVector();

  int i = 5;

  private static final double[][]       sampleSolvableMatrixData                          = sampleSolvableMatrixData();
  private static final Quadruple[][]    sampleSolvableMatrixDataAsQuadruples              = sampleSolvableMatrixDataAsQuadruples();

  private static final double[][]       sampleInconsistentMatrixData                      = sampleInconsistentMatrixData();
  private static final double[][]       sampleUnderdeterminedMatrixData                   = sampleUnderdeterminedMatrixData();

  private static final Quadruple[][]    sampleSpdSolvableMatrixDataAsQuadruples           = sampleSpdSolvableMatrixDataAsQuadruples();
  private static final Quadruple[][]    sampleSpdSolvableMatrixMatrixXAsQuadruples        = sampleSpdMatrixMatrixXAsQuadruples();
  private static final Quadruple[][]    sampleSpdSolvableMatrixMatrixBAsQuadruples        = sampleSpdMatrixMatrixBAsQuadruples();

  private static final Quadruple[][]    sampleSolvableMatrixMatrixXAsQuadruples           = sampleSolvableMatrixMatrixXAsQuadruples();

  private static final double[][]       sampleSolvableMatrixMatrixB                       = sampleSolvableMatrixMatrixB();
  private static final Quadruple[][]    sampleSolvableMatrixMatrixBAsQuadruples           = sampleSolvableMatrixMatrixBAsQuadruples();
  private static final BigDecimal[][]   sampleSolvableMatrixMatrixBAsBigDecimals          = convertToBigDecimals(sampleSolvableMatrixMatrixB());

  private static final Quadruple[][]    sampleSolvableMatrixAnotherMatrixXAsQuadruples    = sampleSolvableMatrixAnotherMatrixXAsQuadruples();

  private static final double[][]       sampleSolvableMatrixAnotherMatrixB                = sampleSolvableMatrixAnotherMatrixB();
  private static final Quadruple[][]    sampleSolvableMatrixAnotherMatrixBAsQuadruples    = sampleSolvableMatrixAnotherMatrixBAsQuadruples();

  private static final double[][]       sampleSolvableMatrixTooLargeMatrix                = sampleSolvableMatrixTooLargeMatrix();
  private static final double[][]       sampleSolvableMatrixNonSquareMatrix               = sampleSolvableMatrixNonSquareMatrix();

  private static final Quadruple[][]    sampleSolvableMatrixTooLargeMatrixAsQuadruples    = convertToQuadruples(sampleSolvableMatrixTooLargeMatrix);
  private static final Quadruple[][]    sampleSolvableMatrixNonSquareMatrixAsQuadruples   = convertToQuadruples(sampleSolvableMatrixNonSquareMatrix);

  private static final BigDecimal[][]   sampleSolvableMatrixTooLargeMatrixAsBigDecimals   = convertToBigDecimals(sampleSolvableMatrixTooLargeMatrix);
  private static final BigDecimal[][]   sampleSolvableMatrixNonSquareMatrixAsBigDecimals  = convertToBigDecimals(sampleSolvableMatrixNonSquareMatrix);


  @BeforeAll
  static void setup() {
    Locale.setDefault(Locale.US);
    QuadrupleMatrix.setDefaultScaling(true);
  }

  @BeforeEach
  void initRandom() {
    random = new Random(RAND_SEED);
  }

  /*####################################################################################
  ## Tested method:
  ##    public Matrix solve(Matrix matrixB)
  ## Behaviors to be tested:
  ## -- solve(Matrix) with null argument throws NullPointerException
  ## -- solve(Matrix) with wrong matrix size throws IllegalArgumentException
  ## -- solve(Matrix) with inconsistent data throws IllegalArgumentException
  ## -- solve(Matrix) with underdetermined data throws IllegalArgumentException
  ## -- solve(Matrix) returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solve(Matrix) with null argument throws NullPointerException")
  void testSolveWithMatrixWithNullArgumentTrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solve((Matrix)null),
        "Attempt to invoke solve() with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Matrix)' and 'is null'").
        contains("solve(matrix)").contains("is null");
  }

  @Test
  @DisplayName("solve(Matrix) with wrong matrix size throws IllegalArgumentException")
  void testSolveWithMatrixWithWrongMatrixSizeThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(new QuadrupleMatrix(sampleSolvableMatrixTooLargeMatrix)),
        "Attempt to solve with wrong matrix size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Matrix)' and 'must have the same size'").
        contains("solve(matrix)").contains("must have the same size");
  }

  @Test
  @DisplayName("solve(Matrix) with inconsistent data throws IllegalArgumentException")
  void testSolveWithMatrixWithInconsistentDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(new QuadrupleMatrix(sampleSolvableMatrixMatrixB)),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(Matrix) with underdetermined data throws IllegalArgumentException")
  void testSolveWithMatrixWithUnderdeterminedDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(new QuadrupleMatrix(sampleSolvableMatrixMatrixB)),
        "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(Matrix) returns a correct solution")
  void testSolveWithMatrixReturnsCorrectSolution() {
    printMethodName();
    final Quadruple[][] expectedSolution = sampleSolvableMatrixMatrixXAsQuadruples;
    final double expectedError = ERROR_TOLERANCE;

    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixDataAsQuadruples);
    final Quadruple[][] actualSolution = matrix.solve(new QuadrupleMatrix(sampleSolvableMatrixMatrixB)).getQuadrupleData();
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:
  ##    public Matrix solve(double[][] matrixB)
  ## Behaviors to be tested:
  ## -- solve(double[][]) with null argument throws NullPointerException
  ## -- solve(double[][]) with wrong matrix size throws IllegalArgumentException
  ## -- solve(double[][]) with a non-square argument throws IllegalArgumentException
  ## -- solve(double[][]) with matrix containing NaN throws IllegalArgumentException
  ## -- solve(double[][]) with matrix containing Infinity throws IllegalArgumentException
  ## -- solve(double[][]) with inconsistent data throws IllegalArgumentException
  ## -- solve(double[][]) with underdetermined data throws IllegalArgumentException
  ## -- solve(double[][]) does not spoil the passed matrixB array
  ## -- solve(double[][]) returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solve(double[][]) with null argument throws NullPointerException")
  void testSolveWithMatrixOfDoublesWithNullArgumentTrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solve((double[][])null),
        "Attempt to invoke solve() with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(double[][])' and 'is null'").
        contains("solve(double[][])").contains("is null");
  }

  @Test
  @DisplayName("solve(double[][]) with wrong matrix size throws IllegalArgumentException")
  void testSolveWithMatrixOfDoublesWithWrongMatrixSizeThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixTooLargeMatrix),
        "Attempt to solve with wrong matrix size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(double[][])' and 'must have the same size'").
        contains("solve(double[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("solve(double[][]) with a non-square argument throws IllegalArgumentException")
  void testSolveWithMatrixOfDoublesWithNonSquareMatrixThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixNonSquareMatrix),
        "Attempt to solve with non-square argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(double[][])' and 'must be square'").
        contains("solve(double[][])").contains("must be square");
  }

  @Test
  @DisplayName("solve(double[][]) with matrix containing NaN throws IllegalArgumentException")
  void testSolveWithMatrixOfDoublesContainingNaNThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final double[][] badMatrix = deepCopyOf(sampleSolvableMatrixData);
    badMatrix[1][1] = Double.NaN;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badMatrix),
        "Attempt to solve with matrix containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(double[][])' and 'must not contain NaN or Infinity'").
        contains("solve(double[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solve(double[][]) with matrix containing Infinity throws IllegalArgumentException")
  void testSolveWithMatrixOfDoublesContainingInfinityThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final double[][] badMatrix = deepCopyOf(sampleSolvableMatrixData);
    badMatrix[1][1] = Double.NEGATIVE_INFINITY;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badMatrix),
        "Attempt to solve with matrix containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(double[][])' and 'must not contain NaN or Infinity'").
        contains("solve(double[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solve(double[][]) with inconsistent data throws IllegalArgumentException")
  void testSolveWithMatrixOfDoublesWithInconsistentDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixMatrixB),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(double[][]) with underdetermined data throws IllegalArgumentException")
  void testSolveWithMatrixOfDoublesWithUnderdeterminedDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixMatrixB),
        "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(double[][]) does not spoil the passed matrixB array")
  void testSolveWithMatrixOfDoublesDoesntSpoilThePassedArray() {
    final double[][] expectedMatrixB = deepCopyOf(sampleSolvableMatrixMatrixB);
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    matrix.solve(sampleSolvableMatrixMatrixB);
    final double[][] actualMatrixB = deepCopyOf(sampleSolvableMatrixMatrixB);

    assertThat(actualMatrixB).
        withFailMessage("matrix.solve(double[][]) spoiled the passed matrix").
        isEqualTo(expectedMatrixB);
  }

  @Test
  @DisplayName("solve(double[][]) returns a correct solution")
  void testSolveWithMatrixOfDoublesReturnsCorrectSolution() {
    printMethodName();
    final Quadruple[][] expectedSolution = sampleSolvableMatrixMatrixXAsQuadruples;
    final double expectedError = ERROR_TOLERANCE;
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixDataAsQuadruples);

    final Quadruple[][] actualSolution = matrix.solve(sampleSolvableMatrixMatrixB).getQuadrupleData();
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:
  ##    public Matrix solve(Number[][] matrixB) (with Quadruple[][])
  ## Behaviors to be tested:
  ## -- solve(Number[][]) with Quadruple[][] with null argument throws NullPointerException
  ## -- solve(Number[][]) with Quadruple[][] with wrong matrix size throws IllegalArgumentException
  ## -- solve(Number[][]) with Quadruple[][] with a non-square argument throws IllegalArgumentException
  ## -- solve(Number[][]) with Quadruple[][] containing NaN throws IllegalArgumentException
  ## -- solve(Number[][]) with Quadruple[][] containing Infinity throws IllegalArgumentException
  ## -- solve(Number[][]) with Quadruple[][] containing null throws IllegalArgumentException
  ## -- solve(Number[][]) with Quadruple[][] with inconsistent data throws IllegalArgumentException
  ## -- solve(Number[][]) with Quadruple[][] with underdetermined data throws IllegalArgumentException
  ## -- solve(Number[][]) with Quadruple[][] does not spoil the passed matrixB array
  ## -- solve(Number[][]) with Quadruple[][] returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solve(Number[][]) with Quadruple[][] with null argument throws NullPointerException")
  void testSolveWithMatrixOfQuadruplesWithNullArgumentTrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solve((Quadruple[][])null),
        "Attempt to invoke solve() with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Number[][])' and 'is null'").
        contains("solve(number[][])").contains("is null");
  }

  @Test
  @DisplayName("solve(Number[][]) with Quadruple[][] with wrong matrix size throws IllegalArgumentException")
  void testSolveWithMatrixOfQuadruplesWithWrongMatrixSizeThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixTooLargeMatrixAsQuadruples),
        "Attempt to solve with wrong matrix size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Number[][])' and 'must have the same size'").
        contains("solve(number[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("solve(Number[][]) with Quadruple[][] with a non-square argument throws IllegalArgumentException")
  void testSolveWithMatrixOfQuadruplesWithNonSquareMatrixThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixNonSquareMatrixAsQuadruples),
        "Attempt to solve with non-square argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Number[][])' and 'must be square'").
        contains("solve(number[][])").contains("must be square");
  }

  @Test
  @DisplayName("solve(Number[][]) with Quadruple[][] containing NaN throws IllegalArgumentException")
  void testSolveWithMatrixOfQuadruplesContainingNaNThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixBAsQuadruples);
    badMatrix[1][1] = Quadruple.nan();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badMatrix),
        "Attempt to solve with an argument containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(number[][])' and 'must not contain NaN or Infinity'").
        contains("solve(number[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solve(Number[][]) with Quadruple[][] containing Infinity throws IllegalArgumentException")
  void testSolveWithMatrixOfQuadruplesContainingInfinityThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixBAsQuadruples);
    badMatrix[1][1] = Quadruple.positiveInfinity();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badMatrix),
        "Attempt to solve with an argument containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(number[][])' and 'must not contain NaN or Infinity'").
        contains("solve(number[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solve(Number[][]) with Quadruple[][] containing null throws IllegalArgumentException")
  void testSolveWithMatrixOfQuadruplesContainingNullThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixBAsQuadruples);
    badMatrix[1][1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badMatrix),
        "Attempt to solve with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(number[][])' and 'must not contain null'").
        contains("solve(number[][])").contains("must not contain null");
  }

  @Test
  @DisplayName("solve(Number[][]) with Quadruple[][] with inconsistent data throws IllegalArgumentException")
  void testSolveWithMatrixOfQuadruplesWithInconsistentDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixMatrixBAsQuadruples),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(Number[][]) with Quadruple[][] with underdetermined data throws IllegalArgumentException")
  void testSolveWithMatrixOfQuadruplesWithUnderdeterminedDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixMatrixBAsQuadruples),
        "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(Number[][]) with Quadruple[][] does not spoil the passed matrixB array")
  void testSolveWithMatrixOfQuadruplesDoesntSpoilThePassedArray() {
    final Quadruple[][] expectedMatrixB = deepCopyOf(sampleSolvableMatrixMatrixBAsQuadruples);

    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixBAsQuadruples);
    final Quadruple[][] actualMatrixB = deepCopyOf(sampleSolvableMatrixMatrixBAsQuadruples);

    assertThat(actualMatrixB).
        withFailMessage("matrix.solve(Number[][]) spoiled the passed matrix").
        isEqualTo(expectedMatrixB);
  }

  @Test
  @DisplayName("solve(Number[][]) with Quadruple[][] returns a correct solution")
  void testSolveWithMatrixOfQuadruplesReturnsCorrectSolution() {
    printMethodName();
    final Quadruple[][] expectedSolution = sampleSolvableMatrixMatrixXAsQuadruples;
    final double expectedError = ERROR_TOLERANCE;
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixDataAsQuadruples);

    final Quadruple[][] actualSolution = matrix.solve(sampleSolvableMatrixMatrixBAsQuadruples).getQuadrupleData();
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError ).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:
  ##    public Matrix solve(Number[][] matrixB) (with BigDecimal[][])
  ## Behaviors to be tested:
  ## -- solve(Number[][]) with BigDecimal[][] with null argument throws NullPointerException
  ## -- solve(Number[][]) with BigDecimal[][] with wrong matrix size throws IllegalArgumentException
  ## -- solve(Number[][]) with BigDecimal[][] with a non-square argument throws IllegalArgumentException
  ## -- solve(Number[][]) with BigDecimal[][] containing null throws IllegalArgumentException
  ## -- solve(Number[][]) with BigDecimal[][] with inconsistent data throws IllegalArgumentException
  ## -- solve(Number[][]) with BigDecimal[][] with underdetermined data throws IllegalArgumentException
  ## -- solve(Number[][]) with BigDecimal[][] does not spoil the passed matrixB array
  ## -- solve(Number[][]) with BigDecimal[][] returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solve(Number[][]) with BigDecimal[][] with null argument throws NullPointerException")
  void testSolveWithMatrixOfBigDecimalsWithNullArgumentTrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solve((BigDecimal[][])null),
        "Attempt to invoke solve() with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Number[][])' and 'is null'").
        contains("solve(number[][])").contains("is null");
  }

  @Test
  @DisplayName("solve(Number[][]) with BigDecimal[][] with wrong matrix size throws IllegalArgumentException")
  void testSolveWithMatrixOfBigDecimalsWithWrongMatrixSizeThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixTooLargeMatrixAsBigDecimals),
        "Attempt to solve with wrong matrix size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Number[][])' and 'must have the same size'").
        contains("solve(number[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("solve(Number[][]) with BigDecimal[][] with a non-square argument throws IllegalArgumentException")
  void testSolveWithMatrixOfBigDecimalsWithNonSquareMatrixThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixNonSquareMatrixAsBigDecimals),
        "Attempt to solve with non-square argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Number[][])' and 'must be square'").
        contains("solve(number[][])").contains("must be square");
  }

  @Test
  @DisplayName("solve(Number[][]) with BigDecimal[][] containing null throws IllegalArgumentException")
  void testSolveWithMatrixOfBigDecimalsContainingNullThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final BigDecimal[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixBAsBigDecimals);
    badMatrix[1][1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(badMatrix),
        "Attempt to solve with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solve(Number[][])' and 'must not contain null'").
        contains("solve(number[][])").contains("must not contain null");
  }

  @Test
  @DisplayName("solve(Number[][]) with BigDecimal[][] with inconsistent data throws IllegalArgumentException")
  void testSolveWithMatrixOfBigDecimalsWithInconsistentDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixMatrixBAsBigDecimals),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(Number[][]) with BigDecimal[][] with underdetermined data throws IllegalArgumentException")
  void testSolveWithMatrixOfBigDecimalsWithUnderdeterminedDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solve(sampleSolvableMatrixMatrixBAsBigDecimals),
        "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solve(Number[][]) with BigDecimal[][] does not spoil the passed matrixB array")
  void testSolveWithMatrixOfBigDecimalsDoesntSpoilThePassedArray() {
    final BigDecimal[][] expectedMatrixB = deepCopyOf(sampleSolvableMatrixMatrixBAsBigDecimals);
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    matrix.solve(sampleSolvableMatrixMatrixBAsBigDecimals);
    final BigDecimal[][] actualMatrixB = deepCopyOf(sampleSolvableMatrixMatrixBAsBigDecimals);

    assertThat(actualMatrixB).
        withFailMessage("matrix.solve(Number[][]) spoiled the passed matrix").
        isEqualTo(expectedMatrixB);
  }

  @Test
  @DisplayName("solve(Number[][]) with BigDecimal[][] returns a correct solution")
  void testSolveWithMatrixOfBigDecimalsReturnsCorrectSolution() {
    printMethodName();
    final Quadruple[][] expectedSolution = sampleSolvableMatrixMatrixXAsQuadruples;
    final double expectedError = ERROR_TOLERANCE;
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixDataAsQuadruples);

    final Quadruple[][] actualSolution = matrix.solve(sampleSolvableMatrixMatrixBAsBigDecimals).getQuadrupleData();
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError ).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:
  ##    public Matrix solveAccurately(Matrix matrixB)
  ## Behaviors to be tested:
  ## -- solveAccurately(Matrix) with null argument throws NullPointerException
  ## -- solveAccurately(Matrix) with wrong matrix size throws IllegalArgumentException
  ## -- solveAccurately(Matrix) with inconsistent data throws IllegalArgumentException
  ## -- solveAccurately(Matrix) with underdetermined data throws IllegalArgumentException
  ## -- solveAccurately(Matrix) returns a correct solution
  ## -- solveAccurately(Matrix) returns a much more accurate solution than solve(Matrix)
  ######################################################################################*/

  @Test
  @DisplayName("solveAccurately(Matrix) with null argument throws NullPointerException")
  void testSolveAccuratelyWithMatrixWithNullArgumentTrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveAccurately((Matrix)null),
        "Attempt to invoke solveAccurately() with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Matrix)' and 'is null'").
        contains("solveaccurately(matrix)").contains("is null");
  }

  @Test
  @DisplayName("solveAccurately(Matrix) with wrong matrix size throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixWithWrongMatrixSizeThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(new QuadrupleMatrix(sampleSolvableMatrixTooLargeMatrix)),
        "Attempt to solve with wrong matrix size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Matrix)' and 'must have the same size'").
        contains("solveaccurately(matrix)").contains("must have the same size");
  }

  @Test
  @DisplayName("solveAccurately(Matrix) with inconsistent data throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixWithInconsistentDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(new QuadrupleMatrix(sampleSolvableMatrixMatrixB)),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solveAccurately(Matrix) with underdetermined data throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixWithUnderdeterminedDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(new QuadrupleMatrix(sampleSolvableMatrixMatrixB)),
        "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solveAccurately(Matrix) returns a correct solution")
  void testSolveAccuratelyWithMatrixReturnsCorrectSolution() {
    printMethodName();
    final Quadruple[][] expectedSolution = sampleSolvableMatrixMatrixXAsQuadruples;
    final double expectedError = ERROR_TOLERANCE;
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixDataAsQuadruples);

    final Quadruple[][] actualSolution = matrix.solveAccurately(new QuadrupleMatrix(sampleSolvableMatrixMatrixB)).getQuadrupleData();
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError ).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /**
   * Generates a series of random matrices (data samples), solves the corresponding systems of linear equations,
   * and estimates the minimum and average ratio of
   * accuracy improvement achieved by the {@link QuadrupleMatrix#solveAccurately(double[])} compared with the
   * {@link QuadrupleMatrix#solve(double[])}. Asserts that the minimum improvement ratio is greater that {@link #EXPECTED_MIN_IMPROVEMENT}
   * and the average improvement ratio is greater than {@link #EXPECTED_AVR_IMPROVEMENT}
   * @see #EXPECTED_MIN_IMPROVEMENT
   */
  @Test
  // @Disabled // TO DO 2023-03-12 17:54:28 Uncomment for higher execution speed while debugging
  @DisplayName("solveAccurately(Matrix) returns a much more accurate solution than solve(Matrix)")
  void testSolveAccuratelyWithMatrixIsAccurate() {
    printMethodName();

    final OperationComparator comparator = makeAccurateVsSimpleMatrixSolutionComparator();

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
  ## Tested method:
  ##    public Matrix solveAccurately(double[][] matrixB)
  ## Behaviors to be tested:
  ## -- solveAccurately(double[][]) with null argument throws NullPointerException
  ## -- solveAccurately(double[][]) with wrong matrix size throws IllegalArgumentException
  ## -- solveAccurately(double[][]) with a non-square argument throws IllegalArgumentException
  ## -- solveAccurately(double[][]) with an argument containing NaN throws IllegalArgumentException
  ## -- solveAccurately(double[][]) with an argument containing Infinity throws IllegalArgumentException
  ## -- solveAccurately(double[][]) with inconsistent data throws IllegalArgumentException
  ## -- solveAccurately(double[][]) with underdetermined data throws IllegalArgumentException
  ## -- solveAccurately(double[][]) does not spoil the passed matrixB array
  ## -- solveAccurately(double[][]) returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solveAccurately(double[][]) with null argument throws NullPointerException")
  void testSolveAccuratelyWithMatrixOfDoublesWithNullArgumentTrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveAccurately((double[][])null),
        "Attempt to invoke solveAccurately() with null argument must throw an exception");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(double[][])' and 'is null'").
        contains("solveaccurately(double[][])").contains("is null");
  }

  @Test
  @DisplayName("solveAccurately(double[][]) with wrong matrix size throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfDoublesWithWrongMatrixSizeThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixTooLargeMatrix),
        "Attempt to solve with wrong matrix size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(double[][])' and 'must have the same size'").
        contains("solveaccurately(double[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("solveAccurately(double[][]) with a non-square argument throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfDoublesWithNonSquareMatrixThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixNonSquareMatrix),
        "Attempt to solve with non-square argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(double[][])' and 'must be square'").
        contains("solveaccurately(double[][])").contains("must be square");
  }

  @Test
  @DisplayName("solveAccurately(double[][]) with an argument containing NaN throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfDoublesContainingNaNThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final double[][] badMatrix = deepCopyOf(sampleSolvableMatrixData);
    badMatrix[1][1] = Double.NaN;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badMatrix),
        "Attempt to call solveAccurately(double[][]) with matrix containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(double[][])' and 'must not contain NaN or Infinity'").
        contains("solveaccurately(double[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveAccurately(double[][]) with an argument containing Infinity throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfDoublesContainingInfinityThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final double[][] badMatrix = deepCopyOf(sampleSolvableMatrixData);
    badMatrix[1][1] = Double.POSITIVE_INFINITY;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badMatrix),
        "Attempt to solve with matrix containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(double[][])' and 'must not contain NaN or Infinity'").
        contains("solveaccurately(double[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveAccurately(double[][]) with inconsistent data throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfDoublesWithInconsistentDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixMatrixB),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solveAccurately(double[][]) with underdetermined data throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfDoublesWithUnderdeterminedDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixMatrixB),
        "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solveAccurately(double[][]) does not spoil the passed matrixB array")
  void testSolveAccuratelyWithMatrixOfDoublesDoesntSpoilThePassedArray() {
    final double[][] expectedMatrixB = deepCopyOf(sampleSolvableMatrixMatrixB);

    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    matrix.solveAccurately(sampleSolvableMatrixMatrixB);
    final double[][] actualMatrixB = deepCopyOf(sampleSolvableMatrixMatrixB);

    assertThat(actualMatrixB).
        withFailMessage("matrix.solveAccurately(double[][]) spoiled the passed matrix").
        isEqualTo(expectedMatrixB);
  }

  @Test
  @DisplayName("solveAccurately(double[][]) returns a correct solution")
  void testSolveAccuratelyWithMatrixOfDoublesReturnsCorrectSolution() {
    printMethodName();
    final Quadruple[][] expectedSolution = sampleSolvableMatrixMatrixXAsQuadruples;
    final double expectedError = ERROR_TOLERANCE;
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] actualSolution = matrix.solveAccurately(sampleSolvableMatrixMatrixB).getQuadrupleData();
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError ).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:
  ##    public Matrix solveAccurately(Number[][] matrixB) (With Quadruple[][])
  ## Behaviors to be tested:
  ## -- solveAccurately(Number[][]) with Quadruple[][] with null argument throws NullPointerException
  ## -- solveAccurately(Number[][]) with Quadruple[][] with wrong matrix size throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with Quadruple[][] with a non-square argument throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with Quadruple[][] containing NaN throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with Quadruple[][] containing Infinity throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with Quadruple[][] containing null throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with Quadruple[][] with inconsistent data throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with Quadruple[][] with underdetermined data throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with Quadruple[][] does not spoil the passed matrixB array
  ## -- solveAccurately(Number[][]) with Quadruple[][] returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solveAccurately(Number[][]) with Quadruple[][] with null argument throws NullPointerException")
  void testSolveAccuratelyWithMatrixOfQuadruplesWithNullArgumentTrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveAccurately((Quadruple[][])null),
        "Attempt to invoke solveAccurately() with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Number[][])' and 'is null'").
        contains("solveaccurately(number[][])").contains("is null");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with Quadruple[][] with wrong matrix size throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfQuadruplesWithWrongMatrixSizeThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixTooLargeMatrixAsQuadruples),
        "Attempt to solve with wrong matrix size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Number[][])' and 'must have the same size'").
        contains("solveaccurately(number[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with Quadruple[][] with a non-square argument throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfQuadruplesWithNonSquareMatrixThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixNonSquareMatrixAsQuadruples),
        "Attempt to solve with non-square argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Number[][])' and 'must be square'").
        contains("solveaccurately(number[][])").contains("must be square");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with Quadruple[][] containing NaN throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfQuadruplesContainingNaNThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixBAsQuadruples);
    badMatrix[1][1] = Quadruple.nan();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badMatrix),
        "Attempt to solve with an argument containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(number[][])' and 'must not contain NaN or Infinity'").
        contains("solveaccurately(number[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with Quadruple[][] containing Infinity throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfQuadruplesContainingInfinityThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixBAsQuadruples);
    badMatrix[1][1] = Quadruple.positiveInfinity();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badMatrix),
        "Attempt to solve with an argument containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(number[][])' and 'must not contain NaN or Infinity'").
        contains("solveaccurately(number[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with Quadruple[][] containing null throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfQuadruplesContainingNullThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixBAsQuadruples);
    badMatrix[1][1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badMatrix),
        "Attempt to solve with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(number[][])' and 'must not contain null'").
        contains("solveaccurately(number[][])").contains("must not contain null");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with Quadruple[][] with inconsistent data throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfQuadruplesWithInconsistentDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixMatrixBAsQuadruples),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with Quadruple[][] with underdetermined data throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfQuadruplesWithUnderdeterminedDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixMatrixBAsQuadruples),
        "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with Quadruple[][] does not spoil the passed matrixB array")
  void testSolveAccuratelyWithMatrixOfQuadruplesDoesntSpoilThePassedArray() {
    final Quadruple[][] expectedMatrixB = deepCopyOf(sampleSolvableMatrixMatrixBAsQuadruples);
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    matrix.solveAccurately(sampleSolvableMatrixMatrixBAsQuadruples);
    final Quadruple[][] actualMatrixB = deepCopyOf(sampleSolvableMatrixMatrixBAsQuadruples);

    assertThat(actualMatrixB).
        withFailMessage("matrix.solveAccurately(double[][]) spoiled the passed matrix").
        isEqualTo(expectedMatrixB);
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with Quadruple[][] returns a correct solution")
  void testSolveAccuratelyWithMatrixOfQuadruplesReturnsCorrectSolution() {
    printMethodName();
    final Quadruple[][] expectedSolution = sampleSolvableMatrixMatrixXAsQuadruples;
    final double expectedError = ERROR_TOLERANCE;
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] actualSolution = matrix.solveAccurately(sampleSolvableMatrixMatrixBAsQuadruples).getQuadrupleData();
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError ).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:
  ##    public Matrix solveAccurately(Number[][] matrixB) (with BigDecimal[][])
  ## Behaviors to be tested:
  ## -- solveAccurately(Number[][]) with BigDecimal[][] with null argument throws NullPointerException
  ## -- solveAccurately(Number[][]) with BigDecimal[][] with wrong matrix size throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with BigDecimal[][] with a non-square argument throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with BigDecimal[][] containing null throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with BigDecimal[][] with inconsistent data throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with BigDecimal[][] with underdetermined data throws IllegalArgumentException
  ## -- solveAccurately(Number[][]) with BigDecimal[][] does not spoil the passed matrixB array
  ## -- solveAccurately(Number[][]) with BigDecimal[][] returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("solveAccurately(Number[][]) with BigDecimal[][] with null argument throws NullPointerException")
  void testSolveAccuratelyWithMatrixOfBigDecimalsWithNullArgumentTrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.solveAccurately((BigDecimal[][])null),
        "Attempt to invoke solveAccurately() with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Number[][])' and 'is null'").
        contains("solveaccurately(number[][])").contains("is null");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with BigDecimal[][] with wrong matrix size throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfBigDecimalsWithWrongMatrixSizeThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixTooLargeMatrixAsBigDecimals),
        "Attempt to solve with wrong matrix size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Number[][])' and 'must have the same size'").
        contains("solveaccurately(number[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with BigDecimal[][] with a non-square argument throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfBigDecimalsWithNonSquareMatrixThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixNonSquareMatrixAsBigDecimals),
        "Attempt to solve with non-square argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Number[][])' and 'must be square'").
        contains("solveaccurately(number[][])").contains("must be square");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with BigDecimal[][] containing null throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfBigDecimalsContainingNullThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final BigDecimal[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixBAsBigDecimals);
    badMatrix[1][1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(badMatrix),
        "Attempt to solve with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'solveAccurately(Number[][])' and 'must not contain null'").
        contains("solveaccurately(number[][])").contains("must not contain null");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with BigDecimal[][] with inconsistent data throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfBigDecimalsWithInconsistentDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixMatrixBAsBigDecimals),
        "Attempt to solve with inconsistent data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with BigDecimal[][] with underdetermined data throws IllegalArgumentException")
  void testSolveAccuratelyWithMatrixOfBigDecimalsWithUnderdeterminedDataThrowsException() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleUnderdeterminedMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.solveAccurately(sampleSolvableMatrixMatrixBAsBigDecimals),
        "Attempt to solve with underdetermined data must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'matrix is unsolvable'").
        contains("matrix is unsolvable");
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with BigDecimal[][] does not spoil the passed matrixB array")
  void testSolveAccuratelyWithMatrixOfBigDecimalsDoesntSpoilThePassedArray() {
    final BigDecimal[][] expectedMatrixB = deepCopyOf(sampleSolvableMatrixMatrixBAsBigDecimals);
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    matrix.solveAccurately(sampleSolvableMatrixMatrixBAsBigDecimals);
    final BigDecimal[][] actualMatrixB = deepCopyOf(sampleSolvableMatrixMatrixBAsBigDecimals);

    assertThat(actualMatrixB).
        withFailMessage("matrix.solveAccurately(double[][]) spoiled the passed matrix").
        isEqualTo(expectedMatrixB);
  }

  @Test
  @DisplayName("solveAccurately(Number[][]) with BigDecimal[][] returns a correct solution")
  void testSolveAccuratelyWithMatrixOfBigDecimalsReturnsCorrectSolution() {
    printMethodName();
    final Quadruple[][] expectedSolution = sampleSolvableMatrixMatrixXAsQuadruples;
    final double expectedError = ERROR_TOLERANCE;
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] actualSolution = matrix.solveAccurately(sampleSolvableMatrixMatrixBAsBigDecimals).getQuadrupleData();
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError ).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:
  ##    public Matrix getMatrixSolution()
  ## Behaviors to be tested:
  ## -- getMatrixSolution() returns null if no solution was found so far
  ## -- getMatrixSolution() returns null if an attempt to solve an equation was unsuccessful
  ## -- getMatrixSolution() after successful solutions returns data equal to those returned by solve()
  ## -- getMatrixSolution() returns the last found solution
  ######################################################################################*/

  @Test
  @DisplayName("getMatrixSolution() returns null if no solution was found so far")
  void testGetMatrixSolutionReturnsNullIfNoSolutionWasFound() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Matrix actualSolution = matrix.getMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getMatrixSolution() returns non-null when no solutions for the matrix was found").
        isNull();
  }

  @Test
  @DisplayName("getMatrixSolution() returns null if an attempt to solve an equation was unsuccessful")
  void testGetMatrixSolutionReturnsNullIfSolutionWasUnsuccessful() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    try {
      matrix.solve(sampleSolvableMatrixMatrixB);
    } catch (final Exception x) {}
    final Matrix actualSolution = matrix.getMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getMatrixSolution() returns non-null after solution has failed").
        isNull();
  }

  @Test
  @DisplayName("getMatrixSolution() after successful solutions returns data equal to those returned by solve()")
  void testGetMatrixSolutionReturnsCorrectData() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    final Matrix expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB);

    final Matrix actualSolution = matrix.getMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getMatrixSolution() returns data that differ from what is returned by solve()").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getMatrixSolution() returns the last found solution")
  void testGetMatrixSolutionReturnsTheLastSolution() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    matrix.solve(sampleSolvableMatrixMatrixB);
    final Matrix expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB);

    final Matrix actualSolution = matrix.getMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getMatrixSolution() returns data that differ from the last found solution").
        isEqualTo(expectedSolution);
  }

  /*####################################################################################
  ## Tested method:
  ##    public Number[][] getNumberMatrixSolution()
  ## Behaviors to be tested:
  ## -- getNumberMatrixSolution() returns null if no solution was found so far
  ## -- getNumberMatrixSolution() returns null if an attempt to solve an equation was unsuccessful
  ## -- getNumberMatrixSolution() after successful solutions returns data equal to those returned by solve()
  ## -- getNumberMatrixSolution() returns a copy of the internally stored data, not a reference of it
  ## -- getNumberMatrixSolution() returns the last found solution
  ######################################################################################*/

  @Test
  @DisplayName("getNumberMatrixSolution() returns null if no solution was found so far")
  void testGetNumberMatrixSolutionReturnsNullIfNoSolutionWasFound() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Number[][] actualSolution = matrix.getNumberMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getNumberMatrixSolution() returns non-null when no solutions for the matrix was found").
        isNull();
  }

  @Test
  @DisplayName("getNumberMatrixSolution() returns null if an attempt to solve an equation was unsuccessful")
  void testGetNumberMatrixSolutionReturnsNullIfSolutionWasUnsuccessful() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    try {
      matrix.solve(sampleSolvableMatrixMatrixB);
    } catch (final Exception x) {}
    final Number[][] actualSolution = matrix.getNumberMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getNumberMatrixSolution() returns non-null after solution has failed").
        isNull();
  }

  @Test
  @DisplayName("getNumberMatrixSolution() after successful solutions returns data equal to those returned by solve()")
  void testGetNumberMatrixSolutionReturnsCorrectData() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    final Number[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getData();

    final Number[][] actualSolution = matrix.getNumberMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getNumberMatrixSolution() returns data that differ from what is returned by solve()").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getNumberMatrixSolution() returns a copy of the internally stored data, not a reference of it")
  void testGetNumberMatrixSolutionReturnsACopyOfInternalData() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    final Number[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getData();

    final Number[][] tmpSolution = matrix.getNumberMatrixSolution();
    ((Quadruple)tmpSolution[0][0]).add(12345);
    final Number[][] actualSolution = matrix.getNumberMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getNumberMatrixSolution() returns a reference of the internally stored data").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getNumberMatrixSolution() returns the last found solution")
  void testGetNumberMatrixSolutionReturnsTheLastSolution() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixB);
    final Number[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getData();

    final Number[][] actualSolution = matrix.getNumberMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getNumberMatrixSolution() returns data that differ from the last found solution").
        isEqualTo(expectedSolution);
  }

  /*####################################################################################
  ## Tested method:
  ##    public double[][] getDoubleMatrixSolution()
  ## Behaviors to be tested:
  ## -- getDoubleMatrixSolution() returns null if no solution was found so far
  ## -- getDoubleMatrixSolution() returns null if an attempt to solve an equation was unsuccessful
  ## -- getDoubleMatrixSolution() after successful solutions returns data equal to those returned by solve()
  ## -- getDoubleMatrixSolution() returns a copy of the internally stored data, not a reference of it
  ## -- getDoubleMatrixSolution() returns the last found solution
  ######################################################################################*/

  @Test
  @DisplayName("getDoubleMatrixSolution() returns null if no solution was found so far")
  void testGetDoubleMatrixSolutionReturnsNullIfNoSolutionWasFound() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final double[][] actualSolution = matrix.getDoubleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getDoubleMatrixSolution() returns non-null when no solutions for the matrix was found").
        isNull();
  }

  @Test
  @DisplayName("getDoubleMatrixSolution() returns null if an attempt to solve an equation was unsuccessful")
  void testGetDoubleMatrixSolutionReturnsNullIfSolutionWasUnsuccessful() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    try {
      matrix.solve(sampleSolvableMatrixMatrixB);
    } catch (final Exception x) {}
    final double[][] actualSolution = matrix.getDoubleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getDoubleMatrixSolution() returns non-null after solution has failed").
        isNull();
  }

  @Test
  @DisplayName("getDoubleMatrixSolution() after successful solutions returns data equal to those returned by solve()")
  void testGetDoubleMatrixSolutionReturnsCorrectData() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    final double[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getDoubleData();

    final double[][] actualSolution = matrix.getDoubleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getDoubleMatrixSolution() returns data that differ from what is returned by solve()").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getDoubleMatrixSolution() returns a copy of the internally stored data, not a reference of it")
  void testGetDoubleMatrixSolutionReturnsACopyOfInternalData() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    final double[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getDoubleData();

    final double[][] tmpSolution = matrix.getDoubleMatrixSolution();
    tmpSolution[0][0] = 12345;
    final double[][] actualSolution = matrix.getDoubleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getDoubleMatrixSolution() returns a reference of the internally stored data").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getDoubleMatrixSolution() returns the last found solution")
  void testGetDoubleMatrixSolutionReturnsTheLastSolution() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixB);
    final double[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getDoubleData();

    final double[][] actualSolution = matrix.getDoubleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getDoubleMatrixSolution() returns data that differ from the last found solution").
        isEqualTo(expectedSolution);
  }

  /*####################################################################################
  ## Tested method:
  ##    public Quadruple[][] getQuadrupleMatrixSolution()
  ## Behaviors to be tested:
  ## -- getQuadrupleMatrixSolution() returns null if no solution was found so far
  ## -- getQuadrupleMatrixSolution() returns null if an attempt to solve an equation was unsuccessful
  ## -- getQuadrupleMatrixSolution() after successful solutions returns data equal to those returned by solve()
  ## -- getQuadrupleMatrixSolution() returns a copy of the internally stored data, not a reference of it
  ## -- getQuadrupleMatrixSolution() returns the last found solution
  ######################################################################################*/

  @Test
  @DisplayName("getQuadrupleMatrixSolution() returns null if no solution was found so far")
  void testGetQuadrupleMatrixSolutionReturnsNullIfNoSolutionWasFound() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] actualSolution = matrix.getQuadrupleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getQuadrupleMatrixSolution() returns non-null when no solutions for the matrix was found").
        isNull();
  }

  @Test
  @DisplayName("getQuadrupleMatrixSolution() returns null if an attempt to solve an equation was unsuccessful")
  void testGetQuadrupleMatrixSolutionReturnsNullIfSolutionWasUnsuccessful() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    try {
      matrix.solve(sampleSolvableMatrixMatrixB);
    } catch (final Exception x) {}
    final Quadruple[][] actualSolution = matrix.getQuadrupleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getQuadrupleMatrixSolution() returns non-null after solution has failed").
        isNull();
  }

  @Test
  @DisplayName("getQuadrupleMatrixSolution() after successful solutions returns data equal to those returned by solve()")
  void testGetQuadrupleMatrixSolutionReturnsCorrectData() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    final Quadruple[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getQuadrupleData();

    final Quadruple[][] actualSolution = matrix.getQuadrupleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getQuadrupleMatrixSolution() returns data that differ from what is returned by solve()").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getQuadrupleMatrixSolution() returns a copy of the internally stored data, not a reference of it")
  void testGetQuadrupleMatrixSolutionReturnsACopyOfInternalData() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    final Quadruple[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getQuadrupleData();

    final Quadruple[][] tmpSolution = matrix.getQuadrupleMatrixSolution();
    tmpSolution[0][0].add(12345);
    final Quadruple[][] actualSolution = matrix.getQuadrupleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getQuadrupleMatrixSolution() returns a reference of the internally stored data").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getQuadrupleMatrixSolution() returns the last found solution")
  void testGetQuadrupleMatrixSolutionReturnsTheLastSolution() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixB);
    final Quadruple[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getQuadrupleData();

    final Quadruple[][] actualSolution = matrix.getQuadrupleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getQuadrupleMatrixSolution() returns data that differ from the last found solution").
        isEqualTo(expectedSolution);
  }

  /*####################################################################################
  ## Tested method:
  ##    BigDecimal[][] getBigDecimalMatrixSolution()
  ## Behaviors to be tested:
  ## -- getBigDecimalMatrixSolution() returns null if no solution was found so far
  ## -- getBigDecimalMatrixSolution() returns null if an attempt to solve an equation was unsuccessful
  ## -- getBigDecimalMatrixSolution() after successful solutions returns data equal to those returned by solve()
  ## -- getBigDecimalMatrixSolution() returns a copy of the internally stored data, not a reference of it
  ## -- getBigDecimalMatrixSolution() returns the last found solution
  ######################################################################################*/

  @Test
  @DisplayName("getBigDecimalMatrixSolution() returns null if no solution was found so far")
  void testGetBigDecimalMatrixSolutionReturnsNullIfNoSolutionWasFound() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    final BigDecimal[][] actualSolution = matrix.getBigDecimalMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getBigDecimalMatrixSolution() returns non-null when no solutions for the matrix was found").
        isNull();
  }

  @Test
  @DisplayName("getBigDecimalMatrixSolution() returns null if an attempt to solve an equation was unsuccessful")
  void testGetBigDecimalMatrixSolutionReturnsNullIfSolutionWasUnsuccessful() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);

    try {
      matrix.solve(sampleSolvableMatrixMatrixB);
    } catch (final Exception x) {}
    final BigDecimal[][] actualSolution = matrix.getBigDecimalMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getBigDecimalMatrixSolution() returns non-null after solution has failed").
        isNull();
  }

  @Test
  @DisplayName("getBigDecimalMatrixSolution() after successful solutions returns data equal to those returned by solve()")
  void testGetBigDecimalMatrixSolutionReturnsCorrectData() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    final BigDecimal[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getBigDecimalData();

    final BigDecimal[][] actualSolution = matrix.getBigDecimalMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getBigDecimalMatrixSolution() returns data that differ from what is returned by solve()").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getBigDecimalMatrixSolution() returns a copy of the internally stored data, not a reference of it")
  void testGetBigDecimalMatrixSolutionReturnsACopyOfInternalData() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    final BigDecimal[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getBigDecimalData();

    final BigDecimal[][] tmpSolution = matrix.getBigDecimalMatrixSolution();
    tmpSolution[0][0] = tmpSolution[0][0].add(new BigDecimal(12345));
    final BigDecimal[][] actualSolution = matrix.getBigDecimalMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getBigDecimalMatrixSolution() returns a reference of the internally stored data").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("getBigDecimalMatrixSolution() returns the last found solution")
  void testGetBigDecimalMatrixSolutionReturnsTheLastSolution() {
    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixB);
    final BigDecimal[][] expectedSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixB).getBigDecimalData();

    final BigDecimal[][] actualSolution = matrix.getBigDecimalMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("getBigDecimalMatrixSolution() returns data that differ from the last found solution").
        isEqualTo(expectedSolution);
  }

  /*####################################################################################
  ## Tested method:
  ##    public String getErrorCode()
  ## Behaviors to be tested:
  ## -- getErrorCode() returns OK after a successful matrix solution
  ## -- getErrorCode() returns appropriate code after an attempt to solve an inconsistent matrix
  ## -- getErrorCode() returns appropriate code after an attempt to solve an underdetermined matrix
  ######################################################################################*/

  @Test
  @DisplayName("getErrorCode() returns OK after a successful matrix solution")
  void testGetErrorCodeReturnsOkIfSolutionWasSuccessful() {
    final String expectedResult = "OK";

    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixB);
    final String actualResult = matrix.getErrorCode();

    assertThat(actualResult).
        withFailMessage("getErrorCode() returns a wrong result after successfull call to solve()").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("getErrorCode() returns appropriate code after an attempt to solve an inconsistent matrix")
  void testGetRightErrorCodeForInconsistentMatrix() {
    final String expectedResult = "NON_INVERTIBLE";

    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleInconsistentMatrixData);
    try {
      matrix.solve(sampleSolvableMatrixMatrixB);
    } catch (final Exception x) {}
    final String actualResult = matrix.getErrorCode();

    assertThat(actualResult).
        withFailMessage("getErrorCode() returns a wrong result after solve() was called for an inconsistent matrix").
        isEqualTo(expectedResult);
  }

  @Test
  @DisplayName("getErrorCode() returns appropriate code after an attempt to solve an underdetermined matrix")
  void testGetRightErrorCodeForUnderdeterminedMatrix() {
    final String expectedResult = "NON_INVERTIBLE";

    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleUnderdeterminedMatrixData);
    try {
      matrix.solve(sampleSolvableMatrixMatrixB);
    } catch (final Exception x) {}
    final String actualResult = matrix.getErrorCode();

    assertThat(actualResult).
        withFailMessage("getErrorCode() returns a wrong result after solve() was called for an underdetermined matrix").
        isEqualTo(expectedResult);
  }

  /*####################################################################################
  ## Tested behaviors:
  ##   Consequent calls to solve() don't interfere with each other, namely:
  ## -- another call to solve(double[][]) with another matrixB also returns a correct solution
  ## -- solve(double[][]) after successful SPD-solution returns a correct solution
  ## -- solve(double[][]) after unsuccessful SPD-solution returns a correct solution
  ######################################################################################*/

  @Test
  @DisplayName("another call to solve(double[][]) with another matrixB also returns a correct solution")
  void testRepeatedSolveWithMatrixOfDoublesReturnsCorrectSolution() {
    printMethodName();
    final Quadruple[][] expectedSolution = sampleSolvableMatrixAnotherMatrixXAsQuadruples;
    final double expectedError = ERROR_TOLERANCE;

    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixDataAsQuadruples);
    matrix.solve(sampleSolvableMatrixMatrixBAsQuadruples);
    final Quadruple[][] actualSolution = matrix.solve(sampleSolvableMatrixAnotherMatrixBAsQuadruples).getQuadrupleData();
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError ).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  @Test
  @DisplayName("solve(double[][]) after successful SPD-solution returns a correct solution")
  void testSolveWithMatrixAfterSolveSPDReturnsCorrectSolution() {
    printMethodName();
    final Quadruple[][] expectedSolution = sampleSpdSolvableMatrixMatrixXAsQuadruples; // It should be a SPD-solvable matrix
    final double expectedError = ERROR_TOLERANCE;

    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSpdSolvableMatrixDataAsQuadruples);
    matrix.solveSPD(sampleSpdSolvableMatrixAnotherVector);
    final Quadruple[][] actualSolution = matrix.solve(sampleSpdSolvableMatrixMatrixBAsQuadruples).getQuadrupleData();
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError ).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  @Test
  @DisplayName("solve(double[][]) after unsuccessful SPD-solution returns a correct solution")
  void testSolveWithMatrixAfterUnsucessfulSolveSPDReturnsCorrectSolution() {
    printMethodName();
    final Quadruple[][] expectedSolution = sampleSolvableMatrixMatrixXAsQuadruples; // It should be non-SPD-solvable matrix
    final double expectedError = ERROR_TOLERANCE;

    final QuadrupleMatrix matrix = new QuadrupleMatrix(sampleSolvableMatrixData);
    try {
      matrix.solveSPD(sampleSpdSolvableMatrixAnotherVector);  // No matter what vector. It's a non-SPD matrix, will throw an exception
    } catch (final Exception x) {                             /* Is expected. Ignore it */ }
    final Quadruple[][] actualSolution = matrix.solve(sampleSolvableMatrixMatrixBAsQuadruples).getQuadrupleData();
    final ErrorSet actualErrors = findErrors(expectedSolution, actualSolution, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError ).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }


  /* ************************************************************************************
  /*** Private methods ******************************************************************
  /**************************************************************************************/

  private static OperationComparator makeAccurateVsSimpleMatrixSolutionComparator() {
    return new OperationComparator(
        () -> MatrixData.makeDataSetForMatrixSolutions(LARGE_MATRIX_SIZE, random),
        MatrixData::quadrupleMatrixSolutionErrors,
        MatrixData::quadrupleAccurateMatrixSolutionErrors
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
