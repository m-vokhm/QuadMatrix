package com.mvohm.quadmatrix;

import java.math.BigDecimal;
import java.util.Locale;
import java.util.Random;

import org.assertj.core.data.Offset;
import org.junit.jupiter.api.*;

import static org.assertj.core.api.Assertions.assertThat;
import static org.junit.jupiter.api.Assertions.*;

import com.mvohm.quadmatrix.test.AuxMethods;
import com.mvohm.quadmatrix.test.MatrixData;
import com.mvohm.quadmatrix.test.OperationComparator;
import com.mvohm.quadmatrix.test.AuxMethods.ErrorSet;
import com.mvohm.quadruple.Quadruple;

import static com.mvohm.quadmatrix.test.MatrixDataSamples.*;
import static com.mvohm.quadmatrix.test.AuxMethods.*;
import static com.mvohm.quadmatrix.test.MatrixDataGenerators.*;

/**
 * Tests for {@link DoubleMatrix} methods for inversion, transposition, matrix arithmetics etc.
 * <br>
 * Tested methods:
 *    public Matrix inverse()
 *    public Matrix inverseAccurately()
 *    public Matrix transpose()
 *    public Matrix unity()
 *    public Matrix multiply(Matrix factor)
 *    public Matrix multiply(double[][] factor)
 *    public Matrix multiply(Number[][] factor)
 *    public Number[] multiply(double[] vector)
 *    public Number[] multiply(Number[] vector)
 *    public Matrix multiply(double scalar)
 *    public Matrix multiply(Number scalar)
 *    public Matrix add(Matrix matrixB)
 *    public Matrix add(double[][] matrixB)
 *    public Matrix add(Number[][] matrixB)
 *    public Matrix subtract(Matrix matrixB)
 *    public Matrix subtract(double[][] matrixB)
 *    public Matrix subtract(Number[][] matrixB)
 *    public Number determinant()
 *    public double determinantAsDouble()
 *    public Quadruple determinantAsQuadruple()
 *    public BigDecimal determinantAsBigDecimal()
 */

public class DoubleMatrixOtherMethodsTests {

  /** Enables progress indication for time-consuming tests and printing values of some computational errors */
  private static final boolean OUTPUT_ENABLED = true;
  private static final boolean DETAILED_OUTPUT = true;

  private static final int    RAND_SEED
//                                          = 123;
//                                          = 12121234;
                                          = 321;

  private static final double ERROR_TOLERANCE             = 1e-15;

  private static final double ANY_DOUBLE_VALUE            = Math.E;
  private static final double DOUBLE_FACTOR               = ANY_DOUBLE_VALUE;

  // To perform iterations to evaluate the accuracy of inverseAccurately()
  // With error estimation based on Quadruple calculations, it takes about 1 minute
  // to run testInverseAccuratelyIsAccurate() with 60 iterations and matrices 157x157 on a i5 3.2 HHz machine
  private static final int    REFINEMENT_TEST_ITERATIONS  = 60; // Approx. 1 min
  private static final int    LARGE_MATRIX_SIZE
//                                                          = 200;
                                                          = 158; // So that accurate inversion test takes approx. 1s on my machine

  private static final double EXPECTED_MIN_IMPROVEMENT    =   20;
  private static final double EXPECTED_AVR_IMPROVEMENT    =   25;

  // Data for test matrices
  private static final double[]       sampleSolvableMatrixVector                = sampleSolvableMatrixVector();
  private static final double[]       sampleSolvableMatrixSolution              = sampleSolvableMatrixSolution();
  private static final Quadruple[]    sampleSolvableMatrixSolutionAsQuadruples  = convertToQuadruples(sampleSolvableMatrixSolution);
  private static final BigDecimal[]   sampleSolvableMatrixSolutionAsBigDecimals = convertToBigDecimals(sampleSolvableMatrixSolution);

  private static final double[]       sampleTooLongMatrixVector                 = sampleTooLongMatrixVector();
  private static final Quadruple[]    sampleTooLongMatrixVectorAsQuadruples     = convertToQuadruples(sampleTooLongMatrixVector);
  private static final BigDecimal[]   sampleTooLongMatrixVectorAsBigDecimals    = convertToBigDecimals(sampleTooLongMatrixVector);

  private static final double[][]
      sampleSolvableMatrixData                          = sampleSolvableMatrixData(),
      sampleInvertibleMatrixData                        = sampleInvertibleMatrixData(),
      sampleInvertibleMatrixInversion                   = sampleInvertibleMatrixInversion(),
      sampleInvertibleMatrixTransposition               = sampleInvertibleMatrixTransposition(),

      sampleSolvableMatrixMatrixX                       = sampleSolvableMatrixMatrixX(),
      sampleSolvableMatrixAnotherMatrixX                = sampleSolvableMatrixAnotherMatrixX(),

      sampleSolvableMatrixMatrixB                       = sampleSolvableMatrixMatrixB(),
      sampleSolvableMatrixTooLargeMatrix                = sampleSolvableMatrixTooLargeMatrix(),
      sampleSolvableMatrixNonSquareMatrix               = sampleSolvableMatrixNonSquareMatrix(),
      sampleInconsistentMatrixData                      = sampleInconsistentMatrixData(),
      sampleUnderdeterminedMatrixData                   = sampleUnderdeterminedMatrixData(),

      sampleRandomMatrixData                            = sampleRandomMatrixData(),
      anyMatrixDataNoMatterWhatItIs                     = randomMatrix(LARGE_MATRIX_SIZE),

      summand1                                          = summand1(),
      summand2                                          = summand2(),
      summand1PlusSummand2                              = summand1PlusSummand2(),

      minuend                                           = minuend(),
      subtrahend                                        = subtrahend(),
      minuendMinusSubtrahend                            = minuendMinusSubtrahend();

  private static final Quadruple[][]
      sampleSolvableMatrixMatrixXAsQuadruples           = convertToQuadruples(sampleSolvableMatrixMatrixX),
      sampleSolvableMatrixTooLargeMatrixAsQuadruples    = convertToQuadruples(sampleSolvableMatrixTooLargeMatrix),
      sampleSolvableMatrixNonSquareMatrixAsQuadruples   = convertToQuadruples(sampleSolvableMatrixNonSquareMatrix),
      summand2AsQuadruples                              = convertToQuadruples(summand2),
      subtrahendAsQuadruples                            = convertToQuadruples(subtrahend);

  private static final BigDecimal[][]
      sampleSolvableMatrixMatrixXAsBigDecimals          = convertToBigDecimals(sampleSolvableMatrixMatrixX),
      sampleSolvableMatrixTooLargeMatrixAsBigDecimals   = convertToBigDecimals(sampleSolvableMatrixTooLargeMatrix),
      sampleSolvableMatrixNonSquareMatrixAsBigDecimals  = convertToBigDecimals(sampleSolvableMatrixNonSquareMatrix),
      summand2AsBigDecimals                             = convertToBigDecimals(summand2),
      subtrahendAsBigDecimals                           = convertToBigDecimals(subtrahend);

  private static double[][]             sampleRandomMultipliedByFactor;
  private static Random                 random;

  @BeforeAll
  static void setup() {
    Locale.setDefault(Locale.US);
    DoubleMatrix.setDefaultScaling(true);

    sampleRandomMultipliedByFactor = deepCopyOf(sampleRandomMatrixData);
    final int size = sampleRandomMultipliedByFactor.length;
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        sampleRandomMultipliedByFactor[i][j] = sampleRandomMatrixData[i][j] * DOUBLE_FACTOR;
      }
    }
  }

  @BeforeEach
  void initRandom() {
    random = new Random(RAND_SEED);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix inverse()
  ## Behaviors to be tested:
  ## -- inverse() does not spoil the source matrix
  ## -- inverse() does not spoil the previously found vector solution
  ## -- inverse() does not spoil the previously found matrix solution
  ## -- inverse() returns the inversion of the matrix
  ######################################################################################*/

  @Test
  @DisplayName("inverse() does not spoil the source matrix")
  void testInverseDoesNotSpoilTheMatrix() {
    final double[][] expectedMatrixData = sampleSolvableMatrixData;

    final DoubleMatrix matrix = new DoubleMatrix(expectedMatrixData);
    matrix.inverse();
    final double[][] actualMatrixData = matrix.getDoubleData();

    assertThat(actualMatrixData).
        withFailMessage("matrix.inverse() spoils the source matrix data").
        isEqualTo(expectedMatrixData);
  }

  @Test
  @DisplayName("inverse() does not spoil the previously found vector solution")
  void testInverseDoesNotSpoilTheVectorSolution() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixVector);
    final double[] expectedSolution = matrix.getDoubleSolution();

    matrix.inverse();
    final double[] actualSolution = matrix.getDoubleSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.inverse() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("inverse() does not spoil the previously found matrix solution")
  void testInverseDoesNotSpoilTheMatrixSolution() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixB);
    final double[][] expectedSolution = matrix.getDoubleMatrixSolution();

    matrix.inverse();
    final double[][] actualSolution = matrix.getDoubleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.inverse() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("inverse() returns the inversion of the matrix")
  void testInverseReturnsTheInversion() {
    printMethodName();
    final double[][] expectedInversion = sampleInvertibleMatrixInversion;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleInvertibleMatrixData);
    final double[][] actualInversion = matrix.inverse().getDoubleData();
    final ErrorSet actualErrors = findErrors(expectedInversion, actualInversion, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /*####################################################################################
  ## Tested method:
  ##   inverseAccurately()
  ## Behaviors to be tested:
  ## -- inverseAccurately() does not spoil the source matrix
  ## -- inverseAccurately() does not spoil the previously found vector solution
  ## -- inverseAccurately() does not spoil the previously found matrix solution
  ## -- inverseAccurately() returns the inversion of the matrix
  ## -- inverseAccurately() returns a more accurate inversion, than inverse()
  ######################################################################################*/

  @Test
  @DisplayName("inverseAccurately() does not spoil the source matrix")
  void testInverseAccuratelyDoesNotSpoilTheMatrix() {
    final double[][] expectedMatrixData = sampleSolvableMatrixData;

    final DoubleMatrix matrix = new DoubleMatrix(expectedMatrixData);
    matrix.inverseAccurately();
    final double[][] actualMatrixData = matrix.getDoubleData();

    assertThat(actualMatrixData).
        withFailMessage("matrix.inverse() spoils the source matrix data").
        isEqualTo(expectedMatrixData);
  }

  @Test
  @DisplayName("inverseAccurately() does not spoil the previously found vector solution")
  void testInverseAccuratelyDoesNotSpoilTheVectorSolution() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixVector);
    final double[] expectedSolution = matrix.getDoubleSolution();

    matrix.inverseAccurately();
    final double[] actualSolution = matrix.getDoubleSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.inverse() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("inverseAccurately() does not spoil the previously found matrix solution")
  void testInverseAccuratelyDoesNotSpoilTheMatrixSolution() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixB);
    final double[][] expectedSolution = matrix.getDoubleMatrixSolution();

    matrix.inverseAccurately();
    final double[][] actualSolution = matrix.getDoubleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.inverse() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("inverseAccurately() returns the inversion of the matrix")
  void testInverseAccuratelyReturnsTheInversion() {
    printMethodName();
    final double[][] expectedInversion = sampleInvertibleMatrixInversion;
    final double expectedError = ERROR_TOLERANCE;

    final DoubleMatrix matrix = new DoubleMatrix(sampleInvertibleMatrixData);
    final double[][] actualInversion = matrix.inverseAccurately().getDoubleData();
    final ErrorSet actualErrors = findErrors(expectedInversion, actualInversion, OUTPUT_ENABLED);

    assertThat(actualErrors.maxError()).
        withFailMessage("Max solution error, %.3e, is greater than expected error, %.3e", actualErrors.maxError(), expectedError).
        isLessThanOrEqualTo(expectedError);
    assertThat(actualErrors.mse()).
        withFailMessage("Solution mse, %.3e, is greater than expected error, %.3e", actualErrors.mse(), expectedError).
        isLessThanOrEqualTo(expectedError);
  }

  /**
   * Generates a series of random matrices (data samples), for each of them finds inversion with
   * inverse() method and another inversion using InverseAccurately(), computes errors for both inversions,
   * and estimates the minimum and average ratio of
   * accuracy improvement achieved by the {@link DoubleMatrix#inverseAccurately()} compared with the
   * {@link DoubleMatrix#inverse()}. Asserts that the minimum improvement ratio is greater that {@link #EXPECTED_MIN_IMPROVEMENT}
   * and the average improvement ratio is greater than {@link #EXPECTED_AVR_IMPROVEMENT}
   * @see #EXPECTED_MIN_IMPROVEMENT
   */
  @Test
  // @Disabled // TO DO 2023-03-12 17:54:28 Uncomment for higher execution speed while debugging
  @DisplayName("inverseAccurately() returns a more accurate inversion, than inverse()")
  void testInverseAccuratelyIsMoreAccurate() {
    printMethodName();
    final OperationComparator comparator = makeAccurateVsSimpleMatrixInversionComparator();

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
  ##   public Matrix transpose()
  ## Behaviors to be tested:
  ## -- transpose() does not spoil the source matrix
  ## -- transpose() does not spoil the previously found vector solution
  ## -- transpose() does not spoil the previously found matrix solution
  ## -- transpose() returns the transposition of the matrix
  ######################################################################################*/

  @Test
  @DisplayName("transpose() does not spoil the source matrix")
  void testTransposeDoesNotSpoilTheMatrix() {
    final double[][] expectedMatrixData = sampleSolvableMatrixData;

    final DoubleMatrix matrix = new DoubleMatrix(expectedMatrixData);
    matrix.transpose();
    final double[][] actualMatrixData = matrix.getDoubleData();

    assertThat(actualMatrixData).
        withFailMessage("matrix.transpose() spoils the source matrix data").
        isEqualTo(expectedMatrixData);
  }

  @Test
  @DisplayName("transpose() does not spoil the previously found vector solution")
  void testTransposeDoesNotSpoilTheVectorSolution() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixVector);
    final double[] expectedSolution = matrix.getDoubleSolution();

    matrix.transpose();
    final double[] actualSolution = matrix.getDoubleSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.transpose() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("transpose() does not spoil the previously found matrix solution")
  void testTransposeDoesNotSpoilTheMatrixSolution() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixB);
    final double[][] expectedSolution = matrix.getDoubleMatrixSolution();

    matrix.transpose();
    final double[][] actualSolution = matrix.getDoubleMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.transpose() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("transpose() returns the transposition of the matrix")
  void testTransposeReturnsTheTransposition() {
    final double[][] expectedTransposition = sampleInvertibleMatrixTransposition;

    final DoubleMatrix matrix = new DoubleMatrix(sampleInvertibleMatrixData);
    final double[][] actualTransposition = matrix.transpose().getDoubleData();

    assertThat(actualTransposition).
        withFailMessage("matrix.transpose() returns a wrong result").
        isEqualTo(expectedTransposition);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix unity()
  ## Behaviors to be tested:
  ## -- unity() returns a unity matrix of the same size
  ######################################################################################*/

  @Test
  @DisplayName("unity() returns a unity matrix of the same size")
  void testUnityReturnsUnityMatrix() {
    final DoubleMatrix matrix = new DoubleMatrix(anyMatrixDataNoMatterWhatItIs);
    final double[][] expectedUnity = unityMatrix(matrix.getSize());

    final double[][] actualUnity = matrix.unity().getDoubleData();

    assertThat(actualUnity).
        withFailMessage("The value returned by DoubleMatrix.unityMatrix() is not a unity matrix").
        isEqualTo(expectedUnity);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix multiply(Matrix factor)
  ## Behaviors to be tested:
  ## -- multiply(Matrix) with null argument throws NullPointerException
  ## -- multiply(Matrix) with wrong argument size throws IllegalArgumentException
  ## -- multiply(Matrix) returns a correct product
  ######################################################################################*/

  @Test
  @DisplayName("multiply(Matrix) with null argument throws NullPointerException")
  void testMultiplyByMatrixWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.multiply((Matrix)null),
        "Attempt to call multiply(Matrix) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Matrix)' and 'is null'").
        contains("multiply(matrix)").contains("is null");
  }

  @Test
  @DisplayName("multiply(Matrix) with wrong argument size throws IllegalArgumentException")
  void testMultiplyByMatrixWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(new DoubleMatrix(sampleSolvableMatrixTooLargeMatrix)),
        "Attempt to call multiply(Matrix) with different size of argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Matrix)' and 'must have the same size'").
        contains("multiply(matrix)").contains("must have the same size");
  }

  @Test
  @DisplayName("multiply(Matrix) returns a correct product")
  void testMultiplyByMatrixReturnsCorrectProduct() {
    final double[][] expectedProduct = sampleSolvableMatrixMatrixB;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[][] actualProduct = matrix.multiply(new DoubleMatrix(sampleSolvableMatrixMatrixX)).getDoubleData();

    assertThat(actualProduct).
        withFailMessage("matrix.multiply(Matrix) returns a wrong result").
        isEqualTo(expectedProduct);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix multiply(double[][] factor)
  ## Behaviors to be tested:
  ## -- multiply(double[][]) with null argument throws NullPointerException
  ## -- multiply(double[][]) with wrong argument size throws IllegalArgumentException
  ## -- multiply(double[][]) with non-square argument throws IllegalArgumentException
  ## -- multiply(double[][]) with an argument containing NaN throws IllegalArgumentException
  ## -- multiply(double[][]) with an argument containing Infinity throws IllegalArgumentException
  ## -- multiply(double[][]) returns a correct product
  ######################################################################################*/

  @Test
  @DisplayName("multiply(double[][]) with null argument throws NullPointerException")
  void testMultiplyByDoubleArrayWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.multiply((double[][])null),
        "Attempt to call multiply(double[][]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(double[][])' and 'is null'").
        contains("multiply(double[][])").contains("is null");
  }

  @Test
  @DisplayName("multiply(double[][]) with wrong argument size throws IllegalArgumentException")
  void testMultiplyByDoubleArrayWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(sampleSolvableMatrixTooLargeMatrix),
        "Attempt to call multiply(double[][]) with different size of argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(double[][])' and 'must have the same size'").
        contains("multiply(double[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("multiply(double[][]) with non-square argument throws IllegalArgumentException")
  void testMultiplyByDoubleArrayWithNonSquareArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(sampleSolvableMatrixNonSquareMatrix),
        "Attempt to call multiply(double[][]) with a non-square array as the argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(double[][])' and 'must be square'").
        contains("multiply(double[][])").contains("must be square");
  }

  @Test
  @DisplayName("multiply(double[][]) with an argument containing NaN throws IllegalArgumentException")
  void testMultiplyByDoubleArrayContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[][] badMatrix = deepCopyOf(sampleSolvableMatrixData);
    badMatrix[1][1] = Double.NaN;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badMatrix),
        "Attempt to call multiply(double[][]) with an argument containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(double[][])' and 'must not contain NaN or Infinity'").
        contains("multiply(double[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("multiply(double[][]) with an argument containing Infinity throws IllegalArgumentException")
  void testMultiplyByDoubleArrayContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[][] badMatrix = deepCopyOf(sampleSolvableMatrixData);
    badMatrix[1][1] = Double.NEGATIVE_INFINITY;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badMatrix),
        "Attempt to call multiply(double[][]) with an argument containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(double[][])' and 'must not contain NaN or Infinity'").
        contains("multiply(double[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("multiply(double[][]) returns a correct product")
  void testMultiplyByDoubleArrayReturnsCorrectProduct() {
    final double[][] expectedProduct = sampleSolvableMatrixMatrixB;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[][] actualProduct = matrix.multiply(sampleSolvableMatrixMatrixX).getDoubleData();

    assertThat(actualProduct).
        withFailMessage("matrix.multiply(double[][]) returns a wrong result").
        isEqualTo(expectedProduct);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix multiply(Number[][] factor) with Quadruple[][] argument
  ## Behaviors to be tested:
  ## -- multiply(Quadruple[][]) with null argument throws NullPointerException
  ## -- multiply(Quadruple[][]) with wrong argument size throws IllegalArgumentException
  ## -- multiply(Quadruple[][]) with non-square argument throws IllegalArgumentException
  ## -- multiply(Quadruple[][]) with an argument containing NaN throws IllegalArgumentException
  ## -- multiply(Quadruple[][]) with an argument containing Infinity throws IllegalArgumentException
  ## -- multiply(Quadruple[][]) with an argument containing null throws IllegalArgumentException
  ## -- multiply(Quadruple[][]) returns a correct product
  ######################################################################################*/

  @Test
  @DisplayName("multiply(Quadruple[][]) with null argument throws NullPointerException")
  void testMultiplyByQuadrupleArrayWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.multiply((Quadruple[][])null),
        "Attempt to call multiply(Quadruple[][]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Number[][])' and 'is null'").
        contains("multiply(number[][])").contains("is null");
  }

  @Test
  @DisplayName("multiply(Quadruple[][]) with wrong argument size throws IllegalArgumentException")
  void testMultiplyByQuadrupleArrayWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(sampleSolvableMatrixTooLargeMatrixAsQuadruples),
        "Attempt to call multiply(Quadruple[][]) with different size of argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Number[][])' and 'must have the same size'").
        contains("multiply(number[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("multiply(Quadruple[][]) with non-square argument throws IllegalArgumentException")
  void testMultiplyByQuadrupleArrayWithNonSquareArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(sampleSolvableMatrixNonSquareMatrixAsQuadruples),
        "Attempt to call multiply(Quadruple[][]) with a non-square array as the argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Number[][])' and 'must be square'").
        contains("multiply(number[][])").contains("must be square");
  }

  @Test
  @DisplayName("multiply(Quadruple[][]) with an argument containing NaN throws IllegalArgumentException")
  void testMultiplyByQuadrupleArrayContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsQuadruples);
    badMatrix[1][1] = Quadruple.nan();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badMatrix),
        "Attempt to call multiply(Quadruple[][]) with an argument containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(number[][])' and 'must not contain NaN or Infinity'").
        contains("multiply(number[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("multiply(Quadruple[][]) with an argument containing Infinity throws IllegalArgumentException")
  void testMultiplyByQuadrupleArrayContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsQuadruples);
    badMatrix[1][1] = Quadruple.negativeInfinity();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badMatrix),
        "Attempt to call multiply(Quadruple[][]) with an argument containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(number[][])' and 'must not contain NaN or Infinity'").
        contains("multiply(number[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("multiply(Quadruple[][]) with an argument containing null throws IllegalArgumentException")
  void testMultiplyByQuadrupleArrayContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsQuadruples);
    badMatrix[1][1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badMatrix),
        "Attempt to call multiply(Quadruple[][]) with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(number[][])' and 'must not contain null'").
        contains("multiply(number[][])").contains("must not contain null");
  }

  @Test
  @DisplayName("multiply(Quadruple[][]) returns a correct product")
  void testMultiplyByQuadrupleArrayReturnsCorrectProduct() {
    final double[][] expectedProduct = sampleSolvableMatrixMatrixB;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[][] actualProduct = matrix.multiply(sampleSolvableMatrixMatrixXAsQuadruples).getDoubleData();

    assertThat(actualProduct).
        withFailMessage("matrix.multiply(Quadruple[][]) returns a wrong result").
        isEqualTo(expectedProduct);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix multiply(Number[][] factor) with BigDecimal[][] argument
  ## Behaviors to be tested:
  ## -- multiply(BigDecimal[][]) with null argument throws NullPointerException
  ## -- multiply(BigDecimal[][]) with wrong argument size throws IllegalArgumentException
  ## -- multiply(BigDecimal[][]) with non-square argument throws IllegalArgumentException
  ## -- multiply(BigDecimal[][]) with an argument containing null throws IllegalArgumentException
  ## -- multiply(BigDecimal[][]) returns a correct product
  ######################################################################################*/

  @Test
  @DisplayName("multiply(BigDecimal[][]) with null argument throws NullPointerException")
  void testMultiplyByBigDecimalArrayWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.multiply((BigDecimal[][])null),
        "Attempt to call multiply(BigDecimal[][]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Number[][])' and 'is null'").
        contains("multiply(number[][])").contains("is null");
  }

  @Test
  @DisplayName("multiply(BigDecimal[][]) with wrong argument size throws IllegalArgumentException")
  void testMultiplyByBigDecimalArrayWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(sampleSolvableMatrixTooLargeMatrixAsBigDecimals),
        "Attempt to call multiply(BigDecimal[][]) with different size of argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Number[][])' and 'must have the same size'").
        contains("multiply(number[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("multiply(BigDecimal[][]) with non-square argument throws IllegalArgumentException")
  void testMultiplyByBigDecimalArrayWithNonSquareArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(sampleSolvableMatrixNonSquareMatrixAsBigDecimals),
        "Attempt to call multiply(BigDecimal[][]) with a non-square array as the argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Number[][])' and 'must be square'").
        contains("multiply(number[][])").contains("must be square");
  }

  @Test
  @DisplayName("multiply(BigDecimal[][]) with an argument containing null throws IllegalArgumentException")
  void testMultiplyByBigDecimalArrayContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final BigDecimal[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsBigDecimals);
    badMatrix[1][1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badMatrix),
        "Attempt to call multiply(BigDecimal[][]) with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Number[][])' and 'must not contain null'").
        contains("multiply(number[][])").contains("must not contain null");
  }

  @Test
  @DisplayName("multiply(BigDecimal[][]) returns a correct product")
  void testMultiplyByBigDecimalArrayReturnsCorrectProduct() {
    final double[][] expectedProduct = sampleSolvableMatrixMatrixB;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[][] actualProduct = matrix.multiply(sampleSolvableMatrixMatrixXAsBigDecimals).getDoubleData();

    assertThat(actualProduct).
        withFailMessage("matrix.multiply(Quadruple[][]) returns a wrong result").
        isEqualTo(expectedProduct);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Number[] multiply(double[] vector)
  ## Behaviors to be tested:
  ## -- multiply(double[]) with null argument throws NullPointerException
  ## -- multiply(double[]) with wrong argument size throws IllegalArgumentException
  ## -- multiply(double[]) with an argument containing NaN throws IllegalArgumentException
  ## -- multiply(double[]) with an argument containing Infinity throws IllegalArgumentException
  ## -- multiply(double[]) returns a correct product
  ######################################################################################*/

  @Test
  @DisplayName("multiply(double[]) with null argument throws NullPointerException")
  void testMultiplyByDoubleVectorWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.multiply((double[])null),
        "Attempt to call multiply(double[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(double[]' and 'is null'").
        contains("multiply(double[])").contains("is null");
  }

  @Test
  @DisplayName("multiply(double[]) with wrong argument size throws IllegalArgumentException")
  void testMultiplyByDoubleVectorWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(sampleTooLongMatrixVector),
        "Attempt to call multiply(double[]) with different size of argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(double[]' and 'the same size'").
        contains("multiply(double[]").contains("the same size");
  }

  @Test
  @DisplayName("multiply(double[]) with an argument containing NaN throws IllegalArgumentException")
  void testMultiplyByDoubleVectorContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[] badVector = sampleSolvableMatrixVector.clone();
    badVector[1] = Double.NaN;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badVector),
        "Attempt to call multiply(double[][]) with an argument containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(double[])' and 'must not contain NaN or Infinity'").
        contains("multiply(double[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("multiply(double[]) with an argument containing Infinity throws IllegalArgumentException")
  void testMultiplyByDoubleVectorContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[] badVector = sampleSolvableMatrixVector.clone();
    badVector[1] = Double.POSITIVE_INFINITY;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badVector),
        "Attempt to call multiply(double[][]) with an argument containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(double[])' and 'must not contain NaN or Infinity'").
        contains("multiply(double[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("multiply(double[]) returns a correct product")
  void testMultiplyByDoubleVectorReturnsCorrectProduct() {
    final double[] expectedProduct = sampleSolvableMatrixVector;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[] actualProduct = convertToDoubles(matrix.multiply(sampleSolvableMatrixSolution));

    assertThat(actualProduct).
        withFailMessage("matrix.multiply(double[]) returns a wrong result").
        isEqualTo(expectedProduct);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Number[] multiply(Number[] vector) with Quadruple[] argument
  ## Behaviors to be tested:
  ## -- multiply(Quadruple[]) with null argument throws NullPointerException
  ## -- multiply(Quadruple[]) with wrong argument size throws IllegalArgumentException
  ## -- multiply(Quadruple[]) with an argument containing NaN throws IllegalArgumentException
  ## -- multiply(Quadruple[]) with an argument containing Infinity throws IllegalArgumentException
  ## -- multiply(Quadruple[]) with an argument containing null throws IllegalArgumentException
  ## -- multiply(Quadruple[]) returns a correct product
  ######################################################################################*/

  @Test
  @DisplayName("multiply(Quadruple[]) with null argument throws NullPointerException")
  void testMultiplyByQuadrupleVectorWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.multiply((Quadruple[])null),
        "Attempt to call multiply(Quadruple[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Number[])' and 'is null'").
        contains("multiply(number[])").contains("is null");
  }

  @Test
  @DisplayName("multiply(Quadruple[]) with wrong argument size throws IllegalArgumentException")
  void testMultiplyByQuadrupleVectorWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(sampleTooLongMatrixVectorAsQuadruples),
        "Attempt to call multiply(Quadruple[]) with different size of argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Number[])' and 'must have the same size'").
        contains("multiply(number[])").contains("must have the same size");
  }

  @Test
  @DisplayName("multiply(Quadruple[]) with an argument containing NaN throws IllegalArgumentException")
  void testMultiplyByQuadrupleVectorContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixSolutionAsQuadruples);
    badVector[1] = Quadruple.nan();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badVector),
        "Attempt to call multiply(Quadruple[]) with an argument containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(number[])' and 'must not contain NaN or Infinity'").
        contains("multiply(number[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("multiply(Quadruple[]) with an argument containing Infinity throws IllegalArgumentException")
  void testMultiplyByQuadrupleVectorContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixSolutionAsQuadruples);
    badVector[1] = Quadruple.negativeInfinity();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badVector),
        "Attempt to call multiply(Quadruple[]) with an argument containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(number[])' and 'must not contain NaN or Infinity'").
        contains("multiply(number[])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("multiply(Quadruple[]) with an argument containing null throws IllegalArgumentException")
  void testMultiplyByQuadrupleVectorContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[] badVector = deepCopyOf(sampleSolvableMatrixSolutionAsQuadruples);
    badVector[1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badVector),
        "Attempt to call multiply(Quadruple[]) with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(number[])' and 'must not contain null'").
        contains("multiply(number[])").contains("must not contain null");
  }

  @Test
  @DisplayName("multiply(Quadruple[]) returns a correct product")
  void testMultiplyByQuadrupleVectorReturnsCorrectProduct() {
    final double[] expectedProduct = sampleSolvableMatrixVector;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[] actualProduct = convertToDoubles(matrix.multiply(sampleSolvableMatrixSolutionAsQuadruples));

    assertThat(actualProduct).
        withFailMessage("matrix.multiply(Quadruple[]) returns a wrong result").
        isEqualTo(expectedProduct);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Number[] multiply(Number[] vector) with BigDecimal[] argument
  ## Behaviors to be tested:
  ## -- multiply(BigDecimal[]) with null argument throws NullPointerException
  ## -- multiply(BigDecimal[]) with wrong argument size throws IllegalArgumentException
  ## -- multiply(BigDecimal[]) with an argument containing null throws IllegalArgumentException
  ## -- multiply(BigDecimal[]) returns a correct product
  ######################################################################################*/

  @Test
  @DisplayName("multiply(BigDecimal[]) with null argument throws NullPointerException")
  void testMultiplyByBigDecimalVectorWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.multiply((BigDecimal[])null),
        "Attempt to call multiply(BigDecimal[]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Number[])' and 'is null'").
        contains("multiply(number[])").contains("is null");
  }

  @Test
  @DisplayName("multiply(BigDecimal[]) with wrong argument size throws IllegalArgumentException")
  void testMultiplyByBigDecimalVectorWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(sampleTooLongMatrixVectorAsBigDecimals),
        "Attempt to call multiply(BigDecimal[]) with different size of argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Number[])' and 'must have the same size'").
        contains("multiply(number[])").contains("must have the same size");
  }

  @Test
  @DisplayName("multiply(BigDecimal[]) with an argument containing null throws IllegalArgumentException")
  void testMultiplyByBigDecimalVectorContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final BigDecimal[] badVector = deepCopyOf(sampleSolvableMatrixSolutionAsBigDecimals);
    badVector[1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(badVector),
        "Attempt to call multiply(BigDecimal[]) with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(number[])' and 'must not contain null'").
        contains("multiply(number[])").contains("must not contain null");
  }

  @Test
  @DisplayName("multiply(BigDecimal[]) returns a correct product")
  void testMultiplyByBigDecimalVectorReturnsCorrectProduct() {
    final double[] expectedProduct = sampleSolvableMatrixVector;

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);
    final double[] actualProduct = convertToDoubles(matrix.multiply(sampleSolvableMatrixSolutionAsBigDecimals));

    assertThat(actualProduct).
        withFailMessage("matrix.multiply(BigDecimal[]) returns a wrong result").
        isEqualTo(expectedProduct);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix multiply(double scalar)
  ## Behaviors to be tested:
  ## -- multiply(double) with NaN throws IllegalArgumentException
  ## -- multiply(double) with Infinity throws IllegalArgumentException
  ## -- multiply(double) returns a correct product
  ######################################################################################*/

  @Test
  @DisplayName("multiply(double) with NaN throws IllegalArgumentException")
  void testMultiplyByDoubleScalarWithNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(Double.NaN),
        "Attempt to call multiply(double) with NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(double)' and 'must not be NaN or Infinity'").
        contains("multiply(double)").contains("must not be nan or infinity");
  }

  @Test
  @DisplayName("multiply(double) with Infinity throws IllegalArgumentException")
  void testMultiplyByDoubleScalarWithInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(Double.POSITIVE_INFINITY),
        "Attempt to call multiply(diuble) with Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(double)' and 'must not be NaN or Infinity'").
        contains("multiply(double)").contains("must not be nan or infinity");
  }

  @Test
  @DisplayName("multiply(double) returns a correct product")
  void testMultiplyByDoubleScalarReturnsCorrectProduct() {
    final double[][] expectedProduct = sampleRandomMultipliedByFactor;

    final DoubleMatrix matrix = new DoubleMatrix(sampleRandomMatrixData);
    final double[][] actualProduct = matrix.multiply(DOUBLE_FACTOR).getDoubleData();

    assertThat(actualProduct).
        withFailMessage("multiply(double) returns a wrong result").
        isEqualTo(expectedProduct);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix multiply(Number scalar) with Quadruple
  ## Behaviors to be tested:
  ## -- multiply(Quadruple) with null argument throws an exception
  ## -- multiply(Quadruple) with NaN throws IllegalArgumentException
  ## -- multiply(Quadruple) with Infinity throws IllegalArgumentException
  ## -- multiply(Quadruple) returns a correct product
  ######################################################################################*/

  @Test
  @DisplayName("multiply(Quadruple) with null argument throws an exception")
  void testMultiplyByQuadrupleScalarWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.multiply((Quadruple)null),
        "Attempt to call multiply(Quadruple) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(number)' and 'is null'").
        contains("multiply(number)").contains("is null");
  }

  @Test
  @DisplayName("multiply(Quadruple) with NaN throws IllegalArgumentException")
  void testMultiplyByQuadrupleScalarWithNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(Quadruple.nan()),
        "Attempt to call multiply(Quadruple) with NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(number)' and 'must not be NaN or Infinity'").
        contains("multiply(number)").contains("must not be nan or infinity");
  }

  @Test
  @DisplayName("multiply(Quadruple) with Infinity throws IllegalArgumentException")
  void testMultiplyByQuadrupleScalarWithInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(Quadruple.negativeInfinity()),
        "Attempt to call multiply(Quadruple) with Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(number)' and 'must not be NaN or Infinity'").
        contains("multiply(number)").contains("must not be nan or infinity");
  }

  @Test
  @DisplayName("multiply(Quadruple) returns a correct product")
  void testMultiplyByQuadrupleScalarReturnsCorrectProduct() {
    final double[][] expectedProduct = sampleRandomMultipliedByFactor;

    final DoubleMatrix matrix = new DoubleMatrix(sampleRandomMatrixData);
    final double[][] actualProduct = matrix.multiply(new Quadruple(DOUBLE_FACTOR)).getDoubleData();

    assertThat(actualProduct).
        withFailMessage("matrix.multiply(Quadruple) returns a wrong result").
        isEqualTo(expectedProduct);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix multiply(Number scalar) with BigDecimal
  ## Behaviors to be tested:
  ## -- multiply(BigDecimal) with null argument throws an exception
  ## -- multiply(BigDecimal) returns a correct product
  ######################################################################################*/

  @Test
  @DisplayName("multiply(BigDecimal) with null argument throws an exception")
  void testMultiplyByBigDecimalScalarWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.multiply((BigDecimal)null),
        "Attempt to call multiply(BigDecimal) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(number)' and 'is null'").
        contains("multiply(number)").contains("is null");
  }

  @Test
  @DisplayName("multiply(BigDecimal) returns a correct product")
  void testMultiplyByBigDecimalScalarReturnsCorrectProduct() {
    final double[][] expectedProduct = sampleRandomMultipliedByFactor;

    final DoubleMatrix matrix = new DoubleMatrix(sampleRandomMatrixData);
    final double[][] actualProduct = matrix.multiply(BigDecimal.valueOf(DOUBLE_FACTOR)).getDoubleData();

    assertThat(actualProduct).
        withFailMessage("matrix.multiply(BigDecimal) returns a wrong result").
        isEqualTo(expectedProduct);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix add(Matrix matrixB)
  ## Behaviors to be tested:
  ## -- add(Matrix) with null argument throws NullPointerException
  ## -- add(Matrix) with wrong argument size throws IllegalArgumentException
  ## -- add(Matrix) returns a correct result
  ######################################################################################*/

  @Test
  @DisplayName("add(Matrix) with null argument throws NullPointerException")
  void testAddMatrixWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(summand1);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.add((Matrix)null),
        "Attempt to call add(Matrix) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add((Matrix)' and 'is null'").
        contains("add(matrix)").contains("is null");
  }

  @Test
  @DisplayName("add(Matrix) with wrong argument size throws IllegalArgumentException")
  void testAddMatrixWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(summand1);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(new DoubleMatrix(sampleSolvableMatrixTooLargeMatrix)),
        "Attempt to call add(Matrix) with argument of wrong size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(Matrix)' and 'must have the same size'").
        contains("add(matrix)").contains("must have the same size");
  }

  @Test
  @DisplayName("add(Matrix) returns a correct result")
  void testAddMatrixReturnsCorrectResult() {
    final double[][] expectedSum = summand1PlusSummand2;

    final DoubleMatrix matrix = new DoubleMatrix(summand1);
    final double[][] actualSum = matrix.add(new DoubleMatrix(summand2)).getDoubleData();

    assertThat(actualSum).
        withFailMessage("matrix.add(Matrix) returns a wrong result").
        isEqualTo(expectedSum);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix add(double[][] matrixB)
  ## Behaviors to be tested:
  ## -- add(double[][]) with null argument throws NullPointerException
  ## -- add(double[][]) with wrong argument size throws IllegalArgumentException
  ## -- add(double[][]) with non-square argument throws IllegalArgumentException
  ## -- add(double[][]) with an argument containing NaN throws IllegalArgumentException
  ## -- add(double[][]) with an argument containing Infinity throws IllegalArgumentException
  ## -- add(double[][]) returns a correct result
  ######################################################################################*/

  @Test
  @DisplayName("add(double[][]) with null argument throws NullPointerException")
  void testAddDoubleArrayWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(summand1);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.add((double[][])null),
        "Attempt to call add(double[][]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(double[][])' and 'is null'").
        contains("add(double[][])").contains("is null");
  }

  @Test
  @DisplayName("add(double[][]) with wrong argument size throws IllegalArgumentException")
  void testAddDoubleArrayWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(summand1);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(sampleSolvableMatrixTooLargeMatrix),
        "Attempt to call add(double[][]) with argument of wrong size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(double[][])' and 'must have the same size'").
        contains("add(double[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("add(double[][]) with non-square argument throws IllegalArgumentException")
  void testAddDoubleArrayWithNonSquareArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(summand1);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(sampleSolvableMatrixNonSquareMatrix),
        "Attempt to call add(double[][]) with a non-square array must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(double[][])' and 'must be square'").
        contains("add(double[][])").contains("must be square");
  }

  @Test
  @DisplayName("add(double[][]) with an argument containing NaN throws IllegalArgumentException")
  void testAddDoubleArrayContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[][] badMatrix = deepCopyOf(sampleSolvableMatrixData);
    badMatrix[1][1] = Double.NaN;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(badMatrix),
        "Attempt to call add(double[][]) with an argument containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(double[][])' and 'must not contain NaN or Infinity'").
        contains("add(double[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("add(double[][]) with an argument containing Infinity throws IllegalArgumentException")
  void testAddDoubleArrayContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[][] badMatrix = deepCopyOf(sampleSolvableMatrixData);
    badMatrix[1][1] = Double.NEGATIVE_INFINITY;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(badMatrix),
        "Attempt to call multiply(double[][]) with an argument containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(double[][])' and 'must not contain NaN or Infinity'").
        contains("add(double[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("add(double[][]) returns a correct result")
  void testAddDoubleArrayReturnsCorrectResult() {
    final double[][] expectedSum = summand1PlusSummand2;

    final DoubleMatrix matrix = new DoubleMatrix(summand1);
    final double[][] actualSum = matrix.add(summand2).getDoubleData();

    assertThat(actualSum).
        withFailMessage("matrix.add(double[][]) returns a wrong result").
        isEqualTo(expectedSum);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix add(Number[][] matrixB) with Quadruple[][] argument
  ## Behaviors to be tested:
  ## -- add(Quadruple[][]) with null argument throws NullPointerException
  ## -- add(Quadruple[][]) with wrong argument size throws IllegalArgumentException
  ## -- add(Quadruple[][]) with non-square argument throws IllegalArgumentException
  ## -- add(Quadruple[][]) with an argument containing NaN throws IllegalArgumentException
  ## -- add(Quadruple[][]) with an argument containing Infinity throws IllegalArgumentException
  ## -- add(Quadruple[][]) with an argument containing null throws IllegalArgumentException
  ## -- add(Quadruple[][]) returns a correct result
  ######################################################################################*/

  @Test
  @DisplayName("add(Quadruple[][]) with null argument throws NullPointerException")
  void testAddQuadrupleArrayWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(summand1);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.add((Quadruple[][])null),
        "Attempt to call add(Quadrupl[][]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(Number[][])' and 'is null'").
        contains("add(number[][])").contains("is null");
  }

  @Test
  @DisplayName("add(Quadruple[][]) with wrong argument size throws IllegalArgumentException")
  void testAddQuadrupleArrayWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(summand1);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(sampleSolvableMatrixTooLargeMatrixAsQuadruples),
        "Attempt to call add(Quadruple[][]) with argument of wrong size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(Number[][])' and 'must have the same size'").
        contains("add(number[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("add(Quadruple[][]) with non-square argument throws IllegalArgumentException")
  void testAddQuadrupleArrayWithNonSquareArrayThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(summand1);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(sampleSolvableMatrixNonSquareMatrixAsQuadruples),
        "Attempt to call add(Quadruple[][]) with a non-square array must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(Number[][])' and 'must be square'").
        contains("add(number[][])").contains("must be square");
  }

  @Test
  @DisplayName("add(Quadruple[][]) with an argument containing NaN throws IllegalArgumentException")
  void testAddQuadrupleArrayContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsQuadruples);
    badMatrix[1][1] = Quadruple.nan();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(badMatrix),
        "Attempt to call add(Quadruple[][]) with an argument containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(number[][])' and 'must not contain NaN or Infinity'").
        contains("add(number[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("add(Quadruple[][]) with an argument containing Infinity throws IllegalArgumentException")
  void testAddQuadrupleArrayContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsQuadruples);
    badMatrix[1][1] = Quadruple.negativeInfinity();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(badMatrix),
        "Attempt to call add(Quadruple[][]) with an argument containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(number[][])' and 'must not contain NaN or Infinity'").
        contains("add(number[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("add(Quadruple[][]) with an argument containing null throws IllegalArgumentException")
  void testAddQuadrupleArrayContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsQuadruples);
    badMatrix[1][1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(badMatrix),
        "Attempt to call add(Quadruple[][]) with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(number[][])' and 'must not contain null'").
        contains("add(number[][])").contains("must not contain null");
  }

  @Test
  @DisplayName("add(Quadruple[][]) returns a correct result")
  void testAddQuadrupleArrayReturnsCorrectResult() {
    final double[][] expectedSum = summand1PlusSummand2;

    final DoubleMatrix matrix = new DoubleMatrix(summand1);
    final double[][] actualSum = matrix.add(summand2AsQuadruples).getDoubleData();

    assertThat(actualSum).
        withFailMessage("matrix.add(Quadruple[][]) returns a wrong result").
        isEqualTo(expectedSum);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix add(Number[][] matrixB) with BigDecimal[][] argument
  ## Behaviors to be tested:
  ## -- add(BigDecimal[][]) with null argument throws NullPointerException
  ## -- add(BigDecimal[][]) with wrong argument size throws IllegalArgumentException
  ## -- add(BigDecimal[][]) with non-square argument throws IllegalArgumentException
  ## -- add(BigDecimal[][]) with an argument containing null throws IllegalArgumentException
  ## -- add(BigDecimal[][]) returns a correct result
  ######################################################################################*/

  @Test
  @DisplayName("add(BigDecimal[][]) with null argument throws NullPointerException")
  void testAddBigDecimalArrayWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(summand1);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.add((BigDecimal[][])null),
        "Attempt to call add(Quadrupl[][]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(Number[][])' and 'is null'").
        contains("add(number[][])").contains("is null");
  }

  @Test
  @DisplayName("add(BigDecimal[][]) with wrong argument size throws IllegalArgumentException")
  void testAddBigDecimalArrayWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(summand1);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(sampleSolvableMatrixTooLargeMatrixAsBigDecimals),
        "Attempt to call add(BigDecimal[][]) with argument of wrong size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(Number[][])' and 'must have the same size'").
        contains("add(number[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("add(BigDecimal[][]) with non-square argument throws IllegalArgumentException")
  void testAddBigDecimalArrayWithNonSquareArrayThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(summand1);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(sampleSolvableMatrixNonSquareMatrixAsBigDecimals),
        "Attempt to call add(BigDecimal[][]) with a non-square array must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(Number[][])' and 'must be square'").
        contains("add(number[][])").contains("must be square");
  }

  @Test
  @DisplayName("add(BigDecimal[][]) with an argument containing null throws IllegalArgumentException")
  void testAddBigDecimalArrayContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final BigDecimal[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsBigDecimals);
    badMatrix[1][1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(badMatrix),
        "Attempt to call add(BigDecimal[][]) with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(Number[][])' and 'must not contain null'").
        contains("add(number[][])").contains("must not contain null");
  }

  @Test
  @DisplayName("add(BigDecimal[][]) returns a correct result")
  void testAddBigDecimalArrayReturnsCorrectResult() {
    final double[][] expectedSum = summand1PlusSummand2;

    final DoubleMatrix matrix = new DoubleMatrix(summand1);
    final double[][] actualSum = matrix.add(summand2AsBigDecimals).getDoubleData();

    assertThat(actualSum).
        withFailMessage("matrix.add(BigDecimal[][]) returns a wrong result").
        isEqualTo(expectedSum);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix subtract(Matrix matrixB)
  ## Behaviors to be tested:
  ## -- subtract(Matrix) with null argument throws NullPointerException
  ## -- subtract(Matrix) with wrong argument size throws IllegalArgumentException
  ## -- subtract(Matrix) returns a correct result
  ######################################################################################*/

  @Test
  @DisplayName("subtract(Matrix) with null argument throws NullPointerException")
  void testSubtractMatrixWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(minuend);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.subtract((Matrix)null),
        "Attempt to call subtract(Matrix) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract((Matrix)' and 'is null'").
        contains("subtract(matrix)").contains("is null");
  }

  @Test
  @DisplayName("subtract(Matrix) with wrong argument size throws IllegalArgumentException")
  void testSubtractMatrixWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(minuend);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(new DoubleMatrix(sampleSolvableMatrixTooLargeMatrix)),
        "Attempt to call subtract(Matrix) with argument of wrong size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(Matrix)' and 'must have the same size'").
        contains("subtract(matrix)").contains("must have the same size");
  }

  @Test
  @DisplayName("subtract(Matrix) returns a correct result")
  void testSubtractMatrixReturnsCorrectResult() {
    final double[][] expectedDifference = minuendMinusSubtrahend;

    final DoubleMatrix matrix = new DoubleMatrix(minuend);
    final double[][] actualDifference = matrix.subtract(new DoubleMatrix(subtrahend)).getDoubleData();

    assertThat(actualDifference).
        withFailMessage("subtract(Matrix) returns a wrong result").
        isEqualTo(expectedDifference);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix subtract(double[][] matrixB)
  ## Behaviors to be tested:
  ## -- subtract(double[][]) with null argument throws NullPointerException
  ## -- subtract(double[][]) with wrong argument size throws IllegalArgumentException
  ## -- subtract(double[][]) with non-square argument throws IllegalArgumentException
  ## -- subtract(double[][])  with an argument containing NaN throws IllegalArgumentException
  ## -- subtract(double[][]) with an argument containing Infinity throws IllegalArgumentException
  ## -- subtract(double[][]) returns a correct result
  ######################################################################################*/

  @Test
  @DisplayName("subtract(double[][]) with null argument throws NullPointerException")
  void testSubtractDoubleArrayWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(minuend);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.subtract((double[][])null),
        "Attempt to call subtract(double[][]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(double[][])' and 'is null'").
        contains("subtract(double[][])").contains("is null");
  }

  @Test
  @DisplayName("subtract(double[][]) with wrong argument size throws IllegalArgumentException")
  void testSubtractDoubleArrayWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(minuend);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(sampleSolvableMatrixTooLargeMatrix),
        "Attempt to call subtract(double[][]) with argument of wrong size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(double[][])' and 'must have the same size'").
        contains("subtract(double[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("subtract(double[][]) with non-square argument throws IllegalArgumentException")
  void testSubtractDoubleArrayWithNonSquareArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(minuend);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(sampleSolvableMatrixNonSquareMatrix),
        "Attempt to call subtract(double[][]) with a non-square array must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(double[][])' and 'must be square'").
        contains("subtract(double[][])").contains("must be square");
  }

  @Test
  @DisplayName("subtract(double[][])  with an argument containing NaN throws IllegalArgumentException")
  void testSubtractDoubleArrayContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[][] badMatrix = deepCopyOf(sampleSolvableMatrixData);
    badMatrix[1][1] = Double.NaN;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(badMatrix),
        "Attempt to call subtract(double[][]) with an argument containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(double[][])' and 'must not contain NaN or Infinity'").
        contains("subtract(double[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("subtract(double[][]) with an argument containing Infinity throws IllegalArgumentException")
  void testSubtractDoubleArrayContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final double[][] badMatrix = deepCopyOf(sampleSolvableMatrixData);
    badMatrix[1][1] = Double.NEGATIVE_INFINITY;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(badMatrix),
        "Attempt to call subtract(double[][]) with an argument containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(double[][])' and 'must not contain NaN or Infinity'").
        contains("subtract(double[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("subtract(double[][]) returns a correct result")
  void testSubtractDoubleArrayReturnsCorrectResult() {
    final double[][] expectedDifference = minuendMinusSubtrahend;

    final DoubleMatrix matrix = new DoubleMatrix(minuend);
    final double[][] actualDifference = matrix.subtract(subtrahend).getDoubleData();

    assertThat(actualDifference).
        withFailMessage("matrix.subtract(double[][]) returns a wrong result").
        isEqualTo(expectedDifference);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix subtract(Number[][] matrixB) with Quadruple[][] argument
  ## Behaviors to be tested:
  ## -- subtract(Quadruple[][]) with null argument throws NullPointerException
  ## -- subtract(Quadruple[][]) with wrong argument size throws IllegalArgumentException
  ## -- subtract(Quadruple[][]) with non-square argument throws IllegalArgumentException
  ## -- subtract(Quadruple[][]) with an argument containing NaN throws IllegalArgumentException
  ## -- subtract(Quadruple[][]) with an argument containing Infinity throws IllegalArgumentException
  ## -- subtract(Quadruple[][]) with an argument containing null throws IllegalArgumentException
  ## -- subtract(Quadruple[][]) returns a correct result
  ######################################################################################*/

  @Test
  @DisplayName("subtract(Quadruple[][]) with null argument throws NullPointerException")
  void testSubtractQuadrupleArrayWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(minuend);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.subtract((Quadruple[][])null),
        "Attempt to call subtract(Quadruple[][]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(Number[][])' and 'is null'").
        contains("subtract(number[][])").contains("is null");
  }

  @Test
  @DisplayName("subtract(Quadruple[][]) with wrong argument size throws IllegalArgumentException")
  void testSubtractQuadrupleArrayWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(minuend);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(sampleSolvableMatrixTooLargeMatrixAsQuadruples),
        "Attempt to call subtract(Quadruple[][]) with argument of wrong size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(Number[][])' and 'must have the same size'").
        contains("subtract(number[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("subtract(Quadruple[][]) with non-square argument throws IllegalArgumentException")
  void testSubtractQuadrupleArrayWithNonSquareArrayThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(minuend);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(sampleSolvableMatrixNonSquareMatrixAsQuadruples),
        "Attempt to call subtract(Quadruple[][]) with a non-square array must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(Number[][])' and 'must be square'").
        contains("subtract(number[][])").contains("must be square");
  }

  @Test
  @DisplayName("subtract(Quadruple[][]) with an argument containing NaN throws IllegalArgumentException")
  void testSubtractQuadrupleArrayContainingNaNThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsQuadruples);
    badMatrix[1][1] = Quadruple.nan();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(badMatrix),
        "Attempt to call subtract(Quadruple[][]) with an argument containing NaN must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(number[][])' and 'must not contain NaN or Infinity'").
        contains("subtract(number[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("subtract(Quadruple[][]) with an argument containing Infinity throws IllegalArgumentException")
  void testSubtractQuadrupleArrayContainingInfinityThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsQuadruples);
    badMatrix[1][1] = Quadruple.negativeInfinity();
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(badMatrix),
        "Attempt to call subtract(Quadruple[][]) with an argument containing Infinity must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(number[][])' and 'must not contain NaN or Infinity'").
        contains("subtract(number[][])").contains("must not contain nan or infinity");
  }

  @Test
  @DisplayName("subtract(Quadruple[][]) with an argument containing null throws IllegalArgumentException")
  void testSubtractQuadrupleArrayContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final Quadruple[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsQuadruples);
    badMatrix[1][1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(badMatrix),
        "Attempt to call subtract(Quadruple[][]) with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(number[][])' and 'must not contain null'").
        contains("subtract(number[][])").contains("must not contain null");
  }

  @Test
  @DisplayName("subtract(Quadruple[][]) returns a correct result")
  void testSubtractQuadrupleArrayReturnsCorrectResult() {
    final double[][] expectedDifference = minuendMinusSubtrahend;

    final DoubleMatrix matrix = new DoubleMatrix(minuend);
    final double[][] actualDifference = matrix.subtract(subtrahendAsQuadruples).getDoubleData();

    assertThat(actualDifference).
        withFailMessage("matrix.subtract(Number[][]) returns a wrong result").
        isEqualTo(expectedDifference);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Matrix subtract(Number[][] matrixB) with BigDecimal[][] argument
  ## Behaviors to be tested:
  ## -- subtract(BigDecimal[][]) with null argument throws NullPointerException
  ## -- subtract(BigDecimal[][]) with wrong argument size throws IllegalArgumentException
  ## -- subtract(BigDecimal[][]) with non-square argument throws IllegalArgumentException
  ## -- subtract(BigDecimal[][]) with an argument containing null throws IllegalArgumentException
  ## -- subtract(BigDecimal[][]) returns a correct result
  ######################################################################################*/

  @Test
  @DisplayName("subtract(BigDecimal[][]) with null argument throws NullPointerException")
  void testSubtractBigDecimalArrayWithNullArgumentThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(minuend);

    final Throwable expectedException = assertThrows(NullPointerException.class,
        () -> matrix.subtract((BigDecimal[][])null),
        "Attempt to call subtract(BigDecimal[][]) with null argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(Number[][])' and 'is null'").
        contains("subtract(number[][])").contains("is null");
  }

  @Test
  @DisplayName("subtract(BigDecimal[][]) with wrong argument size throws IllegalArgumentException")
  void testSubtractBigDecimalArrayWithWrongArgumentSizeThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(minuend);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(sampleSolvableMatrixTooLargeMatrixAsBigDecimals),
        "Attempt to call subtract(BigDecimal[][]) with argument of wrong size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(Number[][])' and 'must have the same size'").
        contains("subtract(number[][])").contains("must have the same size");
  }

  @Test
  @DisplayName("subtract(BigDecimal[][]) with non-square argument throws IllegalArgumentException")
  void testSubtractBigDecimalArrayWithNonSquareArrayThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(minuend);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(sampleSolvableMatrixNonSquareMatrixAsBigDecimals),
        "Attempt to call subtract(BigDecimal[][]) with a non-square array must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(Number[][])' and 'must be square'").
        contains("subtract(number[][])").contains("must be square");
  }

  @Test
  @DisplayName("subtract(BigDecimal[][]) with an argument containing null throws IllegalArgumentException")
  void testSubtractBigDecimalArrayContainingNullThrowsException() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixData);

    final BigDecimal[][] badMatrix = deepCopyOf(sampleSolvableMatrixMatrixXAsBigDecimals);
    badMatrix[1][1] = null;
    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(badMatrix),
        "Attempt to call subtract(BigDecimal[][]) with an argument containing null must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(Number[][])' and 'must not contain null'").
        contains("subtract(number[][])").contains("must not contain null");
  }

  @Test
  @DisplayName("subtract(BigDecimal[][]) returns a correct result")
  void testSubtractBigDecimalArrayReturnsCorrectResult() {
    final double[][] expectedDifference = minuendMinusSubtrahend;

    final DoubleMatrix matrix = new DoubleMatrix(minuend);
    final double[][] actualDifference = matrix.subtract(subtrahendAsBigDecimals).getDoubleData();

    assertThat(actualDifference).
        withFailMessage("matrix.subtract(Number[][]) returns a wrong result").
        isEqualTo(expectedDifference);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Number determinant()
  ## Behaviors to be tested:
  ## -- determinant() for a matrix with a positive determinant returns a correct value
  ## -- determinant() for a matrix with a negative determinant returns a correct value
  ## -- determinant() for an inconsistent matrix returns zero
  ## -- determinant() for an underdetermined matrix returns zero
  ######################################################################################*/

  @Test
  @DisplayName("determinant() for a matrix with a positive determinant returns a correct value")
  void testDeterminantWhenPositiveReturnsCorrectValue() {
    final double expectedDeterminant = sampleSolvableMatrixAnotherMatrixXDeterminant; // 1569.252798052159
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE .. ERROR_TOLERANCE]
    final Offset<Double> delta = Offset.offset(Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant * ERROR_TOLERANCE)));

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixAnotherMatrixX);
    final double actualDeterminant = matrix.determinant().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinant() for a matrix with a positive determinant returns a wrong value").
        isCloseTo(expectedDeterminant, delta);
  }

  @Test
  @DisplayName("determinant() for a matrix with a negative determinant returns a correct value")
  void testDeterminantWhenNegativeReturnsCorrectValue() {
    final double expectedDeterminant = sampleSolvableMatrixMatrixXDeterminant;  // -185.0
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE..ACCEPTED_ERROR]
    final Offset<Double> delta = Offset.offset(Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant * ERROR_TOLERANCE)));

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixMatrixX);
    final double actualDeterminant = matrix.determinant().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinant() for a matrix with a negative determinant returns a wrong value").
        isCloseTo(expectedDeterminant, delta);
  }

  @Test
  @DisplayName("determinant() for an inconsistent matrix returns zero")
  void testDeterminantForInconsistentMatrixReturnsZero() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    final double actualDeterminant = matrix.determinant().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinant() for an inconsistent matrix returns a non-zero result").
        isZero();
  }

  @Test
  @DisplayName("determinant() for an underdetermined matrix returns zero")
  void testDeterminantForUnderdeterminedMatrixReturnsZero() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleUnderdeterminedMatrixData);

    final double actualDeterminant = matrix.determinant().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinant() for an underdetermined matrix returns a non-zero result").
        isZero();
  }

  /*####################################################################################
  ## Tested method:
  ##   public double determinantAsDouble()
  ## Behaviors to be tested:
  ## -- determinantAsDouble() for a matrix with a positive determinant returns a correct value
  ## -- determinantAsDouble() for a matrix with a negative determinant returns a correct value
  ## -- determinantAsDouble() for an inconsistent matrix returns zero
  ## -- determinantAsDouble() for an underdetermined matrix returns zero
  ######################################################################################*/

  @Test
  @DisplayName("determinantAsDouble() for a matrix with a positive determinant returns a correct value")
  void testDeterminantAsDoubleWhenPositiveReturnsCorrectValue() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixAnotherMatrixX);
    final double expectedDeterminant = sampleSolvableMatrixAnotherMatrixXDeterminant; // 1569.252798052159
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE .. ERROR_TOLERANCE]
    final Offset<Double> delta = Offset.offset(Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant * ERROR_TOLERANCE)));

    final double actualDeterminant = matrix.determinantAsDouble();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsDouble() for a matrix with a positive determinant returns a wrong value").
        isCloseTo(expectedDeterminant, delta);
  }

  @Test
  @DisplayName("determinantAsDouble() for a matrix with a negative determinant returns a correct value")
  void testDeterminantAsDoubleWhenNegativeReturnsCorrectValue() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixMatrixX);
    final double expectedDeterminant = sampleSolvableMatrixMatrixXDeterminant;  // -185.0
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE..ACCEPTED_ERROR]
    final Offset<Double> delta = Offset.offset(Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant * ERROR_TOLERANCE)));

    final double actualDeterminant = matrix.determinantAsDouble();
    assertThat(actualDeterminant).
        withFailMessage("determinantAsDouble() for a matrix with a negative determinant returns a wrong value").
        isCloseTo(expectedDeterminant, delta);
  }

  @Test
  @DisplayName("determinantAsDouble() for an inconsistent matrix returns zero")
  void testDeterminantAsDoubleForInconsistentMatrixReturnsZero() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    final double actualDeterminant = matrix.determinantAsDouble();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsDouble() for an inconsistent matrix returns a non-zero result").
        isZero();
  }

  @Test
  @DisplayName("determinantAsDouble() for an underdetermined matrix returns zero")
  void testDeterminantAsDoubleForUnderdeterminedMatrixReturnsZero() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleUnderdeterminedMatrixData);

    final double actualDeterminant = matrix.determinantAsDouble();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsDouble() for an underdetermined matrix returns a non-zero result").
        isZero();
  }

  /*####################################################################################
  ## Tested method:
  ##   public Quadruple determinantAsQuadruple()
  ## Behaviors to be tested:
  ## -- determinantAsQuadruple() for a matrix with a positive determinant returns a correct value
  ## -- determinantAsQuadruple() for a matrix with a negative determinant returns a correct value
  ## -- determinantAsQuadruple() for an inconsistent matrix returns zero
  ## -- determinantAsQuadruple() for an underdetermined matrix returns zero
  ######################################################################################*/

  @Test
  @DisplayName("determinantAsQuadruple() for a matrix with a positive determinant returns a correct value")
  void testdeterminantAsQuadrupleWhenPositiveReturnsCorrectValue() {
    final double expectedDeterminant = sampleSolvableMatrixAnotherMatrixXDeterminant; // 1569.252798052159
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE .. ERROR_TOLERANCE]
    final Offset<Double> delta = Offset.offset(Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant * ERROR_TOLERANCE)));

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixAnotherMatrixX);
    final double actualDeterminant = matrix.determinantAsQuadruple().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsQuadruple() for a matrix with a positive determinant returns a wrong value").
        isCloseTo(expectedDeterminant, delta);
  }

  @Test
  @DisplayName("determinantAsQuadruple() for a matrix with a negative determinant returns a correct value")
  void testdeterminantAsQuadrupleWhenNegativeReturnsCorrectValue() {
    final double expectedDeterminant = sampleSolvableMatrixMatrixXDeterminant;  // -185.0
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE..ACCEPTED_ERROR]
    final Offset<Double> delta = Offset.offset(Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant * ERROR_TOLERANCE)));

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixMatrixX);
    final double actualDeterminant = matrix.determinantAsQuadruple().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsQuadruple() for a matrix with a negative determinant returns a wrong value").
        isCloseTo(expectedDeterminant, delta);
  }

  @Test
  @DisplayName("determinantAsQuadruple() for an inconsistent matrix returns zero")
  void testdeterminantAsQuadrupleForInconsistentMatrixReturnsZero() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    final double actualDeterminant = matrix.determinantAsQuadruple().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsQuadruple() for an inconsistent matrix returns a non-zero result").
        isZero();
  }

  @Test
  @DisplayName("determinantAsQuadruple() for an underdetermined matrix returns zero")
  void testdeterminantAsQuadrupleForUnderdeterminedMatrixReturnsZero() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleUnderdeterminedMatrixData);

    final double actualDeterminant = matrix.determinantAsQuadruple().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsQuadruple() for an underdetermined matrix returns a non-zero result").
        isZero();
  }

  /*####################################################################################
  ## Tested method:
  ##   public BigDecimal determinantAsBigDecimal()
  ## Behaviors to be tested:
  ## -- determinantAsBigDecimal() for a matrix with a positive determinant returns a correct value
  ## -- determinantAsBigDecimal() for a matrix with a negative determinant returns a correct value
  ## -- determinantAsBigDecimal() for an inconsistent matrix returns zero
  ## -- determinantAsBigDecimal() for an underdetermined matrix returns zero
  ######################################################################################*/

  @Test
  @DisplayName("determinantAsBigDecimal() for a matrix with a positive determinant returns a correct value")
  void testdeterminantAsBigDecimalWhenPositiveReturnsCorrectValue() {
    final double expectedDeterminant = sampleSolvableMatrixAnotherMatrixXDeterminant; // 1569.252798052159
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE .. ERROR_TOLERANCE]
    final Offset<Double> delta = Offset.offset(Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant * ERROR_TOLERANCE)));

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixAnotherMatrixX);
    final double actualDeterminant = matrix.determinantAsBigDecimal().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsBigDecimal() for a matrix with a positive determinant returns a wrong value").
        isCloseTo(expectedDeterminant, delta);
  }

  @Test
  @DisplayName("determinantAsBigDecimal() for a matrix with a negative determinant returns a correct value")
  void testdeterminantAsBigDecimalWhenNegativeReturnsCorrectValue() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixMatrixX);
    final double expectedDeterminant = sampleSolvableMatrixMatrixXDeterminant;  // -185.0
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE..ACCEPTED_ERROR]
    final Offset<Double> delta = Offset.offset(Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant * ERROR_TOLERANCE)));

    final double actualDeterminant = matrix.determinantAsBigDecimal().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsBigDecimal() for a matrix with a negative determinant returns a wrong value").
        isCloseTo(expectedDeterminant, delta);
  }

  @Test
  @DisplayName("determinantAsBigDecimal() for an inconsistent matrix returns zero")
  void testdeterminantAsBigDecimalForInconsistentMatrixReturnsZero() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);

    final double actualDeterminant = matrix.determinantAsBigDecimal().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsBigDecimal() for an inconsistent matrix returns a non-zero result").
        isZero();
  }

  @Test
  @DisplayName("determinantAsBigDecimal() for an underdetermined matrix returns zero")
  void testdeterminantAsBigDecimalForUnderdeterminedMatrixReturnsZero() {
    final DoubleMatrix matrix = new DoubleMatrix(sampleUnderdeterminedMatrixData);

    final double actualDeterminant = matrix.determinantAsBigDecimal().doubleValue();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsBigDecimal() for an underdetermined matrix returns a non-zero result").
        isZero();
  }

  /*####################################################################################
  ## Tested method:
  ##   public Number norm()
  ## Behaviors to be tested:
  ## -- norm() for a matrix returns a correct value
  ## -- scaling does not affect the returned norm() value
  ######################################################################################*/

  @Test
  @DisplayName("norm() for a matrix returns a correct value")
  void testNormReturnsCorrectValue() {
    final double expectedNorm = sampleSolvableMatrixAnotherMatrixXNorm; // 1569.252798052159

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixAnotherMatrixX);
    final double actualNorm = matrix.norm().doubleValue();

    assertThat(actualNorm).
        withFailMessage("norm() for a matrix returns a wrong value").
        isEqualTo(expectedNorm);
  }

  @Test
  @DisplayName("scaling does not affect the returned norm() value")
  void testNormDoesNotDependOnScale() {
    final double expectedNorm = sampleSolvableMatrixAnotherMatrixXNorm; // 1569.252798052159

    final DoubleMatrix matrix = new DoubleMatrix(sampleSolvableMatrixAnotherMatrixX, true); // to be scaled
    matrix.solve(sampleSolvableMatrixMatrixB);                                               // force scaling
    final double actualNorm = matrix.norm().doubleValue();

    assertThat(actualNorm).
        withFailMessage("after scaling, norm() returns a wrong value").
        isEqualTo(expectedNorm);
  }

  /*####################################################################################
  ## Tested method:
  ##   public Number cond()
  ## Behaviors to be tested:
  ## -- cond() for an invertible matrix returns a correct value
  ## -- cond() for a non-invertible matrix returns a correct value
  ######################################################################################*/

  @Test
  @DisplayName("cond() for an invertible matrix returns a correct value")
  void testCondInvertibleReturnsCorrectValue() {
    final double expectedCond = sampleInvertibleMatrixConditionNumber; // 82.5
    final Offset<Double> delta = Offset.offset(expectedCond * ERROR_TOLERANCE);

    final DoubleMatrix matrix = new DoubleMatrix(sampleInvertibleMatrixData);
    final double actualCond = matrix.cond();

    assertThat(actualCond).
        withFailMessage("expected cond() value is %s, and actual value is %s", expectedCond, actualCond).
        isCloseTo(expectedCond, delta);
  }

  @Test
  @DisplayName("cond() for a non-invertible matrix returns Double.POSITIVE_INFINITY")
  void testCondNonInvertibleReturnsCorrectValue() {
    final double expectedCond = Double.POSITIVE_INFINITY;

    final DoubleMatrix matrix = new DoubleMatrix(sampleInconsistentMatrixData);
    final double actualCond = matrix.cond();

    assertThat(actualCond).
        withFailMessage("expected cond() value is %s, and actual value is %s", expectedCond, actualCond).
        isEqualTo(expectedCond);
  }


/* ******************************************************************************************************
 *** Private methods ************************************************************************************
 ********************************************************************************************************/

  private static OperationComparator makeAccurateVsSimpleMatrixInversionComparator() {
    return new OperationComparator(
        () -> MatrixData.makeDataSetForInversions(LARGE_MATRIX_SIZE, random),
        MatrixData::matrixInversionErrors,
        MatrixData::accurateMatrixInversionErrors
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
