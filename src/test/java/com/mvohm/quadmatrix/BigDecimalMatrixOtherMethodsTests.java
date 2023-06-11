package com.mvohm.quadmatrix;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
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
 * Tests for {@link BigDecimalMatrix} methods for inversion, transposition, matrix arithmetics etc.
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

public class BigDecimalMatrixOtherMethodsTests {

  /** Enables progress indication for time-consuming tests and printing values of some computational errors */
  private static final boolean OUTPUT_ENABLED = true;
  private static final boolean DETAILED_OUTPUT = true;

  private static final int    RAND_SEED
//                                          = 123;
                                          = 12121234;
//                                          = 321;

  private static final double         ERROR_TOLERANCE                                  = 1e-36;
  private static final MathContext    MC_60                                            = new MathContext(60, RoundingMode.HALF_EVEN);

  private static final double         ANY_DOUBLE_VALUE                                 = Math.E;
  private static final double         DOUBLE_FACTOR                                    = ANY_DOUBLE_VALUE;
  private static final Quadruple      QUADRUPLE_FACTOR                                 = new Quadruple(DOUBLE_FACTOR);

  // To perform iterations to evaluate the accuracy of inverseAccurately()
  // With error estimation based on Quadruple calculations, it takes about 1 minute
  // to run testInverseAccuratelyIsAccurate() with 60 iterations and matrices 157x157 on a i5 3.2 HHz machine
  private static final int            REFINEMENT_TEST_ITERATIONS                        = 30; // Approx. 60 s
  private static final int            LARGE_MATRIX_SIZE                                 = 50; // So that accurate inversion test takes approx. 2 s on my machine

  private static final double         EXPECTED_MIN_IMPROVEMENT                          =   4.0;
  private static final double         EXPECTED_AVR_IMPROVEMENT                          =   7.0;

  // Data for test matrices
  private static final double[][]     sampleRandomMatrixData                            = sampleRandomMatrixData();

  private static final double[][]     sampleSolvableMatrixData                          = sampleSolvableMatrixData();
  private static final BigDecimal[][] sampleSolvableMatrixDataAsBigDecimals             = sampleSolvableMatrixDataAsBigDecimals();

  private static final double[][]     sampleSolvableMatrixMatrixB                       = sampleSolvableMatrixMatrixB();
  private static final Quadruple[][]  sampleSolvableMatrixMatrixBAsQuadruples           = sampleSolvableMatrixMatrixBAsQuadruples();

  private static final double[][]     sampleSolvableMatrixMatrixX                       = sampleSolvableMatrixMatrixX();
  private static final Quadruple[][]  sampleSolvableMatrixMatrixXAsQuadruples           = convertToQuadruples(sampleSolvableMatrixMatrixX);
  private static final BigDecimal[][] sampleSolvableMatrixMatrixXAsBigDecimals          = convertToBigDecimals(sampleSolvableMatrixMatrixX);

  private static final double[]       sampleSolvableMatrixVector                        = sampleSolvableMatrixVector();
  private static final Quadruple[]    sampleSolvableMatrixVectorAsQuadruples            = convertToQuadruples(sampleSolvableMatrixVector);

  private static final double[]       sampleSolvableMatrixSolution                      = sampleSolvableMatrixSolution();
  private static final Quadruple[]    sampleSolvableMatrixSolutionAsQuadruples          = convertToQuadruples(sampleSolvableMatrixSolution);
  private static final BigDecimal[]   sampleSolvableMatrixSolutionAsBigDecimals         = convertToBigDecimals(sampleSolvableMatrixSolution);

  private static final double[]       sampleTooLongMatrixVector                         = sampleTooLongMatrixVector();
  private static final Quadruple[]    sampleTooLongMatrixVectorAsQuadruples             = convertToQuadruples(sampleTooLongMatrixVector);
  private static final BigDecimal[]   sampleTooLongMatrixVectorAsBigDecimals            = convertToBigDecimals(sampleTooLongMatrixVector);

  private static final double[][]     sampleInvertibleMatrixData                        = sampleInvertibleMatrixData();
  private static final BigDecimal[][] sampleInvertibleMatrixInversionAsBigDecimals     = sampleInvertibleMatrixInversionAsBigDecimals();

  private static final double[][]     sampleInvertibleMatrixTransposition               = sampleInvertibleMatrixTransposition();

  private static final double[][]     sampleSolvableMatrixTooLargeMatrix                = sampleSolvableMatrixTooLargeMatrix();
  private static final Quadruple[][]  sampleSolvableMatrixTooLargeMatrixAsQuadruples    = convertToQuadruples(sampleSolvableMatrixTooLargeMatrix);
  private static final BigDecimal[][] sampleSolvableMatrixTooLargeMatrixAsBigDecimals   = convertToBigDecimals(sampleSolvableMatrixTooLargeMatrix);

  private static final double[][]     sampleSolvableMatrixNonSquareMatrix               = sampleSolvableMatrixNonSquareMatrix();
  private static final Quadruple[][]  sampleSolvableMatrixNonSquareMatrixAsQuadruples   = convertToQuadruples(sampleSolvableMatrixNonSquareMatrix);
  private static final BigDecimal[][] sampleSolvableMatrixNonSquareMatrixAsBigDecimals  = convertToBigDecimals(sampleSolvableMatrixNonSquareMatrix);

  private static final double[][]     sampleInconsistentMatrixData                      = sampleInconsistentMatrixData();
  private static final double[][]     sampleUnderdeterminedMatrixData                   = sampleUnderdeterminedMatrixData();

  private static final double[][]     anyMatrixDataNoMatterWhatItIs                     = randomMatrix(LARGE_MATRIX_SIZE);

  private static final double[][]     summand1                                          = summand1();
  private static final Quadruple[][]  summand1AsQuadruples                              = summand1AsQuadruples();

  private static final double[][]     summand2                                          = summand2();
  private static final Quadruple[][]  summand2AsQuadruples                              = summand2AsQuadruples();
  private static final BigDecimal[][] summand2AsBigDecimals                             = convertToBigDecimals(summand2AsQuadruples);

  private static final Quadruple[][]  summand1PlusSummand2AsQuadruples                  = summand1PlusSummand2AsQuadruples();
  private static final BigDecimal[][] summand1PlusSummand2AsBigDecimals                 = summand1PlusSummand2AsBigDecimals();

  private static final double[][]     minuend                                           = minuend();

  private static final Quadruple[][]  sampleSpdSolvableMatrixDataAsQuadruples           = sampleSpdSolvableMatrixDataAsQuadruples();

  private static Random               random;

  private static double[][]           sampleRandomMultipliedByFactor;
  private static Quadruple[][]        sampleRandomMultipliedByFactorAsQuadruples;
  private static BigDecimal[][]       sampleRandomMultipliedByQuadrupleFactorAsBigDecimals;

  private static MathContext          MC = new MathContext(BigDecimalMatrix.getDefaultPrecision(), RoundingMode.HALF_EVEN);

  @BeforeAll
  static void setup() {
    Locale.setDefault(Locale.US);
    BigDecimalMatrix.setDefaultScaling(true);
    BigDecimalMatrix.setDefaultPrecision(40);

    sampleRandomMultipliedByFactor = deepCopyOf(sampleRandomMatrixData);
    final int size = sampleRandomMultipliedByFactor.length;
    sampleRandomMultipliedByFactorAsQuadruples = new Quadruple[size][size];
    sampleRandomMultipliedByQuadrupleFactorAsBigDecimals = new BigDecimal[size][size];

    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        sampleRandomMultipliedByFactor[i][j] =
            BigDecimal.valueOf(sampleRandomMatrixData[i][j]).multiply(BigDecimal.valueOf(DOUBLE_FACTOR), MC).
            doubleValue();
        sampleRandomMultipliedByFactorAsQuadruples[i][j] =
            new Quadruple(BigDecimal.valueOf(sampleRandomMatrixData[i][j]).multiply(QUADRUPLE_FACTOR.bigDecimalValue(), MC));
        sampleRandomMultipliedByQuadrupleFactorAsBigDecimals[i][j] =
            BigDecimal.valueOf(sampleRandomMatrixData[i][j]).multiply(QUADRUPLE_FACTOR.bigDecimalValue(), MC);
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
    final BigDecimal[][] expectedMatrixData = sampleSolvableMatrixDataAsBigDecimals;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(expectedMatrixData);
    matrix.inverse();
    final BigDecimal[][] actualMatrixData = matrix.getBigDecimalData();

    assertThat(actualMatrixData).
        withFailMessage("matrix.inverse() spoils the source matrix data").
        isEqualTo(expectedMatrixData);
  }

  @Test
  @DisplayName("inverse() does not spoil the previously found vector solution")
  void testInverseDoesNotSpoilTheVectorSolution() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixVector);
    final Number[] expectedSolution = matrix.getSolution();

    matrix.inverse();
    final Number[] actualSolution = matrix.getSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.inverse() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("inverse() does not spoil the previously found matrix solution")
  void testInverseDoesNotSpoilTheMatrixSolution() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixB);
    final Number[][] expectedSolution = matrix.getNumberMatrixSolution();

    matrix.inverse();
    final Number[][] actualSolution = matrix.getNumberMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.inverse() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("inverse() returns the inversion of the matrix")
  void testInverseReturnsTheInversion() {
    printMethodName();
    final BigDecimal[][] expectedInversion = sampleInvertibleMatrixInversionAsBigDecimals;
    final double expectedError = ERROR_TOLERANCE;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleInvertibleMatrixData);
    final BigDecimal[][] actualInversion = matrix.inverse().getBigDecimalData();
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
    final BigDecimal[][] expectedMatrixData = sampleSolvableMatrixDataAsBigDecimals;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(expectedMatrixData);
    matrix.inverseAccurately();
    final BigDecimal[][] actualMatrixData = matrix.getBigDecimalData();

    assertThat(actualMatrixData).
        withFailMessage("matrix.inverse() spoils the source matrix data").
        isEqualTo(expectedMatrixData);
  }

  @Test
  @DisplayName("inverseAccurately() does not spoil the previously found vector solution")
  void testInverseAccuratelyDoesNotSpoilTheVectorSolution() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixVector);
    final Number[] expectedSolution = matrix.getSolution();

    matrix.inverseAccurately();
    final Number[] actualSolution = matrix.getSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.inverse() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("inverseAccurately() does not spoil the previously found matrix solution")
  void testInverseAccuratelyDoesNotSpoilTheMatrixSolution() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixB);
    final BigDecimal[][] expectedSolution = matrix.getBigDecimalMatrixSolution();

    matrix.inverseAccurately();
    final BigDecimal[][] actualSolution = matrix.getBigDecimalMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.inverse() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("inverseAccurately() returns the inversion of the matrix")
  void testInverseAccuratelyReturnsTheInversion() {
    printMethodName();
    final BigDecimal[][] expectedInversion = sampleInvertibleMatrixInversionAsBigDecimals;
    final double expectedError = ERROR_TOLERANCE;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleInvertibleMatrixData);
    final BigDecimal[][] actualInversion = matrix.inverseAccurately().getBigDecimalData();
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
   * accuracy improvement achieved by the {@link BigDecimalMatrix#inverseAccurately()} compared with the
   * {@link BigDecimalMatrix#inverse()}. Asserts that the minimum improvement ratio is greater that {@link #EXPECTED_MIN_IMPROVEMENT}
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
    final BigDecimal[][] expectedMatrixData = sampleSolvableMatrixDataAsBigDecimals;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(expectedMatrixData);
    matrix.transpose();
    final BigDecimal[][] actualMatrixData = matrix.getBigDecimalData();

    assertThat(actualMatrixData).
        withFailMessage("matrix.transpose() spoils the source matrix data").
        isEqualTo(expectedMatrixData);
  }

  @Test
  @DisplayName("transpose() does not spoil the previously found vector solution")
  void testTransposeDoesNotSpoilTheVectorSolution() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixVector);
    final BigDecimal[] expectedSolution = matrix.getBigDecimalSolution();

    matrix.transpose();
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.transpose() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("transpose() does not spoil the previously found matrix solution")
  void testTransposeDoesNotSpoilTheMatrixSolution() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    matrix.solve(sampleSolvableMatrixMatrixB);
    final BigDecimal[][] expectedSolution = matrix.getBigDecimalMatrixSolution();

    matrix.transpose();
    final BigDecimal[][] actualSolution = matrix.getBigDecimalMatrixSolution();

    assertThat(actualSolution).
        withFailMessage("matrix.transpose() spoils the previously found solution").
        isEqualTo(expectedSolution);
  }

  @Test
  @DisplayName("transpose() returns the transposition of the matrix")
  void testTransposeReturnsTheTransposition() {
    // It does not pertain to precision, so double[][] data are OK
    final double[][] expectedTransposition = sampleInvertibleMatrixTransposition;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleInvertibleMatrixData);
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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(anyMatrixDataNoMatterWhatItIs);
    final BigDecimal[][] expectedUnity = bigDecimalUnityMatrix(matrix.getSize());

    final BigDecimal[][] actualUnity = matrix.unity().getBigDecimalData();

    assertThat(actualUnity).
        withFailMessage("The value returned by BigDecimalMatrix.unityMatrix() is not a unity matrix").
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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.multiply(new BigDecimalMatrix(sampleSolvableMatrixTooLargeMatrix)),
        "Attempt to call multiply(Matrix) with different size of argument must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'multiply(Matrix)' and 'must have the same size'").
        contains("multiply(matrix)").contains("must have the same size");
  }

  @Test
  @DisplayName("multiply(Matrix) returns a correct product")
  void testMultiplyByMatrixReturnsCorrectProduct() {
    // Use Quadruples to compare products, because BigDecimal.equals(BigDecimal) may return false even when numeric
    // values are equal (they take trailing zeros into account).
    final Quadruple[][] expectedProduct = sampleSolvableMatrixMatrixBAsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    final Quadruple[][] actualProduct = matrix.multiply(new BigDecimalMatrix(sampleSolvableMatrixMatrixX)).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[][] expectedProduct = sampleSolvableMatrixMatrixBAsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    final Quadruple[][] actualProduct = matrix.multiply(sampleSolvableMatrixMatrixX).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[][] expectedProduct = sampleSolvableMatrixMatrixBAsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    final Quadruple[][] actualProduct = matrix.multiply(sampleSolvableMatrixMatrixXAsQuadruples).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[][] expectedProduct = sampleSolvableMatrixMatrixBAsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    final Quadruple[][] actualProduct = matrix.multiply(sampleSolvableMatrixMatrixXAsBigDecimals).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[] expectedProduct = sampleSolvableMatrixVectorAsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    final Quadruple[] actualProduct = convertToQuadruples(matrix.multiply(sampleSolvableMatrixSolution));

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[] expectedProduct = sampleSolvableMatrixVectorAsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    final Quadruple[] actualProduct = convertToQuadruples(matrix.multiply(sampleSolvableMatrixSolutionAsQuadruples));

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[] expectedProduct = sampleSolvableMatrixVectorAsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);
    final Quadruple[] actualProduct = convertToQuadruples(matrix.multiply(sampleSolvableMatrixSolutionAsBigDecimals));

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final double[][] expectedProduct = deepCopyOf(sampleRandomMultipliedByFactor);

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleRandomMatrixData);
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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimal[][] expectedProduct = sampleRandomMultipliedByQuadrupleFactorAsBigDecimals;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleRandomMatrixData);
    final BigDecimal[][] actualProduct = matrix.multiply(QUADRUPLE_FACTOR).getBigDecimalData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[][] expectedProduct = sampleRandomMultipliedByFactorAsQuadruples;
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleRandomMatrixData);

    final Quadruple[][] actualProduct = matrix.multiply(QUADRUPLE_FACTOR.bigDecimalValue()).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.add(new BigDecimalMatrix(sampleSolvableMatrixTooLargeMatrix)),
        "Attempt to call add(Matrix) with argument of wrong size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'add(Matrix)' and 'must have the same size'").
        contains("add(matrix)").contains("must have the same size");
  }

  @Test
  @DisplayName("add(Matrix) returns a correct result")
  void testAddMatrixReturnsCorrectResult() {
    final Quadruple[][] expectedSum = summand1PlusSummand2AsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1AsQuadruples);
    final Quadruple[][] actualSum = matrix.add(new BigDecimalMatrix(summand2AsQuadruples)).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimal[][] expectedSum = summand1PlusSummand2AsBigDecimals;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);
    final BigDecimal[][] actualSum = matrix.add(summand2).getBigDecimalData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[][] expectedSum = summand1PlusSummand2AsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1AsQuadruples);
    final Quadruple[][] actualSum = matrix.add(summand2AsQuadruples).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[][] expectedSum = summand1PlusSummand2AsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1AsQuadruples);
    final Quadruple[][] actualSum = matrix.add(summand2AsBigDecimals).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(minuend);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(minuend);

    final Throwable expectedException = assertThrows(IllegalArgumentException.class,
        () -> matrix.subtract(new BigDecimalMatrix(sampleSolvableMatrixTooLargeMatrix)),
        "Attempt to call subtract(Matrix) with argument of wrong size must throw an exception ");

    assertThat(expectedException.getMessage().toLowerCase()).
        withFailMessage("Exception message should contain 'subtract(Matrix)' and 'must have the same size'").
        contains("subtract(matrix)").contains("must have the same size");
  }

  @Test
  @DisplayName("subtract(Matrix) returns a correct result")
  void testSubtractMatrixReturnsCorrectResult() {
    final Quadruple[][] expectedDifference = summand2AsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1PlusSummand2AsQuadruples);
    final Quadruple[][] actualDifference = matrix.subtract(new BigDecimalMatrix(summand1AsQuadruples)).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(minuend);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(minuend);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(minuend);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[][] expectedDifference = summand1AsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1PlusSummand2AsQuadruples);
    final Quadruple[][] actualDifference = matrix.subtract(summand2AsQuadruples).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(minuend);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(minuend);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(minuend);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[][] expectedDifference = summand1AsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1PlusSummand2AsQuadruples);
    final Quadruple[][] actualDifference = matrix.subtract(summand2AsQuadruples).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(minuend);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(minuend);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(minuend);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixData);

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
    final Quadruple[][] expectedDifference = summand2AsQuadruples;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(summand1PlusSummand2AsQuadruples);
    final Quadruple[][] actualDifference = matrix.subtract(summand1AsQuadruples).getQuadrupleData();

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSpdSolvableMatrixDataAsQuadruples);
    final BigDecimal expectedDeterminant = BigDecimal.valueOf(sampleSpdMatrixDataDeterminant); // 10.86
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE .. ERROR_TOLERANCE]
    final double expectedError = Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant.doubleValue() * ERROR_TOLERANCE));

    final BigDecimal actualDeterminant = (BigDecimal)matrix.determinant();
    final double actualError = Math.abs((actualDeterminant.doubleValue() - expectedDeterminant.doubleValue())
                                       /expectedDeterminant.doubleValue());

    assertThat(actualError).
        withFailMessage("determinant() for a matrix with a positive determinant returns a wrong value").
        isLessThan(expectedError);
  }

  @Test
  @DisplayName("determinant() for a matrix with a negative determinant returns a correct value")
  void testDeterminantWhenNegativeReturnsCorrectValue() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixMatrixXAsQuadruples);
    final BigDecimal expectedDeterminant = BigDecimal.valueOf(sampleSolvableMatrixMatrixXDeterminant);  // -185.0
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE..ACCEPTED_ERROR]
    final double expectedError = Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant.doubleValue() * ERROR_TOLERANCE));

    final BigDecimal actualDeterminant = (BigDecimal)matrix.determinant();
    final double actualError = Math.abs((actualDeterminant.doubleValue() - expectedDeterminant.doubleValue())
                                        /expectedDeterminant.doubleValue());

    assertThat(actualError).
        withFailMessage("determinant() for a matrix with a positive determinant returns a wrong value").
        isLessThan(expectedError);
  }

  @Test
  @DisplayName("determinant() for an inconsistent matrix returns zero")
  void testDeterminantForInconsistentMatrixReturnsZero() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleInconsistentMatrixData);

    final BigDecimal actualDeterminant = (BigDecimal)matrix.determinant();

    assertThat(actualDeterminant.signum() == 0).
        withFailMessage("determinant() for an inconsistent matrix returns a non-zero result").
        isTrue();
  }

  @Test
  @DisplayName("determinant() for an underdetermined matrix returns zero")
  void testDeterminantForUnderdeterminedMatrixReturnsZero() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleUnderdeterminedMatrixData);

    final BigDecimal actualDeterminant = (BigDecimal)matrix.determinant();

    assertThat(actualDeterminant.signum() == 0).
        withFailMessage("determinant() for an underdetermined matrix returns a non-zero result").
        isTrue();
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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSpdSolvableMatrixDataAsQuadruples);
    final double expectedDeterminant = sampleSpdMatrixDataDeterminant; // 10.86
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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixMatrixX);
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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleInconsistentMatrixData);

    final double actualDeterminant = matrix.determinantAsDouble();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsDouble() for an inconsistent matrix returns a non-zero result").
        isZero();
  }

  @Test
  @DisplayName("determinantAsDouble() for an underdetermined matrix returns zero")
  void testDeterminantAsDoubleForUnderdeterminedMatrixReturnsZero() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleUnderdeterminedMatrixData);

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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSpdSolvableMatrixDataAsQuadruples);
    final Quadruple expectedDeterminant = new Quadruple(BigDecimal.valueOf(sampleSpdMatrixDataDeterminant)); // 10.860000000000000000000000000000
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE .. ERROR_TOLERANCE]
    final double expectedError = Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant.doubleValue() * ERROR_TOLERANCE));

    final Quadruple actualDeterminant = matrix.determinantAsQuadruple();
    final double actualError = Quadruple.subtract(actualDeterminant, expectedDeterminant).
        divide(expectedDeterminant).abs().doubleValue();

    assertThat(actualError).
        withFailMessage("determinantAsQuadruple() for a matrix with a positive determinant returns a wrong value").
        isLessThan(expectedError);
  }

  @Test
  @DisplayName("determinantAsQuadruple() for a matrix with a negative determinant returns a correct value")
  void testdeterminantAsQuadrupleWhenNegativeReturnsCorrectValue() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixMatrixXAsQuadruples);
    final Quadruple expectedDeterminant = new Quadruple(sampleSolvableMatrixMatrixXDeterminant);  // -185.0
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE..ACCEPTED_ERROR]
    final double expectedError = Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant.doubleValue() * ERROR_TOLERANCE));

    final Quadruple actualDeterminant = matrix.determinantAsQuadruple();
    final double actualError = Quadruple.subtract(actualDeterminant, expectedDeterminant).
        divide(expectedDeterminant).abs().doubleValue();

    assertThat(actualError).
        withFailMessage("determinantAsQuadruple() for a matrix with a negative determinant returns a wrong value").
        isLessThan(expectedError);
  }

  @Test
  @DisplayName("determinantAsQuadruple() for an inconsistent matrix returns zero")
  void testdeterminantAsQuadrupleForInconsistentMatrixReturnsZero() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleInconsistentMatrixData);

    final Quadruple actualDeterminant = matrix.determinantAsQuadruple();

    assertThat(actualDeterminant.isZero()).
        withFailMessage("determinantAsQuadruple() for an inconsistent matrix returns a non-zero result").
        isTrue();
  }

  @Test
  @DisplayName("determinantAsQuadruple() for an underdetermined matrix returns zero")
  void testdeterminantAsQuadrupleForUnderdeterminedMatrixReturnsZero() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleUnderdeterminedMatrixData);

    final Quadruple actualDeterminant = matrix.determinantAsQuadruple();

    assertThat(actualDeterminant.isZero()).
        withFailMessage("determinantAsQuadruple() for an underdetermined matrix returns a non-zero result").
        isTrue();
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
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSpdSolvableMatrixDataAsQuadruples);
    final BigDecimal expectedDeterminant = BigDecimal.valueOf(sampleSpdMatrixDataDeterminant); // 10.860000000000000000000000000000
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE .. ERROR_TOLERANCE]
    final double expectedError = Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant.doubleValue() * ERROR_TOLERANCE));

    final BigDecimal actualDeterminant = matrix.determinantAsBigDecimal();
    final double actualError = actualDeterminant.subtract(expectedDeterminant).
        divide(expectedDeterminant, MC_60).abs().doubleValue();

    assertThat(actualError).
        withFailMessage("determinantAsBigDecimal() for a matrix with a positive determinant returns a wrong value").
        isLessThan(expectedError);
  }

  @Test
  @DisplayName("determinantAsBigDecimal() for a matrix with a negative determinant returns a correct value")
  void testdeterminantAsBigDecimalWhenNegativeReturnsCorrectValue() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixMatrixX);
    final BigDecimal expectedDeterminant = BigDecimal.valueOf(sampleSolvableMatrixMatrixXDeterminant);  // -185.0
    // Both absolute and relative errors are expected to be within the range of [-ERROR_TOLERANCE..ACCEPTED_ERROR]
    final double expectedError = Math.max(ERROR_TOLERANCE, Math.abs(expectedDeterminant.doubleValue() * ERROR_TOLERANCE));

    final BigDecimal actualDeterminant = matrix.determinantAsBigDecimal();
    final double actualError = actualDeterminant.subtract(expectedDeterminant).
        divide(expectedDeterminant, MC_60).abs().doubleValue();

    assertThat(actualError).
        withFailMessage("determinantAsBigDecimal() for a matrix with a negative determinant returns a wrong value").
        isLessThan(expectedError);
  }

  @Test
  @DisplayName("determinantAsBigDecimal() for an inconsistent matrix returns zero")
  void testdeterminantAsBigDecimalForInconsistentMatrixReturnsZero() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleInconsistentMatrixData);

    final BigDecimal actualDeterminant = matrix.determinantAsBigDecimal();

    assertThat(actualDeterminant).
        withFailMessage("determinantAsBigDecimal() for an inconsistent matrix returns a non-zero result").
        isZero();
  }

  @Test
  @DisplayName("determinantAsBigDecimal() for an underdetermined matrix returns zero")
  void testdeterminantAsBigDecimalForUnderdeterminedMatrixReturnsZero() {
    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleUnderdeterminedMatrixData);

    final BigDecimal actualDeterminant = matrix.determinantAsBigDecimal();

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
    final BigDecimal expectedNorm = BigDecimal.valueOf(sampleSolvableMatrixAnotherMatrixXNorm); // 18.143;
    final double expectedError = Math.max(ERROR_TOLERANCE, Math.abs(expectedNorm.doubleValue() * ERROR_TOLERANCE));

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixAnotherMatrixXAsQuadruples());
    final BigDecimal actualNorm = (BigDecimal)matrix.norm();

    final double actualError = actualNorm.subtract(expectedNorm).
        divide(expectedNorm, MC_60).abs().doubleValue();

    assertThat(actualError).
        withFailMessage("norm() for a matrix returns a wrong value").
        isLessThan(expectedError);
  }

  @Test
  @DisplayName("scaling does not affect the returned norm() value")
  void testNormDoesNotDependOnScale() {
    final BigDecimal expectedNorm = BigDecimal.valueOf(sampleSolvableMatrixAnotherMatrixXNorm); // 18.1430;
    final double expectedError = Math.max(ERROR_TOLERANCE, Math.abs(expectedNorm.doubleValue() * ERROR_TOLERANCE));

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleSolvableMatrixAnotherMatrixXAsQuadruples(), true);
    matrix.solve(sampleSolvableMatrixMatrixB);                                               // force scaling
    final BigDecimal actualNorm = (BigDecimal)matrix.norm();

    final double actualError = actualNorm.subtract(expectedNorm).
        divide(expectedNorm, MC_60).abs().doubleValue();

    assertThat(actualError).
        withFailMessage("after scaling, norm() returns a wrong value").
        isLessThan(expectedError);
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

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleInvertibleMatrixData);
    final double actualCond = matrix.cond();

    assertThat(actualCond).
        withFailMessage("expected cond() value is %s, and actual value is %s", expectedCond, actualCond).
        isCloseTo(expectedCond, delta);
  }

  @Test
  @DisplayName("cond() for a non-invertible matrix returns Double.POSITIVE_INFINITY")
  void testCondNonInvertibleReturnsCorrectValue() {
    final double expectedCond = Double.POSITIVE_INFINITY;

    final BigDecimalMatrix matrix = new BigDecimalMatrix(sampleInconsistentMatrixData);
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
        MatrixData::bigDecimalMatrixInversionErrors,
        MatrixData::bigDecimalAccurateMatrixInversionErrors
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
