package com.mvohm.quadmatrix.test;

import static com.mvohm.quadmatrix.test.AuxMethods.*;
import static com.mvohm.quadmatrix.test.MatrixDataGenerators.*;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Random;

import com.mvohm.quadmatrix.BigDecimalMatrix;
import com.mvohm.quadmatrix.DoubleMatrix;
import com.mvohm.quadmatrix.Matrix;
import com.mvohm.quadmatrix.QuadrupleMatrix;
import com.mvohm.quadmatrix.test.AuxMethods.ErrorSet;
import com.mvohm.quadruple.Quadruple;

/**
 * A helper class to generate random matrix data and examine accuracy of the solutions of the corresponding matrix
 */
public class MatrixData {

  private static final double RANGE_BOTTOM          = -1.0;
  private static final double RANGE_TOP             = 1.0;

  private int size;
  private double[][] matrixData;
  private double[] vector;
  private double[] solution;

  private Quadruple[] quadrupleVector;
  private Quadruple[] quadrupleSolution;

  private BigDecimal[] bigDecimalVector;
  private BigDecimal[] bigDecimalSolution;

  private double[][] matrixB;
  private double[][] matrixX;

  private Quadruple[][] quadrupleMatrixB;
  private Quadruple[][] quadrupleMatrixX;

  private BigDecimal[][] bigDecimalMatrixB;
  private BigDecimal[][] bigDecimalMatrixX;

  private static double[][] unityMatrix = null;
  private static Quadruple[][] quadrupleUnityMatrix = null;
  private static BigDecimal[][] bigDecimalUnityMatrix = null;

  private long time;

  private enum Purpose { VECTOR_SOLUTION, MATRIX_SOLUTION, INVERSION };

  private Purpose purpose;

  private String method_1_name = null;
  private String method_2_name = null;

  private static final MathContext MC = new MathContext(BigDecimalMatrix.getDefaultPrecision(), RoundingMode.HALF_EVEN);

  /* ***************************************************************************
   *********** Creating datasets ***********************************************
   *****************************************************************************/

  /**
   * fills matrixData with random values providing that it is SPD,
   * fills the solution with random values in the range defined by corresponding constants,
   * and finds the vector such that solving the equation system with this vector would give the found solution
   * @param size the size of the matrix to create
   * @param random contains a random seed and ensures reproducibility
   * @return and instance of RandomMatrixData with the generated data
   */
  public static MatrixData makeDataSetForSPDSolutions(int size, Random random) {
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.VECTOR_SOLUTION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomSpdMatrix(size, RANGE_BOTTOM, RANGE_TOP, 1.0));
    data.solution = randomVector(size, RANGE_BOTTOM, RANGE_TOP);
    data.vector = multiply(data.matrixData, data.solution);
    return data;
  }

  /**
   * fills matrixData and solution with random values in the range defined by corresponding constants
   * and finds the vector such that solving the system with this vector would give the found solution
   * @param size the size of the matrix to create
   * @param random contains a random seed and ensures reproducibility
   * @return and instance of RandomMatrixData with the generated data
   */
  public static MatrixData makeDataSetForVectorSolutions(int size, Random random) {
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.VECTOR_SOLUTION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomMatrix(size, RANGE_BOTTOM, RANGE_TOP));
    data.solution = randomVector(size, RANGE_BOTTOM, RANGE_TOP);
    data.vector = multiply(data.matrixData, data.solution);
    return data;
  }

  /**
   * Generates matrix data whose rows differ significantly in magnitude from each other.
   * fills matrixData A and a target solution x with random values in the range -1.0 .. 1.0,
   * then scales the rows of the matrix by a scale gradually growing from 1/sqrt(range) to sqrt(range),
   * and finds the vector b such that solving the system with the scaled matrix
   * and this vector would give the required solution.
   * @param size the size of the matrix to create
   * @param random contains a random seed and ensures reproducibility
   * @return and instance of RandomMatrixData with the generated data
   * @param size the size of the matrix to create
   * @param random contains a random seed and ensures reproducibility
   * @param scaleRange
   * @return
   */
  public static MatrixData makeLargeRangeDataSetForVectorSolutions(int size, Random random, double scaleRange) {
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.VECTOR_SOLUTION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomMatrix(size, RANGE_BOTTOM, RANGE_TOP));

    // Multiply rows by an increasing scale
    final double scaleIncrease = Math.pow(scaleRange, 1.0/(size));  //
    double scale = 1.0 / Math.sqrt(scaleRange);
    for (int i = 0; i < size; i++) {
      data.matrixData[i] = multiply(data.matrixData[i], scale);
      scale *= scaleIncrease;
    }

    data.solution = randomVector(size, RANGE_BOTTOM, RANGE_TOP);
    data.vector = multiply(data.matrixData, data.solution);
    return data;
  }

  /**
   * fills matrixData and solution with non-uniform random values in the range defined by corresponding constants
   * and finds the vector such that solving the system with this vector would give the found solution
   * @param size the size of the matrix to create
   * @param random contains a random seed and ensures reproducibility
   * @return and instance of RandomMatrixData with the generated data
   */
  public static MatrixData makeNonUniformDataSet(int size, Random random, double density, double power, double slope) {
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.VECTOR_SOLUTION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomSparsePowPlusLinearMatrix(size, density, power, slope));
    data.solution = randomVector(size, RANGE_BOTTOM, RANGE_TOP);
    data.vector = multiply(data.matrixData, data.solution);
    return data;
  }

  public static MatrixData makeDataSetForMatrixSolutions(int size, Random random) {
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.MATRIX_SOLUTION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomMatrix(size, RANGE_BOTTOM, RANGE_TOP));
    data.matrixX = randomMatrix(size, RANGE_BOTTOM, RANGE_TOP);
    data.matrixB = multiply(data.matrixData, data.matrixX);
    return data;
  }

  public static MatrixData makeDataSetForInversions(int size, Random random) {
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.INVERSION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomMatrix(size, RANGE_BOTTOM, RANGE_TOP));
    if (unityMatrix == null || unityMatrix.length != size) {
      unityMatrix = unityMatrix(size);
    }
    return data;
  }

  /* ************************************************************************
   ****** Perform operations and return errors ******************************
   **************************************************************************/

  public ErrorSet spdSolutionErrors() {
    setTestedMethodName("DoubleMatrix.solveSPD(double[])");
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    final DoubleMatrix matrix = new DoubleMatrix(matrixData);
    time = -System.nanoTime();
    matrix.solveSPD(vector);
    time += System.nanoTime();
    final double[] actualSolution = matrix.getDoubleSolution();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleSPDSolutionErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("QuadrupleMatrix.solveSPD(Quadruple[])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    quadrupleSolution = makeQuadrupleSolution();
    quadrupleVector = multiply(matrix.getQuadrupleData(), quadrupleSolution);
    time = -System.nanoTime();
    matrix.solveSPD(quadrupleVector);
    time += System.nanoTime();
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();
    return findErrors(quadrupleSolution, actualSolution, false).setTime(time);
  }

  public ErrorSet bigDecimalSPDSolutionErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("BigDecimalMatrix.solveSPD(BigDecimal[])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    bigDecimalSolution = makeBigDecimalSolution();
    bigDecimalVector = multiply(matrix.getBigDecimalData(), bigDecimalSolution);
    time = -System.nanoTime();
    matrix.solveSPD(bigDecimalVector);
    time += System.nanoTime();
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();
    return findErrors(bigDecimalSolution, actualSolution).setTime(time);
  }

  public ErrorSet accurateSPDSolutionErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("DoubleMatrix.solveSPDAccurately(double[])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData);
    time = -System.nanoTime();
    matrix.solveSPDAccurately(vector);
    time += System.nanoTime();
    final double[] actualSolution = matrix.getDoubleSolution();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleAccurateSPDSolutionErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    if (quadrupleSolution == null) {
      throw new RuntimeException("quadrupleSPDSolutionErrors() must be called before quadrupleAccurateSPDSolutionErrors()");
    }
    setTestedMethodName("QuadrupleMatrix.solveSPDAccurately(Quadruple[])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    time = -System.nanoTime();
    matrix.solveSPDAccurately(quadrupleVector);
    time += System.nanoTime();
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();
    return findErrors(quadrupleSolution, actualSolution).setTime(time);
  }

  public ErrorSet bigDecimalAccurateSPDSolutionErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    if (bigDecimalSolution == null) {
      throw new RuntimeException("bifDecimaSPDSolutionErrors() must be called before bigDecimalAccurateSPDSolutionErrors()");
    }
    setTestedMethodName("BigDecimalMatrix.solveSPDAccurately(BigDecimal[])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    time = -System.nanoTime();
    matrix.solveSPDAccurately(bigDecimalVector);
    time += System.nanoTime();
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();
    return findErrors(bigDecimalSolution, actualSolution).setTime(time);
  }

  public ErrorSet luSolutionWithScalingErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(scaled) DoubleMatrix.solve(double[])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, true);
    time = -System.nanoTime();
    matrix.solve(vector);
    time += System.nanoTime();
    final double[] actualSolution = matrix.getDoubleSolution();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleLuSolutionWithScalingErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(scaled) QuadrupleMatrix.solve(Quadruple[])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    quadrupleSolution = makeQuadrupleSolution();
    quadrupleVector = multiply(matrix.getQuadrupleData(), quadrupleSolution);
    time = -System.nanoTime();
    matrix.solve(quadrupleVector);
    time += System.nanoTime();
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();
    return findErrors(quadrupleSolution, actualSolution, false).setTime(time);
  }

  public ErrorSet bigDecimalLuSolutionWithScalingErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(scaled) BigDecimalMatrix.solve(BigDecimal[])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    bigDecimalSolution = makeBigDecimalSolution();
    bigDecimalVector = multiply(matrix.getBigDecimalData(), bigDecimalSolution );
    time = -System.nanoTime();
    matrix.solve(bigDecimalVector);
    time += System.nanoTime();
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();
    return findErrors(bigDecimalSolution, actualSolution, false).setTime(time);
  }

  public ErrorSet luSolutionWithoutScalingErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(unscaled) DoubleMatrix.solve(double[])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, false);
    time = -System.nanoTime();
    matrix.solve(vector);
    time += System.nanoTime();
    final double[] actualSolution = matrix.getDoubleSolution();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleLuSolutionWithoutScalingErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(unscaled) QuadrupleMatrix.solve(Quadruple[])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, false);
    quadrupleSolution = makeQuadrupleSolution();
    quadrupleVector = multiply(matrix.getQuadrupleData(), quadrupleSolution);
    time = -System.nanoTime();
    matrix.solve(quadrupleVector);
    time += System.nanoTime();
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();
    return findErrors(quadrupleSolution, actualSolution).setTime(time);
  }

  public ErrorSet bigDecimalLuSolutionWithoutScalingErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(unscaled) BigDecimalMatrix.solve(BigDecimal[])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, false);
    bigDecimalSolution = makeBigDecimalSolution();
    bigDecimalVector = multiply(matrix.getBigDecimalData(), bigDecimalSolution );
    time = -System.nanoTime();
    matrix.solve(bigDecimalVector);
    time += System.nanoTime();
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();
    return findErrors(bigDecimalSolution, actualSolution, false).setTime(time);
  }

  public ErrorSet accurateLUSolutionWithScalingErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("DoubleMatrix.solveAccurately(double[])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData);
    time = -System.nanoTime();
    matrix.solveAccurately(vector);
    time += System.nanoTime();
    final double[] actualSolution = matrix.getDoubleSolution();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleAccurateLUSolutionWithScalingErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("QuadrupleMatrix.solveAccurately(Quadruple[])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    if (quadrupleSolution == null) {
      throw new RuntimeException("quadrupleLuSolutionWithScalingErrors() must be called before quadrupleAccurateLUSolutionWithScalingErrors()");
    }
    time = -System.nanoTime();
    matrix.solveAccurately(quadrupleVector);
    time += System.nanoTime();
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();
    return findErrors(quadrupleSolution, actualSolution).setTime(time);
  }

  public ErrorSet bigDecimalAccurateLUSolutionWithScalingErrors() {
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("BigDecimalMatrix.solveAccurately(BigDecimal[])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    if (bigDecimalSolution == null) {
      throw new RuntimeException("bigDecimalLuSolutionWithScalingErrors() must be called before bigDecimalAccurateLUSolutionWithScalingErrors()");
    }
    time = -System.nanoTime();
    matrix.solveAccurately(bigDecimalVector);
    time += System.nanoTime();
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();
    return findErrors(bigDecimalSolution, actualSolution, false).setTime(time);
  }

  public ErrorSet matrixSolutionErrors() {
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("DoubleMatrix.solve(double[][])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, true);
    time = -System.nanoTime();
    final double[][] actualSolution = matrix.solve(matrixB).getDoubleData();
    time += System.nanoTime();
    return findErrors(this.matrixX, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleMatrixSolutionErrors() {
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("QuadrupleMatrix.solve(Quadruple[][])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);

    quadrupleMatrixX = makeQuadrupleMatrixX();
    quadrupleMatrixB = multiply(matrix.getQuadrupleData(), quadrupleMatrixX);

    time = -System.nanoTime();
    final Quadruple[][] actualSolution = matrix.solve(quadrupleMatrixB).getQuadrupleData();
    time += System.nanoTime();
    return findErrors(quadrupleMatrixX, actualSolution).setTime(time);
  }

  public ErrorSet bigDecimalMatrixSolutionErrors() {
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("BigDecimalMatrix.solve(BigDecimal[][])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);

    bigDecimalMatrixX = makeBigDecimalMatrixX();
    bigDecimalMatrixB = multiply(matrix.getBigDecimalData(), bigDecimalMatrixX);

    time = -System.nanoTime();
    final BigDecimal[][] actualSolution = matrix.solve(bigDecimalMatrixB).getBigDecimalData();
    time += System.nanoTime();
    return findErrors(bigDecimalMatrixX, actualSolution).setTime(time);
  }

  public ErrorSet accurateMatrixSolutionErrors() {
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("DoubleMatrix.solveAccurately(double[][])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, true);
    time = -System.nanoTime();
    final double[][] actualSolution = matrix.solveAccurately(matrixB).getDoubleData();
    time += System.nanoTime();
    return findErrors(this.matrixX, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleAccurateMatrixSolutionErrors() {
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    if (quadrupleMatrixX == null) {
      throw new RuntimeException("quadrupleMmatrixSolutionErrors() must be called before quadrupleAccurateMatrixSolutionErrors()");
    }
    setTestedMethodName("QuadrupleMatrix.solveAccurately(Quadruple[][])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    time = -System.nanoTime();
    final Quadruple[][] actualSolution = matrix.solveAccurately(quadrupleMatrixB).getQuadrupleData();
    time += System.nanoTime();
    return findErrors(quadrupleMatrixX, actualSolution).setTime(time);
  }

  public ErrorSet bigDecimalAccurateMatrixSolutionErrors() {
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    if (bigDecimalMatrixX == null) {
      throw new RuntimeException("bigDecimalMmatrixSolutionErrors() must be called before bigDecimalAccurateMatrixSolutionErrors()");
    }
    setTestedMethodName("BigDecimalMatrix.solveAccurately(BigDecimal[][])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    time = -System.nanoTime();
    final BigDecimal[][] actualSolution = matrix.solveAccurately(bigDecimalMatrixB).getBigDecimalData();
    time += System.nanoTime();
    return findErrors(bigDecimalMatrixX, actualSolution).setTime(time);
  }

  public ErrorSet matrixInversionErrors() {
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    setTestedMethodName("DoubleMatrix.inverse()");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverse();
    time += System.nanoTime();
    final double[][] product = multiply(matrixData, inverse.getDoubleData());
    return findErrors(unityMatrix, product, 1).setTime(time);
  }

  public ErrorSet quadrupleMatrixInversionErrors() {
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    setTestedMethodName("QuadrupleMatrix.inverse()");
    quadrupleUnityMatrix = convertToQuadruples(unityMatrix);
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverse();
    time += System.nanoTime();
    final Quadruple[][] product = multiply(matrix.getQuadrupleData(), inverse.getQuadrupleData());
    return findErrors(quadrupleUnityMatrix, product).setTime(time);
  }

  public ErrorSet bigDecimalMatrixInversionErrors() {
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    setTestedMethodName("BigDecimalMatrix.inverse()");
    bigDecimalUnityMatrix = convertToBigDecimals(unityMatrix);
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverse();
    time += System.nanoTime();
    final BigDecimal[][] product = multiply(matrix.getBigDecimalData(), inverse.getBigDecimalData());
    return findErrors(bigDecimalUnityMatrix, product).setTime(time);
  }

  public ErrorSet accurateMatrixInversionErrors() {
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    setTestedMethodName("DoubleMatrix.inverseAccurately()");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverseAccurately();
    time += System.nanoTime();
    final double[][] product = multiply(matrixData, inverse.getDoubleData());
    return findErrors(unityMatrix, product, 1).setTime(time);
  }

  public ErrorSet quadrupleAccurateMatrixInversionErrors() {
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    if (quadrupleUnityMatrix == null) {
      throw new RuntimeException("quadrupleMatrixInversionErrors() must be called before quadrupleAccurateMatrixInversionErrors()");
    }
    setTestedMethodName("QuadrupleMatrix.inverseAccurately()");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverseAccurately();
    time += System.nanoTime();
    final Quadruple[][] product = multiply(matrix.getQuadrupleData(), inverse.getQuadrupleData());
    return findErrors(quadrupleUnityMatrix, product).setTime(time);
  }

  public ErrorSet bigDecimalAccurateMatrixInversionErrors() {
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    setTestedMethodName("BigDecimalMatrix.inverseAccurately()");
    if (bigDecimalUnityMatrix == null)
      bigDecimalUnityMatrix = convertToBigDecimals(unityMatrix);
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverseAccurately();
    time += System.nanoTime();
    final BigDecimal[][] product = multiply(matrix.getBigDecimalData(), inverse.getBigDecimalData());
    return findErrors(bigDecimalUnityMatrix, product).setTime(time);
  }

  /* *************************************************************************
  ******** Private methods ***************************************************
  ***************************************************************************/

  private void setMatrixData(double[][] data) {
    matrixData = data;
    size = data.length;
  }

  int getSize() {
    return size;
 }

  private void setTestedMethodName(String methodName) {
    if (method_1_name == null) {
      method_1_name = methodName;
    } else if (method_2_name == null) {
      method_2_name = methodName;
    }
  }

  String testedMethod1name() {
    return method_1_name;
  }

  String testedMethod2name() {
    return method_2_name;
  }

  private static String performerMethodName() {
    final String s = new Exception().getStackTrace()[1].toString();
    return s.replaceFirst("\\(.*\\)", "():").replaceFirst(".*\\.", "");
  }

  private void checkPurpose(Object methodName, Purpose impliedPurpose) {
    if (this.purpose != impliedPurpose) {
      throw new IllegalArgumentException(String.format(
                  "Calling %s on MatrixData implies %s, while this instance was creted for %s",
                  methodName, impliedPurpose, this.purpose));
    }
  }

  public long  getTime() {
    return time;
  }

  private static double[] multiply(double[] vector, double factor) {
    final double[] result = new double[vector.length];
    for (int i = 0; i < vector.length; i++) {
      result[i] = vector[i] * factor;
    }
    return result;
  }

  private Quadruple[] makeQuadrupleSolution() {
    return convertToQuadruples(solution);
  }

  private BigDecimal[] makeBigDecimalSolution() {
    return convertToBigDecimals(solution);
  }

  private Quadruple[][] makeQuadrupleMatrixX() {
    return convertToQuadruples(matrixX);
  }

  private BigDecimal[][] makeBigDecimalMatrixX() {
    return convertToBigDecimals(matrixX);
  }

  private static double[][] multiply(double[][] matrixA, double[][] matrixB) {
    final int size = matrixA.length;
    final double[][] result = new double[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        result[i][j] = multiplyRowByColumn(matrixA, i, matrixB, j);
      }
    }
    return result;
  }

  private static double multiplyRowByColumn(double[][] matrixA, int rowIndex, double[][] matrixB, int colIndex) {
    final int size = matrixA.length;
    final Quadruple productElement = new Quadruple();
    for (int k = 0; k < size; k++) {
      productElement.add(new Quadruple(matrixA[rowIndex][k]).multiply(matrixB[k][colIndex]));
    }
    return productElement.doubleValue();
  }

  private static double[] multiply(double[][] matrix, double[] vector) {
    final int length = matrix.length;
    final double[] result = new double[length];
    for (int i = 0; i < length; i++) {
      final Quadruple productElement = new Quadruple();
      for (int j = 0; j < length; j++) {
        productElement.add(new Quadruple(matrix[i][j]).multiply(vector[j]));
      }
      result[i] = productElement.doubleValue();
    }
    return result;
  }

  private static Quadruple[] multiply(Quadruple[][] matrix, Quadruple[] vector) {
    final int length = matrix.length;
    final Quadruple[] result = new Quadruple[length];
    for (int i = 0; i < length; i++) {
      final Quadruple productElement = new Quadruple();
      for (int j = 0; j < length; j++) {
        productElement.add(new Quadruple(matrix[i][j]).multiply(vector[j]));
      }
      result[i] = productElement;
    }
    return result;
  }

  private static BigDecimal[] multiply(BigDecimal[][] matrix, BigDecimal[] vector) {
    final int length = matrix.length;
    final BigDecimal[] result = new BigDecimal[length];
    for (int i = 0; i < length; i++) {
      BigDecimal productElement = BigDecimal.ZERO;
      for (int j = 0; j < length; j++) {
        productElement = productElement.add(matrix[i][j].multiply(vector[j], MC), MC);
      }
      result[i] = productElement;
    }
    return result;
  }

  private static Quadruple[][] multiply(Quadruple[][] matrixA, Quadruple[][] matrixB) {
    final int length = matrixA.length;
    final Quadruple[][] result = new Quadruple[length][length];
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < length; j++) {
        result[i][j] = multiplyRowByColumn(matrixA, i, matrixB, j);
      }
    }
    return result;
  }

  private static BigDecimal[][] multiply(BigDecimal[][] matrixA, BigDecimal[][] matrixB) {
    final int length = matrixA.length;
    final BigDecimal[][] result = new BigDecimal[length][length];
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < length; j++) {
        result[i][j] = multiplyRowByColumn(matrixA, i, matrixB, j);
      }
    }
    return result;
  }

  private static Quadruple multiplyRowByColumn(Quadruple[][] matrixA, int rowIndex, Quadruple[][] matrixB, int colIndex) {
    final int size = matrixA.length;
    final Quadruple productElement = new Quadruple();
    for (int k = 0; k < size; k++) {
      productElement.add(Quadruple.multiply(matrixA[rowIndex][k], matrixB[k][colIndex]));
    }
    return productElement;
  }

  private static BigDecimal multiplyRowByColumn(BigDecimal[][] matrixA, int rowIndex, BigDecimal[][] matrixB, int colIndex) {
    final int size = matrixA.length;
    BigDecimal productElement = BigDecimal.ZERO;
    for (int k = 0; k < size; k++) {
      productElement = productElement.add(matrixA[rowIndex][k].multiply(matrixB[k][colIndex]));
    }
    return productElement;
  }

} // class MatrixData {