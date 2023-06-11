package com.mvohm.quadmatrix.investigations;

import static com.mvohm.quadmatrix.test.MatrixDataGenerators.*;

import java.util.Random;

import com.mvohm.quadmatrix.Matrix;

class MatrixData {

  /**
   * The types of matrices that can be used for testing
   */
  enum MatrixType {
    RANDOM_UNIFORM_DENSE,
    RANDOM_NONUNIFORM_DENSE,
    RANDOM_UNIFORM_SPARSE,
    RANDOM_NONUNIFORM_SPARSE,
    SPD_DENSE,
    SPD_SPARSE,
  }

  private double[] vector;
  private double[] solution;
  private double[][] matrix;

  private static int    matrixSize;
  private static double matrixDensity;
  private static double rangeBottom;
  private static double rangeTop;
  private static double randPower;
  private static double randRatio;

  static void setMatrixParams(int size, double matrixDensity, double rangeBottom, double rangeTop,
                              double randPower, double randRatio) {
    MatrixData.matrixSize     = size;
    MatrixData.matrixDensity  = matrixDensity;
    MatrixData.rangeBottom    = rangeBottom;
    MatrixData.rangeTop       = rangeTop;
    MatrixData.randPower      = randPower;
    MatrixData.randRatio      = randRatio;
  }

  public double[][] getMatrix() {   return matrix; }
  public double[]   getVector() {   return vector; }
  public double[]   getSolution() { return solution; }


  static MatrixData generate(MatrixType matrixType) {
    MatrixData data = new MatrixData();
    switch(matrixType) {
      case RANDOM_UNIFORM_DENSE:
        return data.generateRandom();
      case RANDOM_UNIFORM_SPARSE:
        return data.generateSparse(matrixDensity);
      case RANDOM_NONUNIFORM_DENSE:
        return data.generatePowPlusLinear(randPower, randRatio);
      case RANDOM_NONUNIFORM_SPARSE:
        return data.generateSparsePowPlusLinear(matrixDensity, randPower, randRatio);
      case SPD_DENSE:
        return data.generateSPD();
      case SPD_SPARSE:
        return data.generateSPDSparse(matrixDensity);
    }
    return null;
  }

  /**
   * Generates a dense random matrix with uniformly-distributed values ranged from -1 to 1
   * @return
   */
  MatrixData generateRandom() {
    solution = randomVector(matrixSize, rangeBottom, rangeTop);
    matrix   = randomMatrix(matrixSize, rangeBottom, rangeTop);
    vector   = findVector(matrix, solution);
    return this;
  }

  MatrixData generateSPD() {
    solution = randomVector(matrixSize, rangeBottom, rangeTop);
    matrix   = randomSpdMatrix(matrixSize, rangeBottom, rangeTop, 1.0);
    vector   = findVector(matrix, solution);
    return this;
  }

  MatrixData generateSPDSparse(double density) {
    solution = randomVector(matrixSize, rangeBottom, rangeTop);
    matrix   = randomSpdMatrix(matrixSize, rangeBottom, rangeTop, density);
    vector   = findVector(matrix, solution);
    return this;
  }

  /**
   * Generates a dense random matrix with a non-unform distribution, where values
   * are computed as<pre>
   *    x = random(-1, 1);
   *    y = x > 0? x ^ power: -(-x^power);
   *    return y * ratio + x * (1 - ratio)
   * </pre>
   * @param power
   * @param ratio
   * @return
   */
  MatrixData generatePowPlusLinear(double power, double ratio) {
    solution = randomVector(matrixSize, rangeBottom, rangeTop);
    matrix   = randomPowPlusLinearMatrix(matrixSize, power, ratio);
    vector   = findVector(matrix, solution);
    return this;
  }

  MatrixData generateSparse(double density) {
    solution = randomVector(matrixSize, rangeBottom, rangeTop);
    matrix   = randomSparseMatrix(matrixSize, density);
    vector   = findVector(matrix, solution);
    return this;
  }

  MatrixData generateSparsePowPlusLinear(double density, double power, double ratio) {
    solution = randomVector(matrixSize, rangeBottom, rangeTop);
    matrix   = randomSparsePowPlusLinearMatrix(matrixSize, density, power, ratio);
    vector   = findVector(matrix, solution);
    return this;
  }

} // static class MatrixData {