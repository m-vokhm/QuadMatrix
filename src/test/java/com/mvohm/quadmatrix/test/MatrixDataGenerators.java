package com.mvohm.quadmatrix.test;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Random;

import com.mvohm.quadruple.Quadruple;

/* Started 2022-01-27 16:56:11
 * An appliance to generate matrices and vectors of various kinds
 *
 */

public class MatrixDataGenerators {

  // Used to generate random SPD-matrices.
  // The value is chosen to produce an SPD matrix with a reasonable condition number
  private static final double SPD_FACTOR = 0.019401;


  private static Random random = new Random(12345);

  /** Set random seed to provide reproducibility */
  public static void setRandomSeed(int seed) {
    random = (seed < 0)? new Random() : new Random(seed);
  }


  public static void setRandomSeed(Random seedContainer) {
    random = seedContainer;
  }

  /**
   * A dense array filled with random values ranged from 0 to 1.0
   * @param length
   * @return
   */
  public static double[] randomVector(int length) {
    final double[] result = new double[length];
    for (int i = 0; i < length; i++)
      result[i] = random.nextDouble();
    return result;
  }

  /**
   * A dense array filled with random values ranged from rangedFrom to rangedTo
   * @param length
   * @return
   */
  public static double[] randomVector(int length, double rangedFrom, double rangedTo) {
    final double[] result = new double[length];
    for (int i = 0; i < length; i++)
      result[i] = randomRanged(rangedFrom, rangedTo);
    return result;
  }


  /**
   * A dense array filled with random values with gaussian distribution ranged from -1.0 to 1.0
   * @param length
   * @return
   */
  public static double[] randomGaussianVector(int length) {
    final double[] result = new double[length];
    for (int i = 0; i < length; i++)
      result[i] = gaussianRandom(random);
    return result;
  }

  public static double[] randomPowPlusLinearVector(int length, double power, double ratio) {
    final double[] result = new double[length];
    for (int i = 0; i < length; i++)
      result[i] = randPowPlusLinear(random, power, ratio);
    return result;
  }

  public static double[] randomSparseVector(int length, double density) {
    final double[] result = new double[length];
    for (int i = 0; i < length; i++) {
      if (random.nextDouble() < density)
        result[i] = randomRanged(-1, 1);
      else
        result[i] = 0;
      }
    return result;
  }

  public static double[] randomSparsePowPlusLinearVector(int length, double density, double power, double slope) {
    final double[] result = new double[length];
    for (int i = 0; i < length; i++) {
      if (random.nextDouble() < density)
        result[i] = randPowPlusLinear(random, power, slope);
      else
        result[i] = 0;
      }
    return result;
  }

  /**
   * x * ratio + x^power * (1 - ratio)
   * @param random
   * @return
   */
  private static double randPowPlusLinear(Random random, double power, double slope) {
    final double x = -1 + 2.0 * random.nextDouble(); // Ftom -1 to +1
    final double y = x < 0? -Math.pow(-x, power) : Math.pow(x, power); // non-linear
    return x * slope + y * (1 - slope);
  }


  public static double gaussianRandom(Random random) {
    double r, x, y;

    do {                               // find a uniform random point (x, y) inside unit circle
       x = 2.0 * random.nextDouble() - 1.0;
       y = 2.0 * random.nextDouble() - 1.0;
       r = x*x + y*y;
    } while (r > 1 || r == 0);         // loop executed 4 / pi = 1.273.. times on average
                                       // http://en.wikipedia.org/wiki/Box-Muller_transform

    // apply the Box-Muller formula to get standard Gaussian z
//    return x * Math.sqrt(-2.0 * Math.log(r) / r) / 5;
    return x * Math.sqrt(-2.0 * Math.log(r) / r) / 5;
  }



  /**
   * @param rangedFrom
   * @param rangedTo
   * @return
   */
  private static double randomRanged(double rangedFrom, double rangedTo) {
    return random.nextDouble() * (rangedTo - rangedFrom) + rangedFrom;
  }

  /**
   * A dense matrix filled with random values in range (0.0, 1.0]
   * @param size
   * @return
   */
  public static double[][] randomMatrix(int size) {
    final double[][] result = new double[size][];
    for (int i = 0; i < size; i++)
      result[i] = randomVector(size);
    return result;
  }

  public static double[][] randomMatrix(int size, double rangedFrom, double rangedTo) {
    final double[][] result = new double[size][];
    for (int i = 0; i < size; i++)
      result[i] = randomVector(size, rangedFrom, rangedTo);
    return result;
  }

  public static double[][] randomGaussianMatrix(int size) {
    final double[][] result = new double[size][];
    for (int i = 0; i < size; i++)
      result[i] = randomGaussianVector(size);
    return result;
  }

  public static double[][] randomPowPlusLinearMatrix(int size, double power, double ratio) {
    final double[][] result = new double[size][];
    for (int i = 0; i < size; i++)
      result[i] = randomPowPlusLinearVector(size, power, ratio);
    return result;
  }

  public static double[][] randomSparseMatrix(int size, double density) {
    final double[][] result = new double[size][];
    for (int i = 0; i < size; i++)
      result[i] = randomSparseVector(size, density);
    return result;
  }

  public static double[][] randomSparsePowPlusLinearMatrix(int size, double density, double power, double slope) {
    final double[][] result = new double[size][];
    for (int i = 0; i < size; i++)
      result[i] = randomSparsePowPlusLinearVector(size, density, power, slope);
    return result;
  }

  public static double[][] randomSpdMatrix(int size, double rangedFrom, double rangedTo, double density) {
    final double[][] matrix = new double[size][size];
    final double[][] lower = new double[size][size];

    // Fill lower triangle with random numbers
    for (int i = 0; i < size; i++) {
      for (int j  = 0; j < i; j++) {
        if (density == 1.0)
          lower[i][j] = randomRanged(rangedFrom, rangedTo);
        else if (random.nextDouble() > density)
          lower[i][j] = 0;
        else
          lower[i][j] = randomRanged(rangedFrom, rangedTo);
      }
      lower[i][i] = random.nextDouble();
    }

    // Cholesky decomposition backwards
    for (int i = 0; i < size; i++) {
      double sum2 = 0;

      for (int j = 0; j < i; j++) {
        double sum = 0;
        for (int k = 0; k < j; k++)
          sum += lower[i][k] * lower[j][k];

        matrix[i][j] = matrix[j][i] = lower[i][j] * lower[j][j] + sum;
        sum2 += lower[i][j] * lower[i][j];
      }

      final double dd = lower[i][i] * lower[i][i];
      matrix[i][i] = sum2 + dd + SPD_FACTOR;
    }
    return matrix;
  }


  /**
   * Element-wise multiplies vectors
   * @param a
   * @param b
   * @return product
   */
  public static double[] multiplyElements(double[] a, double[] b) {
    final int length = a.length;
    final double[] result = new double[length];
    for (int i = 0; i < length; i++)
      result[i] = a[i] * b[i];
    return result;
  }

  public static double[] multiply(double[][] matrix, double[] vector) {
    final int length = matrix.length;
    final double[] result = new double[length];
    final double[] productRow = new double[length];
    for (int i = 0; i < length; i++) {
      if (matrix[i] == null || matrix[i].length != length)
        throw new IllegalArgumentException("the matrix must be square and it's size must be equal to the vector's length");
      for (int j = 0; j < length; j++) {
        productRow[j] = matrix[i][j] * vector[j];
      }
      result[i] = sumOfVector(productRow);
    }
    return result;
  }

  // added 2022-12-24 11:23:55 to test solutions with respect of matrices
  public static double[][] multiply(double[][] matrixA, double[][] matrixB) {
    final int size = matrixA.length;
    final double[][] product = new double[size][size];
    final double[] vectorProduct = new double[size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        for (int k = 0; k < size; k++) {
          vectorProduct[k] = matrixA[i][k] * matrixB[k][j];
        }
        product[i][j] = sumOfVector(vectorProduct);
      }
    }
    return product;
  }

  public static double[][] fastMultiply(double[][] matrixA, double[][] matrixB) {
    final int size = matrixA.length;
    final double[][] product = new double[size][size];
    final double[] vectorProduct = new double[size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        for (int k = 0; k < size; k++) {
          vectorProduct[k] = matrixA[i][k] * matrixB[k][j];
        }
        product[i][j] = fastSumOfVector(vectorProduct); // Without sort
      }
    }
    return product;
  }

  /**
   * Computes the sum of the vector.
   * Kahan summation after sorting by absolute values -- relatively fast while relatively accurate
   * @param vector an array of doubles to sum
   * @return the value of the sum
   */
  public static double sumOfVector(double[] vector) {
    sortByAbs(vector) ;

    double sum = 0.0;
    double c = 0.0;
    for (int i = 0; i < vector.length; i++) {
      final double y = vector[i] - c;
      final double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    return sum;
  }

  public static double fastSumOfVector(double[] vector) {
    sortByAbs(vector) ;

    double sum = 0.0;
    double c = 0.0;
    for (int i = 0; i < vector.length; i++) {
      final double y = vector[i] - c;
      final double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    return sum;
  }

  /**
   * A QuickSort variation comparing items by absolute values
   * Used by {@link #sumOfVector(double[])}
   * @param vector
   */
  private static void sortByAbs(double[] vector) {
    quickSortByAbs(vector, 0, vector.length - 1);
  }

  ///*************************************
  // Taken from F:\.eclipse\Lopt\TestSorts\src\stackoverflow\QSortBenchmarks
  // 2022-11-27 20:06:56
  // Slightly improved 2022-11-28 19:16:39

  private static void quickSortByAbs(double[] array, int lo, int hi) {
    if (hi <= lo) return;
    final int j = partitionByAbs(array, lo, hi);
    quickSortByAbs(array, lo, j - 1);
    quickSortByAbs(array, j + 1, hi);
  }

  private static int partitionByAbs(double[] array, int lo, int hi) {
    int i = lo, j = hi + 1;
    while (true) {
      final double loValueAbs = Math.abs(array[lo]);
      while (Math.abs(array[++i]) < loValueAbs)
        if (i == hi) break;
      while (loValueAbs < Math.abs(array[--j]))
        if (j == lo) break;
      if (i >= j) break;

      final double tmp = array[i];
      array[i] = array[j];
      array[j] = tmp;
    }

    final double tmp = array[lo];
    array[lo] = array[j];
    array[j] = tmp;

    return j;
  }

  //
  //****************************************************


  public static double[][] unityMatrix(int size) {
    final double[][] result = new double[size][size];
    for (int i = 0; i < size; i++) {
      result[i][i] = 1.0;
    }
    return result;
  }

  public static Quadruple[][] quadrupleUnityMatrix(int size) {
    final Quadruple[][] result = new Quadruple[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        if (i == j)
          result[i][j] = Quadruple.one();
        else
          result[i][j] = new Quadruple();
      }
    }
    return result;
  }

  public static BigDecimal[][] bigDecimalUnityMatrix(int size) {
    final BigDecimal[][] result = new BigDecimal[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        if (i == j)
          result[i][j] = BigDecimal.ONE;
        else
          result[i][j] = BigDecimal.ZERO;
      }
    }
    return result;
  }

  public static double[] unityVector(int length) {
    final double[] vector = new double[length];
    for (int i = 0; i < length; i++) {
      vector[i] = 1.0;
    }
    return vector;
  }

  public static double[][] subtractMatrices(double[][] minuend, double[][] subtrahend) {
    final int size = minuend.length;
    final double[][] result = new double[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++)
        result[i][j] = minuend[i][j] - subtrahend[i][j];
    }
    return result;
  }

  public static Quadruple[][] subtractMatrices(Quadruple[][] minuend, Quadruple[][] subtrahend) {
    final int size = minuend.length;
    final Quadruple[][] result = new Quadruple[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++)
        result[i][j] = Quadruple.subtract(minuend[i][j], subtrahend[i][j]);
    }
    return result;
  }

  public static BigDecimal[][] subtractMatrices(BigDecimal[][] minuend, BigDecimal[][] subtrahend, MathContext mc) {
    final int size = minuend.length;
    final BigDecimal[][] result = new BigDecimal[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++)
        result[i][j] = minuend[i][j].subtract(subtrahend[i][j], mc);
    }
    return result;
  }

  /**
   * Finds the vector for the given matrix supposing that all the elements of the solution are 1.0
   * @param matrixData
   * @return
   */
  public static double[] findVector(double[][] matrixData) {
    final int size = matrixData.length;
    final double[] vector  = new double[size];
    for (int i = 0; i < size; i++) {
      vector[i] = sumOfVector(matrixData[i].clone()); // sumOfVector spoils the argument
    }
    return vector;
  }

  /**
   * Finds the vector for the given matrix and the given solution
   * @param matrixData
   * @return
   */
  public static double[] findVector(double[][] matrixData, double[] roots) {
    final int size = matrixData.length;
    final double[] vector  = new double[size];
    for (int i = 0; i < size; i++) {
      final double[] product = multiplyElements(roots, matrixData[i]);
      vector[i] = sumOfVector(product);
    }
    return vector;
  }

}
