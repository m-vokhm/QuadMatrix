package com.mvohm.quadmatrix.test;

import static com.mvohm.quadmatrix.test.AuxMethods.*;

import java.math.BigDecimal;
import java.util.Random;

import com.mvohm.quadruple.Quadruple;

/*
 * Started 2022-01-18 19:58:25
 *
 * Contains data samples for testing matrix calculations
 * */

public class MatrixDataSamples {

  static Random random;
  private static final int RAND_SEED = 432191234;

  static void resetRandom() {
    random = new Random(RAND_SEED);
  }

  /* ***************************************************************************
   ** Data samples (invertible matrix data, vectors and respective solutions) **
   ** to test solving a system of linear equations using LU-decomposition     **
   *****************************************************************************/

  public static double[][] sampleSolvableMatrixData() {
    return new double[][] {
      {1.000,  2.000,  3.000,  1.000,  2.000, },
      {1.000,  2.000,  1.000,  2.000,  1.000, },
      {2.000,  1.000,  0.000,  1.000,  2.000, },
      {2.000,  1.000,  3.000,  2.000,  1.000, },
      {2.000,  2.000,  1.000,  1.000,  5.000, },
    };
  }

  public static Quadruple[][] sampleSolvableMatrixDataAsQuadruples() {
    return new Quadruple[][] {
      {new Quadruple("1.000"),  new Quadruple("2.000"),  new Quadruple(" 3.000"),  new Quadruple(" 1.000"),  new Quadruple(" 2.000"), },
      {new Quadruple("1.000"),  new Quadruple("2.000"),  new Quadruple(" 1.000"),  new Quadruple(" 2.000"),  new Quadruple(" 1.000"), },
      {new Quadruple("2.000"),  new Quadruple("1.000"),  new Quadruple(" 0.000"),  new Quadruple(" 1.000"),  new Quadruple(" 2.000"), },
      {new Quadruple("2.000"),  new Quadruple("1.000"),  new Quadruple(" 3.000"),  new Quadruple(" 2.000"),  new Quadruple(" 1.000"), },
      {new Quadruple("2.000"),  new Quadruple("2.000"),  new Quadruple(" 1.000"),  new Quadruple(" 1.000"),  new Quadruple(" 5.000"), },
    };
  }

  public static BigDecimal[][] sampleSolvableMatrixDataAsBigDecimals() {
    return new BigDecimal[][] {
      { BigDecimal.valueOf(1.000), BigDecimal.valueOf(2.000), BigDecimal.valueOf( 3.000), BigDecimal.valueOf( 1.000), BigDecimal.valueOf( 2.000), },
      { BigDecimal.valueOf(1.000), BigDecimal.valueOf(2.000), BigDecimal.valueOf( 1.000), BigDecimal.valueOf( 2.000), BigDecimal.valueOf( 1.000), },
      { BigDecimal.valueOf(2.000), BigDecimal.valueOf(1.000), BigDecimal.valueOf( 0.000), BigDecimal.valueOf( 1.000), BigDecimal.valueOf( 2.000), },
      { BigDecimal.valueOf(2.000), BigDecimal.valueOf(1.000), BigDecimal.valueOf( 3.000), BigDecimal.valueOf( 2.000), BigDecimal.valueOf( 1.000), },
      { BigDecimal.valueOf(2.000), BigDecimal.valueOf(2.000), BigDecimal.valueOf( 1.000), BigDecimal.valueOf( 1.000), BigDecimal.valueOf( 5.000), },
    };
  }

  public static double[][] sampleInconsistentMatrixData() {
    return new double[][] {
      {1.000,  2.000,  3.000,  1.000,  2.000, },
      {1.000,  2.000,  1.000,  2.000,  1.000, },
      {2.000,  1.000,  0.000,  1.000,  2.000, },
      {2.000,  1.000,  3.000,  2.000,  1.000, },
      {2.000,  2.000,  1.000,  1.000,  3.000, },
    };
  }

  public static final double sampleInconsistentMatrixDataDeterminant = 0;

  public static double[][] sampleUnderdeterminedMatrixData() {
    return new double[][] {
      {1.000,  2.000,  3.000,  1.000,  2.000, },
      {1.000,  2.000,  1.000,  2.000,  1.000, },
      {2.000,  1.000,  0.000,  1.000,  2.000, },
      {2.000,  1.000,  3.000,  2.000,  1.000, },
      {4.000,  2.000,  0.000,  2.000,  4.000, }, // Row 4 is row 2 * 2
    };
  }

  public static double[] sampleSolvableMatrixVector() {
    return new double[] {
      16.000,
       9.000,
       6.000,
      14.000,
      14.000,
    };
  }

  public static double[] sampleSolvableMatrixAnotherVector() {
    return new double[] {
      12.000,
      10.000,
      11.500,
      12.000,
      19.500,
    };
  }

  public static double[] sampleTooLongMatrixVector() {
    return new double[] {
      18.000,
      13.000,
       8.000,
      18.000,
      16.000,
       1.000,
    };
  }

  public static double[] sampleSolvableMatrixSolution() {
    return new double[] {
       1.000,
       2.000,
       3.000,
       0.000,
       1.000,
    };
  }

  public static double[] sampleSolvableMatrixAnotherSolution() {
    return new double[] {
       2.500,
       1.500,
       0.500,
       1.000,
       2.000,
    };
  }

  public static Quadruple[] sampleSolvableMatrixAnotherSolutionAsQuadruples() {
    return new Quadruple[] {
       new Quadruple("2.500"),
       new Quadruple("1.500"),
       new Quadruple("0.500"),
       new Quadruple("1.000"),
       new Quadruple("2.000"),
    };
  }

  public static BigDecimal[] sampleSolvableMatrixAnotherSolutionAsBigDecimals() {
    return new BigDecimal[] {
       new BigDecimal("2.500"),
       new BigDecimal("1.500"),
       new BigDecimal("0.500"),
       new BigDecimal("1.000"),
       new BigDecimal("2.000"),
    };
  }

  public static double[][] sampleSolvableMatrixInversion() {
    return new double[][] {
      { 0.40000000000000000, -0.40000000000000000,  1.10000000000000000, -0.10000000000000000, -0.50000000000000000, },
      { 0.80000000000000000,  0.20000000000000000,  0.70000000000000000, -0.70000000000000000, -0.50000000000000000, },
      { 0.20000000000000000, -0.20000000000000000, -0.20000000000000000,  0.20000000000000000,  0.00000000000000000, },
      {-0.93333333333333333,  0.60000000000000000, -0.90000000000000000,  0.56666666666666667,  0.50000000000000000, },
      {-0.33333333333333333,  0.00000000000000000, -0.50000000000000000,  0.16666666666666667,  0.50000000000000000, },
    };
  }

  public static double[][] sampleInvertibleMatrixData() {
    return new double[][] {
      {1.000,  2.000,  3.000,  1.000,  2.000, },
      {1.000,  2.000,  1.000,  2.000,  1.000, },
      {2.000,  1.000,  3.000,  1.000,  2.000, },
      {2.000,  1.000,  3.000,  2.000,  1.000, },
      {2.000,  2.000,  1.000,  1.000,  5.000, },
    };
  };

  public static double[][] sampleInvertibleMatrixInversion() {
    return new double[][] {
      {-1.25000000000000000,  1.25000000000000000,  2.75000000000000000, -1.75000000000000000, -0.50000000000000000, },
      {-0.25000000000000000,  1.25000000000000000,  1.75000000000000000, -1.75000000000000000, -0.50000000000000000, },
      { 0.50000000000000000, -0.50000000000000000, -0.50000000000000000,  0.50000000000000000,  0.00000000000000000, },
      { 0.41666666666666667, -0.75000000000000000, -2.25000000000000000,  1.91666666666666667,  0.50000000000000000, },
      { 0.41666666666666667, -0.75000000000000000, -1.25000000000000000,  0.91666666666666667,  0.50000000000000000, },
    };
  };

  public static Quadruple[][] sampleInvertibleMatrixInversionAsQuadruples() {
    return new Quadruple[][] {
      { new Quadruple("-5").divide(4 ),  new Quadruple("  5").divide(4),  new Quadruple(" 11").divide(4),  new Quadruple(" -7").divide(4 ),  new Quadruple(" -1").divide(2), },
      { new Quadruple("-1").divide(4 ),  new Quadruple("  5").divide(4),  new Quadruple("  7").divide(4),  new Quadruple(" -7").divide(4 ),  new Quadruple(" -1").divide(2), },
      { new Quadruple(" 1").divide(2 ),  new Quadruple(" -1").divide(2),  new Quadruple(" -1").divide(2),  new Quadruple("  1").divide(2 ),  new Quadruple(0),               },
      { new Quadruple(" 5").divide(12),  new Quadruple(" -3").divide(4),  new Quadruple(" -9").divide(4),  new Quadruple(" 23").divide(12),  new Quadruple("  1").divide(2), },
      { new Quadruple(" 5").divide(12),  new Quadruple(" -3").divide(4),  new Quadruple(" -5").divide(4),  new Quadruple(" 11").divide(12),  new Quadruple("  1").divide(2), },
    };
  };

  public static BigDecimal[][] sampleInvertibleMatrixInversionAsBigDecimals() {
    return convertToBigDecimals(sampleInvertibleMatrixInversionAsQuadruples());
  };

  public static final double sampleInvertibleMatrixConditionNumber = 82.5;

  public static double[][] sampleInvertibleMatrixTransposition() {
    return new double[][] {
      {1.000,  1.000,  2.000,  2.000,  2.000, },
      {2.000,  2.000,  1.000,  1.000,  2.000, },
      {3.000,  1.000,  3.000,  3.000,  1.000, },
      {1.000,  2.000,  1.000,  2.000,  1.000, },
      {2.000,  1.000,  2.000,  1.000,  5.000, },
    };
  };

  // Added 2022-12-22 13:14:22 to test solutions with respect of matrix, A * X = B
  public static double[][] sampleSolvableMatrixMatrixX() {
    return new double[][] {
      {  3.000,  -2.000,   3.000,   6.000,   2.000, },
      {  2.000,   3.000,  -2.000,   1.000,  -1.000, },
      {  2.000,  -2.000,   3.000,   3.000,   1.000, },
      {  2.000,   4.000,  -3.000,   2.000,   1.000, },
      { -2.000,   1.000,   2.000,   1.000,  -2.000, },
    };
  }

  public static Quadruple[][] sampleSolvableMatrixMatrixXAsQuadruples() {
    return new Quadruple[][] {
      { new Quadruple(" 3.000"), new Quadruple("-2.000"), new Quadruple(" 3.000"), new Quadruple(" 6.000"), new Quadruple(" 2.000") },
      { new Quadruple(" 2.000"), new Quadruple(" 3.000"), new Quadruple("-2.000"), new Quadruple(" 1.000"), new Quadruple("-1.000") },
      { new Quadruple(" 2.000"), new Quadruple("-2.000"), new Quadruple(" 3.000"), new Quadruple(" 3.000"), new Quadruple(" 1.000") },
      { new Quadruple(" 2.000"), new Quadruple(" 4.000"), new Quadruple("-3.000"), new Quadruple(" 2.000"), new Quadruple(" 1.000") },
      { new Quadruple("-2.000"), new Quadruple(" 1.000"), new Quadruple(" 2.000"), new Quadruple(" 1.000"), new Quadruple("-2.000") },
    };
  }

  public static BigDecimal[][] sampleSolvableMatrixMatrixXAsBigDecimals() {
    return new BigDecimal[][] {
      { BigDecimal.valueOf( 3.000), BigDecimal.valueOf(-2.000), BigDecimal.valueOf( 3.000), BigDecimal.valueOf( 6.000), BigDecimal.valueOf( 2.000) },
      { BigDecimal.valueOf( 2.000), BigDecimal.valueOf( 3.000), BigDecimal.valueOf(-2.000), BigDecimal.valueOf( 1.000), BigDecimal.valueOf(-1.000) },
      { BigDecimal.valueOf( 2.000), BigDecimal.valueOf(-2.000), BigDecimal.valueOf( 3.000), BigDecimal.valueOf( 3.000), BigDecimal.valueOf( 1.000) },
      { BigDecimal.valueOf( 2.000), BigDecimal.valueOf( 4.000), BigDecimal.valueOf(-3.000), BigDecimal.valueOf( 2.000), BigDecimal.valueOf( 1.000) },
      { BigDecimal.valueOf(-2.000), BigDecimal.valueOf( 1.000), BigDecimal.valueOf( 2.000), BigDecimal.valueOf( 1.000), BigDecimal.valueOf(-2.000) },
    };
  }


  public static final double sampleSolvableMatrixMatrixXDeterminant = -185;

  public static Quadruple[][] sampleSolvableMatrixMatrixBAsQuadruples() {
    return new Quadruple[][] {
      { new Quadruple("11.000"),  new Quadruple(" 4.000"),  new Quadruple(" 9.000"),  new Quadruple(" 21.000"),  new Quadruple(" 0.000"), },
      { new Quadruple("11.000"),  new Quadruple("11.000"),  new Quadruple("-2.000"),  new Quadruple(" 16.000"),  new Quadruple(" 1.000"), },
      { new Quadruple(" 6.000"),  new Quadruple(" 5.000"),  new Quadruple(" 5.000"),  new Quadruple(" 17.000"),  new Quadruple(" 0.000"), },
      { new Quadruple("16.000"),  new Quadruple(" 2.000"),  new Quadruple(" 9.000"),  new Quadruple(" 27.000"),  new Quadruple(" 6.000"), },
      { new Quadruple(" 4.000"),  new Quadruple(" 9.000"),  new Quadruple("12.000"),  new Quadruple(" 24.000"),  new Quadruple("-6.000"), },
    };
  }

  public static BigDecimal[][] sampleSolvableMatrixMatrixBAsBigDecimal() {
    return new BigDecimal[][] {
      { BigDecimal.valueOf(11.000), BigDecimal.valueOf( 4.000), BigDecimal.valueOf( 9.000), BigDecimal.valueOf( 21.000), BigDecimal.valueOf( 0.000), },
      { BigDecimal.valueOf(11.000), BigDecimal.valueOf(11.000), BigDecimal.valueOf(-2.000), BigDecimal.valueOf( 16.000), BigDecimal.valueOf( 1.000), },
      { BigDecimal.valueOf( 6.000), BigDecimal.valueOf( 5.000), BigDecimal.valueOf( 5.000), BigDecimal.valueOf( 17.000), BigDecimal.valueOf( 0.000), },
      { BigDecimal.valueOf(16.000), BigDecimal.valueOf( 2.000), BigDecimal.valueOf( 9.000), BigDecimal.valueOf( 27.000), BigDecimal.valueOf( 6.000), },
      { BigDecimal.valueOf( 4.000), BigDecimal.valueOf( 9.000), BigDecimal.valueOf(12.000), BigDecimal.valueOf( 24.000), BigDecimal.valueOf(-6.000), },
    };
  }

  public static double[][] sampleSolvableMatrixMatrixB() {
    return new double[][] {
      { 11.000,   4.000,   9.000,  21.000,   0.000, },
      { 11.000,  11.000,  -2.000,  16.000,   1.000, },
      {  6.000,   5.000,   5.000,  17.000,   0.000, },
      { 16.000,   2.000,   9.000,  27.000,   6.000, },
      {  4.000,   9.000,  12.000,  24.000,  -6.000, },
    };
  }

  public static double[][] sampleSolvableMatrixAnotherMatrixX() {
    return new double[][] {
      { -4.065,   3.795,   0.099,   0.646,  -1.658, },
      {  1.280,  -0.817,  -1.842,   1.302,  -4.881, },
      {  3.823,  -1.878,   4.287,  -4.179,  -2.667, },
      { -1.503,   1.931,  -4.124,   0.794,   4.063, },
      {  2.060,   3.860,   4.602,  -3.587,  -4.034, },
    };
  }

  public static Quadruple[][] sampleSolvableMatrixAnotherMatrixXAsQuadruples() {
    return new Quadruple[][] {
      { new Quadruple("-4.065"),  new Quadruple(" 3.795"),  new Quadruple(" 0.099"),  new Quadruple(" 0.646"),  new Quadruple("-1.658"), },
      { new Quadruple(" 1.280"),  new Quadruple("-0.817"),  new Quadruple("-1.842"),  new Quadruple(" 1.302"),  new Quadruple("-4.881"), },
      { new Quadruple(" 3.823"),  new Quadruple("-1.878"),  new Quadruple(" 4.287"),  new Quadruple("-4.179"),  new Quadruple("-2.667"), },
      { new Quadruple("-1.503"),  new Quadruple(" 1.931"),  new Quadruple("-4.124"),  new Quadruple(" 0.794"),  new Quadruple(" 4.063"), },
      { new Quadruple(" 2.060"),  new Quadruple(" 3.860"),  new Quadruple(" 4.602"),  new Quadruple("-3.587"),  new Quadruple("-4.034"), },
    };
  }

  public static BigDecimal[][] sampleSolvableMatrixAnotherMatrixXAsBigDecimals() {
    return new BigDecimal[][] {
      { BigDecimal.valueOf(-4.065), BigDecimal.valueOf( 3.795), BigDecimal.valueOf( 0.099), BigDecimal.valueOf( 0.646), BigDecimal.valueOf(-1.658), },
      { BigDecimal.valueOf( 1.280), BigDecimal.valueOf(-0.817), BigDecimal.valueOf(-1.842), BigDecimal.valueOf( 1.302), BigDecimal.valueOf(-4.881), },
      { BigDecimal.valueOf( 3.823), BigDecimal.valueOf(-1.878), BigDecimal.valueOf( 4.287), BigDecimal.valueOf(-4.179), BigDecimal.valueOf(-2.667), },
      { BigDecimal.valueOf(-1.503), BigDecimal.valueOf( 1.931), BigDecimal.valueOf(-4.124), BigDecimal.valueOf( 0.794), BigDecimal.valueOf( 4.063), },
      { BigDecimal.valueOf( 2.060), BigDecimal.valueOf( 3.860), BigDecimal.valueOf( 4.602), BigDecimal.valueOf(-3.587), BigDecimal.valueOf(-4.034), },
    };
  }

  public static final double sampleSolvableMatrixAnotherMatrixXDeterminant = 1569.25279805215911;
  public static final double sampleSolvableMatrixAnotherMatrixXNorm = 18.143;


  public static double[][] sampleSolvableMatrixAnotherMatrixB() {
    return new double[][] {
      {  12.581,   6.178,  14.356, -15.667, -23.426, },
      {   1.372,   8.005,  -2.944,  -2.928,  -9.995, },
      {  -4.233,  16.424,   3.436,  -3.786, -12.202, },
      {   3.673,   8.861,   7.571, -11.942, -12.106, },
      {   7.050,  25.309,  19.687, -17.424, -31.852, },
    };
  }

  public static Quadruple[][] sampleSolvableMatrixAnotherMatrixBAsQuadruples() {
    return new Quadruple[][] {
      { new Quadruple(" 12.581"),  new Quadruple("  6.178"),  new Quadruple(" 14.356"),  new Quadruple("-15.667"),  new Quadruple("-23.426"), },
      { new Quadruple("  1.372"),  new Quadruple("  8.005"),  new Quadruple(" -2.944"),  new Quadruple(" -2.928"),  new Quadruple(" -9.995"), },
      { new Quadruple(" -4.233"),  new Quadruple(" 16.424"),  new Quadruple("  3.436"),  new Quadruple(" -3.786"),  new Quadruple("-12.202"), },
      { new Quadruple("  3.673"),  new Quadruple("  8.861"),  new Quadruple("  7.571"),  new Quadruple("-11.942"),  new Quadruple("-12.106"), },
      { new Quadruple("  7.050"),  new Quadruple(" 25.309"),  new Quadruple(" 19.687"),  new Quadruple("-17.424"),  new Quadruple("-31.852"), },
    };
  }

  public static double[][] sampleSolvableMatrixTooLargeMatrix() {
    return new double[][] {
      {  12.581,   6.178,  14.356, -15.667, -23.426, -23.426, },
      {   1.372,   8.005,  -2.944,  -2.928,  -9.995,  -9.995, },
      {  -4.233,  16.424,   3.436,  -3.786, -12.202, -12.202, },
      {   3.673,   8.861,   7.571, -11.942, -12.106, -12.106, },
      {   7.050,  25.309,  19.687, -17.424, -31.852, -31.852, },
      {   7.050,  25.309,  19.687, -17.424, -31.852, -31.852, },
    };
  }

  public static double[][] sampleSolvableMatrixNonSquareMatrix() {
    return new double[][] {
      {  12.581,   6.178,  14.356, -15.667, -23.426, -23.426, },
      {   1.372,   8.005,  -2.944,  -2.928,  -9.995,  -9.995, },
      {  -4.233,  16.424,   3.436,  -3.786, -12.202, -12.202, },
      {   3.673,   8.861,   7.571, -11.942, -12.106, -12.106, },
      {   7.050,  25.309,  19.687, -17.424, -31.852, -31.852, },
    };

  }

  /* *****************************************************************************
   ** Data samples (invertible matrix data, vectors and respective solutions)   **
   ** to test solving a system of linear equations using Cholesky-decomposition **
   *******************************************************************************/

  public static double[][] sampleSpdSolvableMatrixData() {
    return new double[][] {
      { 1.250, -2.000, -0.250, -1.000, -0.300, },
      {-2.000,  5.000,  1.000,  1.000,  0.500, },
      {-0.250,  1.000,  3.000, -1.500, -1.000, },
      {-1.000,  1.000, -1.500,  3.000,  1.000, },
      {-0.300,  0.500, -1.000,  1.000,  2.000, },
    };
  }

  public static Quadruple[][] sampleSpdSolvableMatrixDataAsQuadruples() {
    return new Quadruple[][] {
      {new Quadruple(" 1.250"), new Quadruple("-2.000"), new Quadruple("-0.250"), new Quadruple("-1.000"), new Quadruple("-0.300"), },
      {new Quadruple("-2.000"), new Quadruple(" 5.000"), new Quadruple(" 1.000"), new Quadruple(" 1.000"), new Quadruple(" 0.500"), },
      {new Quadruple("-0.250"), new Quadruple(" 1.000"), new Quadruple(" 3.000"), new Quadruple("-1.500"), new Quadruple("-1.000"), },
      {new Quadruple("-1.000"), new Quadruple(" 1.000"), new Quadruple("-1.500"), new Quadruple(" 3.000"), new Quadruple(" 1.000"), },
      {new Quadruple("-0.300"), new Quadruple(" 0.500"), new Quadruple("-1.000"), new Quadruple(" 1.000"), new Quadruple(" 2.000"), },
    };
  }
  public static final double sampleSpdMatrixDataDeterminant = 10.86;

  public static double[][] sampleNonPositivelyDefinedMatrixData() {
    return new double[][] {
      { 1.250, -2.000, -0.250, -1.000, -0.300, },
      {-2.000,  5.000,  1.000,  1.000,  0.500, },
      {-0.250,  1.000,  3.000, -1.500, -1.000, },
      {-1.000,  1.000, -1.500,  3.000,  1.000, },
      {-0.300,  0.500, -1.000,  1.000, -1.000, },
    };
  }

  public static double[][] sampleNonSymmetricDMatrixData() {
    return new double[][] {
      { 1.250, -3.000, -0.250, -1.000, -0.300, },
      {-2.000,  5.000,  1.000,  1.000,  0.500, },
      {-0.250,  1.000,  3.000, -1.500, -1.000, },
      {-1.000,  1.000, -1.500,  3.000,  1.000, },
      {-0.300,  0.500, -1.000,  1.000,  2.000, },
    };
  }

  public static double[] sampleSpdSolvableMatrixVector() {
    return new double[] {
      -3.400,
       8.500,
      -1.750,
       6.500,
       7.100,
    };
  }

  public static Quadruple[] sampleSpdSolvableMatrixVectorAsQuadruples() {
    return new Quadruple[] {
      new Quadruple("-3.400"),
      new Quadruple(" 8.500"),
      new Quadruple("-1.750"),
      new Quadruple(" 6.500"),
      new Quadruple(" 7.100"),
    };
  }

  public static double[] sampleSpdTooLongMatrixVector() {
    return new double[] {
      -3.400,
       8.500,
      -1.750,
       6.500,
       7.100,
       1.000,
    };
  }

  public static double[] sampleSpdSolvableMatrixAnotherVector() {
    return new double[] {
      -3.975,
       8.500,
       0.875,
       5.250,
       4.575,
      };
  }

  public static Quadruple[] sampleSpdSolvableMatrixAnotherVectorAsQuadruples() {
    return new Quadruple[] {
      new Quadruple("-3.975"),
      new Quadruple(" 8.500"),
      new Quadruple(" 0.875"),
      new Quadruple(" 5.250"),
      new Quadruple(" 4.575"),
    };
  }

  public static BigDecimal[] sampleSpdSolvableMatrixAnotherVectorAsBigDecimals() {
    return new BigDecimal[] {
      new BigDecimal("-3.975"),
      new BigDecimal("8.500"),
      new BigDecimal("0.875"),
      new BigDecimal("5.250"),
      new BigDecimal("4.575"),
    };
  }

  public static double[] sampleSpdSolvableMatrixSolution() {
    return new double[] {
       3.000,
       2.000,
       1.000,
       2.000,
       3.000,
    };
  }

  public static Quadruple[] sampleSpdSolvableMatrixSolutionAsQuadruples() {
    return new Quadruple[] {
       new Quadruple("3.000"),
       new Quadruple("2.000"),
       new Quadruple("1.000"),
       new Quadruple("2.000"),
       new Quadruple("3.000"),
    };
  }

  public static BigDecimal[] sampleSpdSolvableMatrixSolutionAsBigDecimals() {
    return new BigDecimal[] {
       new BigDecimal("3.000"),
       new BigDecimal("2.000"),
       new BigDecimal("1.000"),
       new BigDecimal("2.000"),
       new BigDecimal("3.000"),
    };
  }

  public static double[] sampleSpdSolvableMatrixAnotherSolution() {
    return new double[] {
       1.000,
       1.250,
       1.500,
       1.750,
       2.000,
    };
  }

  public static Quadruple[] sampleSpdSolvableMatrixAnotherSolutionAsQuadruples() {
    return new Quadruple[] {
       new Quadruple("1.000"),
       new Quadruple("1.250"),
       new Quadruple("1.500"),
       new Quadruple("1.750"),
       new Quadruple("2.000"),
    };
  }

  public static BigDecimal[] sampleSpdSolvableMatrixAnotherSolutionAsBigDecimals() {
    return new BigDecimal[] {
       new BigDecimal("1.000"),
       new BigDecimal("1.250"),
       new BigDecimal("1.500"),
       new BigDecimal("1.750"),
       new BigDecimal("2.000"),
    };
  }

  // matrixX and matrixB satisfy equation sampleSpdMatrix * matrixX = matrixB
  public  static double[][] sampleSpdMatrixMatrixX() {
    // Rounded randoms
    return new double[][] {
      {  2.401,  -1.029,  -0.871,  -0.803,   2.547, },
      { -2.576,   1.227,   2.937,   2.872,   2.489, },
      {  1.570,   0.606,   0.851,   0.619,  -2.513, },
      {  2.953,  -2.636,  -0.484,   0.199,  -1.227, },
      { -2.964,  -0.228,   0.244,   0.867,   0.573, },
    };
  }

  public  static Quadruple[][] sampleSpdMatrixMatrixXAsQuadruples() {
    return new Quadruple[][] {
      { new Quadruple(" 2.401"), new Quadruple(" -1.029"), new Quadruple(" -0.871"), new Quadruple(" -0.803"), new Quadruple("  2.547"), },
      { new Quadruple("-2.576"), new Quadruple("  1.227"), new Quadruple("  2.937"), new Quadruple("  2.872"), new Quadruple("  2.489"), },
      { new Quadruple(" 1.570"), new Quadruple("  0.606"), new Quadruple("  0.851"), new Quadruple("  0.619"), new Quadruple(" -2.513"), },
      { new Quadruple(" 2.953"), new Quadruple(" -2.636"), new Quadruple(" -0.484"), new Quadruple("  0.199"), new Quadruple(" -1.227"), },
      { new Quadruple("-2.964"), new Quadruple(" -0.228"), new Quadruple("  0.244"), new Quadruple("  0.867"), new Quadruple("  0.573"), },
    };
  }

  public  static BigDecimal[][] sampleSpdMatrixMatrixXAsBigDecimals() {
    return new BigDecimal[][] {
      { BigDecimal.valueOf( 2.401), BigDecimal.valueOf( -1.029), BigDecimal.valueOf( -0.871), BigDecimal.valueOf( -0.803), BigDecimal.valueOf(  2.547), },
      { BigDecimal.valueOf(-2.576), BigDecimal.valueOf(  1.227), BigDecimal.valueOf(  2.937), BigDecimal.valueOf(  2.872), BigDecimal.valueOf(  2.489), },
      { BigDecimal.valueOf( 1.570), BigDecimal.valueOf(  0.606), BigDecimal.valueOf(  0.851), BigDecimal.valueOf(  0.619), BigDecimal.valueOf( -2.513), },
      { BigDecimal.valueOf( 2.953), BigDecimal.valueOf( -2.636), BigDecimal.valueOf( -0.484), BigDecimal.valueOf(  0.199), BigDecimal.valueOf( -1.227), },
      { BigDecimal.valueOf(-2.964), BigDecimal.valueOf( -0.228), BigDecimal.valueOf(  0.244), BigDecimal.valueOf(  0.867), BigDecimal.valueOf(  0.573), },
    };
  }

  public  static double[][] sampleSpdMatrixMatrixB() {
    // Calculated by https://matrixcalc.org
    return new double[][] {
      { 113939.0 / 20000,  -23747.0 / 20000,  -67647.0 / 10000,  -4601.0 / 625,   -1109.0 / 10000, },
      { -14641.0 / 1000,     6049.0 / 1000,     4229.0 / 250,     6887.0 / 400,    1559.0 / 400,   },
      {    273.0 / 4000,    29937.0 / 4000,    24759.0 / 4000,   15057.0 / 4000, -17677.0 / 4000,  },
      {  -1437.0 / 1000,    -6789.0 / 1000,     2647.0 / 2000,    8421.0 / 2000,   1207.0 / 2000,  },
      { -65533.0 / 10000,  -13879.0 / 5000,     2207.0 / 2500,   29909.0 / 10000,  7281.0 / 2500,  },
    };
  }

  public  static Quadruple[][] sampleSpdMatrixMatrixBAsQuadruples() {
    // Calculated by https://matrixcalc.org
    return new Quadruple[][] {
      { new Quadruple("113939.0").divide(20000), new Quadruple("-23747.0").divide(20000),  new Quadruple("-67647.0").divide(10000), new Quadruple("  -4601.0").divide(625),  new Quadruple(" -1109.0").divide(10000), },
      { new Quadruple("-14641.0").divide(1000),  new Quadruple("  6049.0").divide(1000),   new Quadruple("  4229.0").divide(250),   new Quadruple("   6887.0").divide(400),  new Quadruple("  1559.0").divide(400),   },
      { new Quadruple("   273.0").divide(4000),  new Quadruple(" 29937.0").divide(4000),   new Quadruple(" 24759.0").divide(4000),  new Quadruple("  15057.0").divide(4000), new Quadruple("-17677.0").divide(4000),  },
      { new Quadruple(" -1437.0").divide(1000),  new Quadruple(" -6789.0").divide(1000),   new Quadruple("  2647.0").divide(2000),  new Quadruple("   8421.0").divide(2000), new Quadruple("  1207.0").divide(2000),  },
      { new Quadruple("-65533.0").divide(10000), new Quadruple("-13879.0").divide(5000),   new Quadruple("  2207.0").divide(2500),  new Quadruple("  29909.0").divide(10000),new Quadruple("  7281.0").divide(2500),  },
    };
  }


  public static double[][] sampleSpdMatrixInverse() {
    // Calculated by https://matrixcalc.org
    return new double[][] {
      { 14725.0 / 4344,    590.0 /  543,   1705.0 / 4344,   8545.0 /  8688,   -85.0 / 1448, },
      {   590.0 /  543,    344.0 /  543,    -46.0 /  543,     79.0 /   543,   -20.0 / 181,  },
      {  1705.0 / 4344,    -46.0 /  543,   2941.0 / 4344,   3733.0 /  8688,   295.0 / 1448, },
      {  8545.0 / 8688,     79.0 /  543,   3733.0 / 8688,  14989.0 / 17376,  -305.0 / 2896, },
      {   -85.0 / 1448,    -20.0 /  181,    295.0 / 1448,   -305.0 /  2896,   975.0 / 1448, },
    };
  }

  public static double[][] sampleHilbertMatrixData() {
    return new double[][] {
      { 1.00000000000000000,  0.50000000000000000,  0.33333333333333333,  0.25000000000000000,  0.20000000000000000, },
      { 0.50000000000000000,  0.33333333333333333,  0.25000000000000000,  0.20000000000000000,  0.16666666666666700, },
      { 0.33333333333333333,  0.25000000000000000,  0.20000000000000000,  0.16666666666666667,  0.14285714285714286, },
      { 0.25000000000000000,  0.20000000000000000,  0.16666666666666667,  0.14285714285714286,  0.12500000000000000, },
      { 0.20000000000000000,  0.16666666666666667,  0.14285714285714286,  0.12500000000000000,  0.11111111111111111, },
    };
  }

  public static double[][] sampleRandomMatrixData() {
    return new double[][] {
      { 0.54302772197448000,  0.85870530420809900,  0.93340589552204500,  0.01788412141919600,  0.67691286202237200, },
      { 0.71922323350908300,  0.36942446355663800,  0.80614928810454900,  0.85087784334787000,  0.09379288873539070, },
      { 0.40775483267354300,  0.05577523178755150,  0.23616878969391600,  0.75973306832430200,  0.08689931871514830, },
      { 0.24752272790281400,  0.93708944205528000,  0.67633394697146000,  0.39618559692542400,  0.60652724269951800, },
      { 0.28063605715027800,  0.73622118035091900,  0.20251964726670000,  0.86691714430168100,  0.93780274231541300, },
    };
  }

  public static final double sampleRandomMatrixDataDeterminant = -0.0027984761776598737;

  public static double[] sampleHilbertMatrixSolution() {
    return new double[] {
       1.000,
       1.000,
       1.000,
       1.000,
       1.000,
    };
  }

  public static double[] sampleRandomMatrixSolution() {
    return new double[] {
       1.000,
       1.000,
       1.000,
       1.000,
       1.000,
    };
  }

  public static double[] sampleHilbertMatrixVector() {
    return new double[] {
      2.28333333333333333,
      1.45000000000000000,
      1.09285714285714286,
      0.884523809523809524,
      0.745634920634920635,
    };
  }

  public static double[] sampleRandomMatrixVector() {
    final double[][] matrix = sampleRandomMatrixData();
    final double[] solution = sampleRandomMatrixSolution();
    final int size = solution.length;
    final double[] vector = new double[size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < 0; j++)
        vector[i] += solution[j] * matrix[i][j];
    }
    return vector;
  }

  public static double[][] summand1() {
    return new double[][] {
      {  -6.230588340,   4.361059587,  -3.153662625,  -8.055240788,  -3.066157090 },
      {  -5.382495544,   5.904873329,  -5.019136199,  -6.363408987,   3.882632664 },
      {  -5.962686302,   8.638076509,   2.331884917,  -5.791651757,   1.745354736 },
      {  -2.264694099,   5.545953749,  -8.763018758,   9.830576561,  -8.823486377 },
      {   5.714912404,  -7.511465851,  -2.812599297,  -2.427320764,   8.373626222 },
    };
  }

  public static Quadruple[][] summand1AsQuadruples() {
    return convertToQuadruples(summand1());
  }

  public static double[][] summand2() {
    return new double[][] {
      {  -3.390424638,  -6.008678105,   1.249308928,  -2.151673213,  -0.085802613 },
      {   1.045805278,   6.237436456,  -4.833196797,  -7.008106090,   2.049592078 },
      {   7.742105335,  -2.355957541,   7.891522917,   5.695333848,   1.437396772 },
      {  -3.521772995,   6.125780420,  -9.094063790,  -4.388552743,  -7.318364965 },
      {  -9.411986940,   7.084248416,  -7.507703330,   3.594828331,   2.354355943 },
    };
  }

  public static Quadruple[][] summand2AsQuadruples() {
    return convertToQuadruples(summand2());
  }

  public static double[][] summand1PlusSummand2() {
    final double[][] a1 = summand1();
    final double[][] a2 = summand2();
    final int size = a1.length;
    final double[][] result = new double[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        result[i][j] = a1[i][j] + a2[i][j];
      }
    }
    return result;
  }

  public static Quadruple[][] summand1PlusSummand2AsQuadruples() {
    final Quadruple[][] a1 = summand1AsQuadruples();
    final Quadruple[][] a2 = summand2AsQuadruples();
    final int size = a1.length;
    final Quadruple[][] result = new Quadruple[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        result[i][j] = Quadruple.add(a1[i][j], a2[i][j]);
      }
    }
    return result;
  }

  public static BigDecimal[][] summand1PlusSummand2AsBigDecimals() {
    final BigDecimal[][] a1 = convertToBigDecimals(summand1());
    final BigDecimal[][] a2 = convertToBigDecimals(summand2());
    final int size = a1.length;
    final BigDecimal[][] result = new BigDecimal[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        result[i][j] = a1[i][j].add(a2[i][j]);
      }
    }
    return result;
  }

  public static double[][] minuend() {
    return new double[][] {
      {  -0.315917634,  -0.032548440,   5.610668225,  -9.379446995,   4.822239459 },
      {  -8.863368286,  -9.484263831,  -3.157013621,   7.928548301,   1.758849970 },
      {  -8.785897178,   6.883451575,   5.847159829,   9.242688400,  -4.857402397 },
      {  -3.293869333,   4.638591953,  -5.820400186,  -6.791636357,  -7.047904448 },
      {   1.268971895,   4.253857385,  -8.400909457,  -8.571061665,   4.831135526 },
    };
  }

  public static double[][] subtrahend() {
    return new double[][] {
      {  -9.716375644,   6.408398918,  -1.393768058,  -6.350510803,   2.187966027 },
      {   2.711658944,  -8.868564573,   7.173928697,   2.790702326,   1.549108143 },
      {   1.324990131,   7.271653728,   8.275688323,  -8.497603557,   4.352750534 },
      {   7.261650389,  -9.037846281,  -4.418400618,  -7.792857417,   4.658355038 },
      {   3.160140868,  -5.294732251,  -4.997648567,   8.767727279,  -5.947035060 },
    };
  }

  public static double[][] minuendMinusSubtrahend() {
    final double[][] a1 = minuend();
    final double[][] a2 = subtrahend();
    final int size = a1.length;
    final double[][] result = new double[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        result[i][j] = a1[i][j] - a2[i][j];
      }
    }
    return result;
  }


  /* *****************************************************************************
   ** Private methods                                                           **
   *******************************************************************************/


}
