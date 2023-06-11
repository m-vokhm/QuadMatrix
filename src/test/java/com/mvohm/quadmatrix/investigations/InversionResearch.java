package com.mvohm.quadmatrix.investigations;

import static com.mvohm.quadmatrix.test.AuxMethods.say;
import static com.mvohm.quadmatrix.test.AuxMethods.say_;

import java.util.Random;

import com.mvohm.quadmatrix.DoubleMatrix;
import com.mvohm.quadmatrix.Matrix;
import com.mvohm.quadmatrix.investigations.MatrixData.MatrixType;
import com.mvohm.quadmatrix.investigations.SolutionRefinementResearch.ErrorSet;
import com.mvohm.quadmatrix.test.MatrixDataGenerators;

// import com.mvohm.quadmatrix.test.SolutionRefinementResearch.Statistics;

//
/**
 * Started 2022-12-18 17:42:07
 * Before writing tests for inversions, see how different types of inversion behave.
 * -- Simple classical inversion via LU-decomposition;
 * -- SPD-inversion based on the Cholesky decomposition;
 * -- Inversion imported from the JAMA library;
 * -- inverseAccurately(yet to be written)
 *
 * @author misa
 *
 */

public class InversionResearch {

  private static final boolean USING_SCALING = false;

  private static Random random = null;
  private static final int RAND_SEED        = 543210; // 123321; // Negative value to ignore constant RAND_SEED
  // private static final int RAND_SEED     = 123321;

  /** Choose the type of matrices to use this time */
  private static final MatrixType MATRIX_TYPE = MatrixType.RANDOM_UNIFORM_DENSE;

  private static final int MATRIX_SIZE      = 200;
  private static final double RANGE_BOTTOM  = -1.0;
  private static final double RANGE_TOP     = 1.0;
  private static final double RANGE         = RANGE_TOP - RANGE_BOTTOM;

  // For non-uniform distribution, we fill the matrix with values computed as
  //    x = ratio * (sign(r) * abs(r)^power) + (1 - ratio) * r, where 'r' is a double random, 'power' and 'ratio' are defined below:
  private static final double RAND_POWER    = 8.0;    // The power to raise a random to obtain a non-linear distribution for generating non-uniform matrices
  private static final double RAND_RATIO    = 0.95;   // A ratio to mix the power of the random and the random value itself,
  static final double MATRIX_DENSITY        = 0.2;    // density for generating sparse matrices

  private static final int MAX_TEST_ITERATIONS = 500;
  private static final int ERRORS_TO_SHOW   = 10;

  private static final long WARMUP_SECONDS      = 3;      // In seconds
  private static final long MILLIS_PER_SECOND   = 1_000;
  private static final long NANOS_PER_SECOND    = 1_000_000_000L;
  private static final long THREE_SECONDS       = 3 * NANOS_PER_SECOND;   // in nanoseconds



  private static long inversionTime;

  private static double[][] unityMatrix;

  public static void main(String[] args) {
    new InversionResearch().run();
  }

  private void run() {
    long t1 = -System.nanoTime();
    DoubleMatrix.setDefaultScaling(USING_SCALING);
    setMatrixParams();
    initUnityMatrix();

    watchInversions(new SimpleInversionTester());
    watchInversions(new JamaInversionTester());
    watchInversions(new AccurateInversionTester());

    t1 = (t1 + System.nanoTime()) / 1_000_000_000L; // seconds
    say("Done! in %2d s (%2d:%02d)", t1, t1 / 60, t1 % 60);
  }

  private static void setMatrixParams() {
    MatrixData.setMatrixParams(MATRIX_SIZE, MATRIX_DENSITY, RANGE_BOTTOM, RANGE_TOP,
        RAND_POWER, RAND_RATIO);
  }

  private static void initUnityMatrix() {
    unityMatrix = new double[MATRIX_SIZE][MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; i++)
      unityMatrix[i][i] = 1.0;
  }

  private static void watchInversions(InversionTester tester) {
    printHeader(tester.getMatrixType(), tester.getInversionMethod());
    try {
      tester.warmup();
      tester.measure();
    } catch (Exception x) {
      // say(x.toString());
      x.printStackTrace();
    }
    say();
  }

  /**
  *
  */
 private static void printHeader(MatrixType matrixType, String inversionMethod) {
   say("RandSeed =    " + RAND_SEED);
   say("Scaling:      " + USING_SCALING);
   say("Matrix type:  " + matrixType);
   say("Inversion:    " + inversionMethod);
   say("Matrix size:  " + MATRIX_SIZE);
   say("Data samples: " + MAX_TEST_ITERATIONS);
   say();
 }

  /**
   * Methods common for all types of inversions
   */
  private abstract class InversionTester {

    protected abstract String getInversionMethod();

    protected abstract MatrixType getMatrixType();

    void warmup() {
      say("Warming up");
      long startTime = System.currentTimeMillis();
      long lastDotTime = startTime;
      initRandom();

      say_(".");
      while (System.currentTimeMillis() < startTime + WARMUP_SECONDS * MILLIS_PER_SECOND) {
        measureNextInversion();         // Generates one of the possible types depending on the MATRIX_TYPE
        if (System.currentTimeMillis() > lastDotTime + MILLIS_PER_SECOND) {
          say_(".");
          lastDotTime = System.currentTimeMillis();
        }
      }
      say();
    }

    void measure() {
      say("Working");
      InversionStatistics stats = new InversionStatistics();
      initRandom();

      @SuppressWarnings("unused")
      long prevStatTime = System.nanoTime();
      for (int iteration = 1; iteration <= MAX_TEST_ITERATIONS; iteration++) {
        InversionShowing showing = measureNextInversion();
        stats.collect(iteration, showing);

        if ((System.nanoTime() > prevStatTime + THREE_SECONDS)
            // || (iteration >= MAX_TEST_ITERATIONS)
            ) {
          prevStatTime = stats.showInterimResults(); // Returns current time in nanoseconds
        }

      }
      stats.showFinalResults();
    }

    /** overriding methods in specific classes generate a matrix, inverse it,
     * measure the inversion time and errors, and create and return a <code>InversionShowing</code>
     * instance containing this information */
    InversionShowing measureNextInversion() {
      // TODO 2022-12-20 11:42:18 Auto-generated method stub
      InversionShowing showing = new InversionShowing();
      DoubleMatrix matrix = new DoubleMatrix(MatrixData.generate(getMatrixType()).getMatrix());
      Matrix inversion = inverseMatrix(matrix);
      showing.inversionTimeMS = inversionTime;

      double[][] product = matrix.multiply(inversion).getDoubleData();

      double mse = 0;;
      double maxErr = 0;
      for (int i = 0; i  < product.length; i++) {
        for (int j = 0; j < product.length; j++) {
          double error = product[i][j] - unityMatrix[i][j];
          maxErr = Math.max(maxErr, Math.abs(error));
          mse += error * error;
        }
      }
      mse = Math.sqrt(mse / (MATRIX_SIZE * MATRIX_SIZE));

      showing.maxRelativeError = maxErr;
      showing.unityMatrixMSE = mse;
      return showing;
    }

    protected abstract Matrix inverseMatrix(DoubleMatrix matrix);

    protected double inversionTime;
  }

  public class SimpleInversionTester extends InversionTester {
    MatrixType matrixType = MatrixType.RANDOM_UNIFORM_DENSE;

    @Override protected String getInversionMethod() { return "Simple inversion"; }
    @Override protected MatrixType getMatrixType() { return matrixType; }

    @Override
    protected Matrix inverseMatrix(DoubleMatrix matrix) {
      long t = -System.nanoTime();
      Matrix inverse = matrix.inverse();
      inversionTime = 1e-6 * (t + System.nanoTime());
      return inverse;
    }
  }

  public class JamaInversionTester extends InversionTester {
    MatrixType matrixType = MatrixType.RANDOM_UNIFORM_DENSE;

    @Override protected String getInversionMethod() { return "JAMA inversion"; }
    @Override protected MatrixType getMatrixType() { return matrixType; }

    @Override
    protected Matrix inverseMatrix(DoubleMatrix matrix) {
      Jama.Matrix jamaMatrix = new Jama.Matrix(matrix.getDoubleData());
      long t = -System.nanoTime();
      Jama.Matrix inverse = jamaMatrix.inverse();
      inversionTime = 1e-6 * (t + System.nanoTime());
      return new DoubleMatrix(inverse.getArray());
    }
  }

  public class AccurateInversionTester extends InversionTester {
    MatrixType matrixType = MatrixType.RANDOM_UNIFORM_DENSE;

    @Override protected String getInversionMethod() { return "Accurate inversion"; }
    @Override protected MatrixType getMatrixType() { return matrixType; }

    @Override
    protected Matrix inverseMatrix(DoubleMatrix matrix) {
      long t = -System.nanoTime();
      Matrix inverse = matrix.inverseAccurately();
      inversionTime = 1e-6 * (t + System.nanoTime());
      return inverse;
    }
  }

//******************************************************************************
//****** class InversionStatistics *********************************************
//******************************************************************************

  private class InversionStatistics {
    private double
      accumilatedInversionMse,
      accumulatedInversionTimeMS,

      maxInversionMSE,
      maxRelativeError,

      lastMaxRelativeError,
      lastUnityMatrixMSE,
      lastInversionTimeMS;

    private int iterationsDone;

    private void collect(int iteration, InversionShowing showing) {
      /**/
      iterationsDone = iteration;

      accumilatedInversionMse += showing.unityMatrixMSE;
      accumulatedInversionTimeMS += showing.inversionTimeMS;

      maxInversionMSE       = Math.max(maxInversionMSE, showing.unityMatrixMSE);
      maxRelativeError      = Math.max(maxRelativeError, showing.maxRelativeError);

      lastMaxRelativeError  = showing.maxRelativeError;  // "max err: %9.3e  "
      lastUnityMatrixMSE    = showing.unityMatrixMSE;    // "mse: %9.3e  "
      lastInversionTimeMS  = showing.inversionTimeMS;   // "time: %7.3f  "

    }

    public void showFinalResults() {
      show();
    }

    public long showInterimResults() {
      show();
      return System.nanoTime();
    }

    private double safelyDivide(double dividend, double divisor) {
      if (divisor == 0)
        return 999.9;
      return dividend / divisor;
    }

    private void show() {
      say(  "#%5d  "
          + "max err: %9.3e  "
          + "mse: %9.3e  "
          + "time: %7.3f  "
          + " \n"

          + "   avr: "
          + "max err: %9.3e  "
          + "mse: %9.3e  "
          + "time: %7.3f  ",

          iterationsDone,             // "#%5d  "
          lastMaxRelativeError,   // "max err: %9.3e  "
          lastUnityMatrixMSE,     // "mse: %9.3e  "
          lastInversionTimeMS,    // "time: %7.3f  "

          // "   avr: "
          maxRelativeError,        // "max err: %9.3e  "
          accumilatedInversionMse / iterationsDone,   // "mse: %9.3e  "
          accumulatedInversionTimeMS / iterationsDone // "time: %7.3f  ",
        );
    } // private void show(ErrorSet errors1, ErrorSet errors2)

  } // private class InversionStatistics {

//******************************************************************************
//****** class InversionShowing  ***********************************************
//******************************************************************************

  /**
   * A container for the measured characteristics shown by inversion of a matrix sample: elapsed by the inversion time,
   * result accuracy etc. Instances of this class are created, filled and returned by the method performing the inversion proper,
   * and used for collecting statistics
   * @author misa
   */
  private class InversionShowing {
    private double  maxRelativeError;
    private double inversionTimeMS;
    private double unityMatrixMSE;
    // TODO Add whatever needed
  }

  //*******************************************************************************************
  //****  Private methods *********************************************************************
  //*******************************************************************************************


  private static void showTime(long inversionTime) {
    say("Average inversion time: %9.3f ms", inversionTime * 1e-6 / MAX_TEST_ITERATIONS );
  }

  @SuppressWarnings("unused")
  private static void initRandom() {
    if (RAND_SEED < 0) {
      random = new Random();
    } else {
      random = new Random(RAND_SEED);
    }
    // say("Random: " +
    random.nextLong(); // Don't remove this in order to keep comaparable results
                        // It affects the matrices data, thea are generated with this random
    MatrixDataGenerators.setRandomSeed(random);
  }


}
