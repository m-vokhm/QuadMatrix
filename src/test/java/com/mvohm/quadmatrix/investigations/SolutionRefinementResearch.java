package com.mvohm.quadmatrix.investigations;

import java.util.Random;
import java.util.stream.DoubleStream;

import com.mvohm.quadmatrix.DoubleMatrix;
import com.mvohm.quadmatrix.investigations.MatrixData.MatrixType;
import com.mvohm.quadmatrix.test.AuxMethods;
import com.mvohm.quadmatrix.test.MatrixDataGenerators;

import static com.mvohm.quadmatrix.test.AuxMethods.*;
import static com.mvohm.quadmatrix.test.MatrixDataGenerators.*;

import static com.mvohm.quadmatrix.investigations.MatrixData.MatrixType;


/**
  Started 2022-02-09 13:05:13
  Accumulate some statistics on refinements of linear system solutions
*/
public class SolutionRefinementResearch {

  private static final boolean USING_SCALING = true;

  private static Random random = null;
  private static final int RAND_SEED  =
                                         // 123456789;
                                         //    543210;
                                             123321; // Negative value to ignore constant RAND_SEED
                                         //       123;
                                         //       456;

  /** Choose the type of matrices to use this time */
  private static final MatrixType MATRIX_TYPE =
                                                    MatrixType.RANDOM_UNIFORM_DENSE;
                                                    // MatrixType.SPD_DENSE;

  /** The solution method to use this time */
  private static final SolutionMethod SOLVE_WITH =
                                                     SolutionMethod.LU_DECOMPOSITION;
                                                    // SolutionMethod.CHOLESKY_DECOMPOSITION; // Chose SPD matrices for this solution method, otherwise it'll throw an exception
                                                    // SolutionMethod.JAMA_SOLUTION;

  private static final int MAX_TEST_ITERATIONS  = 1000; // 10_000;
  private static final int MATRIX_SIZE      = 500;
  private static final double RANGE_BOTTOM  = -1.0;
  private static final double RANGE_TOP     = 1.0;
  private static final double RANGE         = RANGE_TOP - RANGE_BOTTOM;

  // For non-uniform distribution, we fill the matrix with values computed as
  //    x = ratio * (sign(r) * abs(r)^power) + (1 - ratio) * r, where 'r' is a double random, 'power' and 'ratio' are defined below:
  private static final double RAND_POWER    = 8.0;    // The power to raise a random to obtain a non-linear distribution for generating non-uniform matrices
  private static final double RAND_RATIO    = 0.95;   // A ratio to mix the power of the random and the random value itself,
  static final double MATRIX_DENSITY        = 0.2;    // density for generating sparse matrices

  /** Available solution methods */
  enum SolutionMethod {
    LU_DECOMPOSITION,
    CHOLESKY_DECOMPOSITION,
    JAMA_SOLUTION
  }

  private static final int ERRORS_TO_SHOW       = 10;

  private static final long WARMUP_SECONDS      = 3;      // In seconds
  private static final long MILLIS_PER_SECOND   = 1_000;
  private static final long NANOS_PER_SECOND    = 1_000_000_000L;
  private static final long THREE_SECONDS       = 3 * NANOS_PER_SECOND;   // in nanoseconds

  private static long // solveAndRefine() sets the time elapsed by a singular operation here
    solutionTime,     // for solution
    refinementTime;   // for refinement

  public static void main(String[] args) {
    new SolutionRefinementResearch().run();
  }

  private void run() {
    long t1 = -System.nanoTime();
    printHeader();
    DoubleMatrix.setDefaultScaling(USING_SCALING);
    collectStats();
    t1 = (t1 + System.nanoTime()) / 1_000_000_000L; // seconds
    say("Done! in %2d s (%2d:%02d)", t1, t1 / 60, t1 % 60);
  }

  /**
   *
   */
  private static void printHeader() {
    say("RandSeed =    " + RAND_SEED);
    say("Scaling:      " + USING_SCALING);
    say("Matrix type:  " + MATRIX_TYPE.name());
    say("Sol. method:  " + SOLVE_WITH.name());
    say("Matrix size:  " + MATRIX_SIZE);
    say("Data samples: " + MAX_TEST_ITERATIONS);
    say();
  }

  private void collectStats() {
    setMatrixParams();                  // Set the parameters affecting the generated matrices (type, size etc)
    warmUp();
    doWork();
    say();
  }

  private static void setMatrixParams() {
    MatrixData.setMatrixParams(MATRIX_SIZE, MATRIX_DENSITY, RANGE_BOTTOM, RANGE_TOP,
                                RAND_POWER, RAND_RATIO);
  }

  /**
   */
  private void doWork() {
    say("Working");
    final Statistics stats = new Statistics();
    initRandom();
    solutionTime = refinementTime = 0;

    @SuppressWarnings("unused")
    long prevStatTime = System.nanoTime();

    for (int iteration = 1; iteration <= MAX_TEST_ITERATIONS; iteration++) {
      final MatrixData data = MatrixData.generate(MATRIX_TYPE); // Generates one of the possible types depending on the MATRIX_TYPE

      final ErrorSet errors = solveAndRefine(data); // Adds time to static fields solutionTime and refinementTime
      stats.collect(iteration, errors);

      if ((System.nanoTime() > prevStatTime + THREE_SECONDS) || (iteration >= MAX_TEST_ITERATIONS)) {
        prevStatTime = showInterimResults(stats, data, errors); // Returns current time in nanoseconds
      }
    }

    stats.showMinRatios();
    showTimes(solutionTime, refinementTime);
  }

  /**
   * @param stats
   * @param data
   * @param errors
   * @return
   */
  private static long showInterimResults(Statistics stats, MatrixData data, ErrorSet errors) {
    stats.show(errors);
    return System.nanoTime();
  }

  /**
   */
  private static void warmUp() {
    say("Warming up");
    final long startTime = System.currentTimeMillis();
    long lastDotTime = startTime;
    initRandom();

    say_(".");
    while (System.currentTimeMillis() < startTime + WARMUP_SECONDS * MILLIS_PER_SECOND) {
      solveAndRefine(MatrixData.generate(MATRIX_TYPE));         // Generates one of the possible types depending on the MATRIX_TYPE
      if (System.currentTimeMillis() > lastDotTime + MILLIS_PER_SECOND) {
        say_(".");
        lastDotTime = System.currentTimeMillis();
      }
    }
    say();
  }

  private static ErrorSet solveAndRefine(MatrixData data) {
    final DoubleMatrix matrix1 = new DoubleMatrix(data.getMatrix());
    final DoubleMatrix matrix2 = new DoubleMatrix(data.getMatrix());
    final boolean SOLVE_SIMPLY = false, SOLVE_ACCURATELY = true;

    solutionTime += solveMatrix(matrix1, data.getVector(), SOLVE_SIMPLY);
    final double[] actualSolution =
        SOLVE_WITH == SolutionMethod.JAMA_SOLUTION?
            jamaSolution:
            matrix1.getDoubleSolution();

    refinementTime += solveMatrix(matrix2, data.getVector(), SOLVE_ACCURATELY);
    final double[] refinedSolution =
        SOLVE_WITH == SolutionMethod.JAMA_SOLUTION?
            jamaSolution:
            matrix2.getDoubleSolution();

    return ErrorSet.findErrors(data, actualSolution, refinedSolution);
  }

  private static double[] jamaSolution;  // To keep solution found by a Jama matrix, it does not contain the solution

  private static long solveMatrix(DoubleMatrix matrix, double[] vector, boolean accurately) {
    long t1 = -System.nanoTime();

    switch (SOLVE_WITH) {
      case CHOLESKY_DECOMPOSITION: {
        if (accurately)
          matrix.solveSPDAccurately(vector);
        else
          matrix.solveSPD(vector);
        t1 += System.nanoTime();
        return t1;
      }
      case LU_DECOMPOSITION: {
        if (accurately)
          matrix.solveAccurately(vector);
        else
          matrix.solve(vector);
        t1 += System.nanoTime();
        return t1;
      }
      case JAMA_SOLUTION: {
        return solveWithJama(matrix, vector); // It returns actual solution time
      }
      default: { throw new IllegalArgumentException("Unknown solution method"); }
    }
  }

  private static long solveWithJama(DoubleMatrix matrix, double[] vector) {
    final Jama.Matrix jamaMatrix = new Jama.Matrix(matrix.getDoubleData());
    Jama.Matrix matrixB = new Jama.Matrix(vector, 1);
    matrixB = matrixB.transpose();
    long t = -System.nanoTime();
    final Jama.Matrix solution = jamaMatrix.solve(matrixB);
    t += System.nanoTime();
    jamaSolution = solution.getColumnPackedCopy();
    return t;
  }


  private class Statistics {
    double
      accumilatedInitialMse, accumulatedRefinedMse,
      mseImprovementRatio, accumulatedMseImprovemetRatio;

    double vectorImprovementRatio;

    int iterationsDone;

    double
      minMseRatio, minMeanErrorRatio, minMaxErrorRatio;

    {
      minMseRatio = minMeanErrorRatio = minMaxErrorRatio = Double.MAX_VALUE;
    }

    private void collect(int iteration, ErrorSet errors1) {
      /**/
      accumilatedInitialMse += errors1.initialMse;
      accumulatedRefinedMse += errors1.refinedMse;

      mseImprovementRatio = safelyDivide(errors1.initialMse, errors1.refinedMse);
      accumulatedMseImprovemetRatio += mseImprovementRatio;

      vectorImprovementRatio = safelyDivide(errors1.initialVectorMse, errors1.refinedVectorMse);    // How much the reconstructed vector improves after 1 iteration

      iterationsDone = iteration;                        // It's the loop variable of a for loop, 0 for 1st iteration

      minMseRatio       = Math.min(minMseRatio, mseImprovementRatio);
      minMeanErrorRatio = Math.min(minMeanErrorRatio, safelyDivide(errors1.initialMean, errors1.refinedMean));
      minMaxErrorRatio  = Math.min(minMaxErrorRatio, safelyDivide(errors1.initialMax, errors1.refinedMax));
    }

    private double safelyDivide(double dividend, double divisor) {
      if (divisor == 0)
        return 999.9;
      return dividend / divisor;
    }

    private void show(ErrorSet errors) {
      say(  "#%5d  "
          + "e0: %9.3e  "
          + "e1: %9.3e  "
          + "r1: %7.3f    "
          + "v0: %9.3e  "
          + "v1: %9.3e  "
          + "x1: %7.3f  "
          + " \n"

          + "   avr: "
          + "e0: %9.3e  "
          + "e1: %9.3e  "
          + "r1: %7.3f",

          iterationsDone,              // #%5d
          errors.initialMse,           // e0: %9.3e
          errors.refinedMse,           // e1: %9.3e
          mseImprovementRatio,       // r1: %7.3f
          errors.initialVectorMse,     // v0: %9.3e
          errors.refinedVectorMse,     // v1: %9.3e
          vectorImprovementRatio,    // x1: %7.3f  // ratio of initial vector mse to refined vector mse

          // avr:
          accumilatedInitialMse / iterationsDone,             // e0: %9.3e
          accumulatedRefinedMse / iterationsDone,           // e1: %9.3e
          accumulatedMseImprovemetRatio / iterationsDone   // r1: %7.3f
        );
    } // private void show(ErrorSet errors1, ErrorSet errors2)

    private void showMinRatios() {
      say("Min ratios: mse1 = %7.3f, mean1 = %7.3f, max1 = %7.3f",
            minMseRatio, minMeanErrorRatio, minMaxErrorRatio
          );
    }
  }

  static class ErrorSet {
    double
      initialMax,
      initialMean,
      initialMse,
      initialVectorMse,

      refinedMax,
      refinedMean,
      refinedMse,
      refinedVectorMse;


    private ErrorSet (double initialMax, double initialMean, double initialMse, double initialVectorMse) {
      this.initialMax       = initialMax;
      this.initialMean      = initialMean;
      this.initialMse       = initialMse;
      this.initialVectorMse = initialVectorMse;
    }

    static ErrorSet findErrors(MatrixData data, double[] simpleSolution, double[] refinedSolution) {
      final ErrorSet errors = findInitialErrors(data, simpleSolution);
      return errors.setRefinedErrors(data, refinedSolution);
    }

    private static ErrorSet findInitialErrors(MatrixData data, double[] actualSolution) {
      final double[] diff = findDifferences(data, actualSolution);

      final double max  = DoubleStream.of(diff).map(Math::abs).max().getAsDouble();
      final double mean = Math.abs(DoubleStream.of(diff).average().getAsDouble());
      final double mse  = Math.sqrt(DoubleStream.of(diff).map(v -> v * v).average().getAsDouble());

      final double[] computedVector = multVectorByMatrix(actualSolution, data.getMatrix());
      final double[] vectorDiff = AuxMethods.subtractVectors(data.getVector(), computedVector);
      final double vectorMse = Math.sqrt(DoubleStream.of(vectorDiff).map(v -> (v * v) / (RANGE * RANGE)).average().getAsDouble());

      return new ErrorSet(max, mean, mse, vectorMse);
    }

    private ErrorSet setRefinedErrors(MatrixData data, double[] actualSolution) {
      final double[] diff = findDifferences(data, actualSolution);

      refinedMax  = DoubleStream.of(diff).map(Math::abs).max().getAsDouble();
      refinedMean = Math.abs(DoubleStream.of(diff).average().getAsDouble());
      refinedMse  = Math.sqrt(DoubleStream.of(diff).map(v -> v * v).average().getAsDouble());

      final double[] computedVector = multVectorByMatrix(actualSolution, data.getMatrix());
      final double[] vectorDiff = AuxMethods.subtractVectors(data.getVector(), computedVector);
      refinedVectorMse = Math.sqrt(DoubleStream.of(vectorDiff).map(v -> (v * v) / (RANGE * RANGE)).average().getAsDouble());
      return this;
    }

    /**
     * @param data
     * @param actualSolution
     * @return
     */
    private static double[] findDifferences(MatrixData data, double[] actualSolution) {
      double[] diff = AuxMethods.subtractVectors(data.getSolution(), actualSolution);
      diff = DoubleStream.of(diff).map(v -> Math.abs(v/RANGE)).toArray();
      return diff;
    }

    private static double[] multVectorByMatrix(double[] vector, double[][] matrix) {
      final int size = vector.length;
      final double[] result = new double[size];
      for (int i = 0; i < size; i++) {
        final double[] rowProduct = multiplyElements(matrix[i], vector);
        result[i] = sumOfVector(rowProduct);
      }
      return result;
    }

  }

  @SuppressWarnings("unused")
  private static void showDigits(double[] solution, double[] refinedSolution) {
    for (int i = 0; i < ERRORS_TO_SHOW; i++) {
      say_("%21.14e ", solution[i]);
    }
    say();
    for (int i = 0; i < ERRORS_TO_SHOW; i++) {
      say_("%21.14e ", refinedSolution[i]);
    }
    say();
    for (int i = 0; i < ERRORS_TO_SHOW; i++) {
      say_("%21.14e ", solution[i] - refinedSolution[i]);
    }
    say();
  }

  /**
   * Shows average times of the solution, one-pass refinement, and iterative refinement
   * @param solutionTime
   * @param refinement1Time
   * @param refinement2Time
   */
  private static void showTimes(long solutionTime, long refinement1Time) {
    say("Average solution time: %9.3f ms, accurate sol. time: %9.3f ms",
          solutionTime * 1e-6 / MAX_TEST_ITERATIONS,
          refinement1Time * 1e-6 / MAX_TEST_ITERATIONS
        );
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
