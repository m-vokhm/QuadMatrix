package com.mvohm.quadmatrix.test;

import static com.mvohm.quadmatrix.test.AuxMethods.*;
import static org.assertj.core.api.Assertions.assertThat;

import com.mvohm.quadmatrix.test.AuxMethods.ErrorSet;

public class OperationComparator {

  public interface DataGenerator {
    MatrixData generate();
  }

  public interface OperationPerformer {
    ErrorSet perform(MatrixData matrixData);
  }

  private static final int WARMUP_COUNT = 10;

  /** a substitute of infinity, used when we estimate the ratio of e0/e1 and e1 == 0.
   * Used when estimating average refinement effectiveness on a large set of data samples.
   * Will signify a very good ratio for a certain data sample, where the error after refinement for this sample is 0.  */
  private static final int      SUBSTITUTE_INFINITY = 1000;


  private final DataGenerator dataGenerator;
  private final OperationPerformer performerOne, performerTwo;

  private int trialCount, timedTrialCount;

  private ErrorSet errorSetOne, errorSetTwo;
  private long timeOne, timeTwo;
  private double accumulatedMseOne, accumulatedMseTwo;
  private double accumulatedMaxErrOne, accumulatedMaxErrTwo;

  private double accumulatedMseImprovemetRatio;
  private double accumulatedMaxErrImprovemetRatio;
  private double minMseImprovemetRatio = Double.MAX_VALUE;
  private double minMaxErrImprovemetRatio = Double.MAX_VALUE;

  private MatrixData matrixData;

  private double mseImprovementRatio, maxErrImprovementRatio;


  public OperationComparator(DataGenerator generator, OperationPerformer performerOne, OperationPerformer performerTwo) {
    this.dataGenerator = generator;
    this.performerOne = performerOne;
    this.performerTwo = performerTwo;
  }

  /** For a large randomly-generated matrix, compare errors of simple matrix solution
   * (max relative error and MSE) with those of accurate LU solution,
   * and return error ratios */
  public void performOperations() {
    trialCount++;
    matrixData = dataGenerator.generate();

    errorSetOne = performerOne.perform(matrixData);
    accumulatedMseOne += errorSetOne.mse();
    accumulatedMaxErrOne = Math.max(accumulatedMaxErrOne, errorSetOne.maxError());

    errorSetTwo = performerTwo.perform(matrixData);
    accumulatedMseTwo += errorSetTwo.mse();
    accumulatedMaxErrTwo = Math.max(accumulatedMaxErrTwo, errorSetTwo.maxError());

    mseImprovementRatio = safelyDivide(errorSetOne.mse(), errorSetTwo.mse(), SUBSTITUTE_INFINITY);
    maxErrImprovementRatio = safelyDivide(errorSetOne.maxError(), errorSetTwo.maxError(), SUBSTITUTE_INFINITY);

    accumulatedMseImprovemetRatio += mseImprovementRatio;
    accumulatedMaxErrImprovemetRatio += maxErrImprovementRatio;
    minMseImprovemetRatio = Math.min(minMseImprovemetRatio, mseImprovementRatio);
    minMaxErrImprovemetRatio = Math.min(minMaxErrImprovemetRatio, maxErrImprovementRatio);

    if (trialCount > WARMUP_COUNT) {
      timeOne += errorSetOne.getTime();
      timeTwo += errorSetTwo.getTime();
      timedTrialCount++;
    }
  }

  double getAverageTimeOneMs() {
    return 1e-6 * timeOne / timedTrialCount;
  }

  double getAverageTimeTwoMs() {
    return 1e-6 * timeTwo / timedTrialCount;
  }

  double getAverageMseOne() {
    return accumulatedMseOne / trialCount;
  }

  double getAverageMseTwo() {
    return accumulatedMseTwo / trialCount;
  }

  double getAccumulatedMaxErrOne() {
    return accumulatedMaxErrOne;
  }

  double getAccumulatedMaxErrTwo() {
    return accumulatedMaxErrTwo;
  }

  double getAverageMseImprovementRatio() {
    return accumulatedMseImprovemetRatio / trialCount;
  }

  double getAverageMaxErrImprovementRatio() {
    return accumulatedMaxErrImprovemetRatio / trialCount;
  }

  public double getCurrentMseRatio() {
    return errorSetOne.mse() / errorSetTwo.mse();
  }

  public double getCurrentMaxErrRatio() {
    return errorSetOne.maxError() / errorSetTwo.maxError();
  }

  double getMinMaxErrImprovementRatio() {
    return minMaxErrImprovemetRatio;
  }

  double getMinMseImprovementRatio() {
    return minMseImprovemetRatio;
  }

  public void showDetailedReport() {
    say("=== Comparing %s with %s", matrixData.testedMethod1name(), matrixData.testedMethod2name());

    say("avr mse1    = %11.3e, avr mse2    = %11.3e, ratio = %7.3f",
        getAverageMseOne(), getAverageMseTwo(), getAverageMseOne() / getAverageMseTwo());

    say("avr maxErr1 = %11.3e, avr maxErr2 = %11.3e, ratio = %7.3f",
        getAccumulatedMaxErrOne(), getAccumulatedMaxErrTwo(), getAccumulatedMaxErrOne() / getAccumulatedMaxErrTwo());
    say("Matrix size = %7s", matrixData.getSize());
    say("avr time1   = %7.3f  ms, avr time2   = %7.3f  ms, ratio = %7.3f",
        getAverageTimeOneMs(), getAverageTimeTwoMs(), getAverageTimeTwoMs() / getAverageTimeOneMs() );

    say("===================");
  }

  public void assertImprovementAsExpected(double expectedMinImprovement, double expectedAvrImprovement) {
    final double
      actualMinMSEImprovement     = getMinMseImprovementRatio(),
      actualAvrMSEImprovement     = getAverageMseImprovementRatio(),
      actualMinMaxErrImprovement  = getMinMaxErrImprovementRatio(),
      actualAvrMaxErrImprovement  = getAverageMaxErrImprovementRatio();

    say("      Min MSE impr. ratio: %7.3f,        Avr MSE impr. ratio: %7.3f",
        actualMinMSEImprovement, actualAvrMSEImprovement);
    say("Min Max Error impr. ratio: %7.3f,  Avr Max Error impr. ratio: %7.3f",
        actualMinMaxErrImprovement, actualAvrMaxErrImprovement);
    say("Expected min. impr. ratio: %7.3f, expected avr.  impr. ratio: %7.3f\n",
        expectedMinImprovement, expectedAvrImprovement);

    assertImprovementOK("Minimal MSE", actualMinMSEImprovement, expectedMinImprovement);
    assertImprovementOK("Average MSE", actualAvrMSEImprovement, expectedAvrImprovement);
    assertImprovementOK("Minimal max.err.", actualMinMaxErrImprovement, expectedMinImprovement);
    assertImprovementOK("Average max.err.", actualAvrMaxErrImprovement, expectedAvrImprovement);
  }

  private static void assertImprovementOK(String name, double actual, double expected) {
    assertThat(actual)
        .withFailMessage("%s improvement ratio is %.3f that's less than the expected %.3f", name, actual, expected)
        .isGreaterThan(expected);
  }


}



