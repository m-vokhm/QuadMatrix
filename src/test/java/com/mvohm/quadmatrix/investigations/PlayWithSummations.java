package com.mvohm.quadmatrix.investigations;

import static com.mvohm.quadmatrix.test.AuxMethods.*;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Optional;
import java.util.Random;
import java.util.function.Function;
import java.util.stream.Collectors;

// 2022-01-28 18:52:35 Попробовать разные способы суммировать элементы массива
/**
 * Сделано в ходе работы над матрицами, там надо находит суммы массивов часто.
 * Мораль:
 * Collectors.summingDouble(Double::doubleValue) использует алгоритм Кахана,
 * в случае предварительно отсортированного массива ошибка пренебрежимо мала (т.е. неразличима).
 * Тупое суммирование дает ошибку порядка 2e-13 при кол-ве 1_000_000 и диапазоне -100..100
 * Мое собственное избретение, простое сложение предварительно отсортированного массивыа,
 * почему-то работает хуже простого суммирования -- я не понимаю, казалось бы,
 * при сложении сравнимых по магнитуде чисел округляемая часть мантиссы в среднем меньше,
 * должна и ошибка меньше получаться
 */

public class  PlayWithSummations {

  private static final int RAND_SEED    = 4321234; // 1234321;//
  private static final int ARRAY_SIZE   =  1_000_000; // 3_000_000;
  private static final double MIN_VALUE =         -1;
  private static final double MAX_VALUE =          1;
  private static final double POWER_TO_RAISE_RANDOMS = 2.0;

  private static final double[] SPECIFIC_ARRAY = new double[] {
    // *** s1 = -1.00944685539016, s2 = -1.0094468553901603, d = 2.220446049250313E-16
       -0.254168153278599,    -0.19559891804193213,   0.3168261888553615,    0.10183883211436486,  -0.4289729088340693,
        0.3974767255547501,    0.049459199315835194,   0.14699536186156967,  0.4822671207090789,   -0.05435180157820578,
        0.09303880520984541,   0.39692769648611315,    0.40105014590538174,  0.1578511148879448,   -0.04055236759792996,
       -0.021123431805922604, -0.6574868939999959,    -0.20158583501387659, -0.4931111438726348,   -0.13929555982939915,
       -0.07944088346914219,   0.014795384428633228,  -0.46396263692682266,  0.14895996705570394,   0.4490697476708455,
       -0.44159948981915004,   0.18582414474692144,    0.31250494214983254,  0.3961930954286837,   -0.20894641840628247,
        0.1102408616914687,   -0.010683296176456698,  -0.03881496083311635,  0.30048167141501425,  -0.026938360258509157,
       -0.27872133577963515,  -0.00120975763265634,   -0.22073827219287862, -0.3422103979284941,   -0.3851189154535705,
       -0.17963324603324277,  -0.06357202364051234,   -0.5081196599392694,  -0.18536780996822455,   0.21884354567330144,
       -0.0200857072616412,   -0.02849985664909696,    0.2302644217841188,   0.04641984049740752,   0.003134373388930131,
   };


  private static Random random = new Random();
  static double correctSum;
  static double[] vector;
  static Optional<Double> maxValue;
  private static MathContext MC_100 = new MathContext(100, RoundingMode.HALF_EVEN);

  public static void main(String[] args) {
    setRandomSeed(RAND_SEED);
    say("Really working? so fast...");

    vector = randomArray(ARRAY_SIZE, MIN_VALUE, MAX_VALUE);
//    vector = SPECIFIC_ARRAY;

    maxValue  = Arrays.stream(vector).boxed()
        .collect(Collectors.maxBy((a, b) -> Double.compare(Math.abs(a), Math.abs(b))));
    say("Max: " + maxValue.get());

    correctSum = sumAsBigDecimals(vector);
    say("sumAsBigDecimals:    " + correctSum + "\n");

    for (int i = 0; i < 3; i++) {
      testSummation("sumNaive",         PlayWithSummations::sumNaive);
      testSummation("sumNaiveSorted",   PlayWithSummations::sumNaiveSorted);

      // Added 2022-11-27 20:48:09
      testSummation("sumNaiveSorted2",  PlayWithSummations::sumNaiveSorted2);

      testSummation("sumStream",        PlayWithSummations::sumStream);
      testSummation("sumStreamSorted",  PlayWithSummations::sumStreamSorted);
      testSummation("kahanSum",         PlayWithSummations::kahanSum);
      testSummation("kahanSumSorted",   PlayWithSummations::kahanSumSorted);

      // Added 2022-11-27 20:48:09
      testSummation("kahanSumManually", PlayWithSummations::kahanSumManuallySorted);
      testSummation("kahanSumOptimized", PlayWithSummations::kahanSumManuallyOptimized);
      testSummation("kahanSumReverse", PlayWithSummations::kahanSumManuallyReverse);
      say("\n------------------------------\n");
    }
  }

  /** Set random seed to provide reproducibility */
  public static void setRandomSeed(int seed) {
    random = new Random(seed);
  }

  /**
   * A dense array filled with random values ranged from 0 to 1.0
   * @param length
   * @return
   */
  public static double[] randomArray(int length) {
    final double[] result = new double[length];
    for (int i = 0; i < length; i++)
      result[i] = random.nextDouble();
    return result;
  }

  private static void testSummation(String name, Function<double[], Double> tested) {
    long t1 = -System.nanoTime();
    final double sum = tested.apply(vector);
    t1 += System.nanoTime();
    say("%-20s %s", name + ":", sum);
//    say("               " + hexStr(sum));
    say("%-20s %.3e\n%-15s   %9.3f ns\n", "Error:", (sum - correctSum)/maxValue.get(), "time", (double)t1 / vector.length);
  }

  /**
   * A dense array filled with random values ranged from 0 to 1.0
   * @param length
   * @return
   */
  public static double[] randomArray(int length, double rangedFrom, double rangedTo) {

    final double power = POWER_TO_RAISE_RANDOMS;

    // Find roots to limit range when generating randoms
    rangedFrom = rangedFrom >= 0 ? Math.pow(rangedFrom, 1.0/power) : -Math.pow(-rangedFrom, 1.0/power);
    rangedTo = rangedTo >= 0? Math.pow(rangedTo, 1.0/power) : -Math.pow(-rangedTo, 1.0/power);

    final double[] result = new double[length];
    for (int i = 0; i < length; i++) {
      final double d = random.nextDouble() * (rangedTo - rangedFrom) + rangedFrom;
      // Each random raise to the given power
      result[i] = d >= 0? Math.pow(d, power) : -Math.pow(-d, power);
    }
    return result;
  }

  public static double sumAsBigDecimals(double[] vector) {
    BigDecimal sum = BigDecimal.ZERO;
    final int size = vector.length;
    for (int i = 0; i < size; i++) {
//      sum = sum.add(BigDecimal.valueOf(vector[i]), MC_100);
      sum = sum.add(new BigDecimal(vector[i]), MC_100);
    }
    return sum.doubleValue();
  }

  //**********************************************************************
  //********** Tested methods ********************************************
  //**********************************************************************

  public static double sumNaive(double[] vector) {
    double sum = 0;
    for (final double d: vector) {
      sum += d;
    }
    return sum;
  }

  public static double sumNaiveSorted(double[] vector) {
    final Double[] negatives = Arrays.stream(vector).filter(d -> d < 0).boxed()
                        .sorted(Comparator.reverseOrder())
                        .toArray(Double[]::new);

    final Double[] positives = Arrays.stream(vector).filter(d -> d > 0).boxed()
                        .sorted().toArray(Double[]::new);


    double positiveSum = 0;
    double negativeSum = 0;
    for (final double d: positives) {
      positiveSum += d;
    }
    for (final double d: negatives) { // 2022-11-28 11:25:59 Жопа!!!! Надо в обратном порядке
      negativeSum += d;
    }
    return positiveSum + negativeSum;
  }

  public static double sumNaiveSorted2(double[] vector) {
    vector = vector.clone();
    Arrays.sort(vector);

    // Найдем минимальное неотрицательное число
    int lo = 0, hi = vector.length - 1;
    while (hi >= lo) {
      final int m = (lo + hi) / 2;
      if (vector[m] < 0) lo = m + 1;
      else hi = m - 1;
    }

    double positiveSum = 0;
    double negativeSum = 0;
    for (int i = lo; i < vector.length; i++) {
      positiveSum += vector[i];
    }
    for (int i = lo - 1; i >=0; i--) {
      negativeSum += vector[i];
    }

    return positiveSum + negativeSum;
  }

  public static double sumStream(double[] vector) {
    return Arrays.stream(vector).boxed()
                      .collect(Collectors.summingDouble(Double::doubleValue));
  }

  public static double sumStreamSorted(double[] vector) {
    return Arrays.stream(vector).boxed()
                      .sorted((a, b) -> Double.compare(Math.abs(a), Math.abs(b)))
                      .collect(Collectors.summingDouble(Double::doubleValue));
  }

  public static double kahanSum(double[] vector) {
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

  public static double kahanSumSorted(double[] vector) {
    final Double[] sorted = Arrays.stream(vector).boxed()
        .sorted((a, b) -> Double.compare(Math.abs(a), Math.abs(b)))
        .toArray(Double[]::new);

    double sum = 0.0;
    double c = 0.0;
    for (int i = 0; i < sorted.length; i++) {
      final double y = sorted[i] - c;
      final double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    return sum;
  }

  //*******************************************************
  // Added 2022-11-27 20:48:38 -- other additions are extremely ineffective
  public static double kahanSumManuallySorted(double[] vector) {
    vector = vector.clone();
    sortByAbs(vector);

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
   * @param vector
   */
  private static void sortByAbs(double[] vector) {
    quickSort(vector, 0, vector.length - 1);
  }

  ///*************************************
  // Taken from F:\.eclipse\Lopt\TestSorts\src\stackoverflow\QSortBenchmarks
  // 2022-11-27 20:06:56
  private static void quickSort(double[] array, int lo, int hi) {
    if (hi <= lo) return;
    final int j = partition(array, lo, hi);
    quickSort(array, lo, j - 1);
    quickSort(array, j + 1, hi);
  }

  private static int partition(double[] array, int lo, int hi) {
    int i = lo, j = hi + 1;
    while (true) {
      while (compare(array[++i], array[lo]) < 0) // while (array[++i] < array[lo])
        if (i == hi) break;
      while (compare(array[lo], array[--j]) < 0) // while (array[lo] < array[--j])
        if (j == lo) break;
      if (i >= j) break;
      swapItems(array, i, j);
    }
    swapItems(array, lo, j);
    return j;
  }

  private static int compare(double v1, double v2) {
    return Double.compare(Math.abs(v1), Math.abs(v2));
  }

  private static void swapItems(double[] array, int i, int j) {
    final double tmp = array[i];
    array[i] = array[j];
    array[j] = tmp;
  }
  //
  //****************************************************

  //*******************************************************
  // Added 2022-11-28 18:55:43 -- some more variations of kahan summation after manual sorting
  public static double kahanSumManuallyOptimized(double[] vector) {
    vector = vector.clone();
    sortByAbsOptimized(vector);

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
   * A QuickSort variation comparing items by absolute values, slightly improved
   * @param vector
   */
  private static void sortByAbsOptimized(double[] vector) {
    quickSortOpt(vector, 0, vector.length - 1);
  }

  ///*************************************
  // Taken from F:\.eclipse\Lopt\TestSorts\src\stackoverflow\QSortBenchmarks
  // 2022-11-27 20:06:56
  private static void quickSortOpt(double[] array, int lo, int hi) {
    if (hi <= lo) return;
    final int j = partitionOpt(array, lo, hi);
    quickSortOpt(array, lo, j - 1);
    quickSortOpt(array, j + 1, hi);
  }

  private static int partitionOpt(double[] array, int lo, int hi) {
    int i = lo, j = hi + 1;
    while (true) {
      double loValueAbs = Math.abs(array[lo]);
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

  public static double kahanSumManuallyReverse(double[] vector) {
    vector = vector.clone();
    sortByAbsReverse(vector);

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
   * A QuickSort variation comparing items by absolute values, slightly improved
   * @param vector
   */
  private static void sortByAbsReverse(double[] vector) {
    quickSortReverse(vector, 0, vector.length - 1);
  }

  ///*************************************
  // Taken from F:\.eclipse\Lopt\TestSorts\src\stackoverflow\QSortBenchmarks
  // 2022-11-27 20:06:56
  private static void quickSortReverse(double[] array, int lo, int hi) {
    if (hi <= lo) return;
    final int j = partitionReverse(array, lo, hi);
    quickSortOpt(array, lo, j - 1);
    quickSortOpt(array, j + 1, hi);
  }

  private static int partitionReverse(double[] array, int lo, int hi) {
    int i = lo, j = hi + 1;
    while (true) {
      double loValueAbs = Math.abs(array[lo]);
      while (Math.abs(array[++i]) >= loValueAbs)
        if (i == hi) break;
      while (loValueAbs >= Math.abs(array[--j]))
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
}
