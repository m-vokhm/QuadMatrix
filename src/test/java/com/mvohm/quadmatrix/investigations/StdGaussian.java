package com.mvohm.quadmatrix.investigations;

import static com.mvohm.quadmatrix.test.AuxMethods.*;

// 2022-12-05 19:30:17
// Поиграемся с разными неравномерно распределнными случайными числами

import java.util.Arrays;


/*************************************************************************
 *  Author:       Kevin Wayne
 *  Date:         8/20/04
 *  Compilation:  javac StdGaussian.java
 *  Execution:    java StdGaussian
 *
 *  Computes a standard Gaussian random deviate using the cartesian
 *  form of the Box-Muller transform. The method is to compute a
 *  random point (x, y) inside the circle centered at (0, 0) with
 *  radius 1. Then
 *
 *     x * Math.sqrt(-2.0 * Math.log(r) / r)
 *
 *  and
 *
 *     y * Math.sqrt(-2.0 * Math.log(r) / r)
 *
 *  are each Gaussian random variables with mean 0 and standard deviation 1.
 *  This formula appears in
 *
 *     Knuth, The Art of Computer Programming, Volume 2, p. 122.
 *
 *
 *  Sample executions
 *  ---------------------
 *  % java StdGaussian
 *  -1.2927277357189828
 *
 *  % java StdGaussian
 *  0.32433796089430544
 *
 *  % java StdGaussian
 *  -0.1174251833833895
 *
 *  % java StdGaussian
 *  0.053192581194524566
 *
 *************************************************************************/

public class StdGaussian {

  /*
   * Classic Kevin Waine code:
      double r, x, y;

      // find a uniform random point (x, y) inside unit circle
      do {
         x = 2.0 * Math.random() - 1.0;
         y = 2.0 * Math.random() - 1.0;
         r = x*x + y*y;
      } while (r > 1 || r == 0);    // loop executed 4 / pi = 1.273.. times on average
                                    // http://en.wikipedia.org/wiki/Box-Muller_transform

      // apply the Box-Muller formula to get standard Gaussian z
      double z = x * Math.sqrt(-2.0 * Math.log(r) / r);

      // print it to standard output
      System.out.println(z);
   */

  enum RandType {
      GAUSS,
      POW_1_5,
      POW_2,
      POW_4,
      POW_4_LIN_1,
      TAN,
      };



  RandType usedType = RandType.POW_4_LIN_1;
  // RandType usedType = RandType.TAN;

  int[] counts = new int[51];

  public static void main(String[] args) {
    new StdGaussian().run();
  }

   private void run() {
     fillArray();
     printArray();
  }

  private void printArray() {
    double startValue = -1.02;
    for (int i = 0; i < counts.length; i++) {
      // Выведем один символ на каждые 100 счетчика
      final char[] str = new char[(int)Math.round(counts[i] / 10000.0)];
      Arrays.fill(str, '#');
      final String s = new String(str);
      say("%5.2f - %5.2f: %,8d %s", startValue, startValue += 0.04, counts[i], s);
    }

  }

  private void fillArray() {
    double max = 1;
    for (int i = 0; i < 5_000_000; i++) {
      double x = 0;
      do {
        x = thisRandom();
//        if (x < -1.0 || 1.0 < x)
//          say("%s", x);
        max = Math.max(Math.abs(x), max);
      }
      while (x < -1.0 || 1.0 < x);

      if (x < -10) counts[0]++;
      else if (x > 10) counts[50]++;
      else {
        final int idx = (int)Math.round(x * 25) + 25;
        counts[idx]++;
      }
    }
    say("Max: " + max);
  }


  private static double RAND_POWER = 12.0;
  private static double RAND_RATIO = 0.98;

  private double thisRandom() {
    switch (usedType) {
      case GAUSS:   return gaussian();
      case POW_1_5: return randPow1_5();
      case POW_2:   return randPow2();
      case POW_4:   return randPow4();
      case POW_4_LIN_1: return randPow4PlusLinear(RAND_POWER, RAND_RATIO);
      case TAN: return randTan();
      default:      return 0;
    }

//     return (gaussian() /* *+ randPow4() */ ) ;
  }

  private double randPow1_5() {
    final double x = 2.0 * Math.random() - 1.0;
    return x < 0? -Math.pow(-x, 1.5) : Math.pow(x, 1.5);
  }

  private double randPow2() {
    final double x = 2.0 * Math.random() - 1.0;
    return 0.5 * x + 0.5 * (x < 0? -Math.pow(-x, 2) : Math.pow(x, 2));
  }

  private double randPow4() {
    final double x = 2.0 * Math.random() - 1.0;
    return x < 0? -Math.pow(-x, 4) : Math.pow(x, 4);
  }

  private double randPow4PlusLinear(double power, double ratio) {
    final double x = 2.0 * Math.random() - 1.0;
    final double y = x < 0? -Math.pow(-x, power) : Math.pow(x, power);
    return ratio * y + (1 - ratio) * x;
  }

  private double randTan() {
    final double x = 2.0 * Math.random() - 1.0;
    return Math.tan(x);
  }

  public static double gaussian() {
     double r, x, y;

     do {                               // find a uniform random point (x, y) inside unit circle
        x = 2.0 * Math.random() - 1.0;
        y = 2.0 * Math.random() - 1.0;
        r = x*x + y*y;
     } while (r > 1 || r == 0);         // loop executed 4 / pi = 1.273.. times on average
                                        // http://en.wikipedia.org/wiki/Box-Muller_transform

     // apply the Box-Muller formula to get standard Gaussian z
     return x * Math.sqrt(-2.0 * Math.log(r) / r) / 5;
   }

}
