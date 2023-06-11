package com.mvohm.quadmatrix.investigations;

import com.mvohm.quadmatrix.DoubleMatrix;

import static com.mvohm.quadmatrix.test.AuxMethods.*;

import java.util.stream.DoubleStream;


/**
 * TODO
 * Nota bene!
 * При показателе степени -8
 * этот классический алгоритм решения дает совершенно неверный результат:

Delta: 1.0E-8
Matrix:
  1,0000000000000000, 0,9999999900000000
  0,9999999900000000, 0,9999999800000000
Vector:
  1,9999999900000000, 1,9999999700000000

First solution:
  -1,9999999700000000, 4,0000000000000000  (mse = 3,000e+00)
Refined solution:
  -1,9999999700000000, 4,0000000000000000  (mse = 3,000e+00) *

 * Ему следовало бы явно заявить, что эта система нерешаемая.
 * Там, где в разложении делается проверка на решаемость,
 *
      final double[] row_i = decompositionLU[i];
      if (row_i[i] != 0) {
        ... // Решаем
      } else {
        throwNonInvertibleError(i);
      }

 *  получается row_i[i] = 5.551115123125783E-17

 // 2022-12-19 17:40:49
   Исправлено в DoubleMatrixSolver.decomposeLU заменой
     if (row_i[i] != 0) {
   на
     if (Math.abs(row_i[i]) > 1e-16) {

 * @author misa
 *
 */

public class IllConditionedMatrixSolution {


  public static void main(String[] args) {

    // say(1 - Double.longBitsToDouble(0x3FEF_FFFF_FFFF_FFFFL));

    for (int exp = 2; exp < 10; exp++) {

      double delta = Math.pow(10, -exp);
      say("\nDelta: " + delta);
      double k = 1;

      double[][] badMatrix = {
          {1.00 * k,       (1 - delta) * k},
          {(1 - delta) * k,  (1 - 2 * delta) * k},
        };

      say("Matrix: \n\t%.16f, %.16f\n\t%.16f, %.16f",
          badMatrix[0][0], badMatrix[0][1], badMatrix[1][0], badMatrix[1][1]);

      double[] vector = {(2 - delta) * k, (2 - 3 * delta) * k};
      say("Vector: \n\t%.16f, %.16f\n",
          vector[0], vector[1]);


      DoubleMatrix m = new DoubleMatrix(badMatrix);
      m.solve(vector);
      double[] x0 = m.getDoubleSolution();
      m.solveAccurately(vector);
      double[] x1 = m.getDoubleSolution();
      say("First solution:  \n\t%.16f, %.16f  (mse = %.3e)", x0[0], x0[1], findMse(x0));
      say("Refined solution:\n\t%.16f, %.16f  (mse = %.3e)", x1[0], x1[1], findMse(x1));
    }
  }
  private static double findMse(double[] solution) {
    return Math.sqrt(DoubleStream.of(solution).map(v -> (1 - v) * (1 - v)).average().getAsDouble());
  }

}
