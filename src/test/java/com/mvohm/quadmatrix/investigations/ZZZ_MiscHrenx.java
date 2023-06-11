package com.mvohm.quadmatrix.investigations;

import static com.mvohm.quadmatrix.test.AuxMethods.*;

import java.math.BigDecimal;

import com.mvohm.quadruple.Quadruple;

public class ZZZ_MiscHrenx {

  public static void main(String[] args) {
//    testBigDecimalThrowsNullpointerException();
//    showHexFormOfDouble(1.3);
    convertFromBigDecimalToQuadruple(1.3);
  }

  private static void convertFromBigDecimalToQuadruple(double d) {
    showHexFormOfDouble(d);
//    final BigDecimal bd = BigDecimal.valueOf(d);
    final BigDecimal bd = new BigDecimal(d);
    say("BigDec: " + bd);
    final Quadruple q = new Quadruple(bd);
    say(q);
    say(q.toHexString());

  }

  private static void showHexFormOfDouble(double d) {
    final long bits = Double.doubleToLongBits(d);
    say("Hex form of %s: %s", d, hexStr(bits));
  }

  private static void testBigDecimalThrowsNullpointerException() {
    final BigDecimal a = new BigDecimal(123);
    final BigDecimal b = null;
    final BigDecimal c = a.add(b); //
    say("Sum = " + c);
  }

}
