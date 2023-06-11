package com.mvohm.quadmatrix;

import java.math.BigDecimal;

import com.mvohm.quadruple.Quadruple;

class Util {

  /****************************************************************************************
   *** Array Conversions double[] to BigDecimal[] and Quadruple[] and vice versa
   ****************************************************************************************/

  public static Quadruple[] convertToQuadruples(double[] source) {
    if (source == null)
      return null;
    final Quadruple[] result = new Quadruple[source.length];
    for (int i = 0; i < source.length; i++)
      result[i] = new Quadruple(source[i]);
    return result;
  }

  public static Quadruple[][] convertToQuadruples(double[][] source) {
    if (source == null)
      return null;
    final Quadruple[][] result = new Quadruple[source.length][];
    for (int i = 0; i < source.length; i++) {
      result[i] = convertToQuadruples(source[i]);
    }
    return result;
  }

  public static Quadruple[] convertToQuadruples(Number[] source) {
    if (source == null)
      return null;
    final Quadruple[] result = new Quadruple[source.length];
    if (source.getClass() == Double[].class) {
      for (int i = 0; i < source.length; i++)
        result[i] = new Quadruple((Double)source[i]);
    } else if (source.getClass() == Quadruple[].class) {
      for (int i = 0; i < source.length; i++)
        result[i] = new Quadruple((Quadruple)source[i]);
    } else if (source.getClass() == BigDecimal[].class) {
      for (int i = 0; i < source.length; i++)
        result[i] = new Quadruple((BigDecimal)source[i]);
    } else {
      throw new IllegalArgumentException("Invalid argument type: BigDecimal[], Double[] and Quadruple[] are only allowed");
    }
    return result;
  }

  public static Quadruple[][] convertToQuadruples(Number[][] source) {
    if (source == null)
      return null;
    final Quadruple[][] result = new Quadruple[source.length][];
    for (int i = 0; i < source.length; i++) {
      result[i] = convertToQuadruples(source[i]);
    }
    return result;
  }

  public static BigDecimal[] convertToBigDecimals(double[] source) {
    if (source == null)
      return null;
    final BigDecimal[] result = new BigDecimal[source.length];
    for (int i = 0; i < source.length; i++)
      result[i] = BigDecimal.valueOf(source[i]);
    return result;
  }

  public static BigDecimal[][] convertToBigDecimals(double[][] source) {
    if (source == null)
      return null;
    final BigDecimal[][] result = new BigDecimal[source.length][];
    for (int i = 0; i < source.length; i++) {
      result[i] = convertToBigDecimals(source[i]);
    }
    return result;
  }

  public static BigDecimal[] convertToBigDecimals(Number[] source) {
    if (source == null)
      return null;
    final BigDecimal[] result = new BigDecimal[source.length];
    if (source.getClass() == Double[].class) {
      for (int i = 0; i < source.length; i++)
        result[i] = BigDecimal.valueOf((Double)source[i]);
    } else if (source.getClass() == Quadruple[].class) {
      for (int i = 0; i < source.length; i++)
        result[i] = ((Quadruple)source[i]).bigDecimalValue();
    } else if (source.getClass() == BigDecimal[].class) {
      for (int i = 0; i < source.length; i++)
        result[i] = (BigDecimal)source[i];
    } else {
      throw new RuntimeException("AuxMethods.convertToQuadruples(Number[]) got a "
          + source.getClass().getSimpleName() + " argument");
    }
    return result;
  }

  public static BigDecimal[][] convertToBigDecimals(Number[][] source) {
    if (source == null)
      return null;
    final BigDecimal[][] result = new BigDecimal[source.length][];
    for (int i = 0; i < source.length; i++) {
      result[i] = convertToBigDecimals(source[i]);
    }
    return result;
  }

  public static Number[] convertToDoublesAsNumbers(double[] source) {
    if (source == null)
      return null;
    final Number[] result = new Number[source.length];
    for (int i = 0; i < source.length; i++) {
      result[i] = Double.valueOf(source[i]);
    }
    return result;
  }

  public static Number[][] convertToDoublesAsNumbers(double[][] source) {
    if (source == null)
      return null;
    final Number[][] result = new Number[source.length][];
    for (int i = 0; i < source.length; i++) {
      result[i] = convertToDoublesAsNumbers(source[i]);
    }
    return result;
  }

  public static double[] convertToDoubles(Number[] source) {
    if (source == null)
      return null;
    final double[] result = new double[source.length];
    for (int i = 0; i < source.length; i++)
      result[i] = source[i].doubleValue();
    return result;
  }

  public static double[][] convertToDoubles(Number[][] source) {
    if (source == null)
      return null;
    final double[][] result = new double[source.length][];
    for (int i = 0; i < source.length; i++) {
      result[i] = convertToDoubles(source[i]);
    }
    return result;
  }

  public static double[][] deepCopyOf(double[][] source) {
    if (source == null)
      return null;
    final double[][] result = new double[source.length][];
    for (int i = 0; i < source.length; i++)
      result[i] = source[i].clone();
    return result;
  }

  public static Quadruple[] deepCopyOf(Quadruple[] source) {
    if (source == null)
      return null;
    final Quadruple[] result = new Quadruple[source.length];
    for (int i = 0; i < source.length; i++)
      result[i] = new Quadruple(source[i]);
    return result;
  }

  public static Quadruple[][] deepCopyOf(Quadruple[][] source) {
    if (source == null)
      return null;
    final Quadruple[][] result = new Quadruple[source.length][];
    for (int i = 0; i < source.length; i++)
      result[i] = deepCopyOf(source[i]);
    return result;
  }

  public static BigDecimal[] deepCopyOf(BigDecimal[] source) {
    if (source == null)
      return null;
    final BigDecimal[] result = new BigDecimal[source.length];
    for (int i = 0; i < source.length; i++)
      result[i] = source[i];
    return result;
  }

  public static BigDecimal[][] deepCopyOf(BigDecimal[][] source) {
    if (source == null)
      return null;
    final BigDecimal[][] result = new BigDecimal[source.length][];
    for (int i = 0; i < source.length; i++)
      result[i] = deepCopyOf(source[i]);
    return result;
  }

  /****************************************************************************************
   *** Output to the console **************************************************************
   ****************************************************************************************/
   // Was used for debugging, may still come in handy

   /** == System.out.println(); */
   public static void say()  { System.out.println(); }

   /** == System.out.println(Object o);
    * @param o {@code Object} to print */
   public static void say(Object o)  { System.out.println(o); }

   /** == System.out.print(Object o);
    * @param o {@code Object} to print */
   public static void say_(Object o)   { System.out.print(o); }

   /** == System.out.println(String.format(String format, Object... args)
    * @param format a format string to format the {@code args}
    * @param args arguments to format
    * @see String#format(String, Object...)
    */
   public static void say(String format, Object... args) { System.out.println(String.format(format, args)); }

   /** == System.out.print(String.format(String format, Object... args)
    * @param format a format string to format the {@code args}
    * @param args arguments to format
    * @see String#format(String, Object...)
    */
   public static void say_(String format, Object... args) { System.out.print(String.format(format, args)); }


}
