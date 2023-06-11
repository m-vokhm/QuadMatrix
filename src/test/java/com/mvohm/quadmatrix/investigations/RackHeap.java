package com.mvohm.quadmatrix.investigations;

import static com.mvohm.quadmatrix.test.AuxMethods.*;

import java.util.Random;

import com.mvohm.quadruple.Quadruple;


/**
 * Memory allocation test.
 *
 * To create 200_000_000 Quadruples,
 * the heap must be at least -Xms9216M -Xmx9216M
 *
 * Then the output will be something like the following:
 *
heap:   9 261 023 232 free:   9 164 386 440 used:      96 636 792, max:   9 261 023 232
       5,000 ms:               0 instances, heap:   9 261 023 232, used:     896 636 808
     179,000 ms:       2 000 000 instances, heap:   9 261 023 232, used:     944 955 224
     .....
     165,000 ms:      36 000 000 instances, heap:   9 261 023 232, used:   2 346 189 352
   2 007,000 ms:      38 000 000 instances, heap:   9 261 023 232, used:   2 370 872 792
     153,000 ms:      40 000 000 instances, heap:   9 261 023 232, used:   2 418 224 832
     .....
     146,000 ms:      96 000 000 instances, heap:   9 261 023 232, used:   4 715 281 960
  28 204,000 ms:      98 000 000 instances, heap:   9 261 023 232, used:   4 805 742 440
     146,000 ms:     100 000 000 instances, heap:   9 261 023 232, used:   4 853 263 584
     .....
     149,000 ms:     156 000 000 instances, heap:   9 261 023 232, used:   7 125 821 968
  19 899,000 ms:     158 000 000 instances, heap:   9 261 023 232, used:   7 156 196 792
     146,000 ms:     160 000 000 instances, heap:   9 261 023 232, used:   7 203 548 808
     .....
     144,000 ms:     198 000 000 instances, heap:   9 261 023 232, used:   8 776 372 616
heap:   9 261 023 232 free:     437 298 600 used:   8 823 724 632, max:   9 261 023 232
43,635 bytes per object
Elapsed 64,642 s
 *
 * @author M.Vokhmentev
 *
 */
public class RackHeap {

  private final static int SIZE = 200_000_000;

  public static void main(String[] args) {
    final Random rand = new Random(123);
    final long max1 = Runtime.getRuntime().maxMemory();
    final long heapSize1 = Runtime.getRuntime().totalMemory();
    final long heapFreeSize1 = Runtime.getRuntime().freeMemory();
    final long used1 = heapSize1 - heapFreeSize1;
    say("heap: %,15d free: %,15d used: %,15d, max: %,15d", heapSize1, heapFreeSize1, used1, max1);

    final Quadruple[] data = new Quadruple[SIZE];
    final double startTime = System.currentTimeMillis();
    double time = startTime;
    for (int i = 0; i < SIZE; i++) {
      data[i] = Quadruple.nextRandom(rand);
      if (i % (SIZE/100) == 0) {
        final long heapSize = Runtime.getRuntime().totalMemory();
        final long heapFreeSize = Runtime.getRuntime().freeMemory();
        final long used = heapSize - heapFreeSize;

        time = (System.currentTimeMillis() - time);
        say("%,12.3f ms: %,15d instances, heap: %,15d, used: %,15d", time, i, heapSize, used);
        time = (System.currentTimeMillis());
      }
    }

    final long max2 = Runtime.getRuntime().maxMemory();
    final long heapSize2 = Runtime.getRuntime().totalMemory();
    final long heapFreeSize2 = Runtime.getRuntime().freeMemory();
    final long used2 =heapSize2 - heapFreeSize2;
    say("heap: %,15d free: %,15d used: %,15d, max: %,15d", heapSize2, heapFreeSize2, used2, max2);
    say("%.3f bytes per object", (double)(used2 - used1) / SIZE);
    time = (System.currentTimeMillis());
    say("Elapsed %.3f s", (time - startTime)/1000);


  }

}
