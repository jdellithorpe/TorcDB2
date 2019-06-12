/* Copyright (c) 2015-2019 Stanford University
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR(S) DISCLAIM ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL AUTHORS BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
package net.ellitron.torcdb2;

import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.UUID;

/**
 * Container for a 128 bit unsigned integer stored in big-endian format. This
 * class is not intended to be used for doing arithmetic, comparison, or
 * logical operations and hence has no such methods for doing so. For that see
 * {@link java.math.BigInteger}. Its primary use is to convert other number
 * representations to their 128 bit unsigned representation, including Strings
 * and BigIntegers, and return a byte array of length 16 that contains the
 * unsigned representation.
 *
 * @author Jonathan Ellithorpe (jde@cs.stanford.edu)
 */
public class UInt128 implements Comparable<UInt128> {

  public static final int SIZE = 128;
  public static final int BYTES = SIZE / Byte.SIZE;

  private final long upperLong;
  private final long lowerLong;
  private final byte[] val;

  /**
   * Constructs a UInt128 from a byte array value.
   *
   * @param val Byte array in big-endian format.
   */
  public UInt128(byte[] val) {
    ByteBuffer buf = ByteBuffer.allocate(BYTES).order(ByteOrder.LITTLE_ENDIAN);
    buf.put(val);
    buf.flip();
    this.upperLong = buf.getLong();
    this.lowerLong = buf.getLong();
    this.val = buf.array();
  }

  /**
   * Constructs a UInt128 from two longs, one representing the upper 64 bits
   * and the other representing the lower 64 bits.
   *
   * @param upperLong Upper 64 bits.
   * @param lowerLong Lower 64 bits.
   */
  public UInt128(final long upperLong, final long lowerLong) {
    ByteBuffer buf = ByteBuffer.allocate(BYTES).order(ByteOrder.LITTLE_ENDIAN);
    buf.putLong(upperLong);
    buf.putLong(lowerLong);
    this.upperLong = upperLong;
    this.lowerLong = lowerLong;
    this.val = buf.array();
  }

  /**
   * Returns the upper 64 bits in a long.
   *
   * @return Upper 64 bits.
   */
  public long getUpperLong() {
    return upperLong;
  }

  /**
   * Returns the lower 64 bits in a long.
   *
   * @return Lower 64 bits.
   */
  public long getLowerLong() {
    return lowerLong;
  }

  /**
   * Returns a byte array containing this 128 bit unsigned integer in
   * big-endian format.
   *
   * @return Byte array containing this number in big-endian format.
   */
  public byte[] toByteArray() {
    return val;
  }

  /**
   * Returns a hexadecimal String representing this 128 bit unsigned integer
   * with the minimum number of digits.
   *
   * @return Formatted string representing this number.
   */
  @Override
  public String toString() {
    if (upperLong == 0) {
      return String.format("0x%X", lowerLong);
    } else {
      return String.format("0x%X%016X", upperLong, lowerLong);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int compareTo(UInt128 that) {
    int ret = Long.compareUnsigned(this.upperLong, that.upperLong);
    return ret != 0 ? ret : Long.compareUnsigned(this.lowerLong, that.lowerLong);
  }

  /**
   * Tests for bit-wise equality.
   */
  @Override
  public boolean equals(Object that) {
    if (that instanceof UInt128) {
      return ((UInt128) that).upperLong == this.upperLong
          && ((UInt128) that).lowerLong == this.lowerLong;
    } else {
      return false;
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int hashCode() {
    int hash = 7;
    hash = 83 * hash + (int) (this.upperLong ^ (this.upperLong >>> 32));
    hash = 83 * hash + (int) (this.lowerLong ^ (this.lowerLong >>> 32));
    return hash;
  }
}
