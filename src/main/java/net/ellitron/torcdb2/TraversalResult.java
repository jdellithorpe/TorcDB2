/* Copyright (c) 2019-2019 Stanford University
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

import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Represents the result of a traversal.
 *
 * @author Jonathan Ellithorpe (jde@cs.stanford.edu)
 */
public class TraversalResult {
  public Map<Vertex, List<Vertex>> vMap;
  public Map<Vertex, List<Map<Object, Object>>> pMap;
  public Set<Vertex> vSet;

  public TraversalResult(
      Map<Vertex, List<Vertex>> vMap, 
      Map<Vertex, List<Map<Object, Object>>> pMap,
      Set<Vertex> vSet) {
    this.vMap = vMap;
    this.pMap = pMap;
    this.vSet = vSet;
  }

  public String toString() {
    StringBuilder sb = new StringBuilder();
    if (vMap != null)
      sb.append("TraversalResult[vMap: " + vMap.toString() + ", ");
    else
      sb.append("[vMap: null, ");

    if (pMap != null)
      sb.append("pMap: " + pMap.toString() + ", ");
    else
      sb.append("pMap: null, ");

    if (vSet != null)
      sb.append("vSet: " + vSet.toString() + "]");
    else
      sb.append("vSet: null]");


    return sb.toString();
  }
}
