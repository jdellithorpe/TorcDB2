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

import java.util.Map;

/**
 * Represents a parsed edge from a RAMCloud edge list.
 *
 * @author Jonathan Ellithorpe (jde@cs.stanford.edu)
 */
public class ParsedEdge {
  public Map<Object,Object> properties;
  public UInt128 vertexId;

  public ParsedEdge(Map<Object,Object> properties, UInt128 vertexId) {
    this.properties = properties;
    this.vertexId = vertexId;
  }
}
