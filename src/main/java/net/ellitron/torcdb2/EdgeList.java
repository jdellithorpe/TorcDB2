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

import edu.stanford.ramcloud.*;
import edu.stanford.ramcloud.ClientException.*;
import edu.stanford.ramcloud.multiop.*;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.*;

/**
 * A collection of static methods for reading and writing edge lists to
 * RAMCloud in TorcDB2.
 *
 * Here is an illustration of the layout of edge lists in RAMCloud and what
 * they look like before and after a split:
 *
 * +-------------------------------------------------------------------------+
 * |                                                                         |
 * |          RAMCloud Key           RAMCloud Value                          |
 * |                                                                         |
 * |                                   +numTailSegments                      |
 * |          +-----------+---+      +-v-+-------+-------+-------+---------+ |
 * | HeadSeg  |Key Prefix | 0 | +--> | 3 | edge0 | edge1 |  ...  | edge324 | |
 * |    +     +-----------+---+      +---+-------+-------+-------+---------+ |
 * |    |                                                                    |
 * |    v     +-----------+---+      +---------+---------+-------+---------+ |
 * | TailSeg0 |Key Prefix | 3 | +--> | edge325 | edge326 |  ...  | edge683 | |
 * |    +     +-----------+---+      +---------+---------+-------+---------+ |
 * |    |                                                                    |
 * |    v     +-----------+---+      +---------+---------+------+----------+ |
 * | TailSeg1 |Key Prefix | 2 | +--> | edge684 | edge685 |  ... | edge1245 | |
 * |    +     +-----------+---+      +---------+---------+------+----------+ |
 * |    |                                                                    |
 * |    v     +-----------+---+      +----------+----------+----+----------+ |
 * | TailSeg2 |Key Prefix | 1 | +--> | edge1246 | edge1247 | .. | edge1545 | |
 * |          +-----------+---+      +----------+----------+----+----------+ |
 * |                                                                         |
 * +----------------------------------+--------------------------------------+
 *                                    |
 *                       Split after edge160 creates
 *                       new tail seg w/ edges161-324
 *                                    |
 *                                    v
 * +----------------------------------+--------------------------------------+
 * |                                                                         |
 * |          RAMCloud Key           RAMCloud Value                          |
 * |                                                                         |
 * |                                   +numTailSegments                      |
 * |          +-----------+---+      +-v-+-------+-------+-------+---------+ |
 * | HeadSeg  |Key Prefix | 0 | +--> | 4 | edge0 | edge1 |  ...  | edge160 | |
 * |    +     +-----------+---+      +---+-------+-------+-------+---------+ |
 * |    |                                                                    |
 * |    v     +-----------+---+      +---------+---------+-------+---------+ |
 * | TailSeg0 |Key Prefix | 4 | +--> | edge161 | edge162 |  ...  | edge324 | |
 * |    +     +-----------+---+      +---------+---------+-------+---------+ |
 * |    |                                                                    |
 * |    v     +-----------+---+      +---------+---------+-------+---------+ |
 * | TailSeg1 |Key Prefix | 3 | +--> | edge325 | edge326 |  ...  | edge683 | |
 * |    +     +-----------+---+      +---------+---------+-------+---------+ |
 * |    |                                                                    |
 * |    v     +-----------+---+      +---------+---------+------+----------+ |
 * | TailSeg2 |Key Prefix | 2 | +--> | edge684 | edge685 |  ... | edge1245 | |
 * |    +     +-----------+---+      +---------+---------+------+----------+ |
 * |    |                                                                    |
 * |    v     +-----------+---+      +----------+----------+----+----------+ |
 * | TailSeg3 |Key Prefix | 1 | +--> | edge1246 | edge1247 | .. | edge1545 | |
 * |          +-----------+---+      +----------+----------+----+----------+ |
 * |                                                                         |
 * +-------------------------------------------------------------------------+
 *
 *
 * @author Jonathan Ellithorpe (jde@cs.stanford.edu)
 */
public class EdgeList {

  /*
   * Limit placed on the number of bytes allowed to be stored in a single RAMCloud object. For
   * large edge lists, if this limit is too high then simple operations like prepend will need to
   * read a lot of data only to add a relatively small number of bytes to the list. If this limit
   * is too small, then operations like reading all of the edges in the list will require reading
   * many RAMCloud objects and incur high read overhead.
   */
  private static final int DEFAULT_SEGMENT_SIZE_LIMIT = 1 << 13;

  /*
   * When a RAMCloud object exceeds its size limit (DEFAULT_SEGMENT_SIZE_LIMIT), the object is
   * split into two parts. This parameter specifies the byte offset into the segment where the
   * split should occur. Most of the time this will not land exactly between two edges in the list,
   * and in this case the nearest boundary to the split point is selected, unless that happens to
   * be past the size limit, in which case the lower boundary is selected.
   */
  private static final int DEFAULT_SEGMENT_TARGET_SPLIT_POINT = 0;

  /*
   * Limit placed on the number of asynchronous reads that can be outstanding at any one time.
   */
  private static final int DEFAULT_MAX_ASYNC_READS = 1 << 7;

  /*
   * Limit placed on the maximum size of multireads.
   */
  private static final int DEFAULT_MAX_MULTIREAD_SIZE = 1 << 11;

  public static boolean prepend(
      RAMCloudTransaction rctx,
      RAMCloud client,
      long rcTableId,
      byte[] keyPrefix,
      UInt128 neighborId, 
      byte[] serializedProperties) {
    return prepend(rctx, client, rcTableId, keyPrefix, neighborId, serializedProperties, 
        DEFAULT_SEGMENT_SIZE_LIMIT, DEFAULT_SEGMENT_TARGET_SPLIT_POINT);
  }

  /**
   * Prepends the edge represented by the given neighbor vertex and serialized properties to this
   * edge list. If the edge list does not exist, then this method will create a new one and return
   * true to indicate that a new edge list has been created (otherwise the method returns false).
   *
   * Note that this method allows multiple edges with the same neighbor vertex to exist in the list
   * (it will not check for duplicates).
   *
   * This method will perform the operation in a transaction context if it is given (rctx) and
   * non-transactionally if not.
   *
   * @param rctx RAMCloud transaction in which to perform the operation if not null.
   * @param client RAMCloud client to use to perform the operation if rctx is null.
   * @param rcTableId The table in which the edge list is (to be) stored.
   * @param keyPrefix Key prefix for the edge list.
   * @param neighborId Remote vertex Id for this edge.
   * @param serializedProperties Pre-serialized properties for this edge.
   * @param segment_size_limit Limit on the max size of segments.
   * @param segment_target_split_point Where to split when splitting is needed.
   *
   * @return True if a new edge list was created, false otherwise.
   */
  public static boolean prepend(
      RAMCloudTransaction rctx,
      RAMCloud client,
      long rcTableId,
      byte[] keyPrefix,
      UInt128 neighborId, 
      byte[] serializedProperties,
      int segment_size_limit,
      int segment_target_split_point) {
    /* Read out the head segment. */
    ByteBuffer headSeg;
    byte[] headSegKey = getSegmentKey(keyPrefix, 0);
    boolean newList = false;
    try {
      RAMCloudObject headSegObj;
      if (rctx != null)
        headSegObj  = rctx.read(rcTableId, headSegKey);
      else
        headSegObj  = client.read(rcTableId, headSegKey);

      if (headSegObj != null) {
        headSeg = ByteBuffer.allocate(headSegObj.getValueBytes().length)
            .order(ByteOrder.LITTLE_ENDIAN)
            .put(headSegObj.getValueBytes());
        headSeg.flip();
      } else {
        headSeg = ByteBuffer.allocate(Integer.BYTES)
            .order(ByteOrder.LITTLE_ENDIAN).putInt(0);
        headSeg.flip();
        newList = true;
      }

      //System.out.println(String.format("Head Seg Size: %d", headSeg.capacity()));
    } catch (ClientException e) {
      throw new RuntimeException(e);
    }

    int serializedEdgeLength =
        UInt128.BYTES + Short.BYTES + serializedProperties.length;
      
    //System.out.println(String.format("Serialized Edge Length: %d", serializedEdgeLength));

    ByteBuffer serializedEdge = ByteBuffer.allocate(serializedEdgeLength)
        .order(ByteOrder.LITTLE_ENDIAN);
    serializedEdge.put(neighborId.toByteArray());
    serializedEdge.putShort((short) serializedProperties.length);
    serializedEdge.put(serializedProperties);
    serializedEdge.flip();

    /* Prepend edge to head segment. */
    ByteBuffer prependedSeg =
        ByteBuffer.allocate(serializedEdge.capacity() + headSeg.capacity())
        .order(ByteOrder.LITTLE_ENDIAN);
    int majorSegments = headSeg.getInt();
    prependedSeg.putInt(majorSegments);
    prependedSeg.put(serializedEdge);
    prependedSeg.put(headSeg);
    prependedSeg.flip();

    /* Check if we need to split the head segment. */
    if (prependedSeg.capacity() <= segment_size_limit) {
      /* Common case, don't need to split. */
      if (rctx != null)
        rctx.write(rcTableId, headSegKey, prependedSeg.array());
      else
        client.write(rcTableId, headSegKey, prependedSeg.array(), null);
    } else {
      /* Head segment is too big, we need to find a good split point. In some
       * special cases we won't be able to split, like when the segment is just
       * one enormous edge. The following code sets splitIndex to the right
       * point to split the head segment. */
      int splitIndex = prependedSeg.capacity();
      int currentNumTailSegments = prependedSeg.getInt();
      while (prependedSeg.hasRemaining()) {
        int edgeStartPos = prependedSeg.position();
        int nextEdgeStartPos = edgeStartPos + UInt128.BYTES + Short.BYTES 
            + prependedSeg.getShort(edgeStartPos + UInt128.BYTES);

        if (nextEdgeStartPos >= segment_target_split_point) {
          /*
           * The current edge either stradles the split point, or is right up
           * against it.
           *
           *                                       nextEdgeStartPos
           *            <--left-->          <--right-->   V
           * ------|--------------------|-----------------|--------
           *       ^                    ^
           * edgeStartPos     DEFAULT_SEGMENT_TARGET_SPLIT_POINT
           */
          int left = segment_target_split_point - edgeStartPos;
          int right = nextEdgeStartPos - segment_target_split_point;

          if (right < left) {
            /* Target split point is closer to the start of the next edge in
             * the list than the start of this edge. In this case we generally
             * want to split at the start of the next edge, except for a
             * special case handled here. */
            if (nextEdgeStartPos > segment_size_limit) {
              /* Special case, the current edge extends beyond the size limit.
               * To still enforce the size limit policy we choose not to keep
               * this edge in the head segment. */
              splitIndex = edgeStartPos;
              break;
            } else {
              splitIndex = nextEdgeStartPos;
              break;
            }
          } else {
            /* Target split point is closer to the start of this edge than the
             * next. In this case we choose to make this edge part of the newly
             * created segment. */
            splitIndex = edgeStartPos;
            break;
          }
        }

        prependedSeg.position(nextEdgeStartPos);
      }

      prependedSeg.rewind();

      if (splitIndex == prependedSeg.capacity()) {
        /* We have chosen not to split this segment. */
        if (rctx != null)
          rctx.write(rcTableId, headSegKey, prependedSeg.array());
        else
          client.write(rcTableId, headSegKey, prependedSeg.array(), null);
      } else {
        /* Split based on splitIndex. */
        ByteBuffer newHeadSeg = ByteBuffer.allocate(splitIndex)
            .order(ByteOrder.LITTLE_ENDIAN);
        ByteBuffer newTailSeg = ByteBuffer.allocate(prependedSeg.capacity() 
            - splitIndex).order(ByteOrder.LITTLE_ENDIAN);

        //System.out.println(String.format("New Head Segment: %d", newHeadSeg.capacity()));
        //System.out.println(String.format("New Tail Segment: %d", newTailSeg.capacity()));

        int newNumTailSegments = currentNumTailSegments + 1;

        newHeadSeg.put(prependedSeg.array(), 0, splitIndex);
        newHeadSeg.rewind();
        newHeadSeg.putInt(newNumTailSegments);

        newTailSeg.put(prependedSeg.array(), splitIndex,
            prependedSeg.capacity() - splitIndex);

        byte[] newTailSegKey = getSegmentKey(keyPrefix, newNumTailSegments);

        if (rctx != null) {
          rctx.write(rcTableId, headSegKey, newHeadSeg.array());
          rctx.write(rcTableId, newTailSegKey, newTailSeg.array());
        } else {
          client.write(rcTableId, headSegKey, newHeadSeg.array(), null);
          client.write(rcTableId, newTailSegKey, newTailSeg.array(), null);
        }
      }
    }

    return newList;
  }

  /**
   * This method takes an array of edges and creates the same set of RAMCloud key / value pairs
   * that would be created had the edges been added one by one, starting with the first edge in the
   * array and ending with the last edge in the array, and writes the resulting key/value pairs
   * into the given RAMCloud image file. The resulting key / value pairs written into the image
   * should exactly match the key / value pairs that would be in RAMCloud had the
   * EdgeList.prepend() method been called for each of the edges, starting with the first in the
   * array and ending with the last in the array.
   *
   * @param edgeListTableOS The image file to write to.
   * @param keyPrefix Key prefix for the edge list.
   * @param neighborIds Remote vertex Ids for this edge list. List is in the order these edges
   * would have been added in (0th edge is the first edge added).
   * @param propMaps Property maps for the edges. Same ordering as neighborIds. Can be an empty
   * list, which signals that these edges do not have properties.
   */
  public static void writeListToFile(
      OutputStream edgeListTableOS,
      byte[] keyPrefix,
      List<UInt128> neighborIds, 
      List<byte[]> serializedPropList) {
    /* General strategy here is to simulate the prepending of edges by prepending edge lengths
     * instead of actual edges and split by the sum of the edge lengths in the list, and thus
     * calculate how many edges should go in each segment had they been prepended one by one in
     * that order. Using this information, we can then take the list of edges and directly pack the
     * correct number into the correct segments. 
     */

    // We only need to keep track of the series of edge lengths in the head segment, because once a
    // tail segment is "pinched" off after a split of the head segment, it will remain unchanged,
    // and the only information we need to keep around is the number of edges that made it into the
    // segment.
    LinkedList<Integer> headSegEdgeLengths = new LinkedList<>();
    
    // As we split off tail segments from the head, we record the number of edges that made it into
    // the resulting tail segment in this list. Elements are added to the end of this list as the
    // simulation proceeds, therefore the first element of the list represents the number of edges
    // in the last segment of the edge list, and the last element represents the number of edges in
    // the head segment.
    ArrayList<Integer> edgesPerSegment = new ArrayList<>();
    
    // Here we record the sizes, in bytes, of segments created during the simulation (in the same
    // ordering as the edgesPerSegment list). These data are used after the simulation to allocate
    // appropriately sized ByteBuffers that represent the edge list segments. Although these data
    // could be derived from edgesPerSegment and the argument list of edges post-simulation, this
    // information is calculated already during the simulation and so it is more efficient to
    // simply save it for later use.
    ArrayList<Integer> segmentSizes = new ArrayList<>();
   
    // Head segment starts with an integer field containing the total number of tail segments for
    // this edge list, so this is our starting length for the head segment.
    int headSegLen = Integer.BYTES;

    // Simulate prepending the edges, starting with the first in the argument list and ending with
    // the last in the argument list.
    for (int i = 0; i < neighborIds.size(); i++) {
      int edgeLength;
      if (serializedPropList.size() > 0) {
        edgeLength = UInt128.BYTES + Short.BYTES + serializedPropList.get(i).length;
      } else {
        edgeLength = UInt128.BYTES + Short.BYTES;
      }
      headSegLen += edgeLength;
      headSegEdgeLengths.addFirst(edgeLength);

      if (headSegLen >= DEFAULT_SEGMENT_SIZE_LIMIT) {
        int edgesInNewTailSeg = 0;
        // In the head segment, edges start after the integer field that stores
        // the total number of tail segments for the edge list.
        int edgeStartPos = Integer.BYTES;
        int nextEdgeStartPos = Integer.BYTES;
        for (int j = 0; j < headSegEdgeLengths.size(); j++) {
          edgeStartPos = nextEdgeStartPos;
          nextEdgeStartPos = edgeStartPos + headSegEdgeLengths.get(j);

          if (nextEdgeStartPos >= DEFAULT_SEGMENT_TARGET_SPLIT_POINT) {
            /*
             * The current edge either stradles the split point, or is right up against it.
             *
             *                                       nextEdgeStartPos
             *            <--left-->          <--right-->   V
             * ------|--------------------|-----------------|--------
             *       ^                    ^
             * edgeStartPos     DEFAULT_SEGMENT_TARGET_SPLIT_POINT
             */
            int left = DEFAULT_SEGMENT_TARGET_SPLIT_POINT - edgeStartPos;
            int right = nextEdgeStartPos - DEFAULT_SEGMENT_TARGET_SPLIT_POINT;

            if (right < left) {
              /* Target split point is closer to the start of the next edge in the list than the
               * start of this edge. In this case we generally want to split at the start of the
               * next edge, except for a special case handled here. */
              if (nextEdgeStartPos > DEFAULT_SEGMENT_SIZE_LIMIT) {
                /* Special case, the current edge extends beyond the size limit.  To still enforce
                 * the size limit policy we choose not to keep this edge in the head segment. */
                edgesInNewTailSeg = headSegEdgeLengths.size() - j;
                break;
              } else {
                edgesInNewTailSeg = headSegEdgeLengths.size() - (j + 1);
                break;
              }
            } else {
              /* Target split point is closer to the start of this edge than the next. In this case
               * we choose to make this edge part of the newly created segment. */
              edgesInNewTailSeg = headSegEdgeLengths.size() - j;
              break;
            }
          }
        }

        // At this point we have figured out how many edges go in the new tail segment (which could
        // potentially be zero, which means the edge is actually NOT split. In this case just move
        // on).
        
        if (edgesInNewTailSeg > 0) {
          edgesPerSegment.add(edgesInNewTailSeg);

          int segmentSize = 0;
          for (int j = 0; j < edgesInNewTailSeg; j++) {
            segmentSize += headSegEdgeLengths.getLast();
            headSegEdgeLengths.removeLast();
          }
          headSegLen -= segmentSize;
          
          segmentSizes.add(segmentSize);
        }
      } // if (headSegLen >= DEFAULT_SEGMENT_SIZE_LIMIT) 
    } // for (int i = 0; i < neighborIds.size(); i++) 

    // Whatever is left in headSegEdgeLengths after the simulation is over represents the final
    // state of the head segment.
    edgesPerSegment.add(headSegEdgeLengths.size());
    segmentSizes.add(headSegLen);
  
    // Now edgesPerSegment and segmentSizes contain the metadata for all the segments that
    // represent this edge list in RAMCloud. Time to pack the edges into ByteBuffers and write them
    // out to the edge image file.

    int neighborListSegOffset = 0;
    ByteBuffer keyLen = ByteBuffer.allocate(Integer.BYTES).order(ByteOrder.LITTLE_ENDIAN);
    ByteBuffer valLen = ByteBuffer.allocate(Integer.BYTES).order(ByteOrder.LITTLE_ENDIAN);
    for (int i = 0; i < edgesPerSegment.size(); i++) {
      int edgesInSegment = edgesPerSegment.get(i);
      int segmentSize = segmentSizes.get(i);
      ByteBuffer segment = ByteBuffer.allocate(segmentSize).order(ByteOrder.LITTLE_ENDIAN);
      
      byte[] segKey;
      if (i == edgesPerSegment.size() - 1) {
        // This is the head segment.
        segKey = getSegmentKey(keyPrefix, 0);
        // Special field in head segment for total number of tail segments.
        segment.putInt(edgesPerSegment.size() - 1);
      } else {
        // This is a tail segment.
        segKey = getSegmentKey(keyPrefix, i + 1);
      }

      // Remember that the given edges were prepended, so a given segment actually starts with the
      // edges in the end of the range and finishes with the first edge in the range.
      for (int j = edgesInSegment - 1; j >= 0; j--) {
        UInt128 neighborId = neighborIds.get(neighborListSegOffset + j);
        if (serializedPropList.size() > 0) {
          byte[] serializedProps = serializedPropList.get(neighborListSegOffset + j);
          segment.put(neighborId.toByteArray());
          segment.putShort((short) serializedProps.length);
          segment.put(serializedProps);
        } else {
          segment.put(neighborId.toByteArray());
          segment.putShort((short) 0);
        }
      }

      byte[] segVal = segment.array();

      keyLen.rewind();
      keyLen.putInt(segKey.length);
      valLen.rewind();
      valLen.putInt(segVal.length);

      try {
        edgeListTableOS.write(keyLen.array());
        edgeListTableOS.write(segKey);
        edgeListTableOS.write(valLen.array());
        edgeListTableOS.write(segVal);
      } catch (IOException e) {
        throw new RuntimeException(e);
      }

      neighborListSegOffset += edgesInSegment;
    }
  }

  /* Metadata we want to keep track of for MutliReadObjects. */
  private static class MultiReadSpec {
    public byte[] keyPrefix;
    public Vertex baseVertex;
    public String neighborLabel;
    public boolean isHeadSeg;

    public MultiReadSpec(byte[] keyPrefix, 
        Vertex baseVertex, 
        String neighborLabel, 
        boolean isHeadSeg) {
      this.keyPrefix = keyPrefix;
      this.baseVertex = baseVertex;
      this.neighborLabel = neighborLabel;
      this.isHeadSeg = isHeadSeg;
    }
  }

  public static TraversalResult batchReadSingleThreaded(
      RAMCloudTransaction rctx,
      RAMCloud client,
      long rcTableId,
      Collection<Vertex> vCol,
      String eLabel, 
      Direction dir, 
      boolean parseProps,
      String ... nLabels) {
    return batchReadThread(rctx, client, null, rcTableId, vCol, 1, 0, eLabel, dir, parseProps, nLabels);
  }

  public static TraversalResult batchReadMultiThreaded(
      ExecutorService threadPool,
      RAMCloudTransaction rctx,
      RAMCloud client,
      Lock rclock,
      long rcTableId,
      Collection<Vertex> vCol,
      String eLabel, 
      Direction dir, 
      boolean parseProps,
      String ... nLabels) {
    long startTime = System.nanoTime();
    long invokeTime = 0;
    long threadExecutionTime = 0;
    long reduceTime = 0;

    final int nThreads = 7;

    long invokeStartTime = System.nanoTime();

		List<Callable<TraversalResult>> callables = new ArrayList<>();
    for (int i = 0; i < nThreads; i++) {
      final int threadId = i;
      callables.add(() -> {
        return batchReadThread(rctx, client, rclock, rcTableId, vCol, nThreads, threadId, eLabel, dir, parseProps, nLabels);
      });
    }

    final Map<Vertex, List<Vertex>> vMap = new HashMap<>();
    final Map<Vertex, List<Map<Object, Object>>> pMap;
    if (parseProps)
      pMap = new HashMap<>();
    else
      pMap = null;
    final Set<Vertex> vSet = new HashSet<>();
    
    try {
      List<Future<TraversalResult>> futureResults = threadPool.invokeAll(callables);
      invokeTime = System.nanoTime() - invokeStartTime;

      long threadStartTime = System.nanoTime(); 
      List<TraversalResult> results = new ArrayList<>(nThreads);
      for (Future<TraversalResult> futureResult : futureResults) {
        results.add(futureResult.get());
      }
      threadExecutionTime = System.nanoTime() - threadStartTime;

      long reduceStartTime = System.nanoTime(); 
      for (TraversalResult result : results) {
        vMap.putAll(result.vMap);
        if (parseProps)
          pMap.putAll(result.pMap);
        vSet.addAll(result.vSet);
      }
      reduceTime = System.nanoTime() - reduceStartTime;
    } catch (Exception e) {
      throw new RuntimeException(e);
    }

    long endTime = System.nanoTime();

//    System.out.println(String.format(
//          "{\"tag\": \"EdgeList.batchReadMultiThreaded()\", "
//          + "\"invokeTime\": %d, "
//          + "\"threadExecutionTime\": %d, "
//          + "\"reduceTime\": %d, "
//          + "\"time\": %d}",
//          invokeTime/1000,
//          threadExecutionTime/1000,
//          reduceTime/1000,
//          (endTime - startTime)/1000));

    return new TraversalResult(vMap, pMap, vSet);
  }

  /**
   * Batch reads in parallel all of the edges for all the given vertices.
   *
   * This method will perform the operation in a transaction context if it is given (rctx) and
   * non-transactionally if not.
   *
   * @param rctx RAMCloud transaction in which to perform the operation if not null.
   * @param client RAMCloud client to use to perform the operation if rctx is null.
   * @param rcTableId The table in which the edge list is (to be) stored.
   * @param keyPrefix List of key prefixes for the edge lists.
   * @param parseProps Whether or not to parse property information from edges.
   *
   * @return List of all the Edges contained in the edge lists.
   */ 
  public static TraversalResult batchReadThread(
      RAMCloudTransaction rctx,
      RAMCloud client,
      Lock rclock,
      long rcTableId,
      Collection<Vertex> vCol,
      int numThreads,
      int threadId,
      String eLabel, 
      Direction dir,
      boolean parseProps,
      String ... nLabels) {
    long startTime = System.nanoTime();
    long totalReadTime = 0;
    long totalLockWaitTime = 0;
    long totalRequests = 0;
    long totalDeserializationTime = 0;
    long totalEdgeParseTime = 0;
    long totalDedupTime = 0;
    long totalEdges = 0;
    long totalHeadSegments = 0;
    long totalTailSegments = 0;
    long samplerCounter = 0;

    List<byte[]> keyPrefixes = GraphHelper.getEdgeListKeyPrefixes(vCol, eLabel, dir, nLabels);

    // Future read requests are appended to this queue as we figure out what we need to read. We
    // store either RAMCloudTransactionReadOps or MultiReadObjets in this queue, depending on if we
    // are executing in a transaction context or not.
    LinkedList<Object> requestQ = new LinkedList<>();
    LinkedList<MultiReadSpec> specQ = new LinkedList<>();
    
    Map<Vertex, List<Vertex>> vMap = new HashMap<>();
    Map<Vertex, List<Map<Object, Object>>> pMap = null;
    if (parseProps)
      pMap = new HashMap<>();

    Map<UInt128, Vertex> vDedupMap = new HashMap<>();

    Set<Vertex> vSet = new HashSet<>();

    int index = 0;
    for (String nLabel : nLabels) {
      for (Vertex vertex : vCol) {
        if (index % numThreads == threadId)
          specQ.addLast(new MultiReadSpec(keyPrefixes.get(index), vertex, nLabel, true));
        index++;
      }
    }

    /* Add head segments to queue and prepare edgeMap. */
    for (int i = 0; i < keyPrefixes.size(); i++) {
      if (i % numThreads == threadId) {
        byte[] headSegKey = getSegmentKey(keyPrefixes.get(i), 0);
        if (rctx != null)
          requestQ.addLast(new RAMCloudTransactionReadOp(rctx, rcTableId, headSegKey, true));
        else
          requestQ.addLast(new MultiReadObject(rcTableId, headSegKey));
      }
    }

    /* Go through request queue and read at most DEFAULT_MAX_MULTIREAD_SIZE at a time. */
    ByteBuffer seg = ByteBuffer.allocate(DEFAULT_SEGMENT_SIZE_LIMIT*2)
      .order(ByteOrder.LITTLE_ENDIAN);
    while (requestQ.size() > 0) {
      int batchSize = Math.min(requestQ.size(), DEFAULT_MAX_MULTIREAD_SIZE);

      RAMCloudObject[] rcobjs = new RAMCloudObject[batchSize];
      if (rctx != null) {
        for (int i = 0; i < batchSize; i++) {
          RAMCloudTransactionReadOp readOp = (RAMCloudTransactionReadOp)requestQ.removeFirst();
          try {
            rcobjs[i] = readOp.getValue();
          } catch (ClientException e) {
            throw new RuntimeException(e);
          } finally {
            readOp.close();
          }
        }
      } else {
        MultiReadObject[] mrobjs = new MultiReadObject[batchSize];
        for (int i = 0; i < batchSize; i++) {
          mrobjs[i] = (MultiReadObject)requestQ.removeFirst();
          rcobjs[i] = (RAMCloudObject)mrobjs[i];
        }

        long multireadStartTime;
        long multireadEndTime;
        long lockStartTime = 0;
        long lockEndTime = 0;
        if (rclock != null) {
          lockStartTime = System.nanoTime();
          rclock.lock();
          multireadStartTime = System.nanoTime();
          client.read(mrobjs);
          multireadEndTime = System.nanoTime();
          rclock.unlock();
          lockEndTime = System.nanoTime();
        } else {
          multireadStartTime = System.nanoTime();
          client.read(mrobjs);
          multireadEndTime = System.nanoTime();
        }
        totalReadTime += multireadEndTime - multireadStartTime;
        totalLockWaitTime += lockEndTime - lockStartTime;
        totalRequests += mrobjs.length;
      }

      long deserStartTime = System.nanoTime();

      /* Process this batch, adding more requests to the queue if needed. */
      for (int i = 0; i < batchSize; i++) {
        MultiReadSpec spec = specQ.removeFirst();

        // See if this object exists. If not, short circuit this loop iteration.
        if (rctx != null) {
          if (rcobjs[i] == null)
            continue;
        } else {
          if (((MultiReadObject)rcobjs[i]).getStatus() != Status.STATUS_OK) {
            if (((MultiReadObject)rcobjs[i]).getStatus() == Status.STATUS_OBJECT_DOESNT_EXIST) {
              continue;
            } else {
              throw new RuntimeException("Segment had status " + ((MultiReadObject)rcobjs[i]).getStatus());
            }
          }
        }

        List<Vertex> vList;
        List<Map<Object, Object>> pList = null;
        if (vMap.containsKey(spec.baseVertex)) {
          vList = vMap.get(spec.baseVertex);
          if (parseProps)
            pList = pMap.get(spec.baseVertex);
        } else {
          vList = new ArrayList<>();
          vMap.put(spec.baseVertex, vList);
          if (parseProps) {
            pList = new ArrayList<>();
            pMap.put(spec.baseVertex, pList);
          }
        }

        seg.clear();
        seg.put(rcobjs[i].getValueBytes());
        seg.flip();

        if (spec.isHeadSeg) {
          /* Queue up async. reads for tail segments. */
          int numTailSegments = seg.getInt();
          for (int j = numTailSegments; j > 0; --j) {
            byte[] tailSegKey = getSegmentKey(spec.keyPrefix, j);
            if (rctx != null)
              requestQ.addLast(new RAMCloudTransactionReadOp(rctx, rcTableId, tailSegKey, true));
            else
              requestQ.addLast(new MultiReadObject(rcTableId, tailSegKey));
            spec.isHeadSeg = false;
            specQ.addLast(spec);
          }

          totalHeadSegments++;
        } else
          totalTailSegments++;

        long edgeParseStartTime = System.nanoTime();

        byte[] neighborIdBytes = new byte[UInt128.BYTES];
        while (seg.hasRemaining()) {
          totalEdges++;

          seg.get(neighborIdBytes);

          
          UInt128 neighborId = new UInt128(neighborIdBytes);

//          long dedupStartTime = System.nanoTime();

//          long dedupGetStartTime = System.nanoTime();
          Vertex nv = vDedupMap.get(neighborId);
//          long dedupGetEndTime = System.nanoTime();

//          if (samplerCounter % 999 == 0)
//            System.out.println(dedupGetEndTime - dedupGetStartTime);
//          samplerCounter++;

          if (nv != null) {
            vList.add(nv);
          } else {
            Vertex v = new Vertex(neighborId, spec.neighborLabel);
//            if (samplerCounter % 999 == 0) {
//              long vListAddStartTime = System.nanoTime();
//              vList.add(v);
//              long vListAddEndTime = System.nanoTime();
//              System.out.println(vListAddEndTime - vListAddStartTime);
//            } else {
//              vList.add(v);
//            }
//            samplerCounter++;
            vList.add(v);
            vDedupMap.put(neighborId, v);
            vSet.add(v);
          }

//          totalDedupTime += System.nanoTime() - dedupStartTime;

          short propLen = seg.getShort();

          if (parseProps) {
            if (propLen > 0) {
              pList.add((Map<Object, Object>)GraphHelper.deserializeObject(seg));
            } else {
              pList.add(new HashMap<Object, Object>());
            }
          } else {
            seg.position(seg.position() + propLen);
          }
        }

        totalEdgeParseTime += System.nanoTime() - edgeParseStartTime;
      }
      long deserEndTime = System.nanoTime();
      totalDeserializationTime += deserEndTime - deserStartTime;
    }

    long endTime = System.nanoTime();

    if (totalLockWaitTime == 0)
      totalLockWaitTime = totalReadTime;

//    System.out.println(String.format(
//          "{\"tag\": \"EdgeList.batchReadThread()\", "
//          + "\"numThreads\": %d, "
//          + "\"threadId\": %d, "
//          + "\"totalRequests\": %d, "
//          + "\"totalHeadSegments\": %d, "
//          + "\"totalTailSegments\": %d, "
//          + "\"totalReadTime\": %d, "
//          + "\"totalLockWaitTime\": %d, "
//          + "\"totalDeserializationTime\": %d, "
//          + "\"totalEdgeParseTime\": %d, "
//          + "\"totalDedupTime\": %d, "
//          + "\"totalEdges\": %d, "
//          + "\"parseProps\": %s, "
//          + "\"avgDeserTimePerEdge\": %.3f, "
//          + "\"avgDeserTimePerSeg\": %.3f, "
//          + "\"totalTime\": %d, "
//          + "\"percentTimeRead\": %.2f, "
//          + "\"percentTimeLockWait\": %.2f, "
//          + "\"percentTimeDeser\": %.2f}",
//          numThreads,
//          threadId,
//          totalRequests,
//          totalHeadSegments,
//          totalTailSegments,
//          totalReadTime/1000,
//          (totalLockWaitTime/1000) - (totalReadTime/1000),
//          totalDeserializationTime/1000,
//          totalEdgeParseTime/1000,
//          totalDedupTime/1000,
//          totalEdges,
//          parseProps,
//          (double)totalDeserializationTime/(double)totalEdges/1000.0,
//          (double)totalDeserializationTime/(double)totalRequests/1000.0,
//          (endTime - startTime)/1000,
//          (double)totalReadTime/(double)(endTime - startTime),
//          (double)(totalLockWaitTime - totalReadTime)/(double)(endTime - startTime),
//          (double)totalDeserializationTime/(double)(endTime - startTime)));

    return new TraversalResult(vMap, pMap, vSet);
  }

  /**
   * Creates a RAMCloud key for the given edge list segment.
   *
   * @param keyPrefix RAMCloud key prefix for this list.
   * @param segmentNumber Number of the segment.
   *
   * @return Byte array representing the RAMCloud key.
   */
  private static byte[] getSegmentKey(byte[] keyPrefix, int segmentNumber) {
    ByteBuffer buffer =
        ByteBuffer.allocate(keyPrefix.length + Integer.BYTES)
        .order(ByteOrder.LITTLE_ENDIAN);
    buffer.put(keyPrefix);
    buffer.putInt(segmentNumber);
    return buffer.array();
  }
}
