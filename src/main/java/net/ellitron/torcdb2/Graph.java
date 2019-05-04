/* Copyright (c) 2015-2019 Stanford University
 *
 * Permission to use, copy, modify, and distribute this software for any purpose with or without
 * fee is hereby granted, provided that the above copyright notice and this permission notice
 * appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS
 * SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL
 * AUTHORS BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
 * OF THIS SOFTWARE.
 */
package net.ellitron.torcdb2;

import edu.stanford.ramcloud.*;
import edu.stanford.ramcloud.ClientException.*;
import edu.stanford.ramcloud.multiop.*;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Graph {

  /* Configured parameters. */
  private final Map<String, String> config;
  private final String graphName;
  private final String coordLocator;
  private final int dpdkPort;
  private final int numServers;

  /* Internal state. */
  private final RAMCloud client;                        // RAMCloud client interface.
  private RAMCloudTransaction tx;                       // Current RAMCloud tx context.
  private final long vertexTableId, edgeListTableId;    // RAMCloud tableIds.
  private OutputStream vertexTableOS, edgeListTableOS;  // For writing RAMCloud image files.

  /* **************************************************************************
   *
   * Construction and Destruction
   *
   * *************************************************************************/

  public Graph(final Map<String, String> config) {
    /* Set parameters according to the configuration. */
    this.config = config;

    if (config.containsKey("graphName"))
      this.graphName = config.get("graphName");
    else
      this.graphName = "default";

    if (config.containsKey("coordinatorLocator"))
      this.coordLocator = config.get("coordinatorLocator");
    else
      this.coordLocator = "tcp:host=127.0.0.1,port=12247";

    if (config.containsKey("dpdkPort"))
      this.dpdkPort = Integer.parseInt(config.get("dpdkPort"));
    else
      this.dpdkPort = -1; // Default is unset

    if (config.containsKey("numServers"))
      this.numServers = Integer.parseInt(config.get("numServers"));
    else
      this.numServers = 1;

    // Configure database for RAMCloud image file generation.
    if (config.containsKey("ramcloudImageGenerationOutputDir")) {
      String rcImageDir = config.get("ramcloudImageGenerationOutputDir");

      try {
        vertexTableOS = new BufferedOutputStream(new FileOutputStream(
              rcImageDir + "/" + graphName + "_vertexTable.img"));

        edgeListTableOS = new BufferedOutputStream(new FileOutputStream(
              rcImageDir + "/" + graphName + "_edgeListTable.img"));
      } catch (FileNotFoundException e) {
        throw new RuntimeException(e);
      }

      this.client = null;
      this.tx = null;
      this.vertexTableId = -1;
      this.edgeListTableId = -1;
    } else {
      this.client = new RAMCloud(coordLocator, "main", dpdkPort);
      this.tx = null;
      this.vertexTableId = client.createTable(graphName + "_vertexTable", numServers);
      this.edgeListTableId = client.createTable(graphName + "_edgeListTable", numServers);
      this.vertexTableOS = null;
      this.edgeListTableOS = null;
    }
  }

  /**
   * Close our connection to the database. 
   *
   * Need to free up anything that was allocated in native libraries.
   *
   * If we were generating RAMCloud image files, then we need to close those file handles.
   */
  public void close() {
    if (tx != null)
      abortTx();

    if (client != null)
      client.disconnect();

    try {
      if (vertexTableOS != null)
        vertexTableOS.close();

      if (edgeListTableOS != null)
        edgeListTableOS.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  /* **************************************************************************
   *
   * Transactions Control
   *
   * *************************************************************************/

  /** 
   * Start a new transaction. 
   *
   * If there is currently an open transaction, it is aborted. 
   */
  public void beginTx() {
    if (tx != null)
      abortTx();

    tx = new RAMCloudTransaction(client);
  }

  /** 
   * Attempt to commit the transaction, waiting for all participant servers to confirm that the
   * transaction has been applied.
   * 
   * @return True if the transaction committed successfully, false otherwise. 
   */
  public boolean commitAndSyncTx() {
    // An empty transaction always commit successfully.
    if (tx == null)
      return true;

    boolean committed = tx.commitAndSync();
    tx.close();
    tx = null;
    return committed;
  }

  /**
   * Attempt to commit the transaction, waiting only until a commit decision is reached, but not
   * waiting for the participant servers to confirm that the transaction has been applied. 
   *
   * Although this method may return true, it may still be the case that the transaction has not
   * yet been applied, and future reads may not see the writes of this transaction, thus violating
   * "read your own writes" semantics.
   *
   * @return True if the transaction committed successfully, false otherwise. 
   */
  public boolean commitTx() {
    // An empty transaction always commit successfully.
    if (tx == null)
      return true;

    boolean committed = tx.commit();
    tx.close();
    tx = null;
    return committed;
  }

  /**
   * Abort the current transaction. 
   */
  public void abortTx() {
    if (tx != null)
      tx.close();

    tx = null;
  }

  /* **************************************************************************
   *
   * Graph Traversal and Reading Operations.
   *
   * *************************************************************************/

  public TraversalResult traverse(
      Vertex v, 
      String eLabel, 
      Direction dir, 
      boolean fillEdge,
      String ... neighborLabels) {
    return traverse(Collections.singleton(v), eLabel, dir, fillEdge, neighborLabels);
  }

  public TraversalResult traverse(
      TraversalResult r, 
      String eLabel, 
      Direction dir, 
      boolean fillEdge,
      String ... neighborLabels) {
    return traverse(r.vSet, eLabel, dir, fillEdge, neighborLabels);
  }

  /** 
   * Traverses an edge type for a set of vertices.
   *
   * @param vCol Collection of vertices to start from.
   * @param eLabel Label of edge to traverse.
   * @param dir Direction of edge.
   * @param fillEdge Whether or not to fill in edge properties in the return result
   * @param nLabels Labels of neighbor vertices.
   *
   * @return TraversalResult describing the result of the traversal.
   */
  public TraversalResult traverse(
      Collection<Vertex> vCol,
      String eLabel, 
      Direction dir, 
      boolean fillEdge,
      String ... nLabels) {
    List<byte[]> keyPrefixes = GraphHelper.getEdgeListKeyPrefixes(vCol, eLabel, dir, nLabels);

    Map<byte[], List<SerializedEdge>> serEdgeLists;
    serEdgeLists = EdgeList.batchRead(tx, client, edgeListTableId, keyPrefixes);

    Map<Vertex, List<Vertex>> nbrListMap = new HashMap<>();

    Map<Vertex, List<Map<Object, Object>>> ePropListMap = null;
    if (fillEdge)
      ePropListMap = new HashMap<>();

    Map<UInt128, Vertex> nbrDedupMap = new HashMap<>();
    Set<Vertex> uniqNbrSet = new HashSet<>();

    int i = 0;
    for (String nLabel : nLabels) {
      for (Vertex vertex : vCol) {
        byte[] keyPrefix = keyPrefixes.get(i);

        if (serEdgeLists.containsKey(keyPrefix)) {
          List<SerializedEdge> serEdgeList = serEdgeLists.get(keyPrefix);

          List<Vertex> nList;
          List<Map<Object, Object>> ePropList = null;
          if (nbrListMap.containsKey(vertex)) {
            nList = nbrListMap.get(vertex);
            if (fillEdge)
              ePropList = ePropListMap.get(vertex);
          } else {
            nList = new ArrayList<>(serEdgeList.size());
            nbrListMap.put(vertex, nList);
            if (fillEdge) {
              ePropList = new ArrayList<>(serEdgeList.size());
              ePropListMap.put(vertex, ePropList);
            }
          }

          for (SerializedEdge serEdge : serEdgeList) {
            if (nbrDedupMap.containsKey(serEdge.vertexId)) {
              nList.add(nbrDedupMap.get(serEdge.vertexId));
            } else {
              Vertex v = new Vertex(serEdge.vertexId, nLabel);
              nList.add(v);
              nbrDedupMap.put(serEdge.vertexId, v);
              uniqNbrSet.add(v);
            }

            if (fillEdge)
              ePropList.add((Map<Object, Object>)GraphHelper.deserializeObject(
                    serEdge.serializedProperties));
          }
        }

        i++;
      }
    }    

    return new TraversalResult(nbrListMap, ePropListMap, uniqNbrSet);
  }

  public void fillProperties(Vertex ... vertices) {
    List<Vertex> vList = new ArrayList<>(vertices.length);
    for (Vertex v : vertices)
      vList.add(v);
    fillProperties(vList);
  }

  public void fillProperties(TraversalResult ... rs) {
    if (rs.length == 1) {
      fillProperties(rs[0].vSet);
    } else {
      int totalSize = 0;
      for (TraversalResult r : rs)
        totalSize += r.vSet.size();

      List<Vertex> vList = new ArrayList<>(totalSize);
      for (TraversalResult r : rs)
        vList.addAll(r.vSet);

      fillProperties(vList);
    }
  }

  /**
   * Read in the properties of all the given vertices. If specific keys are specified, may perform
   * space saving optimizations to store only that key or keys.
   *
   * @param keys (Optional) set of keys to fetch.
   */
  public void fillProperties(Iterable<Vertex> vertices, String ... keys) {
    // Max number of reads to issue in a multiread / batch
    int DEFAULT_MAX_MULTIREAD_SIZE = 1 << 11; 

    Iterator<Vertex> it = vertices.iterator();
    LinkedList<Object> requestQ = new LinkedList<>();
    LinkedList<Vertex> vertexQ = new LinkedList<>();
    while (it.hasNext()) {
      Vertex v = it.next();
      if (tx != null) {
        requestQ.addLast(new RAMCloudTransactionReadOp(tx, vertexTableId, 
              GraphHelper.getVertexPropertiesKey(v.id()), true));
      } else {
        requestQ.addLast(new MultiReadObject(vertexTableId, 
              GraphHelper.getVertexPropertiesKey(v.id())));
      }

      vertexQ.addLast(v);

      // If we've reached the multiread size limit or we're at the end, then issue the reads.
      if (requestQ.size() == DEFAULT_MAX_MULTIREAD_SIZE || !it.hasNext()) {
        if (tx != null) {
          while(requestQ.size() > 0) {
            RAMCloudTransactionReadOp readOp = (RAMCloudTransactionReadOp)requestQ.removeFirst();
            v = vertexQ.removeFirst();

            RAMCloudObject obj;
            try {
              obj = readOp.getValue();
              if (obj == null) {
                // This vertex has no properties set.
                v.setProperties(new HashMap<>());
                continue;
              }
            } catch (ClientException e) {
              throw new RuntimeException(e);
            } finally {
              readOp.close();
            }

            Map<Object, Object> properties = (Map<Object, Object>)
              GraphHelper.deserializeObject(obj.getValueBytes());
            if (keys.length == 1) {
              Map<Object, Object> minimap = new ArrayMap<>(1);
              minimap.put(keys[0], properties.get(keys[0]));
              properties = minimap;
            }

            v.setProperties(properties);
          }
        } else {
          MultiReadObject[] requests = new MultiReadObject[requestQ.size()];
          for (int i = 0; i < requests.length; i++) {
            requests[i] = (MultiReadObject)requestQ.removeFirst();
          }

          client.read(requests);
        
          for (int i = 0; i < requests.length; i++) {
            v = vertexQ.removeFirst();

            if (requests[i].getStatus() != Status.STATUS_OK) {
              if (requests[i].getStatus() == Status.STATUS_OBJECT_DOESNT_EXIST) {
                // This vertex has no properties set.
                v.setProperties(new HashMap<>());
                continue;
              } else {
                throw new RuntimeException("Vertex properties RAMCloud object had status " + 
                    requests[i].getStatus());
              }
            }

            Map<Object, Object> properties = (Map<Object, Object>)
              GraphHelper.deserializeObject(requests[i].getValueBytes());
            if (keys.length == 1) {
              Map<Object, Object> minimap = new ArrayMap<>(1);
              minimap.put(keys[0], properties.get(keys[0]));
              properties = minimap;
            }

            v.setProperties(properties);
          } 
        }
      } 
    }
  }

  /* **************************************************************************
   *
   * Graph Update Operations.
   *
   * *************************************************************************/

  /**
   * Add a vertex into the graph, updating the vertex if it exists, and creating the vertex if it
   * does not exist in the graph.
   *
   * @param v Vertex to add (create or update).
   */
  public void addVertex(Vertex v) {
    if (v.getProperties() == null)
      v.setProperties(new HashMap<>(0));

    byte[] serializedProps = GraphHelper.serializeObject(v.getProperties());
    byte[] propKeyByteArray = GraphHelper.getVertexPropertiesKey(v.id());

    if (tx != null)
      tx.write(vertexTableId, propKeyByteArray, serializedProps);
    else
      client.write(vertexTableId, propKeyByteArray, serializedProps, null);
  }

  /**
   * Add an edge in the graph from the source vertex to the destination vertex.
   *
   * Notes: 
   *   - Multiple edges of the same label between the same vertices is allowed.
   *   - Vertex labels must be set, as they are used to index edge lists.
   *
   * @param src Source vertex.
   * @param edgeLabel Label for the edge.
   * @param dst Destination vertex.
   * @param props Properties for the edge.
   */
  public void addEdge(Vertex src, String edgeLabel, Vertex dst, Map<Object, Object> props) {
    byte[] serializedProps;
    if (props == null || props.size() == 0)
      serializedProps = new byte[0];
    else
      serializedProps = GraphHelper.serializeObject(props);

    /* Add one vertex to the other's edge list, and vice versa, choosing one as the base and the
     * other as the neighbor. */
    for (int i = 0; i < 2; ++i) {
      Vertex baseVertex;
      Vertex neighborVertex;
      Direction direction;
      if (i == 0) {
        baseVertex = src;
        neighborVertex = dst;
        direction = Direction.OUT;
      } else {
        baseVertex = dst;
        neighborVertex = src;
        direction = Direction.IN;
      }

      byte[] keyPrefix = GraphHelper.getEdgeListKeyPrefix(baseVertex.id(), edgeLabel, direction,
              neighborVertex.label());

      EdgeList.prepend(tx, client, edgeListTableId, keyPrefix, neighborVertex.id(), serializedProps);
    }
  }

  /* **************************************************************************
   *
   * RAMCloud Image Generation Methods
   *
   * *************************************************************************/

  /**
   * Writes a vertex to the vertex table RAMCloud image file.
   * 
   * @param v Vertex to write to RAMCloud image file.
   */
  public void loadVertex(final Vertex v) {
    byte[] key = GraphHelper.getVertexPropertiesKey(v.id());
    byte[] value = GraphHelper.serializeObject(v.getProperties());

    ByteBuffer buffer = ByteBuffer.allocate(
        Integer.BYTES +
        key.length +
        Integer.BYTES +
        value.length)
        .order(ByteOrder.LITTLE_ENDIAN);

    buffer.putInt(key.length);
    buffer.put(key);
    buffer.putInt(value.length);
    buffer.put(value);

    try {
      vertexTableOS.write(buffer.array());
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Writes edge lists to the edge list table RAMCloud image file.
   *
   * @param baseVertexId ID of the vertex that owns this edge list.
   * @param edgeLabel Edge label for the edges in this list.
   * @param direction Direction of the edges in this list relative to base.
   * @param neighborIds IDs of the neighbors in the edge list.
   * @param propMaps Properties on edges in the list.
   */
  public void loadEdges(
      final UInt128 baseVertexId, 
      final String edgeLabel,
      final Direction direction, 
      final String neighborLabel, 
      final List<UInt128> neighborIds, 
      final List<Map<Object, Object>> propMaps) {
    byte[] keyPrefix = GraphHelper.getEdgeListKeyPrefix(baseVertexId, edgeLabel, direction,
            neighborLabel);

    List<byte[]> serializedPropList = new ArrayList<>(propMaps.size());
    for (int i = 0; i < propMaps.size(); i++)
      serializedPropList.add(GraphHelper.serializeObject(propMaps.get(i)));

    EdgeList.writeListToFile(edgeListTableOS, keyPrefix, neighborIds, serializedPropList);
  }
}
