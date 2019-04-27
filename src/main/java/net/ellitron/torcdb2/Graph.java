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
  private final RAMCloud client;
  private RAMCloudTransaction tx;
  private final long vertexTableId, edgeListTableId;

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

    if (config.containsKey("coordLocator"))
      this.coordLocator = config.get("coordLocator");
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

    this.client = new RAMCloud(coordLocator, "main", dpdkPort);
    this.tx = null;
    this.vertexTableId = client.createTable(graphName + "_vertexTable", numServers);
    this.edgeListTableId = client.createTable(graphName + "_edgeListTable", numServers);
  }

  /**
   * Close our connection to the database. 
   *
   * Need to free up anything that was allocated in native libraries.
   */
  public void close() {
    abortTx();
    client.disconnect();
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
   * Traversal
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
    if (tx != null) {
      serEdgeLists = EdgeList.batchRead(tx, edgeListTableId, keyPrefixes);
    } else {
      serEdgeLists = EdgeList.batchRead(client, edgeListTableId, keyPrefixes);
    }

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

  public void fillProperties(Vertex v) {
    fillProperties(Collections.singletonList(v));
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
   * Read in the properties of all the given vertices.
   */
  public void fillProperties(Iterable<Vertex> vertices, String ... keys) {
    // Max number of reads to issue in a multiread / batch
    int DEFAULT_MAX_MULTIREAD_SIZE = 1 << 11; 

    if (tx != null) {
//      RAMCloudTransaction rctx = torcGraphTx.getThreadLocalRAMCloudTx();
//
//      // Keeps track of where we are in the vList.
//      int vListMarker = 0;
//
//      while (vListMarker < vList.size()) {
//        // Queue up a batchSize of asynchronous ReadOps.
//        int batchSize = Math.min(vList.size() - vListMarker, 
//            DEFAULT_MAX_MULTIREAD_SIZE);
//
//        RAMCloudTransactionReadOp[] readOps = 
//          new RAMCloudTransactionReadOp[batchSize];
//        for (int i = 0; i < batchSize; i++) {
//          readOps[i] = new RAMCloudTransactionReadOp(rctx, vertexTableId,
//              GraphHelper.getVertexPropertiesKey(vList.get(vListMarker + i).id()), 
//              true);
//        }
//
//        for (int i = 0; i < batchSize; i++) {
//          Vertex v = vList.get(vListMarker + i);
//
//          RAMCloudObject obj;
//          try {
//            obj = readOps[i].getValue();
//            if (obj == null) {
//              // This vertex has no properties set.
//              v.setProperties(new HashMap<>());
//              continue;
//            }
//          } catch (ClientException e) {
//            throw new RuntimeException(e);
//          } finally {
//            readOps[i].close();
//          }
//
//          Map<Object, Object> properties = 
//            (Map<Object, Object>)GraphHelper.deserializeObject(obj.getValueBytes());
//          v.setProperties(properties);
//        }
//
//        vListMarker += batchSize;
//      }
    } else {
      Iterator<Vertex> it = vertices.iterator();
      LinkedList<MultiReadObject> requestQ = new LinkedList<>();
      LinkedList<Vertex> vertexQ = new LinkedList<>();
      while (it.hasNext()) {
        Vertex v = it.next();
        requestQ.addLast(new MultiReadObject(vertexTableId, 
              GraphHelper.getVertexPropertiesKey(v.id())));
        vertexQ.addLast(v);
        if (requestQ.size() == DEFAULT_MAX_MULTIREAD_SIZE) {
          MultiReadObject[] requests = new MultiReadObject[requestQ.size()];
          for (int i = 0; i < requests.length; i++) {
            requests[i] = requestQ.removeFirst();
          }

          client.read(requests);
        
          for (int i = 0; i < requests.length; i++) {
            v = vertexQ.removeFirst();

            if (requests[i].getStatus() != Status.STATUS_OK) {
              if (requests[i].getStatus() == Status.STATUS_OBJECT_DOESNT_EXIST) {
                // This vertex has no properties set.
                continue;
              } else {
                throw new RuntimeException(
                    "Vertex properties RAMCloud object had status " + 
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

      if (requestQ.size() > 0) {
        MultiReadObject[] requests = new MultiReadObject[requestQ.size()];
        for (int i = 0; i < requests.length; i++) {
          requests[i] = requestQ.removeFirst();
        }

        client.read(requests);
      
        for (int i = 0; i < requests.length; i++) {
          Vertex v = vertexQ.removeFirst();

          if (requests[i].getStatus() != Status.STATUS_OK) {
            if (requests[i].getStatus() == Status.STATUS_OBJECT_DOESNT_EXIST) {
              // This vertex has no properties set.
              continue;
            } else {
              throw new RuntimeException(
                  "Vertex properties RAMCloud object had status " + 
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
