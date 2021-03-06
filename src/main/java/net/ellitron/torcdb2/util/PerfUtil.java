/* 
 * Copyright (C) 2019 Stanford University
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package net.ellitron.torcdb2.util;

import net.ellitron.torcdb2.*;

import org.docopt.Docopt;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.lang.Math;
import java.lang.reflect.InvocationTargetException;
import java.lang.Thread;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import java.text.SimpleDateFormat;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Random;
import java.util.Set;

/**
 * A utility for running microbenchmarks on TorcDB2.
 *
 * @author Jonathan Ellithorpe (jde@cs.stanford.edu)
 */
public class PerfUtil {
  private static final String doc =
      "PerfUtil: A utility for running microbenchmarks on TorcDB2"
      + "\n"
      + "Usage:\n"
      + "  PerfUtil [options] edge_rdwr\n"
      + "  PerfUtil [options] addVertex\n"
      + "  PerfUtil [options] addEdge\n"
      + "  PerfUtil [options] getProperties\n"
      + "  PerfUtil [options] traverse_unique_neighbors\n"
      + "  PerfUtil [options] traverse_one_neighbor\n"
      + "  PerfUtil [options] hashmapuint128keys\n"
      + "  PerfUtil [options] traverse_hashmap_put\n"
      + "  PerfUtil (-h | --help)\n"
      + "  PerfUtil --version\n"
      + "\n"
      + "Options:\n"
      + "  --config=<file>      PerfUtil configuration file\n"
      + "                       [default: ./config/perfutil.properties].\n"
      + "  --script=<script>    File to use as command script. Commands will\n"
      + "                       be executed from this script as if the args\n"
      + "                       were supplied at the command line, one per\n"
      + "                       line.\n"
      + "  --output=<file>      Output file.\n"
      + "  -h --help            Show this screen.\n"
      + "  --version            Show version.\n"
      + "\n";

  /** 
   * Returns a range of numbers based on the given parameters.
   *
   * @param start - Start of the range
   * @param end - End of the range
   * @param points - Number of points in the range
   * @param mode - Progression type of the range (linear or geometric)
   *
   * @return List of numbers representing the range.
   */
  public static List<Long> range(long start, long end, long points, String mode) {
    List<Long> range = new ArrayList<>();

    if (points > 1)
      if (mode.equals("linear"))
        for ( long i = start; 
              i <= end; 
              i += (end - start) / (points - 1) )
          range.add(i);    
      else if (mode.equals("geometric")) {
        double c = Math.pow(10, Math.log10((double)end/(double)start) 
                                            / (double)(points - 1));
        for (long i = start; i <= end; i = (long)Math.ceil(c * i))
          range.add(i);
      } else
        throw new RuntimeException(
            "Error: Unrecognized value for parameter traverse.nodes.mode");
    else 
      range.add(start);

    return range;
  }

  public static void main(String[] args) throws Exception {
    if (args.length == 1)
     args = args[0].split("\\s+");

    Map<String, Object> opts = new Docopt(doc).withVersion("PerfUtil 1.0").parse(args);

    System.out.println(opts);

    // Load properties from the configuration file.
    String configFilename = (String) opts.get("--config");
    Properties prop = new Properties();
    prop.load(new FileInputStream(configFilename));

    Map<String, String> config = new HashMap<>();
    prop.stringPropertyNames().stream()
        .forEach((propName) -> {
          config.put(propName, prop.getProperty(propName));
        });

    String rccoordlocator = prop.getProperty("rccoordlocator");

    BufferedReader scriptFile = null;
    if (opts.get("--script") != null) {
      String scriptFilename = (String) opts.get("--script");

      Path path = Paths.get(scriptFilename);
      scriptFile = Files.newBufferedReader(path, StandardCharsets.UTF_8);
    }

    boolean done = false;
    while (!done) {
      if (scriptFile != null) {
        String line;
        while (((line = scriptFile.readLine()) != null) && line.startsWith("#")) {
          // Skip over commented out lines
          continue;
        }

        if (line != null) {
          args = line.split(" (?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)");
          for (int i = 0; i < args.length; i++)
            args[i] = args[i].replaceAll("^\"|\"$", "");

          Map<String, Object> scriptOpts =
              new Docopt(doc).withVersion("QueryTester 1.0").parse(args);
          scriptOpts.put("cmdstring", line);
          opts = scriptOpts;
        } else {
          scriptFile.close();
          break;
        }
      } else {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < args.length; i++) {
          sb.append(args[i]);
          if (i < args.length - 1)
            sb.append(" ");
        }
        opts.put("cmdstring", sb.toString());
        done = true;
      }

      String cmdstring = (String) opts.get("cmdstring");

      if ((Boolean) opts.get("edge_rdwr")) {
        Graph graph = new Graph(config);
        String edgeLabel = "knows";
        Direction direction = Direction.OUT;
        int sizes[] = {8192,4096,2048,1024+512,1024,512+256,512,256+128,256,128+64,128,64+32,64};

        // WarmUp
        {
          Vertex baseVertex = new Vertex(new UInt128(0,0), "Warmup");
          Vertex neigVertex = new Vertex(new UInt128(0,0), "Warmup");
          byte[] keyPrefix = GraphHelper.getEdgeListKeyPrefix(baseVertex.id(), edgeLabel, direction,
                  neigVertex.label());

          int numElements = 10000;
          for (int j = 0; j < numElements; j++) {
            EdgeList.prepend(null, graph.getClient(), graph.getEdgeListTableId(), keyPrefix, 
                neigVertex.id(), new byte[0], 1024, 0);
          }

          for (int j = 0; j < 10000; j++) {
            EdgeList.batchReadSingleThreaded(null, graph.getClient(), graph.getEdgeListTableId(),
                Collections.singleton(baseVertex), edgeLabel, direction, false, "Warmup");
          }
        }

        System.out.print(String.format("numElements,segSize,prepend50th,read50th\n"));
        for (int i = 0; i < sizes.length; i++) {
          Vertex baseVertex = new Vertex(new UInt128(0,i), "Person");
          Vertex neigVertex = new Vertex(new UInt128(0,0), "Person");
          byte[] keyPrefix = GraphHelper.getEdgeListKeyPrefix(baseVertex.id(), edgeLabel, direction,
                  neigVertex.label());

          int size = sizes[i];
          int numElements = 1000;
          Long[] writeLatency = new Long[numElements];
          for (int j = 0; j < numElements; j++) {
            long startTime = System.nanoTime();
            EdgeList.prepend(null, graph.getClient(), graph.getEdgeListTableId(), keyPrefix, 
                neigVertex.id(), new byte[0], size, 0);
            writeLatency[j] = System.nanoTime() - startTime;
          }

          Arrays.sort(writeLatency);

          Long[] readLatency = new Long[1000];
          for (int j = 0; j < 1000; j++) {
            long startTime = System.nanoTime();
            EdgeList.batchReadSingleThreaded(null, graph.getClient(), graph.getEdgeListTableId(),
                Collections.singleton(baseVertex), edgeLabel, direction, false, "Person");
            readLatency[j] = System.nanoTime() - startTime;
          }

          Arrays.sort(readLatency);

          long sum = 0;
          long min = Long.MAX_VALUE;
          long max = 0;
          for (int k = 0; k < readLatency.length; k++) {
            sum += readLatency[k];

            if (readLatency[k] < min) {
              min = readLatency[k];
            }

            if (readLatency[k] > max) {
              max = readLatency[k];
            }
          }

          long mean = sum / readLatency.length;

          int p25 = (int) (0.25 * (float) readLatency.length);
          int p50 = (int) (0.50 * (float) readLatency.length);
          int p75 = (int) (0.75 * (float) readLatency.length);
          int p90 = (int) (0.90 * (float) readLatency.length);
          int p95 = (int) (0.95 * (float) readLatency.length);
          int p99 = (int) (0.99 * (float) readLatency.length);

          System.out.print(String.format("%d, %d, %d, %d\n", 
                numElements,
                size,
                writeLatency[(int) (0.50 * (float) writeLatency.length)]/1000, 
                readLatency[(int) (0.50 * (float) readLatency.length)]/1000));
        }
      } else if ((Boolean) opts.get("addVertex")) {
        Graph graph = new Graph(config);
       
        long props_start = Long.decode(config.get("addvertex.props.start"));
        long props_end = Long.decode(config.get("addvertex.props.end"));
        long props_points = Long.decode(config.get("addvertex.props.points"));
        String props_mode = config.get("addvertex.props.mode");

        long vertices_start = Long.decode(config.get("addvertex.vertices.start"));
        long vertices_end = Long.decode(config.get("addvertex.vertices.end"));
        long vertices_points = Long.decode(config.get("addvertex.vertices.points"));
        String vertices_mode = config.get("addvertex.vertices.mode");

        long warmup_maxtime = Long.decode(config.get("addvertex.warmup.maxtime"));

        long samples = Long.decode(config.get("addvertex.samples"));
       
        FileWriter outfile = null;
        if (opts.get("--output") !=  null) {
          try {
            outfile = new FileWriter((String)opts.get("--output"));
          } catch (Exception e) {
            throw new RuntimeException(e);
          }
        }

        List<Long> props_params = range(props_start, props_end, props_points, props_mode);
        List<Long> vertices_params = range(vertices_start, vertices_end, vertices_points, vertices_mode);

        // WarmUp
        long warmupStartTime = System.currentTimeMillis();
        System.out.println(String.format("Beginning warmup (max %d minutes)...", warmup_maxtime));
        {
          for (long numVertices : vertices_params) {
            for (long numProps : props_params) {
              // Setup the vertices and properties
              Vertex[] v = new Vertex[(int)numVertices];
              for (long i = 0; i < numVertices; i++) {
                v[(int)i] = new Vertex(new UInt128(0, i), "Person");
              }

              Map<Object, Object> props = new HashMap<>();
              for (long i = 0; i < numProps; i++) {
                String str = String.format("%010d", i);
                props.put(str, str);
              }

              // Execution times are recorded in nanoseconds.
              Long[] execTimes = new Long[(int)samples];
              for (long i = 0; i < samples; i++) {
                long startTime = System.nanoTime();
                graph.beginTx();
                for (long j = 0; j < numVertices; j++)
                  graph.addVertex(v[(int)j], props);
                graph.commitAndSyncTx();
                long endTime = System.nanoTime();
                execTimes[(int)i] = endTime - startTime;
              }

              if ((System.currentTimeMillis() - warmupStartTime)/(1000*60) >= warmup_maxtime)
                break;
            }

            if ((System.currentTimeMillis() - warmupStartTime)/(1000*60) >= warmup_maxtime)
              break;
          }
        }

        System.out.println(String.format(
              "Warmup phase complete. Total time: %d seconds",
              (System.currentTimeMillis() - warmupStartTime)/1000));

        // Output the header
        // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
        if (outfile == null)
          System.out.println(
              "props, vertices, min, mean, max, samples, 50th, 90th, 95th, 99th, 99.9th");
        else {
          outfile.append(
              "props, vertices, min, mean, max, samples, 50th, 90th, 95th, 99th, 99.9th\n");
          outfile.flush();
        }

        for (long numVertices : vertices_params) {
          for (long numProps : props_params) {
            // Setup the vertices and properties
            Vertex[] v = new Vertex[(int)numVertices];
            for (long i = 0; i < numVertices; i++) {
              v[(int)i] = new Vertex(new UInt128(0, i), "Person");
            }

            Map<Object, Object> props = new HashMap<>();
            for (long i = 0; i < numProps; i++) {
              String str = String.format("%010d", i);
              props.put(str, str);
            }

            // Execution times are recorded in nanoseconds.
            Long[] execTimes = new Long[(int)samples];
            for (long i = 0; i < samples; i++) {
              long startTime = System.nanoTime();
              graph.beginTx();
              for (long j = 0; j < numVertices; j++)
                graph.addVertex(v[(int)j], props);
              graph.commitAndSyncTx();
              long endTime = System.nanoTime();
              execTimes[(int)i] = endTime - startTime;
            }

            // Calculate latency statistics
            Arrays.sort(execTimes);

            long sum = 0;
            long min = Long.MAX_VALUE;
            long max = 0;
            for (int k = 0; k < execTimes.length; k++) {
              sum += execTimes[k];

              if (execTimes[k] < min)
                min = execTimes[k];

              if (execTimes[k] > max)
                max = execTimes[k];
            }

            long mean = sum / execTimes.length;

            long p25  = (int) (0.250 * (float) execTimes.length);
            long p50  = (int) (0.500 * (float) execTimes.length);
            long p75  = (int) (0.750 * (float) execTimes.length);
            long p90  = (int) (0.900 * (float) execTimes.length);
            long p95  = (int) (0.950 * (float) execTimes.length);
            long p99  = (int) (0.990 * (float) execTimes.length);
            long p999 = (int) (0.999 * (float) execTimes.length);

            // Print out the results (in microseconds)
            // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
            if (outfile == null)
              System.out.println(String.format("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d", 
                                                numProps,
                                                numVertices,
                                                min/1000,
                                                mean/1000,
                                                max/1000,
                                                samples,
                                                execTimes[(int)p50]/1000,
                                                execTimes[(int)p90]/1000,
                                                execTimes[(int)p95]/1000,
                                                execTimes[(int)p99]/1000,
                                                execTimes[(int)p999]/1000));
            else {
              System.out.println(String.format(
                    "\t{NumProps: %d, NumVertices: %d, Time: %d}",
                    numProps,
                    numVertices,
                    min/1000));

              outfile.append(String.format("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", 
                                            numProps,
                                            numVertices,
                                            min/1000,
                                            mean/1000,
                                            max/1000,
                                            samples,
                                            execTimes[(int)p50]/1000,
                                            execTimes[(int)p90]/1000,
                                            execTimes[(int)p95]/1000,
                                            execTimes[(int)p99]/1000,
                                            execTimes[(int)p999]/1000));
              outfile.flush();
            }
          }
        }

        if (outfile != null)
          outfile.close();

        graph.delete();
        graph.close();
      } else if ((Boolean) opts.get("addEdge")) {
        Graph graph = new Graph(config);
       
        long props_start = Long.decode(config.get("addedge.props.start"));
        long props_end = Long.decode(config.get("addedge.props.end"));
        long props_points = Long.decode(config.get("addedge.props.points"));
        String props_mode = config.get("addedge.props.mode");

        long edges = Long.decode(config.get("addedge.edges"));

        long samples = Long.decode(config.get("addedge.samples"));

        long warmup_maxtime = Long.decode(config.get("addedge.warmup.maxtime"));
       
        FileWriter outfile = null;
        if (opts.get("--output") !=  null) {
          try {
            outfile = new FileWriter((String)opts.get("--output"));
          } catch (Exception e) {
            throw new RuntimeException(e);
          }
        }

        List<Long> props_params = range(props_start, props_end, props_points, props_mode);

        // WarmUp
        long warmupStartTime = System.currentTimeMillis();
        System.out.println(String.format("Beginning warmup (max %d minutes)...", warmup_maxtime));
        {
          for (long numProps = 0; numProps < 16; numProps++) {
            Map<Object, Object> props = new HashMap<>();
            for (long i = 0; i < numProps; i++) {
              String str = String.format("%010d", i);
              props.put(str, str);
            }

            Vertex v1 = new Vertex(new UInt128(1, numProps), "Warmup");
            Vertex v2 = new Vertex(new UInt128(1, 1), "Warmup");
            // Execution times are recorded in nanoseconds.
            Long[] execTimes = new Long[(int)edges];
            for (long i = 0; i < edges; i++) {
              long startTime = System.nanoTime();
              graph.beginTx();
              graph.addEdge(v1, "posts", v2, props);
              graph.commitAndSyncTx();
              long endTime = System.nanoTime();
              execTimes[(int)i] = endTime - startTime;
            }
          }
        }

        System.out.println(String.format(
              "Warmup phase complete. Total time: %d seconds",
              (System.currentTimeMillis() - warmupStartTime)/1000));

        // Output the header
        // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
        if (outfile == null)
          System.out.println(
              "props, edge, time");
        else {
          outfile.append(
              "props, edge, time\n");
          outfile.flush();
        }

        for (long numProps : props_params) {
          Map<Object, Object> props = new HashMap<>();
          for (long i = 0; i < numProps; i++) {
            String str = String.format("%010d", i);
            props.put(str, str);
          }

          // Execution times are recorded in nanoseconds.
          Long[][] execTimes = new Long[(int)samples][(int)edges];
          for (long j = 0; j < samples; j++) {
            Vertex v1 = new Vertex(new UInt128(j, numProps), "Person");
            Vertex v2 = new Vertex(new UInt128(j, numProps), "Person");
            for (long i = 0; i < edges; i++) {
              long startTime = System.nanoTime();
              graph.beginTx();
              graph.addEdge(v1, "knows", v2, props);
              graph.commitAndSyncTx();
              long endTime = System.nanoTime();
              execTimes[(int)j][(int)i] = endTime - startTime;
              //TraversalResult result = graph.traverse(v1, "knows", Direction.OUT, false, "Person");
              //System.out.println(String.format("%d", result.vMap.get(v1).size()));
            }

            if (outfile != null)
              System.out.println(String.format(
                    "\t{NumProps: %d, Sample: %d}",
                    numProps,
                    j));
          }

          for (long i = 0; i < edges; i++) {
            Long[] sampleArray = new Long[(int)samples];
            for (long j = 0; j < samples; j++)
              sampleArray[(int)j] = execTimes[(int)j][(int)i];

            Arrays.sort(sampleArray);

            long sum = 0;
            long min = Long.MAX_VALUE;
            long max = 0;
            for (int k = 0; k < sampleArray.length; k++) {
              sum += sampleArray[k];

              if (sampleArray[k] < min)
                min = sampleArray[k];

              if (sampleArray[k] > max)
                max = sampleArray[k];
            }

            long mean = sum / sampleArray.length;

            long p25  = (int) (0.250 * (float) sampleArray.length);
            long p50  = (int) (0.500 * (float) sampleArray.length);
            long p75  = (int) (0.750 * (float) sampleArray.length);
            long p90  = (int) (0.900 * (float) sampleArray.length);
            long p95  = (int) (0.950 * (float) sampleArray.length);
            long p99  = (int) (0.990 * (float) sampleArray.length);
            long p999 = (int) (0.999 * (float) sampleArray.length);

            // Print out the results (in microseconds)
            // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
            if (outfile == null)
              System.out.println(String.format("%d, %d, %d", 
                                                numProps,
                                                i+1,
                                                sampleArray[(int)p50]/1000));
            else {
              outfile.append(String.format("%d, %d, %d\n", 
                                            numProps,
                                            i+1,
                                            sampleArray[(int)p50]/1000));
              outfile.flush();
            }
          }

          if (outfile != null)
            System.out.println(String.format(
                  "\t{NumProps: %d}",
                  numProps));
        }

        if (outfile != null)
          outfile.close();

        graph.delete();
        graph.close();
      } else if ((Boolean) opts.get("getProperties")) {
        Graph graph = new Graph(config);
       
        long props_start = Long.decode(config.get("getproperties.props.start"));
        long props_end = Long.decode(config.get("getproperties.props.end"));
        long props_points = Long.decode(config.get("getproperties.props.points"));
        String props_mode = config.get("getproperties.props.mode");

        long vertices_start = Long.decode(config.get("getproperties.vertices.start"));
        long vertices_end = Long.decode(config.get("getproperties.vertices.end"));
        long vertices_points = Long.decode(config.get("getproperties.vertices.points"));
        String vertices_mode = config.get("getproperties.vertices.mode");

        long warmup_maxtime = Long.decode(config.get("getproperties.warmup.maxtime"));

        long samples = Long.decode(config.get("getproperties.samples"));
       
        FileWriter outfile = null;
        if (opts.get("--output") !=  null) {
          try {
            outfile = new FileWriter((String)opts.get("--output"));
          } catch (Exception e) {
            throw new RuntimeException(e);
          }
        }

        List<Long> props_params = range(props_start, props_end, props_points, props_mode);
        List<Long> vertices_params = range(vertices_start, vertices_end, vertices_points, vertices_mode);

        // WarmUp
        long warmupStartTime = System.currentTimeMillis();
        System.out.println(String.format("Beginning warmup (max %d minutes)...", warmup_maxtime));
        {
          for (long numProps : props_params) {
            Map<Object, Object> props = new HashMap<>();
            for (long i = 0; i < numProps; i++) {
              String str = String.format("%010d", i);
              props.put(str, str);
            }

            long numVertices_prev = 0;
            Set<Vertex> vSet = new HashSet<>();
            for (long numVertices : vertices_params) {
              // Setup the vertices and properties
              for (long i = numVertices_prev; i < numVertices; i++) {
                Vertex v = new Vertex(new UInt128(1, i), "Person");
                graph.addVertex(v, props);
                vSet.add(v);
              }

              // Execution times are recorded in nanoseconds.
              Long[] execTimes = new Long[(int)samples];
              for (long i = 0; i < samples; i++) {
                long startTime = System.nanoTime();
                graph.beginTx();
                Map<Vertex, Map<Object, Object>> vPropMap = new HashMap<>();
                graph.getProperties(vPropMap, vSet);
                graph.commitAndSyncTx();
                long endTime = System.nanoTime();
                execTimes[(int)i] = endTime - startTime;
              }

              if ((System.currentTimeMillis() - warmupStartTime)/(1000*60) >= warmup_maxtime)
                break;

              numVertices_prev = numVertices;
            }

            if ((System.currentTimeMillis() - warmupStartTime)/(1000*60) >= warmup_maxtime)
              break;
          }
        }

        System.out.println(String.format(
              "Warmup phase complete. Total time: %d seconds",
              (System.currentTimeMillis() - warmupStartTime)/1000));

        // Output the header
        // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
        if (outfile == null)
          System.out.println(
              "props, vertices, min, mean, max, samples, 50th, 90th, 95th, 99th, 99.9th");
        else {
          outfile.append(
              "props, vertices, min, mean, max, samples, 50th, 90th, 95th, 99th, 99.9th\n");
          outfile.flush();
        }

        for (long numProps : props_params) {
          Map<Object, Object> props = new HashMap<>();
          for (long i = 0; i < numProps; i++) {
            String str = String.format("%010d", i);
            props.put(str, str);
          }

          long numVertices_prev = 0;
          Set<Vertex> vSet = new HashSet<>();
          for (long numVertices : vertices_params) {
            // Setup the vertices and properties
            for (long i = numVertices_prev; i < numVertices; i++) {
              Vertex v = new Vertex(new UInt128(0, i), "Person");
              graph.addVertex(v, props);
              vSet.add(v);
            }

            // Execution times are recorded in nanoseconds.
            Long[] execTimes = new Long[(int)samples];
            for (long i = 0; i < samples; i++) {
              long startTime = System.nanoTime();
              graph.beginTx();
              Map<Vertex, Map<Object, Object>> vPropMap = new HashMap<>();
              graph.getProperties(vPropMap, vSet);
              graph.commitAndSyncTx();
              long endTime = System.nanoTime();
              execTimes[(int)i] = endTime - startTime;
            }

            // Calculate latency statistics
            Arrays.sort(execTimes);

            long sum = 0;
            long min = Long.MAX_VALUE;
            long max = 0;
            for (int k = 0; k < execTimes.length; k++) {
              sum += execTimes[k];

              if (execTimes[k] < min)
                min = execTimes[k];

              if (execTimes[k] > max)
                max = execTimes[k];
            }

            long mean = sum / execTimes.length;

            long p25  = (int) (0.250 * (float) execTimes.length);
            long p50  = (int) (0.500 * (float) execTimes.length);
            long p75  = (int) (0.750 * (float) execTimes.length);
            long p90  = (int) (0.900 * (float) execTimes.length);
            long p95  = (int) (0.950 * (float) execTimes.length);
            long p99  = (int) (0.990 * (float) execTimes.length);
            long p999 = (int) (0.999 * (float) execTimes.length);

            // Print out the results (in microseconds)
            // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
            if (outfile == null)
              System.out.println(String.format("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d", 
                                                numProps,
                                                numVertices,
                                                min/1000,
                                                mean/1000,
                                                max/1000,
                                                samples,
                                                execTimes[(int)p50]/1000,
                                                execTimes[(int)p90]/1000,
                                                execTimes[(int)p95]/1000,
                                                execTimes[(int)p99]/1000,
                                                execTimes[(int)p999]/1000));
            else {
              System.out.println(String.format(
                    "\t{NumProps: %d, NumVertices: %d, Time: %d}",
                    numProps,
                    numVertices,
                    min/1000));

              outfile.append(String.format("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", 
                                            numProps,
                                            numVertices,
                                            min/1000,
                                            mean/1000,
                                            max/1000,
                                            samples,
                                            execTimes[(int)p50]/1000,
                                            execTimes[(int)p90]/1000,
                                            execTimes[(int)p95]/1000,
                                            execTimes[(int)p99]/1000,
                                            execTimes[(int)p999]/1000));
              outfile.flush();
            }

            numVertices_prev = numVertices;
          }
        }

        if (outfile != null)
          outfile.close();

        graph.delete();
        graph.close();
      } else if ((Boolean) opts.get("traverse_unique_neighbors")) {
        Graph graph = new Graph(config);

        long warmup_nodes_start = Long.decode(config.get("traverse_unique_neighbors.warmup.nodes.start"));
        long warmup_nodes_end = Long.decode(config.get("traverse_unique_neighbors.warmup.nodes.end"));
        long warmup_nodes_points = Long.decode(config.get("traverse_unique_neighbors.warmup.nodes.points"));
        String warmup_nodes_mode = config.get("traverse_unique_neighbors.warmup.nodes.mode");

        long warmup_degree_start = Long.decode(config.get("traverse_unique_neighbors.warmup.degree.start"));
        long warmup_degree_end = Long.decode(config.get("traverse_unique_neighbors.warmup.degree.end"));
        long warmup_degree_points = Long.decode(config.get("traverse_unique_neighbors.warmup.degree.points"));
        String warmup_degree_mode = config.get("traverse_unique_neighbors.warmup.degree.mode");

        long warmup_samples = Long.decode(config.get("traverse_unique_neighbors.warmup.samples"));

        long warmup_maxtime = Long.decode(config.get("traverse_unique_neighbors.warmup.maxtime"));
       
        long nodes_start = Long.decode(config.get("traverse_unique_neighbors.nodes.start"));
        long nodes_end = Long.decode(config.get("traverse_unique_neighbors.nodes.end"));
        long nodes_points = Long.decode(config.get("traverse_unique_neighbors.nodes.points"));
        String nodes_mode = config.get("traverse_unique_neighbors.nodes.mode");

        long degree_start = Long.decode(config.get("traverse_unique_neighbors.degree.start"));
        long degree_end = Long.decode(config.get("traverse_unique_neighbors.degree.end"));
        long degree_points = Long.decode(config.get("traverse_unique_neighbors.degree.points"));
        String degree_mode = config.get("traverse_unique_neighbors.degree.mode");

        long samples = Long.decode(config.get("traverse_unique_neighbors.samples"));
       
        FileWriter outfile = null;
        if (opts.get("--output") !=  null) {
          try {
            outfile = new FileWriter((String)opts.get("--output"));
          } catch (Exception e) {
            throw new RuntimeException(e);
          }
        }

        List<Long> warmup_nodes_params = range(warmup_nodes_start, warmup_nodes_end, 
            warmup_nodes_points, warmup_nodes_mode);
        List<Long> warmup_degree_params = range(warmup_degree_start, warmup_degree_end, 
            warmup_degree_points, warmup_degree_mode);

        List<Long> nodes_params = range(nodes_start, nodes_end, nodes_points, nodes_mode);
        List<Long> degree_params = range(degree_start, degree_end, degree_points, degree_mode);

        // WarmUp
        long warmupStartTime = System.currentTimeMillis();
        System.out.println(String.format("Beginning warmup (max %d minutes)...", warmup_maxtime));
        {
          String edgeLabel = "knows";
          Direction direction = Direction.OUT;
          Vertex baseVertex = new Vertex(new UInt128(1,0), "Warmup");
          Vertex neigVertex = new Vertex(new UInt128(1,0), "Warmup");
          byte[] keyPrefix = GraphHelper.getEdgeListKeyPrefix(baseVertex.id(), edgeLabel, direction,
                  neigVertex.label());

          int numElements = 10000;
          for (int j = 0; j < numElements; j++) {
            EdgeList.prepend(null, graph.getClient(), graph.getEdgeListTableId(), keyPrefix, 
                neigVertex.id(), new byte[0], 1024, 0);
          }

          for (int j = 0; j < 10000; j++) {
            EdgeList.batchReadSingleThreaded(null, graph.getClient(), graph.getEdgeListTableId(),
                Collections.singleton(baseVertex), edgeLabel, direction, false, "Warmup");
          }

          long degree_prev = 0;
          long neighborIdCounter = 0;
          for (long degree : warmup_degree_params) {
            // Experiment setup: For each degree we setup all the nodes to have that degree. Each time
            // degree is increased, we add the additional edges needed to the existing nodes to bring
            // the degree up to the specified amount.
            for (long i = 0; i < warmup_nodes_end; i++) {
              Vertex startNode = new Vertex(new UInt128(1,i), "Person");
              for (long j = degree_prev; j < degree; j++) {
                // Vertices need labels, need to create vertex objects
                graph.addEdge(startNode, 
                              "knows", 
                              new Vertex(new UInt128(1, neighborIdCounter), 
                                "Person"), 
                              null);
                neighborIdCounter++;
              }
            }

            long nodes_prev = 0;
            Set<Vertex> startNodeSet = new HashSet<>();
            for (long nodes : warmup_nodes_params) {
              for (long i = nodes_prev; i < nodes; i++)
                startNodeSet.add(new Vertex(new UInt128(1,i), "Person"));

              for (long i = 0; i < warmup_samples; i++) {
                TraversalResult result = graph.traverse(startNodeSet, 
                                                        "knows", 
                                                        Direction.OUT, 
                                                        false, 
                                                        "Person");
              }

              if ((System.currentTimeMillis() - warmupStartTime)/(1000*60) >= warmup_maxtime)
                break;

              nodes_prev = nodes;
            }

            if ((System.currentTimeMillis() - warmupStartTime)/(1000*60) >= warmup_maxtime)
              break;

            degree_prev = degree;
          }
        }

        System.out.println(String.format(
              "Warmup phase complete. Total time: %d seconds",
              (System.currentTimeMillis() - warmupStartTime)/1000));

        // Output the header
        // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
        if (outfile == null)
          System.out.println(
              "degree, nodes, min, mean, max, samples, 50th, 90th, 95th, 99th, 99.9th");
        else {
          outfile.append(
              "degree, nodes, min, mean, max, samples, 50th, 90th, 95th, 99th, 99.9th\n");
          outfile.flush();
        }

        long degree_prev = 0;
//        long neighborIdCounter = 0;
        Random random = new Random(0);
        for (long degree : degree_params) {
          // Experiment setup: For each degree we setup all the nodes to have that degree. Each time
          // degree is increased, we add the additional edges needed to the existing nodes to bring
          // the degree up to the specified amount.
          for (long i = 0; i < nodes_end; i++) {
            Vertex startNode = new Vertex(new UInt128(0,i), "Person");
            for (long j = degree_prev; j < degree; j++) {
              // Vertices need labels, need to create vertex objects
              graph.addEdge(startNode, 
                            "knows", 
                            new Vertex(new UInt128(0, random.nextLong()), "Person"), 
                            null);
//              neighborIdCounter++;
            }
          }

          long nodes_prev = 0;
          Set<Vertex> startNodeSet = new HashSet<>();
          for (long nodes : nodes_params) {
            for (long i = nodes_prev; i < nodes; i++)
              startNodeSet.add(new Vertex(new UInt128(0,i), "Person"));

            // Execution times are recorded in nanoseconds.
            Long[] execTimes = new Long[(int)samples];
            long numNeighbors = 0;
            for (long i = 0; i < samples; i++) {
              long startTime = System.nanoTime();
              TraversalResult result = graph.traverse(startNodeSet, 
                                                      "knows", 
                                                      Direction.OUT, 
                                                      false, 
                                                      "Person");
              long endTime = System.nanoTime();
              execTimes[(int)i] = endTime - startTime;
              numNeighbors = result.vSet.size();
            }

            // Calculate latency statistics
            Arrays.sort(execTimes);

            long sum = 0;
            long min = Long.MAX_VALUE;
            long max = 0;
            for (int k = 0; k < execTimes.length; k++) {
              sum += execTimes[k];

              if (execTimes[k] < min)
                min = execTimes[k];

              if (execTimes[k] > max)
                max = execTimes[k];
            }

            long mean = sum / execTimes.length;

            long p25  = (int) (0.250 * (float) execTimes.length);
            long p50  = (int) (0.500 * (float) execTimes.length);
            long p75  = (int) (0.750 * (float) execTimes.length);
            long p90  = (int) (0.900 * (float) execTimes.length);
            long p95  = (int) (0.950 * (float) execTimes.length);
            long p99  = (int) (0.990 * (float) execTimes.length);
            long p999 = (int) (0.999 * (float) execTimes.length);

            // Print out the results (in microseconds)
            // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
            if (outfile == null)
              System.out.println(String.format("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d", 
                                                degree,
                                                nodes,
                                                min/1000,
                                                mean/1000,
                                                max/1000,
                                                samples,
                                                execTimes[(int)p50]/1000,
                                                execTimes[(int)p90]/1000,
                                                execTimes[(int)p95]/1000,
                                                execTimes[(int)p99]/1000,
                                                execTimes[(int)p999]/1000));
            else {
              System.out.println(String.format(
                    "\t{Degree: %d, Nodes: %d, Nbrs: %d, Time: %d, EPS: %01.2f}",
                    degree,
                    nodes,
                    numNeighbors,
                    min/1000,
                    1000.0 * (double)(degree * nodes) / (double)min));

              outfile.append(String.format("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", 
                                                degree,
                                                nodes,
                                                min/1000,
                                                mean/1000,
                                                max/1000,
                                                samples,
                                                execTimes[(int)p50]/1000,
                                                execTimes[(int)p90]/1000,
                                                execTimes[(int)p95]/1000,
                                                execTimes[(int)p99]/1000,
                                                execTimes[(int)p999]/1000));
              outfile.flush();
            }
              
            nodes_prev = nodes;
          }

          degree_prev = degree;
        }

        if (outfile != null)
          outfile.close();

        graph.delete();
        graph.close();
      } else if ((Boolean) opts.get("hashmapuint128keys")) {
        /*
         * A test that measures the performance of hashmap operations when used with UInt128 objects
         * as elements. The reason this test exists in the first place is because I found that the
         * performance of get/put operations was affected by the values of UInt128 objects. The
         * performance of hashCode() was measured, and found not to be the source of the performance
         * variations.
         */

        // Warmup 
        System.out.println("Beginning warmup phase...");
        { 
          {
            HashMap<UInt128, Object> map = new HashMap<>();

            for (int i = 0; i < 1000; i++) {
              for (int j = 0; j < 1000; j++) {
                UInt128 idValue = new UInt128(0, i*1000 + j);
                map.put(idValue, null);
              }
            }
          }

          {
            HashMap<UInt128, Object> map = new HashMap<>();

            for (int i = 0; i < 1000; i++) {
              for (int j = 0; j < 1000; j++) {
                UInt128 idValue = new UInt128(0, j*1000 + i);
                map.put(idValue, null);
              }
            }
          }
        }

        System.out.println("Completed warmup phase.");

        // Sequential
        System.out.println("Performing sequential id test...");
        long nodes = 24;
        long degree = 128;
        long numRuns = 1000;
        Long[] latency = new Long[(int)(nodes * degree * numRuns)];
        for (int r = 0; r < numRuns; r++) {
          HashMap<UInt128, Object> map = new HashMap<>();

          for (int i = 0; i < nodes; i++) {
            for (int j = 0; j < degree; j++) {
              UInt128 idValue = new UInt128(0, i*degree + j);
              long startTime = System.nanoTime();
              map.put(idValue, null);
              long endTime = System.nanoTime();
              latency[(int)(r*nodes*degree + i*degree + j)] = endTime - startTime;
            }
          }
        }

        {
          // Calculate latency statistics
          Arrays.sort(latency);

          long sum = 0;
          long min = Long.MAX_VALUE;
          long max = 0;
          for (int k = 0; k < latency.length; k++) {
            sum += latency[k];

            if (latency[k] < min)
              min = latency[k];

            if (latency[k] > max)
              max = latency[k];
          }

          long mean = sum / latency.length;

          long p25  = (int) (0.250 * (float) latency.length);
          long p50  = (int) (0.500 * (float) latency.length);
          long p75  = (int) (0.750 * (float) latency.length);
          long p90  = (int) (0.900 * (float) latency.length);
          long p95  = (int) (0.950 * (float) latency.length);
          long p99  = (int) (0.990 * (float) latency.length);
          long p999 = (int) (0.999 * (float) latency.length);

          System.out.println(String.format("%d, %d, %d, %d",
                min,
                mean,
                max,
                latency[(int)p50]));
        }

        // Random
        System.out.println("Performing random id test...");
        Random random = new Random();
        for (int r = 0; r < numRuns; r++) {
          HashMap<UInt128, Object> map = new HashMap<>();

          for (int i = 0; i < nodes; i++) {
            for (int j = 0; j < degree; j++) {
              UInt128 idValue = new UInt128(0, random.nextLong());
              long startTime = System.nanoTime();
              map.put(idValue, null);
              long endTime = System.nanoTime();
              latency[(int)(r*nodes*degree + i*degree + j)] = endTime - startTime;
            }
          }
        }

        {
          // Calculate latency statistics
          Arrays.sort(latency);

          long sum = 0;
          long min = Long.MAX_VALUE;
          long max = 0;
          for (int k = 0; k < latency.length; k++) {
            sum += latency[k];

            if (latency[k] < min)
              min = latency[k];

            if (latency[k] > max)
              max = latency[k];
          }

          long mean = sum / latency.length;

          long p25  = (int) (0.250 * (float) latency.length);
          long p50  = (int) (0.500 * (float) latency.length);
          long p75  = (int) (0.750 * (float) latency.length);
          long p90  = (int) (0.900 * (float) latency.length);
          long p95  = (int) (0.950 * (float) latency.length);
          long p99  = (int) (0.990 * (float) latency.length);
          long p999 = (int) (0.999 * (float) latency.length);

          System.out.println(String.format("%d, %d, %d, %d",
                min,
                mean,
                max,
                latency[(int)p50]));
        }

        // Big block
        System.out.println("Performing block id test...");
        for (int r = 0; r < numRuns; r++) {
          HashMap<UInt128, Object> map = new HashMap<>();

          for (int i = 0; i < nodes; i++) {
            for (int j = 0; j < degree; j++) {
              UInt128 idValue = new UInt128(0, j*1024 + i);
              long startTime = System.nanoTime();
              map.put(idValue, null);
              long endTime = System.nanoTime();
              latency[(int)(r*nodes*degree + i*degree + j)] = endTime - startTime;
            }
          }
        }

        {
          // Calculate latency statistics
          Arrays.sort(latency);

          long sum = 0;
          long min = Long.MAX_VALUE;
          long max = 0;
          for (int k = 0; k < latency.length; k++) {
            sum += latency[k];

            if (latency[k] < min)
              min = latency[k];

            if (latency[k] > max)
              max = latency[k];
          }

          long mean = sum / latency.length;

          long p25  = (int) (0.250 * (float) latency.length);
          long p50  = (int) (0.500 * (float) latency.length);
          long p75  = (int) (0.750 * (float) latency.length);
          long p90  = (int) (0.900 * (float) latency.length);
          long p95  = (int) (0.950 * (float) latency.length);
          long p99  = (int) (0.990 * (float) latency.length);
          long p999 = (int) (0.999 * (float) latency.length);

          System.out.println(String.format("%d, %d, %d, %d",
                min,
                mean,
                max,
                latency[(int)p50]));
        }
      } else if ((Boolean) opts.get("traverse_hashmap_put")) {
        long nodes_start = Long.decode(config.get("traverse_hashmap_put.nodes.start"));
        long nodes_end = Long.decode(config.get("traverse_hashmap_put.nodes.end"));
        long nodes_points = Long.decode(config.get("traverse_hashmap_put.nodes.points"));
        String nodes_mode = config.get("traverse_hashmap_put.nodes.mode");

        long degree_start = Long.decode(config.get("traverse_hashmap_put.degree.start"));
        long degree_end = Long.decode(config.get("traverse_hashmap_put.degree.end"));
        long degree_points = Long.decode(config.get("traverse_hashmap_put.degree.points"));
        String degree_mode = config.get("traverse_hashmap_put.degree.mode");

        long samples = Long.decode(config.get("traverse_hashmap_put.samples"));
       
        FileWriter outfile = null;
        if (opts.get("--output") !=  null) {
          try {
            outfile = new FileWriter((String)opts.get("--output"));
          } catch (Exception e) {
            throw new RuntimeException(e);
          }
        }

        List<Long> nodes_params = range(nodes_start, nodes_end, nodes_points, nodes_mode);
        List<Long> degree_params = range(degree_start, degree_end, degree_points, degree_mode);

        // Warmup 
        System.out.println("Beginning warmup phase...");

        { 
//          {
//            HashMap<UInt128, Object> map = new HashMap<>();
//
//            for (int i = 0; i < 1000; i++) {
//              for (int j = 0; j < 1000; j++) {
//                UInt128 idValue = new UInt128(0, i*1000 + j);
//                map.put(idValue, null);
//              }
//            }
//          }
//
//          {
//            HashMap<UInt128, Object> map = new HashMap<>();
//
//            for (int i = 0; i < 1000; i++) {
//              for (int j = 0; j < 1000; j++) {
//                UInt128 idValue = new UInt128(0, j*1000 + i);
//                map.put(idValue, null);
//              }
//            }
//          }


          long degree_prev = 0;
          Random random = new Random(0);
         
          // Use array of UInt128s to simulate edge lists.
          List<List<UInt128>> edgeLists = new ArrayList<>((int)nodes_end);
          for (long i = 0; i < nodes_end; i++)
            edgeLists.add(new ArrayList<>((int)degree_end));

          for (long degree : degree_params) {
            for (long i = 0; i < nodes_end; i++)
              for (long j = degree_prev; j < degree; j++)
                edgeLists.get((int)i).add(new UInt128(0, random.nextLong()));

            for (long nodes : nodes_params) {
              // Execution times are recorded in nanoseconds.
              Long[] execTimes = new Long[(int)samples];
              long numNeighbors = 0;
              for (long i = 0; i < samples; i++) {
                Set<UInt128> vSet = new HashSet<>();
                long startTime = System.nanoTime();
                for (long j = 0; j < nodes; j++) {
                  for (UInt128 neighborId : edgeLists.get((int)j)) {
                    vSet.add(neighborId);
                  }
                }
                long endTime = System.nanoTime();
                execTimes[(int)i] = endTime - startTime;
                numNeighbors = vSet.size();
              }
            }

            degree_prev = degree;
          }
        }

        System.out.println("Completed warmup phase.");

        Thread.sleep(5000);

        // Output the header
        // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
        if (outfile == null)
          System.out.println(
              "degree, nodes, min, mean, max, samples, 50th, 90th, 95th, 99th, 99.9th");
        else {
          outfile.append(
              "degree, nodes, min, mean, max, samples, 50th, 90th, 95th, 99th, 99.9th\n");
          outfile.flush();
        }

        long degree_prev = 0;
        Random random = new Random(0);
       
        // Use array of UInt128s to simulate edge lists.
        List<List<UInt128>> edgeLists = new ArrayList<>((int)nodes_end);
        for (long i = 0; i < nodes_end; i++)
          edgeLists.add(new ArrayList<>((int)degree_end));

        for (long degree : degree_params) {
          for (long i = 0; i < nodes_end; i++)
            for (long j = degree_prev; j < degree; j++)
              edgeLists.get((int)i).add(new UInt128(0, random.nextLong()));

          for (long nodes : nodes_params) {
            // Execution times are recorded in nanoseconds.
            Long[] execTimes = new Long[(int)samples];
            long numNeighbors = 0;
            for (long i = 0; i < samples; i++) {
              Set<UInt128> vSet = new HashSet<>();
              long startTime = System.nanoTime();
              for (long j = 0; j < nodes; j++) {
                for (UInt128 neighborId : edgeLists.get((int)j)) {
                  vSet.add(neighborId);
                }
              }
              long endTime = System.nanoTime();
              execTimes[(int)i] = endTime - startTime;
              numNeighbors = vSet.size();
            }

            // Calculate latency statistics
            Arrays.sort(execTimes);

            long sum = 0;
            long min = Long.MAX_VALUE;
            long max = 0;
            for (int k = 0; k < execTimes.length; k++) {
              sum += execTimes[k];

              if (execTimes[k] < min)
                min = execTimes[k];

              if (execTimes[k] > max)
                max = execTimes[k];
            }

            long mean = sum / execTimes.length;

            long p25  = (int) (0.250 * (float) execTimes.length);
            long p50  = (int) (0.500 * (float) execTimes.length);
            long p75  = (int) (0.750 * (float) execTimes.length);
            long p90  = (int) (0.900 * (float) execTimes.length);
            long p95  = (int) (0.950 * (float) execTimes.length);
            long p99  = (int) (0.990 * (float) execTimes.length);
            long p999 = (int) (0.999 * (float) execTimes.length);

            // Print out the results (in nanoseconds)
            // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
            if (outfile == null)
              System.out.println(String.format("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d", 
                                                degree,
                                                nodes,
                                                min,
                                                mean,
                                                max,
                                                samples,
                                                execTimes[(int)p50],
                                                execTimes[(int)p90],
                                                execTimes[(int)p95],
                                                execTimes[(int)p99],
                                                execTimes[(int)p999]));
            else {
              System.out.println(String.format(
                    "\t{Degree: %d, Nodes: %d, Nbrs: %d, Time: %d, EPS: %01.2f}",
                    degree,
                    nodes,
                    numNeighbors,
                    min,
                    1000.0 * (double)(degree * nodes) / (double)min));

              outfile.append(String.format("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", 
                                                degree,
                                                nodes,
                                                min,
                                                mean,
                                                max,
                                                samples,
                                                execTimes[(int)p50],
                                                execTimes[(int)p90],
                                                execTimes[(int)p95],
                                                execTimes[(int)p99],
                                                execTimes[(int)p999]));
              outfile.flush();
            }

            degree_prev = degree;
          }
        }

        if (outfile != null)
          outfile.close();
      } else if ((Boolean) opts.get("traverse_one_neighbor")) {
        Graph graph = new Graph(config);

        long warmup_nodes_start = Long.decode(config.get("traverse_one_neighbor.warmup.nodes.start"));
        long warmup_nodes_end = Long.decode(config.get("traverse_one_neighbor.warmup.nodes.end"));
        long warmup_nodes_points = Long.decode(config.get("traverse_one_neighbor.warmup.nodes.points"));
        String warmup_nodes_mode = config.get("traverse_one_neighbor.warmup.nodes.mode");

        long warmup_degree_start = Long.decode(config.get("traverse_one_neighbor.warmup.degree.start"));
        long warmup_degree_end = Long.decode(config.get("traverse_one_neighbor.warmup.degree.end"));
        long warmup_degree_points = Long.decode(config.get("traverse_one_neighbor.warmup.degree.points"));
        String warmup_degree_mode = config.get("traverse_one_neighbor.warmup.degree.mode");

        long warmup_samples = Long.decode(config.get("traverse_one_neighbor.warmup.samples"));

        long warmup_maxtime = Long.decode(config.get("traverse_one_neighbor.warmup.maxtime"));
       
        long nodes_start = Long.decode(config.get("traverse_one_neighbor.nodes.start"));
        long nodes_end = Long.decode(config.get("traverse_one_neighbor.nodes.end"));
        long nodes_points = Long.decode(config.get("traverse_one_neighbor.nodes.points"));
        String nodes_mode = config.get("traverse_one_neighbor.nodes.mode");

        long degree_start = Long.decode(config.get("traverse_one_neighbor.degree.start"));
        long degree_end = Long.decode(config.get("traverse_one_neighbor.degree.end"));
        long degree_points = Long.decode(config.get("traverse_one_neighbor.degree.points"));
        String degree_mode = config.get("traverse_one_neighbor.degree.mode");

        long samples = Long.decode(config.get("traverse_one_neighbor.samples"));
       
        FileWriter outfile = null;
        if (opts.get("--output") !=  null) {
          try {
            outfile = new FileWriter((String)opts.get("--output"));
          } catch (Exception e) {
            throw new RuntimeException(e);
          }
        }

        List<Long> warmup_nodes_params = range(warmup_nodes_start, warmup_nodes_end, 
            warmup_nodes_points, warmup_nodes_mode);
        List<Long> warmup_degree_params = range(warmup_degree_start, warmup_degree_end, 
            warmup_degree_points, warmup_degree_mode);

        List<Long> nodes_params = range(nodes_start, nodes_end, nodes_points, nodes_mode);
        List<Long> degree_params = range(degree_start, degree_end, degree_points, degree_mode);

        // WarmUp
        long warmupStartTime = System.currentTimeMillis();
        System.out.println(String.format("Beginning warmup (max %d minutes)...", warmup_maxtime));
        {
          String edgeLabel = "knows";
          Direction direction = Direction.OUT;
          Vertex baseVertex = new Vertex(new UInt128(1,0), "Warmup");
          Vertex neigVertex = new Vertex(new UInt128(1,0), "Warmup");
          byte[] keyPrefix = GraphHelper.getEdgeListKeyPrefix(baseVertex.id(), edgeLabel, direction,
                  neigVertex.label());

          int numElements = 10000;
          for (int j = 0; j < numElements; j++) {
            EdgeList.prepend(null, graph.getClient(), graph.getEdgeListTableId(), keyPrefix, 
                neigVertex.id(), new byte[0], 1024, 0);
          }

          for (int j = 0; j < 10000; j++) {
            EdgeList.batchReadSingleThreaded(null, graph.getClient(), graph.getEdgeListTableId(),
                Collections.singleton(baseVertex), edgeLabel, direction, false, "Warmup");
          }

          long degree_prev = 0;
          for (long degree : warmup_degree_params) {
            // Experiment setup: For each degree we setup all the nodes to have that degree. Each time
            // degree is increased, we add the additional edges needed to the existing nodes to bring
            // the degree up to the specified amount.
            for (long i = 0; i < warmup_nodes_end; i++) {
              Vertex startNode = new Vertex(new UInt128(1,i), "Person");
              for (long j = degree_prev; j < degree; j++) {
                // Vertices need labels, need to create vertex objects
                graph.addEdge(startNode, 
                              "knows", 
                              new Vertex(new UInt128(1, 0), 
                                "Person"), 
                              null);
              }
            }

            long nodes_prev = 0;
            Set<Vertex> startNodeSet = new HashSet<>();
            for (long nodes : warmup_nodes_params) {
              for (long i = nodes_prev; i < nodes; i++)
                startNodeSet.add(new Vertex(new UInt128(1,i), "Person"));

              for (long i = 0; i < warmup_samples; i++) {
                TraversalResult result = graph.traverse(startNodeSet, 
                                                        "knows", 
                                                        Direction.OUT, 
                                                        false, 
                                                        "Person");
              }

              if ((System.currentTimeMillis() - warmupStartTime)/(1000*60) >= warmup_maxtime)
                break;

              nodes_prev = nodes;
            }

            if ((System.currentTimeMillis() - warmupStartTime)/(1000*60) >= warmup_maxtime)
              break;

            degree_prev = degree;
          }
        }

        System.out.println(String.format(
              "Warmup phase complete. Total time: %d seconds",
              (System.currentTimeMillis() - warmupStartTime)/1000));

        // Output the header
        // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
        if (outfile == null)
          System.out.println(
              "degree, nodes, min, mean, max, samples, 50th, 90th, 95th, 99th, 99.9th");
        else {
          outfile.append(
              "degree, nodes, min, mean, max, samples, 50th, 90th, 95th, 99th, 99.9th\n");
          outfile.flush();
        }


        long degree_prev = 0;
        for (long degree : degree_params) {
          // Experiment setup: For each degree we setup all the nodes to have that degree. Each time
          // degree is increased, we add the additional edges needed to the existing nodes to bring
          // the degree up to the specified amount.
          for (long i = 0; i < nodes_end; i++) {
            Vertex startNode = new Vertex(new UInt128(0,i), "Person");
            for (long j = degree_prev; j < degree; j++) {
              // Vertices need labels, need to create vertex objects
              graph.addEdge(startNode, 
                            "knows", 
                            new Vertex(new UInt128(0, 0), "Person"), 
                            null);
            }
          }

          long nodes_prev = 0;
          Set<Vertex> startNodeSet = new HashSet<>();
          for (long nodes : nodes_params) {
            for (long i = nodes_prev; i < nodes; i++)
              startNodeSet.add(new Vertex(new UInt128(0,i), "Person"));

            // Execution times are recorded in nanoseconds.
            Long[] execTimes = new Long[(int)samples];
            long numNeighbors = 0;
            for (long i = 0; i < samples; i++) {
              long startTime = System.nanoTime();
              TraversalResult result = graph.traverse(startNodeSet, 
                                                      "knows", 
                                                      Direction.OUT, 
                                                      false, 
                                                      "Person");
              long endTime = System.nanoTime();
              execTimes[(int)i] = endTime - startTime;
              numNeighbors = result.vSet.size();
            }

            // Calculate latency statistics
            Arrays.sort(execTimes);

            long sum = 0;
            long min = Long.MAX_VALUE;
            long max = 0;
            for (int k = 0; k < execTimes.length; k++) {
              sum += execTimes[k];

              if (execTimes[k] < min)
                min = execTimes[k];

              if (execTimes[k] > max)
                max = execTimes[k];
            }

            long mean = sum / execTimes.length;

            long p25  = (int) (0.250 * (float) execTimes.length);
            long p50  = (int) (0.500 * (float) execTimes.length);
            long p75  = (int) (0.750 * (float) execTimes.length);
            long p90  = (int) (0.900 * (float) execTimes.length);
            long p95  = (int) (0.950 * (float) execTimes.length);
            long p99  = (int) (0.990 * (float) execTimes.length);
            long p999 = (int) (0.999 * (float) execTimes.length);

            // Print out the results (in microseconds)
            // degree, nodes, min, mean, max, 50th, 90th, 95th, 99th, 99.9th
            if (outfile == null)
              System.out.println(String.format("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d", 
                                                degree,
                                                nodes,
                                                min/1000,
                                                mean/1000,
                                                max/1000,
                                                samples,
                                                execTimes[(int)p50]/1000,
                                                execTimes[(int)p90]/1000,
                                                execTimes[(int)p95]/1000,
                                                execTimes[(int)p99]/1000,
                                                execTimes[(int)p999]/1000));
            else {
              System.out.println(String.format(
                    "\t{Degree: %d, Nodes: %d, Nbrs: %d, Time: %d, EPS: %01.2f}",
                    degree,
                    nodes,
                    numNeighbors,
                    min/1000,
                    1000.0 * (double)(degree * nodes) / (double)min));

              outfile.append(String.format("%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", 
                                                degree,
                                                nodes,
                                                min/1000,
                                                mean/1000,
                                                max/1000,
                                                samples,
                                                execTimes[(int)p50]/1000,
                                                execTimes[(int)p90]/1000,
                                                execTimes[(int)p95]/1000,
                                                execTimes[(int)p99]/1000,
                                                execTimes[(int)p999]/1000));
              outfile.flush();
            }
              
            nodes_prev = nodes;
          }

          degree_prev = degree;
        }

        if (outfile != null)
          outfile.close();

        graph.delete();
        graph.close();
      }
    }
  }
}
