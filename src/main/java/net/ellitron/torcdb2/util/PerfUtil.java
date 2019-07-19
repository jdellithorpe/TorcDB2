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
import java.lang.reflect.InvocationTargetException;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
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
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Properties;

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
      + "  PerfUtil [options] edgetest\n"
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
      + "  -h --help            Show this screen.\n"
      + "  --version            Show version.\n"
      + "\n";

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

      if ((Boolean) opts.get("edgetest")) {
        Graph graph = new Graph(config);
        String edgeLabel = "knows";
        Direction direction = Direction.OUT;
        int sizes[] = {8192,4096,2048,1024,512,256,128,64};

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

        for (int i = 0; i < sizes.length; i++) {
          Vertex baseVertex = new Vertex(new UInt128(0,i), "Person");
          Vertex neigVertex = new Vertex(new UInt128(0,0), "Person");
          byte[] keyPrefix = GraphHelper.getEdgeListKeyPrefix(baseVertex.id(), edgeLabel, direction,
                  neigVertex.label());

          int size = sizes[i];
          int numElements = 10000;
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

          System.out.println(String.format("%d, %d\n", 
                writeLatency[(int) (0.50 * (float) writeLatency.length)]/1000, 
                readLatency[(int) (0.50 * (float) readLatency.length)]/1000));
        }
      }
    }
  }
}
