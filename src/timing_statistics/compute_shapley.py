'''
  File:   compute_shapley.py
  Author: Nicholas Mattei (nicholas.mattei@gmail.com)
  Date: May 31, 2016

  * Copyright (c) 2016, Nicholas Mattei and NICTA and CSIRO
  * All rights reserved.
  *
  * Developed by: Nicholas Mattei
  *               Data61/CSIRO and UNSW
  *               http://www.nickmattei.net
  *               http://www.preflib.org
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions are met:
  *     * Redistributions of source code must retain the above copyright
  *       notice, this list of conditions and the following disclaimer.
  *     * Redistributions in binary form must reproduce the above copyright
  *       notice, this list of conditions and the following disclaimer in the
  *       documentation and/or other materials provided with the distribution.
  *     * Neither the name of NICTA nor the
  *       names of its contributors may be used to endorse or promote products
  *       derived from this software without specific prior written permission.
  *
  * THIS SOFTWARE IS PROVIDED BY NICTA ''AS IS'' AND ANY
  * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  * DISCLAIMED. IN NO EVENT SHALL NICTA BE LIABLE FOR ANY
  * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


About
--------------------
  wrapper around bellman.py to do all the file handling, distance matrix making, and file writing.


'''
import Bellman
import itertools
import sys
import glob
import numpy as np

# Shamelessly stolen from: http://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python

import re

def sorted_nicely( l ):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


if __name__ == '__main__':

  # Box (EX from paper)
  dist = np.array([[0,    3,    4.25, 3   ],
                   [3,    0,    3,    4.25],
                   [4.25, 3,    0,    3   ],
                   [3,   4.25, 3,    0   ]])

  #Symmetric (no triangle)
  '''
  dist = np.array([[0, 2, 9, 10],
                   [2, 0, 6, 4],
                   [9, 6, 0, 8],
                   [10, 4, 8, 0]])
  '''

  #Asymetric
  '''
  dist = np.array([[0, 2, 9, 10],
                   [1, 0, 6, 4],
                   [15, 7, 0, 8],
                   [6, 3, 12, 0]])
  '''

  # dist = np.random.rand(11,11)

  #f = "/Users/Nick/repo/vrg.git-csirogithub/Collection/euctsp_size=5_sample=100.tsp"
  '''
    File Format:
    NAME: shapley.tsp
    TYPE: TSP
    COMMENT: Shapley Experiment
    DIMENSION: 5
    EDGE_WEIGHT_TYPE: EUC_2D
    DISPLAY_DATA_TYPE: COORD_DISPLAY
    NODE_COORD_SECTION
    1 4.16303 45.4492
    2 83.4817 33.5986
    3 56.5489 0.176691
    4 18.759 99.0434
    5 75.0497 36.6274
    EOF
  '''
  OUT_PATH = "/Users/Nick/repo/vrg.git-csirogithub/real-world/Final_Real_World/scaled_and_resolved/"

  locations = [10]

  for num_loc in locations:
    #files = glob.glob("/Users/Nick/repo/vrg.git-csirogithub/Collection/euctsp_size="+str(num_loc)+"*.tsp")
    files = glob.glob("/Users/Nick/repo/vrg.git-csirogithub/real-world/Final_Real_World/original/Size_10/*.tsp")
    out = ""
    for f in sorted_nicely(files):
      print(f)
      with open(f, "r") as infile:
        lines = infile.readlines()
      ### This block reads the TSP files that are points, not distance matricies....

      # lines = lines[7:len(lines)-1]
      # #print(lines)
      # points = {}
      # for i,l in enumerate(lines):
      #   bits = l.strip().split(" ")
      #   points[i] = (float(bits[1].strip()), float(bits[2].strip()))
      # #print(points)
      # # Make distance matrix...
      # dist = np.zeros((len(points), len(points)))
      # for i,j in itertools.combinations(points.keys(), 2):
      #   d = np.linalg.norm(np.array(points[i]) - np.array(points[j]))
      #   dist[i,j] = d
      #   dist[j,i] = d

      ### This block handles the TSP files that are distance matricies.
      print(lines)
      headder = lines[:7]
      headder[0] = "NAME: " + f.strip().split("/")[-1] + "\n"
      headder[2] = "COMMENT: " + str(num_loc) + " location problem fro C2S Optimization for " + f.strip().split("/")[-1].split("-")[0] + "\n"
      print(headder[3])
      lines = lines[7:len(lines)-1]
      print(lines)
      print(len(lines))
      if len(lines) == num_loc:
        dist = np.zeros((num_loc, num_loc))
        for i, l in enumerate(lines):
          for j, n in enumerate(l.strip().split(" ")):
            dist[i,j] = float(n)

        np.set_printoptions(suppress=True, linewidth=200)
        print("Distance Matrix is: \n" + str(dist))

        # Fix matricies, multipley by 100.0 and resolve to the largest...
        for i,j in itertools.combinations(range(num_loc), 2):
          # dist[i,j] *= 10.0
          # dist[j,i] *= 10.0
          if dist[j,i] != dist[i,j]:
            larger = max(dist[j,i], dist[i,j])
            dist[i,j] = larger
            dist[j,i] = larger

        np.set_printoptions(suppress=True, linewidth=200)
        print("Distance Matrix is: \n" + str(dist))

        # Write out a new file...
        wf = f.strip().split("/")[-1].split(".")[0]

        with open(OUT_PATH + wf + "." + str(num_loc) + ".tsp", "w") as outfile:
          for i in headder:
            outfile.write(i)
          for row in dist:
            out = ""
            for e in row:
              out += str(e) + " "
            outfile.write(out.strip() + "\n")
          outfile.write("EOF\n")


    #   sv, gt = Bellman.compute_shapley_bellman(dist)
    #   #print("Grand Tour is: " + str(gt))
    #   #print("Shapley Values : ")
    #   wf = f.strip().split("/")[-1]
    #   out += wf + ",0"
    #   for i in sorted(sv.keys()):
    #     out += "," + str(sv[i])
    #   out += "\n"
    # with open(OUT_PATH + str(num_loc) + ".csv", "w") as outfile:
    #   outfile.write(out)




