'''
  File:   Bellman.py
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
  This file implements the Bellman-Held-Karp DP solution, drawing primarily
  from https://people.eecs.berkeley.edu/~vazirani/algorithms/chap6.pdf and
  then adds the special sauce from our paper to compute Shapley allocations.

  We never charge the 0th location anything -- this is the depot.


'''
import itertools
import math
import numpy as np


'''
  Compute and return the Shapley value and the min cost
  tour of the TSP.

  Parameters
  -----------
  dist: numpy array
    symmetric, metric TSP distance matrix.  Entry i,j represents
    the distance from i to j.

  Returns
  -----------
  sv: dict
    A dict mapping from agent_i --> Shapley value for i.

  gt: Cost of the grand tour

  Notes
  -----------
'''
def compute_shapley_bellman(dist):
  # Note that all these start at 0 implicity because we are only doing
  # the DP update from the initial set that all goes through 0.
  C = {}
  n = dist.shape[0]
  players = set(range(1,n))

  # All tours start and end at 0.
  tours = {}

  # Shapley Value vector
  sv = {i:0.0 for i in players}

  '''
    C[(S,k)] where S is a tuple of sorted locations
    and k is the location at the end.  C[(S,k)] is the cost
    to travel through all points in S and end at k.

    tours[S] is the tour starting and ending at 0 (depot)
    and visiting all locations in S where S is sorted.

    sv[i] is the Shapley value of location i.

    Seed the first few otherwise we have an empty sequence.
  '''
  #Using n implicitly omits the depot.
  bias = (math.factorial(len(players)-1)) / math.factorial(len(players))
  for i in range(1, n):
    C[(tuple([0,i]), i)] = dist[0, i]
    tours[(i,)] = dist[0,i] + dist[i,0]
    #Update Shapley:
    sv[i] += bias * tours[(i,)]

  # For all subsets larger than 2
  for s in range(2, n):
    #Iterate over all possible subsets
    for T in itertools.combinations(list(range(1,n)), s):
      S = set(list(T) + [0])
      #Sort the keys so we lookup the correct thing.
      kS = tuple(sorted(S))
      for j in T:
        # Pick the minimum length route that starts at 0, hits every
        # element of S, and ends at j.
        kSj = tuple(sorted(S - set([j])))
        C[(kS,j)] = min([C[(kSj,i)] + dist[i,j] for i in T if i != j])

      # For Shapley, we need S back to 0 -- IE shortest tour to here.
      kT = tuple(sorted(T))
      tours[kT] = min([C[(kS,j)] + dist[j,0] for j in T])
      # Now we have the shortest tour for T, we can update all T \ i
      #print(T)
      ### TODO: Something is wrong here with the biasing factor....
      # Bias is w.r.t. T-1 because we are *REMOVING* each i in T.
      bias = (math.factorial(len(T)-1) * math.factorial(len(players) - (len(T) - 1) - 1)) / math.factorial(len(players))
      for i in T:
        without = tuple(sorted(set(T) - set([i])))
        sv[i] += bias * (tours[kT] - tours[without])
  '''
  # Debug output.
  print("DP Table ending going through S and ending at j")
  for i in sorted(C.keys()):
    print(str(i) + " :: " + str(C[i]))
  print("Cost for tour of S")
  for i in sorted(tours.keys()):
    print(str(i) + " :: " + str(tours[i]))
  '''
  return sv, tours[tuple(players)]


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


  print("Distance Matrix is: \n" + str(dist))

  sv, gt = compute_shapley_bellman(dist)

  print("Grand Tour is: " + str(gt))
  print("Shapley Values : ")
  for i in sorted(sv.keys()):
    print(str(i) + " :: " + str(sv[i]))
