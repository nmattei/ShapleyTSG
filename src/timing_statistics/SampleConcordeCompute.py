'''
Author: Nicholas Mattei
Date:	25-3-2015

Test out Reading/Writing/Calling Concord and other fun things.

Running Concorde for a specified number of iterations with ApproShapley and
returning the result.

Create the RAMDISK using the RAMDISK creator APP in the src directory. Then
feed this thing a number of iterations and it kicks out some Shapley values.

This is working for now to create and cache the files.

NOTE: This script uses a depot location at 0 which is not included in the inital
call.  So if you run it with 5 players, it is really using 5 + depot = 6.

'''

import random
import sys
import os
import glob
import time
import subprocess
import re
import itertools
import math
import copy
import collections
import numpy as np
import argparse

'''
		Read/Write Paths
'''
#Server Concorde:
CONCORDE_EXE = "./concorde/concorde"
#Mac Concorde
#CONCORDE_EXE = "./concorde/TSP/concorde"

#Server RAMDISK Path
RAMDISK_PATH = "/mnt/ramdisk/"
#Mac Path..
#RAMDISK_PATH = "/Volumes/RAMDISK"

# Charles made the datafiles: euctsp_size=n_sample=0.tsp
#DATA_PATH = ""
#FNAME_PREFIX = ""
#FNAME_SUFFIX = ""


'''
	Generate num_points in the 2D plane with 
	the size specified.  This number of cities
	is inclusive of the "depot"
	
	map is a list of tuples(x,y) rounded to 3 decimals.
'''
def generate_map(size_x, size_y, num_points):
	map = []
	for i in range(num_points):
		point_x = round(random.random()*size_x, 3)
		point_y = round(random.random()*size_y, 3)
		map.append((point_x, point_y))
	return map
	
def write_for_concord(map, fname):
	'''
	Write out a map in the Concord sense.
	Concord File Spec.
	
	NAME: ulysses22.tsp
	TYPE: TSP
	COMMENT: Shapley Experiment
	DIMENSION: 10
	EDGE_WEIGHT_TYPE: EUC_2D
	DISPLAY_DATA_TYPE: COORD_DISPLAY
	NODE_COORD_SECTION
	1 6433.62 3903.88
	2 9472.1 9819.56
	3 864.555 700.644
	4 6458.31 3541.85
	5 4810.06 3740.25
	6 2025.45 3233.4
	7 3583 3442.15
	8 9440.08 1942.7
	9 302.076 546.492
	10 3103.93 7758.99
	EOF
	'''
	#Open file
	outstr = ""
	
	outstr +="NAME: " + fname + "\n"
	outstr+="TYPE: TSP\n"
	outstr+="COMMENT: Shapley Experiment\n"
	outstr+="DIMENSION: " + str(len(map)) + "\n"
	outstr+="EDGE_WEIGHT_TYPE: EUC_2D\nDISPLAY_DATA_TYPE: COORD_DISPLAY\nNODE_COORD_SECTION\n"
	
	for i in range(len(map)):
		outstr+= str(i+1) + " " + str(map[i][0]) + " " + str(map[i][1]) + "\n"
	outstr += "EOF\n"
	
	#write the string.	
	outf = open(fname, 'w')
	outf.write(outstr)	
	outf.close()

#
# Simple interface, give it a list of [(x1,y1), (x2,y2) ... ] and let it fly.
#
def eval_with_concord(points):
	## Concorde processing...
	#Extract the optimal output...
	regex = re.compile("Optimal Solution: \d+.\d+$",re.MULTILINE)

	if len(points) >= 4:
		#Write lmap out to a file.
		concord_file_name = RAMDISK_PATH + "testmap.tsp"
		write_for_concord(points, concord_file_name)
				
		#Capture
		out = subprocess.check_output([CONCORDE_EXE, concord_file_name], stderr=subprocess.STDOUT)
		#print(out)
		#Extract the optimal output...
		opt = regex.findall(out.decode())
		if opt == []:
			print("No optimal")
			exit()
		else:
			v = float(opt[0].split(" ")[2].strip())	
	else:
		v = compute_complete(points)

	return v


## Should be 2 or 3. compute the distances...
def compute_complete(pts):
	if len(pts) < 2 or len(pts) > 3:
		print("Error computing Char function -- need to call Concord")
		exit()
	dist = 0.0
	if len(pts) == 2:
		dist = 2.0 * math.sqrt( math.pow(pts[0][0] - pts[1][0], 2) + math.pow(pts[0][1] - pts[1][1], 2))
	elif len(pts) == 3:
		#distance is p0 -> p1 -> p2 -> p0.
		dist = math.sqrt( math.pow(pts[0][0] - pts[1][0], 2) + math.pow(pts[0][1] - pts[1][1], 2))
		dist += math.sqrt( math.pow(pts[1][0] - pts[2][0], 2) + math.pow(pts[1][1] - pts[2][1], 2))
		dist += math.sqrt( math.pow(pts[2][0] - pts[0][0], 2) + math.pow(pts[2][1] - pts[0][1], 2))
	return dist

# Take the set of points and compute one appro shapley update. 
# Return the cache size and record the unormalized Shapley
# values at each "Frequency" step
#

def ApproConcorde(players, points, samples, frequency):
	sv_hash = {}
	time_hash = {}

	## Set each Player's start to 0.0.
	sv = {key:0.0 for key in players}

	## Eval Cache (has the charistic function evaluations we already did).
	## Players ---> Char Fc'n
	request_hash = {}
	hits_hash = {}
	cache_size_hash = {}
	cache = {}
	requests = 0
	hits = 0

	start_time = time.clock()
	## We know that the depot has cost 0 so pre-cache that
	cache["0"] = 0.0
	## We can precache the grand tour as well.
	grand_c = "0 " + "".join(" " + str(x) for x in players).strip()
	cache[grand_c] = eval_with_concord(points)

	temp_list = copy.copy(players)
	for cSample in range(0, samples+1):
		# Update Bookeeping.
		if (cSample) % frequency == 0 and cSample != 0:
			now = time.clock()
			#Hashes are mutable..
			sv_hash[cSample] = copy.copy(sv)
			time_hash[cSample] = now - start_time
			request_hash[cSample] = requests
			hits_hash[cSample] = hits
			cache_size_hash[cSample] = len(cache.keys())

		random.shuffle(temp_list)
		#Always include the depot location...
		key = [0]
		c_points = [points[0]]
		for cE in temp_list:
			## Run with and without cE.
			without_key = "".join(" " + str(x) for x in key).strip()
			key.append(int(cE))
			key.sort()
			with_key = "".join(" " + str(x) for x in key).strip()

			#only need to track points with cE because N\cE is cached.
			c_points.append(points[cE])
			#print(str(temp_list) + " :: " + without_key + " ---> " + with_key)
			
			#Check for prev run (we know withOut exists in the cache by def).
			
			requests += 2
			v_with = -1000.
			if with_key in cache.keys():
				v_with = cache[with_key]
				hits += 1
			else:
				v_with = eval_with_concord(c_points)
				#Cache the value...
				cache[with_key] = v_with

			hits += 1
			sv[cE] += (v_with - cache[without_key]) / cache[grand_c]


	# Go through and normalize the results...
	#[print(k,v) for k,v in sorted(cache.items())]
	#divide by the number of iterations.
	for c_sample in sv_hash.keys():
		for k, v in sv_hash[c_sample].items():
			sv_hash[c_sample][k] = v / c_sample

	return sv_hash, time_hash, request_hash, hits_hash, cache_size_hash

def compute_charistic_function(players, points):
	### Big ass function to compute
	### the charistic function key --> value
	### for all subsets of the points.

	char_func = {}
	
	## Concorde processing...
	#Extract the optimal output...
	regex = re.compile("Optimal Solution: \d+.\d+$",re.MULTILINE)
	
	for subset in itertools.chain.from_iterable([itertools.combinations(players, x) for x in range(1, len(players)+1)]):
		subset_key = "".join(" " + str(x) for x in subset).strip()
		c_points = []
		# Location 0 is the depot and in every set...
		c_points.append(points[0])
		
		for i in subset:
			c_points.append(points[i])
		
		if len(c_points) >= 4:
			#Write lmap out to a file.
			concord_file_name = RAMDISK_PATH + "testmap.tsp"
			write_for_concord(c_points, concord_file_name)
			#print("Sleeping to make sure we wrote the file correctly....")
			#time.sleep(2)
		
			#Capture
			out = subprocess.check_output([CONCORDE_EXE, concord_file_name], stderr=subprocess.STDOUT)
			#Extract the optimal output...
			opt = regex.findall(out.decode())
			if opt == []:
				print("No optimal")
				exit()
			else:
				opt = float(opt[0].split(" ")[2].strip())
				#print(str(opt))	
				char_func[subset_key] = opt
		else:
			d = compute_complete(c_points)
			#print(d)
			char_func[subset_key] = d
	
	#for i in sorted(char_func.keys()):
	#	print(str(i) + " :: " + str(char_func[i]))	
	return char_func
	
# computeShapley(cF, play)
#
# Computes and returns the normalized Shapley value for a set of players (play)
# cF must be a map where play --> value and play is a sorted string of the 
# elements of play (i.e. "0 1 2 3" is the coalition of players 0 1 2 and 3)
def computeShapley(play, cF):
	
	#Set each player's start value to 0..
	sV = {key:0.0 for key in play}

	#Determine the Grand Tour value...
	gt_value = cF["".join(" " + str(x) for x in play).strip()]

	#Iterate over all permutations....
	for perm in itertools.permutations(play):
		#For each Permutation... break it apart...
		allElets = []
		for cE in perm:
			#Build the appropriate current key and compute...
			if len(allElets) > 0:
				withOut = "".join(" " + str(x) for x in allElets).strip()
				allElets.append(int(cE))
				allElets.sort()
				withEle = "".join(" " + str(x) for x in allElets).strip()
				sV[cE] += (cF[str(withEle)] - cF[str(withOut)]) / gt_value
			else:
				#add it to the list of elements so far and update....
				allElets.append(int(cE))
				sV[cE] += cF[str(cE)] / gt_value
	
	# At each step we updated by a fraction 0--1 of the Grand Tour.
	# Divide it back out (value / samples)
	for k, v in sV.items():
		sV[k] = v / math.factorial(len(play))

	return sV	
	

# Write a game out in Char form.  For use with old version 
# of Shapley tester
def write_char_function(char_func, out):
	for k in sorted(char_func.keys()):
		out.write(str(k) + " = " + str(char_func[k]) + "\n")

# Helper function to randomly return a subset.
def genSubSet(subset):
	tset = []
	for i in range(len(subset)):
		if random.randint(0, 1):
			tset.append(int(subset[i]))
	return tset

#
# Modification of Subset to work with Concorde
# also note that this is modified so that each call
# updates a sample for EACH player so we can compare
# the number of samples as each sample here and in ApproShapley
# update the SV estimate for each player exactly once.
#

def SubsetConcorde(players, points, samples, frequency):
	sv_hash = {}
	time_hash = {}

	#Set each player's start value to 0..
	sV = {key:0.0 for key in players}
	

	## Eval Cache (has the charistic function evaluations we already did).
	## Players ---> Char Fc'n
	request_hash = {}
	hits_hash = {}
	cache_size_hash = {}
	cache = {}
	requests = 0
	hits = 0
	
	start_time = time.clock()
	## We know that the depot has cost 0 so pre-cache that
	cache["0"] = 0.0
	## We can precache the grand tour as well.
	grand_c = "0 " + "".join(" " + str(x) for x in players).strip()
	cache[grand_c] = eval_with_concord(points)
	
	#Iterate for samples...
	for cSample in range(0, samples+1):
		# Update Bookeeping.
		if cSample % frequency == 0 and cSample != 0:
			now = time.clock()
			sv_hash[cSample] = copy.copy(sV)
			time_hash[cSample] = now - start_time
			request_hash[cSample] = requests
			hits_hash[cSample] = hits
			cache_size_hash[cSample] = len(cache.keys())

		#Note that we do 1 update PER PLAYER so comparable...
		for cplayer in players:
			subset = copy.copy(players)
			subset.remove(cplayer)
			tset = genSubSet(subset)
			
			#Why are we biasing again? Check with Haris?
			bias = math.factorial(len(tset)) * math.factorial(len(players) - len(tset) - 1)
			# Removing this, it empirically makes no sense as we don't do it for permutations
			# either (# of permutations where that set of agents is added in front of i).

			#Get without_v
			without_v = -1000.
			tset.sort()
			without_key = "0 " + "".join(" " + str(x) for x in tset).strip()
			c_points = [points[0]] + [points[x] for x in tset]
			
			#Check cache and compute if necessary...
			requests += 1
			if len(c_points) > 1:
				if without_key in cache.keys():
					without_v = cache[without_key]
					hits+=1
				else:
					without_v = eval_with_concord(c_points)
					cache[without_key] = without_v
			else:
				hits += 1
				without_v = 0

			#get with_v
			with_v = -1000.
			tset.append(int(cplayer))
			tset.sort()
			with_key = "0 " + "".join(" " + str(x) for x in tset).strip()
			c_points = [points[0]] + [points[x] for x in tset]
			requests += 1
			if with_key in cache.keys():
				with_v = cache[with_key]
				hits+=1
			else:
				with_v = eval_with_concord(c_points)
				cache[with_key] = with_v

			#print(with_key + " :: " + without_key)
			sV[cplayer] += ((with_v - without_v)*bias) / cache[grand_c]
			
	# Go through and normalize the results...
	#[print(k,v) for k,v in sorted(cache.items())]
	#print(cache[grand_c])
	#divide by the number of iterations.
	for c_sample in sv_hash.keys():
		for k, v in sV.items():
			sV[k] = v / samples

	return sv_hash, time_hash, request_hash, hits_hash, cache_size_hash


#
# Compute some simple error functions..
#
def compute_error(shapley_value, appro_shapley_value):
	error_vector = []
	rmse = []
	for i in shapley_value.keys():
		error_vector.append(abs(shapley_value[i] - appro_shapley_value[i]))
		rmse.append(math.pow(shapley_value[i] - appro_shapley_value[i], 2.0))

	mae = sum(error_vector) / float(len(shapley_value.keys()))
	point_min = min(error_vector)
	point_max = max(error_vector)
	rmse = math.sqrt(sum(rmse) / float(len(shapley_value.keys())))

	return rmse, mae, point_min, point_max

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Concorde Based Shapley \
		Value Testing Suite. Run a variety of timing experiments with \
		Concorde as the backend. The number of agents IS NOT \
		INCLUSIVE OF THE DEPOT!')

	# Experiment Number
	parser.add_argument('-e', '--experiment', type=int, dest='experiment', \
		required=True, help="Experiment number to run. \n \
		(1) Time Samples ApproConcorde \n \
		(2) Time Samples SubsetConcorde \n \
		(3) Time Shapley Computation \n \
		(4) Error For Sample Number ApproConcorde\n \
		(5) Error For Sample Number SubsetConcorde\n ")

	# Number of locations to sweep.
	parser.add_argument('-ml', '--min_locations', type=int, dest='min_locations', default=2, help='Minimum number of agents in settings.')
	parser.add_argument('-Ml', '--max_locations', required=True, type=int, dest='max_locations', help='Maximum number of agents in settings.')

	# Seed value
	parser.add_argument('-sd', '--seed', type=int, default=-1, dest='random_seed', help='Random Seed, Default is time.')

	# Tracking results.
	parser.add_argument('-sa', '--samples', default=1000, type=int, dest='num_samples', help='Number of samples to run.')	
	parser.add_argument('-f', '--frequency', default=1, type=int, dest='write_frequency', help='Frequency with which to record timing/error samples.')
	parser.add_argument('-r', '--runs', default=1, type=int, dest='num_runs', help='Numer of times to repeat experiment step for jitter.')

	results = parser.parse_args()
	

	print("* Running experiment " + str(results.experiment) + " for " + \
		str(results.num_runs) + " run(s) per output. \n" + "* " + \
		str(results.num_samples) + " samples with write frequency " + str(results.write_frequency) + \
			" for agent range [" + str(results.min_locations) + "," + str(results.max_locations) + 
			"] and random seed " + str(results.random_seed))


	## Call the Particular Experiment ##

	if results.experiment == 1 or results.experiment == 2:
		#
		# Track timing and run just the sampler
		# for a certain number of steps with ApproConcorde. 
		# Write the time stamp out at the specified iteration number.
		#
		print("* ApproConcorde") if results.experiment == 1 else print("* SubsetConcorde")
		print("num_runs,num_players,sample_number,mean_requests,std_requests,mean_cache_size,std_cache_size,mean_hits,std_hits,mean_time,std_time")
		
		#Set the random seed.
		random.seed(results.random_seed)
		
		for c_players in range(results.min_locations, results.max_locations+1):
			
			#Hash from {SAMPLE NUMBER} ---> List of values for the run..
			total_requests_hash = collections.defaultdict(list)
			total_cache_size_hash = collections.defaultdict(list)
			total_hits_hash = collections.defaultdict(list)
			total_time_hash = collections.defaultdict(list)

			for c_run in range(results.num_runs):
				#Generate a random map to run on.
				lmap = generate_map(1000.0, 1000.0, c_players+1)
				p_keys = list(range(1,c_players+1))
				if results.experiment == 1:
					sv_hash, time_hash, request_hash, hits_hash, cache_size_hash = \
						ApproConcorde(p_keys, lmap, results.num_samples, results.write_frequency)	
				else:
					sv_hash, time_hash, request_hash, hits_hash, cache_size_hash = \
						SubsetConcorde(p_keys, lmap, results.num_samples, results.write_frequency)

				#Copy over.
				for c_sample in sorted(sv_hash.keys()):
					total_requests_hash[c_sample].append(request_hash[c_sample])
					total_cache_size_hash[c_sample].append(cache_size_hash[c_sample])
					total_hits_hash[c_sample].append(hits_hash[c_sample])
					total_time_hash[c_sample].append(time_hash[c_sample])

			for c_sample in sorted(total_requests_hash.keys()):
				print(",".join([str(item) for item in \
					[results.num_runs, c_players, c_sample, \
					np.mean(total_requests_hash[c_sample]), np.std(total_requests_hash[c_sample]), \
					np.mean(total_cache_size_hash[c_sample]), np.std(total_cache_size_hash[c_sample]), \
					np.mean(total_hits_hash[c_sample]), np.std(total_hits_hash[c_sample]), \
					np.mean(total_time_hash[c_sample]), np.std(total_time_hash[c_sample])]]))
					
	elif results.experiment == 3:
		#
		# Compute the shapley value for the agent 
		# numbers and record timings.
		#
		#
		# Attempt to do this with the input datafiles..
		#
		print("num_runs,num_players,mean_time,std_time")
		for c_players in range(results.min_locations, \
				results.max_locations+1):
			
			total_time = []
			for c_run in range(results.num_runs):
				random.seed(results.random_seed)
			
				lmap = generate_map(1000.0, 1000.0, c_players+1)
				p_keys = list(range(1,c_players+1))
			
				start_time = time.time()
				char_function = compute_charistic_function(p_keys, lmap)
				shapley_value = computeShapley(p_keys, char_function)
				end_time = time.time()

				total_time.append(end_time-start_time)

			print(",".join([str(item) for item in \
				[results.num_runs,c_players, \
				np.mean(total_time), np.std(total_time)]]))

	elif results.experiment == 4 or results.experiment == 5:
		#
		# Compute the shapley value for the agents 
		# and then run the sampler, computing the error
		# at the specified time step intervals.
		#
		print("* ApproConcorde") if results.experiment == 4 else print("* SubsetConcorde")
		print("num_runs,num_players,sample_number,mean_requests,std_requests," + \
			"mean_cache_size,std_cache_size,mean_hits,std_hits,mean_time,std_time," +\
			"mean_rmse,std_rmse,mean_mae,std_mae,mean_point_max,std_point_max," + \
			"mean_point_min,std_point_min")
		

		#Reset random Seed each set of iterations
		random.seed(results.random_seed)

		for c_players in range(results.min_locations, \
				results.max_locations+1):
			
			##Hash from {SAMPLE NUMBER} ---> List of values for the run..
			total_time_hash = collections.defaultdict(list)
			total_rmse_hash = collections.defaultdict(list)
			total_mae_hash = collections.defaultdict(list)
			total_point_min_hash = collections.defaultdict(list)
			total_point_max_hash = collections.defaultdict(list)
			total_requests_hash = collections.defaultdict(list)
			total_cache_hits_hash = collections.defaultdict(list)
			total_cache_size_hash = collections.defaultdict(list)


			#Re-run a specified number of times...
			for c_run in range(results.num_runs):
				lmap = generate_map(1000.0, 1000.0, c_players+1)
				p_keys = list(range(1,c_players+1))

				#Compute exact SV for all points.
				char_function = compute_charistic_function(p_keys, lmap)
				shapley_value = computeShapley(p_keys, char_function)
				#Normalize the SV for points
				total = sum(shapley_value.values())
				for key, value in shapley_value.items():
					shapley_value[key] = value / total

				if results.experiment == 4:
					sv_hash, time_hash, request_hash, hits_hash, cache_size_hash = \
						ApproConcorde(p_keys, lmap, results.num_samples, results.write_frequency)		
				else:
					sv_hash, time_hash, request_hash, hits_hash, cache_size_hash = \
						SubsetConcorde(p_keys, lmap, results.num_samples, results.write_frequency)

				# Need to extract per sample data:
				#Normalize each sample step, compute error, and append..
				for c_sample in sorted(sv_hash.keys()):
					total = sum(sv_hash[c_sample].values())
					for key, value in sv_hash[c_sample].items():
						sv_hash[c_sample][key] = value / total
				

					#compute statistics
					rmse, mae, point_min, point_max =  \
						compute_error(shapley_value, sv_hash[c_sample])

					#Aggregate...
					total_time_hash[c_sample].append(time_hash[c_sample])
					total_rmse_hash[c_sample].append(rmse)
					total_mae_hash[c_sample].append(mae)
					total_point_min_hash[c_sample].append(point_min)
					total_point_max_hash[c_sample].append(point_max)
					total_requests_hash[c_sample].append(request_hash[c_sample])
					total_cache_hits_hash[c_sample].append(hits_hash[c_sample])
					total_cache_size_hash[c_sample].append(cache_size_hash[c_sample])

			for c_sample in sorted(total_requests_hash.keys()):
				print(",".join([str(item) for item in \
					[results.num_runs, c_players, c_sample, \
					np.mean(total_requests_hash[c_sample]), np.std(total_requests_hash[c_sample]), \
					np.mean(total_cache_size_hash[c_sample]), np.std(total_cache_size_hash[c_sample]), \
					np.mean(total_cache_hits_hash[c_sample]), np.std(total_cache_hits_hash[c_sample]), \
					np.mean(total_time_hash[c_sample]), np.std(total_time_hash[c_sample]), \
					np.mean(total_rmse_hash[c_sample]), np.std(total_rmse_hash[c_sample]), \
					np.mean(total_mae_hash[c_sample]), np.std(total_mae_hash[c_sample]), \
					np.mean(total_point_min_hash[c_sample]), np.std(total_point_min_hash[c_sample]), \
					np.mean(total_point_max_hash[c_sample]), np.std(total_point_max_hash[c_sample])]]))

	else:
		print(" *** You have selected an invalid experiment. *** ")
