import random
import sys
import os
import glob
import copy
import argparse


'''
	Recs are:
	4 0 0 0 0 0 0
	4 612.674 449.554 184.446 196.4 1621.71 98.8617
	4 670.715 480.612 1273.9 280.7 1979.6 178.947
	4 696.214 1109.33 119.401 1502.5 1979.6 819.503
	
	The fields are:
	(0) The number of customers (including depot)
	(1) ApproShapley Approx (1K iters) with Concorde for TSP solution
	(2) ApproShapley Approx (1K iters) with Christofedias Tree approx for TSP solution
	(3) Routing Margin Allocation (keeping the TSP solution for the grand tour, delete a customer and charge them the delta in the tour cost while keeping all the routes the same (simply skipping agent i in the tour).)
	(4) Moat Packing allocation
	(5) Reroute Marginal Allocation (delete customer i and then re-compute the optimal route.  Charging the deleted customer the delta in cost from the optimal route with i to the optimal route without i.)
	(6) Depot distance

	Quick and dirty script to merge the two different data files 
	Charles is geneating.  We need to take the Moat information
	out of the moat data and merge it with the rest of the columns 
	from the other raw_samples data. So from the --moat directory 
	we take element (5) and merge it with all the elements from 
	the --other directory.  writing merged
	files to the --out directory.

'''

# Parse a sample file and return it's contents as a 
# 2D array..
def parse_file(inf):
    ''' Parse info for the RAW Files... 
    Recs are:
    4 0 0 0 0 0 0
    4 612.674 449.554 184.446 196.4 1621.71 98.8617
    4 670.715 480.612 1273.9 280.7 1979.6 178.947
    4 696.214 1109.33 119.401 1502.5 1979.6 819.503
    
    The fields are:
    (1) The number of customers (including depot)
    (2) ApproShapley Approx (1K iters) with Concorde for TSP solution
    (3) ApproShapley Approx (1K iters) with Christofedias Tree approx for TSP solution
    (4) Routing Margin Allocation (keeping the TSP solution for the grand tour, delete a customer and charge them the delta in the tour cost while keeping all the routes the same (simply skipping agent i in the tour).)
    (5) Moat Packing allocation
    (6) Reroute Marginal Allocation (delete customer i and then re-compute the optimal route.  Charging the deleted customer the delta in cost from the optimal route with i to the optimal route without i.)
    (7) Depot distance
    
    '''
    
    #Read it and turn it into a 2d matrix for each line..
    return_data = []
    for line in inf:
        c_rec = line.strip().split(" ")
        c_rec = [float(x) for x in c_rec]
        return_data.append(c_rec)
    return return_data



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Quick and dirty script \
        to merge the two different data files Charles is geneating.      \
        We need to take the Moat information out of the moat data and    \
        merge it with the rest of the columns from the other raw_samples \
        data. So from the --moat directory we take element (5) and merge \
        it with all the elements from the --other directory.  Writing    \
        merged files to the --out directory.")

	# Experiment Number
    parser.add_argument('-m', '--moat', type=str, dest='moat_dir', \
        required=True, help='Directory containing moat files.')

    parser.add_argument('-s', '--other', type=str, dest='sample_dir', \
        required=True, help='Directory containing moat files.')

    parser.add_argument('-o', '--out', type=str, dest='out_dir', \
        required=True, help='Directory containing moat files.')

    results = parser.parse_args()

    print("\n *** Taking MOAT files from " + results.moat_dir + " and other \
        files from " + results.sample_dir + " and writing to " + results.out_dir + "\n\n")

    
    __SAMPLE_SEARCH_NAME = "/raw_samples-s"
    __MOAT_SEARCH_NAME = "/raw_samples-s"
    __N_LOCATIONS = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]


    for c_sample_num in __N_LOCATIONS:
        sample_files = glob.glob(results.sample_dir + __SAMPLE_SEARCH_NAME + \
            str(c_sample_num) + "*")
        moat_files = glob.glob(results.moat_dir + __MOAT_SEARCH_NAME + \
            str(c_sample_num) + "*")
        all_data = []
        for s,m in zip(sorted(sample_files), sorted(moat_files)):
            print("Matched: " + s + " :: " + m)
            with open(s, 'r') as input_file:
                sample_data = parse_file(input_file)
            with open(m, 'r') as input_file:
                moat_data = parse_file(input_file)

            for sample_row, moat_row in zip(sample_data, moat_data):
                # DEPOT should match (throw error if not) ELEMENT 7..
                if (sample_row[6] != moat_row[6]):
                    print("Error, Depot Numbers Do Not Match!")
                    exit()
                if (moat_row[4] == 0.0 and moat_row[6] != 0.0):
                    print("error..")
                    exit()

                # Replace the raw_sample with the info from raw_moat.
                sample_row[4] = moat_row[4]
            all_data += sample_data

        #Write back out...
        with open(results.out_dir + "merged_"+os.path.basename(s), "w") \
        as out_file:
            for line in all_data:
                out_file.write(" ".join(['{0:g}'.format(x) for x in line]) + "\n")








