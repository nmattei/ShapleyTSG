'''
	File		: 	OwaExp-MakeDisplay
	Author		:	Nicholas Mattei
	
	Parse and display some data from Experimental runs.
'''
import datetime
import itertools
import math
import copy
import time
import random
import sys

import glob

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats




# Create a class for the datainstances... just makes life easier.
class Metrics:
	SHAPLEY = "Shapley"
	CHRIS = "Christofides"
	ROUTEMARGIN = "Shortcut Distance"
	MOAT = "Moat Packing"
	REROUTE = "Re-routed Margin"
	DIST = "Depot Distance"
	MIX = "60/40 Moat/Depot"
	MIX2 = "25/75 Moat/Depot"
	TUP = (SHAPLEY, CHRIS, ROUTEMARGIN, MOAT, REROUTE, DIST, MIX, MIX2)
	GRAPH = (ROUTEMARGIN, REROUTE, DIST, MOAT, CHRIS, MIX, MIX2)
	STAT = (ROUTEMARGIN, REROUTE, DIST, MOAT, CHRIS, MIX, MIX2)
	
class route_test:
	def __init__(self):
		self.number_of_location = 0
		# Values is a dict from Location# --> Dictionary of Metrics with values...
		self.values = {}
		self.normalized = {}
		
	def print_values(self):
		print("Number of Locations: " + str(self.number_of_locations))
		outstr = ""
		for i in self.values.keys():
			for k in self.values[i]:
				print("Location " + str(i) + " : " + k + " = " + str(self.values[i][k]))
	
	def print_normalized(self):
		print("Number of Locations: " + str(self.number_of_locations))
		for i in self.normalized.keys():
			for k in self.normalized[i]:
				print("Location " + str(i) + " : " + k + " = " + str(self.normalized[i][k]))
				
	def compute_normalized(self):
		# For each type, sum every rec of that type, then divide each by the total.
		if self.values == {}:
			print("Big Error... Values aren't computed")
			exit()
		
		blank = {c:-1.00 for c in Metrics.TUP}
		self.normalized = {i:copy.copy(blank) for i in range(1, self.number_of_locations)}
		
		for c_metric in Metrics.TUP:
			running_sum = 0
			for c_index in self.values.keys():
				running_sum += self.values[c_index][c_metric]
			if running_sum == 0:
				print("GOT NO VALUES!")
				print(self.values)

			for k in self.values.keys():
				self.normalized[k][c_metric] = self.values[k][c_metric] / running_sum			
		


#parse off the data...
def parseDataFile(inf):
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
	
	#Print the first line of whatever we got...
	fin = inf.readlines()
	allrecs = []
	
	while len(fin) > 0:
		#Parse number of rec's
		current_example = route_test()
		current_example.number_of_locations = int(fin[0].split(" ")[0].strip())
		fin.pop(0)
		
		smallest_negitive_reroute = (-1000.0, 1)
		smallest_negitive_routemargin = (-1000.0, 1)
		for i in range(current_example.number_of_locations - 1):
			#Break Up an Entry...
			clocation = fin.pop(0)
			bits = clocation.strip().split(" ")
			bits = bits[1:len(bits)]
			values = {}
			for c in range(len(bits)):
				values[Metrics.TUP[c]] = float(bits[c])
			
			# 60/40 Blend...
			values[Metrics.MIX] = 0.60*values[Metrics.MOAT] + 0.40*values[Metrics.DIST]
			# 25/75 Blend...
			values[Metrics.MIX2] = 0.25*values[Metrics.MOAT] + 0.75*values[Metrics.DIST]

			## Check and 0 for negitive Shortcut and Re-route...
			if values[Metrics.REROUTE] < 0.0:
				#Capture the least negitive guy.. (Fix Kendall Tau error for unranked list)
				if values[Metrics.REROUTE] > smallest_negitive_reroute[0]:
					smallest_negitive_reroute = (values[Metrics.REROUTE], i+1)
				values[Metrics.REROUTE] = 0.0
			if values[Metrics.ROUTEMARGIN] < 0.0:
				if values[Metrics.ROUTEMARGIN] > smallest_negitive_routemargin[0]:
					smallest_negitive_routemargin = (values[Metrics.ROUTEMARGIN], i+1)
				values[Metrics.ROUTEMARGIN] = 0.0
			current_example.values[i+1] = values
		#current_example.print_values()
		
		# Check for potential all 0's for REROUTE or SHORTCUT...
		if all(current_example.values[i][Metrics.REROUTE] == 0.0000 for i in current_example.values.keys()):
			for i in current_example.values.keys():
				current_example.values[i][Metrics.REROUTE] = 1.0
			current_example.values[smallest_negitive_reroute[1]][Metrics.REROUTE] = 1.1
		
		if all(current_example.values[i][Metrics.ROUTEMARGIN] == 0.0000 for i in current_example.values.keys()):
			for i in current_example.values.keys():
				current_example.values[i][Metrics.ROUTEMARGIN] = 1.0
			print(str(smallest_negitive_routemargin))
			current_example.print_values()
			current_example.values[smallest_negitive_routemargin[1]][Metrics.ROUTEMARGIN] = 1.1
			print(str(smallest_negitive_routemargin))
			exit()
		

		allrecs.append(current_example)
	
	return allrecs

#Make the error box plots....	
def make_error_boxplot(pp, Big_Plot, x_labels, y_max):
	for c_metric in Metrics.STAT:
		#Create Figure Object
		fig = plt.figure(figsize=(10, 10))
		# Create an Axes object.
		ax1 = fig.add_subplot(1,1,1) # one row, two column, first plot
		
		#make pointers for labels and boxes...
		box_plots = []
		for i in sorted(Big_Plot[c_metric].keys()):
			box_plots.append(Big_Plot[c_metric][i])
		
		ax1.boxplot(box_plots)
		xtickNames = plt.setp(ax1, xticklabels=x_labels)
		#plt.setp(xtickNames)
		ax1.set_title("RMSE Box Plot -- " + str(c_metric))
		ax1.set_ylabel("Root Mean Squared Error")
		ax1.set_xlabel("Number of Locations")
		ax1.set_ylim(0, y_max)
		#draw it to the pdf.
		plt.draw()
		plt.tight_layout()
		plt.savefig(pp, format='pdf')

		
#Make the Snakes...
def make_error_snakes(pp, Big_Plot, x_labels, y_max):
	Mins = {c:[] for c in Metrics.STAT}
	Maxes = {c:[] for c in Metrics.STAT}
	Avg = {c:[] for c in Metrics.STAT}
	StDev = {c:[] for c in Metrics.STAT}
	ErMax = 0.7
	
	for c_metric in Metrics.STAT:
		for c_loc in sorted(Big_Plot[c_metric]):
			Mins[c_metric].append(min(Big_Plot[c_metric][c_loc]))
			Maxes[c_metric].append(max(Big_Plot[c_metric][c_loc]))
			Avg[c_metric].append(sum(Big_Plot[c_metric][c_loc]) / len(Big_Plot[c_metric][c_loc]))
			StDev[c_metric].append(np.std(Big_Plot[c_metric][c_loc]))

	print("Creating RMSE Graphs")
	
	#
	### THIS SECTION MAKES THE LINE GRAPH WITH ERROR BARS!!
	#
	
	# Create a Figure object.
	fig2 = plt.figure(figsize=(10, 10))
	# Create an Axes object.
	ax2 = fig2.add_subplot(1,1,1) # one row, one column, first plot
	colors = [ "green", "purple","blue","red", "yellow",  "cyan", "magenta"]
	markers = ["*", "h", "^", "o", "D", "s", "h"]

	c_element = 0
	for c_metric in Metrics.GRAPH:
		#ax2.plot(labels, Avg[c_metric], color=colors[c_element], marker=markers[c_element], label=c_metric)
		if len(Avg[c_metric]) != len(StDev[c_metric]):
			print(len(StDev[c_metric]))
		print(x_labels)
		print(Avg[c_metric])
		print(StDev[c_metric])
		ax2.errorbar(x_labels, Avg[c_metric], yerr=StDev[c_metric], color=colors[c_element], marker=markers[c_element], label=c_metric, alpha=0.9)
		#ax2.errorbar(labels, Avg[c_metric], yerr=[Mins[c_metric], Maxes[c_metric]], color=colors[c_element], marker=markers[c_element], label=c_metric, alpha=0.9)
		c_element += 1
		
	# Add a title.
	ax2.set_title("Average Root Mean Squared Error by Location")
	# Add some axis labels.
	ax2.set_xlabel("Number of Locations")
	ax2.set_ylabel("Average Root Mean Squared Error")
	#ax2.set_xlim(x_labels[0]-1, x_labels[len(x_labels)-1]+1)
	ax2.set_ylim(0, y_max)
	ax2.legend()
	plt.draw()
	plt.tight_layout()
	plt.savefig(pp, format='pdf')
	
	#
	### THIS SECTION MAKES SOME SNAKES!!
	#
	
	#Attempt to plot the Average Error for them all...
	# Create a Figure object.
	fig2 = plt.figure(figsize=(10, 10))
	# Create an Axes object.
	ax2 = fig2.add_subplot(1,1,1) # one row, one column, first plot
	colors = [ "green", "purple","blue","red", "yellow",  "cyan", "magenta"]
	markers = ["*", "h", "^", "o", "D", "s", "h"]

	c_element = 0
	for c_metric in Metrics.GRAPH:
		thigh = []
		tlow = []
		for i in range(len(Avg[c_metric])):
			thigh.append(Avg[c_metric][i]+StDev[c_metric][i])
			tlow.append(Avg[c_metric][i]-StDev[c_metric][i])
	
		ax2.plot(x_labels, Avg[c_metric], color=colors[c_element], marker=markers[c_element], label=c_metric, mec="black", linestyle="None")
		ax2.fill_between(x_labels, tlow, thigh, color=colors[c_element], alpha=0.2)
		c_element += 1
		
	# Add a title.
	ax2.set_title("Average Root Mean Squared Error by Location")
	# Add some axis labels.
	ax2.set_xlabel("Number of Locations")
	ax2.set_ylabel("Average Root Mean Squared Error")
	ax2.set_xlim(x_labels[0]-1, x_labels[len(x_labels)-1]+1)
	ax2.set_ylim(0, y_max)
	ax2.legend()
	plt.draw()
	plt.tight_layout()
	plt.savefig(pp, format='pdf')
	
	## Dump the stuff to the screen for tables:
	print("Average and SD RMSE by Metric and Size")
	outline = "{},".format('Metric')
	for i in sorted(Big_Plot[Metrics.CHRIS].keys()):
		outline += "{},".format(str(i)) + "{:^10},".format("SD")
	print(outline)
	for c_metric in Metrics.STAT:
		outline = "{},".format(c_metric)
		for i in range(len(Avg[c_metric])):
			outline += "{},".format(str(round(Avg[c_metric][i], 4))) + "{},".format(str(round(StDev[c_metric][i], 4)))
		print(outline)


	# print("Creating Min-Max RMSE Snake")
	
	# #
	# ### Snake fills with Min and Max as top and bottom elements.
	# #
	
	# #Attempt to plot the Average Error for them all...
	# # Create a Figure object.
	# fig2 = plt.figure(figsize=(10, 10))
	# # Create an Axes object.
	# ax2 = fig2.add_subplot(1,1,1) # one row, one column, first plot
	# colors = [ "green", "purple","blue","red", "yellow",  "cyan", "magenta"]
	# markers = ["*", "h", "^", "o", "D", "s", "h"]

	# c_element = 0
	# for c_metric in Metrics.GRAPH:
	# 	thigh = []
	# 	tlow = []
	# 	for i in range(len(Mins[c_metric])):
	# 		thigh.append(Mins[c_metric][i])
	# 		tlow.append(Maxes[c_metric][i])
	
	# 	ax2.plot(x_labels, Mins[c_metric], color=colors[c_element], marker=markers[c_element], label=c_metric, mec="black", linestyle="None")
	# 	ax2.plot(x_labels, Maxes[c_metric], color=colors[c_element], marker=markers[c_element], label=c_metric, mec="black", linestyle="None")
	# 	ax2.fill_between(x_labels, tlow, thigh, color=colors[c_element], alpha=0.2)
	# 	c_element += 1
		
	# # Add a title.
	# ax2.set_title("Min and Max Root Mean Squared Error by Number of Locations")
	# # Add some axis labels.
	# ax2.set_xlabel("Number of Locations")
	# ax2.set_ylabel("Min and Max Root Mean Squared Error")
	# ax2.set_xlim(x_labels[0]-1, x_labels[len(x_labels)-1]+1)
	# ax2.set_ylim(0, y_max)
	# ax2.legend()
	# plt.draw()
	# plt.tight_layout()
	# plt.savefig(pp, format='pdf')
	
	
#Make the boxplots for the KT-Distances for each # of locations...
def make_kt_boxplots(pp, kt_dist, x_labels):
	
	Avg = {c:[] for c in Metrics.STAT}
	StDev = {c:[] for c in Metrics.STAT}
	
	#agree_top = = {c:[] for c in Metrics.STAT}
	#agree_bottom = {c:[] for c in Metrics.STAT}
	median_p = {c:[] for c in Metrics.STAT}
	
	for c_metric in Metrics.STAT:
		#Create Figure Object
		fig = plt.figure(figsize=(10, 10))
		# Create an Axes object.
		ax1 = fig.add_subplot(1,1,1) # one row, two column, first plot
		
		#make pointers for labels and boxes...
		box_plots = []
		p_values = []
		for c_size in sorted(kt_dist[c_metric].keys()):
			bp = []
			#peal out all the tau's
			for c_t in sorted(kt_dist[c_metric][c_size].keys()):
				bp.append(kt_dist[c_metric][c_size][c_t][0])
				p_values.append(kt_dist[c_metric][c_size][c_t][1])
				
			Avg[c_metric].append(sum(bp) / float(len(bp)))
			StDev[c_metric].append(np.std(bp))
			box_plots.append(bp)
			median_p[c_metric].append(np.median(p_values))
		
		
		ax1.boxplot(box_plots)
		xtickNames = plt.setp(ax1, xticklabels=x_labels)
		#plt.setp(xtickNames)
		ax1.set_title("KT Distance Tau Box Plot -- " + str(c_metric))
		ax1.set_ylabel("KT-Distance")
		ax1.set_xlabel("Number of Locations")
		ax1.set_ylim(-1, 1)
		#draw it to the pdf.
		plt.draw()
		plt.tight_layout()
		plt.savefig(pp, format='pdf')
		
	## Make Charts for all the KT_errors...
	#
	### THIS SECTION MAKES THE LINE GRAPH WITH ERROR BARS!!
	#
	
	# Create a Figure object.
	fig2 = plt.figure(figsize=(10, 10))
	# Create an Axes object.
	ax2 = fig2.add_subplot(1,1,1) # one row, one column, first plot
	colors = [ "green", "purple","blue","red", "yellow",  "cyan", "magenta"]
	markers = ["*", "h", "^", "o", "D", "s", "h"]

	c_element = 0
	for c_metric in Metrics.GRAPH:
		#ax2.plot(labels, Avg[c_metric], color=colors[c_element], marker=markers[c_element], label=c_metric)
		ax2.errorbar(x_labels, Avg[c_metric], yerr=StDev[c_metric], color=colors[c_element], marker=markers[c_element], label=c_metric, alpha=0.9)
		#ax2.errorbar(labels, Avg[c_metric], yerr=[Mins[c_metric], Maxes[c_metric]], color=colors[c_element], marker=markers[c_element], label=c_metric, alpha=0.9)
		c_element += 1
		
	# Add a title.
	ax2.set_title("Average Kendall's Tau by Number of Location")
	# Add some axis labels.
	ax2.set_xlabel("Number of Locations")
	ax2.set_ylabel("Average Kendall's Tau")
	ax2.set_xlim(x_labels[0]-1, x_labels[len(x_labels)-1]+1)
	ax2.set_ylim(-1, 1)
	ax2.legend(loc=4)
	plt.draw()
	plt.tight_layout()
	plt.savefig(pp, format='pdf')
	
	#
	### THIS SECTION MAKES SOME SNAKES!!
	#
	
	#Attempt to plot the Average Error for them all...
	# Create a Figure object.
	fig2 = plt.figure(figsize=(10, 10))
	# Create an Axes object.
	ax2 = fig2.add_subplot(1,1,1) # one row, one column, first plot
	colors = [ "green", "purple","blue","red", "yellow",  "cyan", "magenta"]
	markers = ["*", "h", "^", "o", "D", "s", "h"]

	c_element = 0
	for c_metric in Metrics.GRAPH:
		thigh = []
		tlow = []
		for i in range(len(Avg[c_metric])):
			thigh.append(Avg[c_metric][i]+StDev[c_metric][i])
			tlow.append(Avg[c_metric][i]-StDev[c_metric][i])
	
		ax2.plot(x_labels, Avg[c_metric], color=colors[c_element], marker=markers[c_element], label=c_metric, mec="black", linestyle="None")
		ax2.fill_between(x_labels, tlow, thigh, color=colors[c_element], alpha=0.2)
		c_element += 1
		
	# Add a title.
	ax2.set_title("Average Kendall's Tau by Number of Location")
	# Add some axis labels.
	ax2.set_xlabel("Number of Locations")
	ax2.set_ylabel("Average Kendall's Tau")
	ax2.set_xlim(x_labels[0]-1, x_labels[len(x_labels)-1]+1)
	ax2.set_ylim(-1, 1)
	ax2.legend()
	ax2.invert_yaxis()
	plt.draw()
	plt.tight_layout()
	plt.savefig(pp, format='pdf')
	
	
	print("Average and SD KT Correlation Values by Metric and Size")
	outline = "{},".format('Metric')
	for i in sorted(kt_dist[Metrics.CHRIS].keys()):
		outline += "{},".format(str(i)) + "{},".format("SD")
	print(outline)
	for c_metric in Metrics.STAT:
		outline = "{},".format(c_metric)
		for i in range(len(Avg[c_metric])):
			outline += "{},".format(str(round(Avg[c_metric][i], 4))) + "{},".format(str(round(StDev[c_metric][i], 4)))
		print(outline)
	

# For each of the Metrics, across all sizes, we want to 
# make a scatter plot that is (SV Allocation %, Average RMSE).
# This may give us a visual on where the errors are occuring.
#
#
def make_alloc_error_plot(pp, x_val, y_val):
	for c_metric in Metrics.STAT:
				
		#One for each size and one for the union of all the sizes...
		x_union = []
		y_union = []
		
		for c_size in sorted(x_val[c_metric].keys()):
			x_clip = copy.copy(x_val[c_metric][c_size])
			y_clip = copy.copy(y_val[c_metric][c_size])
			
			#Clip off negitive values as they are the result of some numerical instablility
			for c in range(len(x_clip)):
				if x_clip[c] <= 0.000000000:
					x_clip[c] = 0.00000001
				if y_clip[c] <= 0.000000000:
					y_clip[c] = 0.00000001
			# Scale up y union to %'s
			#y_clip = [x*100 for x in y_clip]
			
			x_union += x_clip
			y_union += y_clip
		
			#Create Figure
			fig = plt.figure(figsize=(10,10))
			#create an axis
			plt.minorticks_on()
			sc_ax = fig.add_subplot(1,1,1)			
			sc_ax.scatter(x_clip, y_clip)
			sc_ax.set_title(str(c_size) + " Locations, Proxy: " + c_metric)
			sc_ax.set_ylabel("Absolute Difference between Fractional Shapley and Proxy")
			sc_ax.set_xlabel("Fraction of Total Shapley Allocation")
			sc_ax.set_xlim(0, 100)
			sc_ax.set_ylim(0, 20*100)
			#Logging...
			plt.xscale('symlog')
			sc_ax.xaxis.grid(True, which='major') 
			plt.yscale('symlog')
			sc_ax.yaxis.grid(True, which='major') 
			
			plt.draw()
			plt.tight_layout()
			plt.savefig(pp, format='pdf')
		
		#Create Figure
		fig = plt.figure(figsize=(10,10))
		plt.minorticks_on()
		#create an axis
		sc_ax = fig.add_subplot(1,1,1)
		#One more for the union...
		
		sc_ax.scatter(x_union, y_union)
		sc_ax.set_title("All Locations, Proxy: " + c_metric)
		sc_ax.set_ylabel("Absolute Difference between Fractional Shapley and Proxy")
		sc_ax.set_xlabel("Fraction of Total Shapley Allocation")
		sc_ax.set_xlim(0, 100)
		sc_ax.set_ylim(0, 20*100)
		#Logging...
		plt.xscale('symlog')
		sc_ax.xaxis.grid(True, which='major') 
		plt.yscale('symlog')
		sc_ax.yaxis.grid(True, which='major') 
		
		plt.draw()
		plt.tight_layout()
		plt.savefig(pp, format='pdf')
	
	
	
	
if __name__ == '__main__':

	if(len(sys.argv) != 3):
		print("Use: " + str(sys.argv[0]) + " <Input Directory> <OutputFile>\n")
		sys.exit()

	# Glob off all files in the input directory
	files = glob.glob(sys.argv[1] + "/merged_raw_samples*.data")
	all_rmse = {c:{} for c in Metrics.STAT}
	all_max_pointwise = {c:{} for c in Metrics.STAT}
	all_min_pointwise = {c:{} for c in Metrics.STAT}
	KT_Tests = {c:{} for c in Metrics.STAT}
	error_scatter_x = {c:{} for c in Metrics.STAT}
	error_scatter_y = {c:{} for c in Metrics.STAT}
	
	kt_top = {c:{} for c in Metrics.STAT}
	kt_bottom = {c:{} for c in Metrics.STAT}
	kt_nsig_p = {c:{} for c in Metrics.STAT}
	kt_median_p = {c:{} for c in Metrics.STAT}
	kt_worst_p = {c:{} for c in Metrics.STAT}
	
	
	x_labels = []
	
	YMAX = 0.0
	
	n_trials = []
	
	print(files)
	
	# Peal off each of the files and try to parse it...
	for cfile in sorted(files):
		print("Parsing: " + cfile)
		inf = open(cfile, 'r')
		recs = parseDataFile(inf)
		inf.close()
		
		#Normalize the records.
		[i.compute_normalized() for i in recs]
		n_trials.append(len(recs))

		# Print the Normalized Rec's...
		#[i.print_normalized() for i in recs]
		#exit()

		#Keep track of our x_labels...
		# Remove the Depot from n locs 
		x_labels.append(recs[0].number_of_locations - 1)
		x_labels = sorted(x_labels)
	
		#ease of use.
		n_locations = recs[0].number_of_locations

		#Create empty space in the KT_Test Map
		for c_metric in KT_Tests.keys():
			KT_Tests[c_metric][n_locations] = {}
		
		for c_metric in Metrics.STAT:
			all_rmse[c_metric][n_locations] = []
			all_max_pointwise[c_metric][n_locations] = []
			all_min_pointwise[c_metric][n_locations] = []

			error_scatter_x[c_metric][n_locations] = []
			error_scatter_y[c_metric][n_locations] = []
			
			kt_top[c_metric][n_locations] = 0
			kt_bottom[c_metric][n_locations] = 0
			m_p = []
			
			c_trial = 0
			for c_rec in recs:
				RMSE = 0
				sv = []
				metricrank = []
				t_error = []
				for c_location in sorted(c_rec.normalized.keys()):
					# RSME Calc
					point_error = c_rec.normalized[c_location][c_metric] - c_rec.normalized[c_location][Metrics.SHAPLEY]
					RMSE += point_error * point_error
					
					#KT stat Calc
					sv.append(c_rec.values[c_location][Metrics.SHAPLEY])
					metricrank.append(c_rec.values[c_location][c_metric])
					
					#Error Size Calc
					error_scatter_x[c_metric][n_locations].append(round(c_rec.normalized[c_location][Metrics.SHAPLEY] * 100.0, 4))
					#error_scatter_y[c_metric][recs[0].number_of_locations].append(round(abs(c_rec.normalized[c_location][c_metric] - c_rec.normalized[c_location][Metrics.SHAPLEY]) / c_rec.normalized[c_location][Metrics.SHAPLEY], 4))
					
					t_error.append(round(abs(c_rec.normalized[c_location][c_metric] - c_rec.normalized[c_location][Metrics.SHAPLEY]) * 100, 4))
					
					# error_scatter_y[c_metric][n_locations].append( \
					# 	round(abs(c_rec.normalized[c_location][c_metric] - c_rec.normalized[c_location][Metrics.SHAPLEY]) * 100, 4))
					
				
				all_max_pointwise[c_metric][n_locations].append(max(t_error))
				all_min_pointwise[c_metric][n_locations].append(min(t_error))
				error_scatter_y[c_metric][n_locations].append(t_error)	
				all_rmse[c_metric][n_locations].append(math.sqrt(RMSE))
				if math.sqrt(RMSE) > YMAX:
					YMAX = math.sqrt(RMSE)
					
				#KT test stat testing...
				KT_Tests[c_metric][n_locations][c_trial] = stats.mstats.kendalltau(sv, metricrank,use_ties=True)
				#KT_Tests[c_metric][recs[0].number_of_locations][c_trial] = stats.kendalltau(sv, metricrank)
				#KT_Tests[c_metric][recs[0].number_of_locations][c_trial] = -1, -1
				
				m_p.append(KT_Tests[c_metric][n_locations][c_trial][1])
				c_trial += 1
				if sv.index(max(sv)) == metricrank.index(max(metricrank)):
					kt_top[c_metric][n_locations] += 1
				if sv.index(min(sv)) == metricrank.index(min(metricrank)):
					kt_bottom[c_metric][n_locations] += 1
			
			kt_nsig_p[c_metric][n_locations] = sum(x < 0.05 for x in m_p)
			kt_median_p[c_metric][n_locations] = np.median(m_p)
			kt_worst_p[c_metric][n_locations] = max(m_p)
	
	## Attempt to set the font...
	font = {'size'   : 18}

	matplotlib.rc('font', **font)
				
	
	### Make a graph of the errors...
	#Try to make a box plot with the labels 'n' such...
	pp = PdfPages(sys.argv[2] + "AggError.pdf")
	
	print("Generating RMSE Error Boxes")
	make_error_boxplot(pp, all_rmse, x_labels, YMAX)
	
	#Let's graph the Average, Max, and Min points all together...
	print("Generating RMSE Snakes")
	make_error_snakes(pp, all_rmse, x_labels, YMAX)

	# Graph out the minimum and maximum pointwise errors for the locations.
	print("Generating Pointwise Max/Min Plots")

	
	
	# Setup some box plits...
	print("KT Graphs")
	make_kt_boxplots(pp, KT_Tests, x_labels)
	#print out some KT stats...
	print("\n")
	print("\nMedian KT p Values by Metric and Size")
	outline = "{},".format('Metric')
	for i in sorted(kt_median_p[Metrics.CHRIS].keys()):
		outline += "{},".format(str(i))
	print(outline)
	for c_metric in Metrics.STAT:
		outline = "{},".format(c_metric)
		for i in sorted(kt_median_p[c_metric].keys()):
			outline += "{},".format(str(round(kt_median_p[c_metric][i], 4)))
		print(outline)
	
	print("\nMax KT p Values by Metric and Size")
	outline = "{},".format('Metric')
	for i in sorted(kt_worst_p[Metrics.CHRIS].keys()):
		outline += "{},".format(str(i))
	print(outline)
	for c_metric in Metrics.STAT:
		outline = "{},".format(c_metric)
		for i in sorted(kt_worst_p[c_metric].keys()):
			outline += "{},".format(str(round(kt_worst_p[c_metric][i], 4)))
		print(outline)
	
	print("\nN KT p Values p < 0.05 by Metric and Size")
	outline = "{},".format('Metric')
	for i in sorted(kt_worst_p[Metrics.CHRIS].keys()):
		outline += "{},".format(str(i))
	print(outline)
	for c_metric in Metrics.STAT:
		outline = "{},".format(c_metric)
		for i in sorted(kt_worst_p[c_metric].keys()):
			outline += "{},".format(str(kt_nsig_p[c_metric][i]))
		print(outline)
	
	#Pull out possibly different number of trials per game size.
	print("\nNumber out of " + str(n_trials) + " with matching top ranked (costliest) elemnent.")
	outline = "{},".format('Metric')
	for i in sorted(kt_top[Metrics.CHRIS].keys()):
		outline += "{},".format(str(i))
	print(outline)
	for c_metric in Metrics.STAT:
		outline = "{},".format(c_metric)
		for i in sorted(kt_top[c_metric].keys()):
			outline += "{},".format(str(kt_top[c_metric][i]))
		print(outline)
	
	
	
	
	pp.close()
	
	# pd = PdfPages(sys.argv[2] + "Scatter_Error.pdf")
	# #Display the Average Normed RMSE as a function of allocation % in SV per metric.
	# make_alloc_error_plot(pd, error_scatter_x, error_scatter_y)
	# pd.close()
	