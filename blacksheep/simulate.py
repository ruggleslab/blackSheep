'''
File: simulate.py
Emily Kawaler, last updated 10/01/19

# Do we want to change the rounding to 2 decimal places? It's nicer but provides less information.
# Automatically replace existing pvalue table file?
# Better error detection, especially with the input arguments
# Future, if I'm having a slow week sometime: threading simulations
# Make number of points sampled from the KDE an input parameter (?)
# Make KDE vs true distribution an input parameter (?) 
'''

import random
import numpy as np
from scipy.stats import ttest_1samp
from scipy.stats import percentileofscore
from scipy.stats import gaussian_kde


def get_full_gene_list(ind_sep, infile):
	# Gets a list of all molecules in the file
	f = open(infile, 'r')
	genes = set([])
	line = f.readline() # Header line
	line = f.readline()
	while line:
		genes.add(line.split()[0].split(ind_sep)[0])
		line = f.readline()
	return genes


def get_values(gene, ind_sep, infile):
	# Digs around in the input file for the specified gene
	f = open(infile,'r')
	values = {}
	missings = {}
	line = f.readline()
	total = len(line.split())-1
	while line:
		if line.split(ind_sep)[0] == gene and len(line.split()) > 1:
			values[line.split()[0]] = [float(num) for num in line.split()[1:]]
			missings[line.split()[0]] = 1-(len(values[line.split()[0]])/float(total))
		line = f.readline()
		
	return values, missings


def outlier_thresholds(values, thresh):
	# Fixes the value that is (threshold) IQR above the median for each phospho
	o_thresh = {}
	for val in values.keys():
		dist = values[val]
		iqr = np.percentile(dist,75)-np.percentile(dist,25)
		outlier_threshold = np.percentile(dist,50)+(iqr*thresh)
		o_thresh[val] = outlier_threshold
	return o_thresh


def simulate(values, missings, o_thresh, reps):
	# Generates a random value for each p-site and counts outliers
	# Does this (reps) number of times
	r = random.Random()
	ol_list = []
	for i in range(reps):
		outliers = 0
		for val in values.keys():
			m = r.uniform(0.0,1.0)
			if m < missings[val]:
				# Makes sure to allow for missing psites.
				continue
			u = r.choice(values[val])
			if u > o_thresh[val]:
				outliers +=1
		ol_list.append(outliers)
	return ol_list


def simulate_kde(values, missings, o_thresh, reps):
	# Generates a KDE for each p-site from the existing values,
	# generates a random value for each p-site from that KDE, and counts outliers
	# Does this (reps) number of times
	r = random.Random()
	ol_list = []
	stored_kdes = {}
	for val in values.keys():
		# Make the size an editable parameter?
		stored_kdes[val] = gaussian_kde(values[val]).resample(size=reps)[0]
	for i in range(reps):
		outliers = 0
		for val in values.keys():
			m = r.uniform(0.0,1.0)
			if m < missings[val]:
				# Makes sure to allow for missing psites.
				continue
			# I used to pick randomly from the stored KDEs but realized that would
			# skew the distribution toward the mode, probably?
			u = stored_kdes[val][i]
			if u > o_thresh[val]:
				outliers +=1
		ol_list.append(outliers)
	return ol_list


def alpha_thresh(ol_dist):
	# Spits out outlier (p<0.05) number of outliers
	return np.percentile(ol_dist,95)
	#return np.mean(ol_dist)+1.64*np.std(ol_dist) # I think this is the 0.05, if you have a normal distribution, which you probably won't. 


def run_simulations(infile, ind_sep, thresh, reps, outfile, genes):

	w = open(outfile+"_pvals.tsv", 'w')
	
	# Write header line - not actually a thing
	#f = open(infile, 'r')
	#line = f.readline()
	#w.write(line)
	#f.close()
	
	if len(genes) == 0:
		genes = get_full_gene_list(ind_sep, infile)
	
	for gene in genes: # Can make this much more efficient
		# Get values for your gene from the input file
		values, missings = get_values(gene, ind_sep, infile)
		# Figure out the outlier threshold for each phosphosite
		o_thresh = outlier_thresholds(values, thresh)
		# Do the actual simulation
		ol_dist = simulate_kde(values, missings, o_thresh, reps)
		# Grab out the significance threshold for outliers
		out_line = alpha_thresh(ol_dist)
	
		# Write to the output file
		w.write(gene)
		if int(out_line) == len(values):
			pass
		else:
			for i in range(int(out_line), len(values)):
			#for i in range(int(out_line)+1,len(values)+1): This is how it used to be. I think I was wrong? If something downstream is broken, though, look here for a possible cause.
				w.write('\t'+str(i+1)+'\t'+str(round(((100-percentileofscore(ol_dist,i,kind='weak'))/100.0),3)))
		w.write('\n')
	w.close()
