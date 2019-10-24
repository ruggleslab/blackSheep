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
import sys
import argparse
import numpy as np
from scipy.stats import ttest_1samp
from scipy.stats import percentileofscore
from scipy.stats import gaussian_kde

# Argparser, when testing on its own
def _make_parser():
    parser = argparse.ArgumentParser(prog="blacksheep", description="")
    parser.add_argument("--version", "-v", action="version", version="%(prog)s 0.0.1")

    parser.add_argument(
        "values",
        type=str,
        help="File path to input values. Samples are columns and genes/sites are rows. Only .tsv "
             "and .csv accepted.",
    )
    parser.add_argument(
        "--ind_sep",
        type=str,
        default="-",
        help="Delimiter between the parent molecule (e.g. a gene name such "
             "as ATM) and a site identifier (e.g. S365). Default is -",
    )
    parser.add_argument(
        "--iqrs",
        type=float,
        default=1.5,
        help="Number of inter-quartile ranges (IQRs) above or below the "
             "median to consider a value an outlier. Default is 1.5.",
    )
    parser.add_argument(
        "--reps",
        type=int,
        default=1000000,
        help="Number of repetitions for the simulation to perform. "
             "Default is 1,000,000.",
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        default="simulated_pvals",
        help="Output prefix for writing files. Default is 'simulated_pvals'.",
    )
    parser.add_argument(
        "--molecules",
        nargs='+',
        default=[],
        help="List of parent molecules of interest. Empty list or absence of "
             "argument defaults to all parent molecules in input file.",
    )
    parser.add_argument(
        "--pval",
        type=float,
        default=0.05,
        help="p-value threshold for significant results. Must be between 0 and 1."
        "Default is 0.05.",
    )

    return parser

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
	values = {} # only non-missing values
	missings = {} # missing values
	all_values = {} # all values, including missing and present
	line = f.readline()
	total = len(line.split())-1
	while line:
		if line.split(ind_sep)[0] == gene and len(line.split()) > 1:
			temp_vals = [float(num) for num in line.split()[1:]]
			if len(temp_vals) > 1:
				values[line.split()[0]] = [float(num) for num in line.split()[1:]]
				missings[line.split()[0]] = 1-(len(values[line.split()[0]])/float(total))
				all_values[line.split()[0]] = line.rstrip('\n').split('\t')[1:]
		line = f.readline()
		
	return values, missings, all_values


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


def alpha_thresh(ol_dist, pval):
	# Spits out outlier number of outliers
	pv = 100-(100*pval)
	return np.percentile(ol_dist,pv)
	#return np.mean(ol_dist)+1.64*np.std(ol_dist) # I think this is the 0.05, if you have a normal distribution, which you probably won't. 

def generate_output_line(o_thresh, values, ol_dist, pval):
	# Determines the p-value for each sample being significantly hyperphosphorylated
	# for the given gene.
		
	output_list = []

	for s in range(len(list(values.values())[0])): # for each sample
		tot_outliers = 0
		for o in o_thresh.keys():
			if values[o][s]:
				if float(values[o][s]) > o_thresh[o]:
					tot_outliers += 1
		pv = round(((100-percentileofscore(ol_dist,tot_outliers-1,kind='weak'))/100.0),3) # pval for i+1 outliers
		if pv <= pval:
			output_list.append(str(pv))
		else:
			output_list.append('NS')
		
	return output_list
	
def run_simulations(infile, ind_sep, thresh, reps, outfile, genes, pval):
	w = open(outfile+"_pvals.tsv", 'w')
	
	# writing header to file
	f = open(infile, 'r')
	w.write(f.readline())
	f.close()
	
	if len(genes) == 0:
		genes = get_full_gene_list(ind_sep, infile)
	
	for gene in genes: # Can make this much more efficient
		# Get values for your gene from the input file
		values, missings, all_values = get_values(gene, ind_sep, infile)
		# Figure out the outlier threshold for each phosphosite
		o_thresh = outlier_thresholds(values, thresh)
		# Do the actual simulation
		ol_dist = simulate_kde(values, missings, o_thresh, reps)
		# Grab out the significance threshold for outliers
		out_line = alpha_thresh(ol_dist, pval)
		
		# Write to the output file
		if len(all_values) > 0: # won't write it out if the gene has no phosphosites with more than one value
			w.write(gene+'\t')
			to_write = generate_output_line(o_thresh, all_values, ol_dist, pval)
			w.write('\t'.join(to_write)+'\n')
		
	w.close()

if __name__=='__main__':

	args = sys.argv[1:]
	args = _make_parser().parse_args(args) 
    
	run_simulations(
            args.values,
            args.ind_sep,
            args.iqrs,
            int(args.reps),
            args.output_prefix,
            args.molecules,
            args.pval,
        )