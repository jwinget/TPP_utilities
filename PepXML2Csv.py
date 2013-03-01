#-- PepXML2Csv
# Filter PepXML based on FDR
# Write results to CSV file
# Author: Jason Winget
# Version: 0.1
#--

import argparse
import csv
from PepXMLFilter import *

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('pepXML',
						help = 'An iprophet pepXML file')
	parser.add_argument('-f', '--fdr', type=float,
						help = 'An FDR cutoff value, e.g. 0.01.  Default is 0.01')
	args = parser.parse_args()
	pepxml = args.pepXML
	fdr = 0.01
	if args.fdr: # If unset, the default will be used
		fdr = args.fdr
	
	outfn = pepxml[:-3]+'csv'
	rh = RunHits(pepxml, fdr)
	with open(outfn, 'wb') as f:
		# Manual sorting, more prone to breaking than pulling the fieldnames from the data
		fieldnames = ['rt', 'scan', 'mass', 'charge', 'peptide', 'iprob', 'protein', 'desc', 'alt_prots', 'spectrum']
		# fieldnames = vars(rh[rh.keys()[0]][0]) # Grab attributes from first class instance
		writer = csv.DictWriter(f, delimiter=',', fieldnames=fieldnames)
		writer.writerow(dict((fn, fn) for fn in fieldnames)) # Write headers
		for run, hits in rh.iteritems():
			for hit in hits:
				writer.writerow(vars(hit))
	
	print('Output written to '+outfn)
