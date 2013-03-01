#-- ProtXML2Csv
# Filter ProtXML based on FDR
# Write results to CSV file
# Author: Jason Winget
# Version: 0.1
#--

import argparse
import csv
from ProtXMLFilter import *

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('protXML',
						help = 'A protXML file')
	parser.add_argument('-f', '--fdr', type=float, default=0.01,
						help = 'An FDR cutoff value, e.g. 0.01.  Default is 0.01')
	args = parser.parse_args()
	protxml = args.protXML
	fdr = args.fdr
	
	outfn = protxml[:-3]+'csv'
	rh = RunHits(protxml, fdr)
	with open(outfn, 'wb') as f:
		# Manual sorting, more prone to breaking than pulling the fieldnames from the data
		fieldnames = ['name', 'desc', 'probability', 'unique_peps', 'total_peps', 'coverage', 'perc_spec_ids']
		# fieldnames = vars(rh[rh.keys()[0]][0]) # Grab attributes from first class instance
		writer = csv.DictWriter(f, delimiter=',', fieldnames=fieldnames)
		writer.writerow(dict((fn, fn) for fn in fieldnames)) # Write headers
		for hit in rh:
			writer.writerow(vars(hit))
	
	print('Output written to '+outfn)
