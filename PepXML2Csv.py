#-- PepXML2Csv
# Filter PepXML based on FDR
# Write results to CSV file
# Author: Jason Winget
# Version: 0.1
#--

import argparse
import csv
from PepXMLFilter import *

def GetFieldnames(hits):
	''' Find all field names '''
	print('Parsing field names')
	fn = ['run']
	for run, hitlist in hits.iteritems():
		for hit in hitlist:
			for k in hit.keys():
				if k not in fn:
					fn.append(k)
	return fn

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Filter PepXML for a given FDR and write results to a csv file'
	)
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

	ip = ProbCutoff(pepxml, fdr)
	at = GetAnalysisType(pepxml)
	hits = HitParse(pepxml, ip, at)
	fn = GetFieldnames(hits)

	with open(outfn, 'wb') as f:
		writer = csv.DictWriter(f, delimiter=',', fieldnames=fn)
		writer.writerow(dict((f, f) for f in fn)) # Write headers
		for run, hitlist in hits.iteritems():
			for hit in hitlist:
				hit['run'] = run
				writer.writerow(hit)
	print('Output written to '+outfn)
