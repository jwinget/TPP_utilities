#-- PepXMLFilter
#
# Filters PepXML based on FDR
# Maps resulting hits to class
# Author: Jason Winget
# Version: 0.2
#
# For speed, most XML parsing uses XPath expressions.
# These are likely to fail if the XML structure changes.
#--

from numpy import interp
from lxml import etree
from utils import Timer

# Just for testing
PEPXML = '1a_comet.interact.iproph.pep.xml'
FDR = 0.01

# PepXML namespace, required by parser
NS = { 'n': 'http://regis-web.systemsbiology.net/pepXML' }

class PepHit(object):
	''' A class to hold hits from PepXML '''
	def __init__(self, spectrum, scan, mass, charge, rt, peptide, protein, desc, iprob):
		self.spectrum = spectrum
		self.scan = scan
		self.mass = mass
		self.charge = charge
		self.rt = rt
		self.peptide = peptide
		self.protein = protein
		self.desc = desc
		self.iprob = iprob
		self.alt_prots = {}
	def add_alt_prot(self, altprot, desc):
		self.alt_prots[altprot] = desc

@Timer
def ReadXML(xmlfile):
	''' Read an XML file into memory '''
	with open(xmlfile, 'rb') as f:
		tree = etree.parse(xmlfile)
	return tree

def ProbCutoff(tree, fdr):
	''' Convert FDR to iprobability cutoff '''
	root = tree.getroot()
	nodes = root.xpath('n:analysis_summary/n:roc_data_point',
						namespaces = NS)
	prob = []
	error = []
	for node in nodes:
		prob.append(float(node.get('min_prob')))
		error.append(float(node.get('error')))
	ip = interp(fdr, error, prob) # This is the whole reason for requiring numpy
	print('iProbability of '+str(ip)+' gives FDR of '+str(fdr))
	return ip

def SepAnalyses(tree):
	''' Split PepXML into individual msms runs '''
	runs = tree.xpath('n:msms_run_summary', namespaces = NS)
	rundict = {}
	for run in runs:
		basename = run.get('base_name').split('/')[-1]
		rundict[basename] = run
	return rundict

@Timer
def HitParse(pepclass, subtree, ip):
	''' Find hits above the iprobability cutoff
	Parse to PepHit class '''
	# Get all spectrum queries
	sqs = subtree.xpath('n:spectrum_query', namespaces = NS)
	hits = []
	for sq in sqs:
		expr = 'n:search_result/n:search_hit/n:analysis_result[@analysis="interprophet"]/n:interprophet_result'
		iprob = float(sq.xpath(expr, namespaces = NS)[0].get('probability'))
		if iprob > ip: # Parse this spectrum query and add it to the output
			spectrum = sq.get('spectrum').split('.')[0]
			scan = int(sq.get('start_scan'))
			mass = float(sq.get('precursor_neutral_mass'))
			charge = int(sq.get('assumed_charge'))
			rt = float(sq.get('retention_time_sec')) / 60.0 # Convert to minutes
			bh = sq.xpath('n:search_result/n:search_hit', namespaces = NS)[0]
			peptide = bh.get('peptide')
			protein = bh.get('protein')
			try:
				desc = bh.get('protein_descr')
			except:
				desc = ''
			hit = pepclass(spectrum, scan, mass, charge, rt, peptide, protein, desc, iprob)
			try:
				altprots = bh.xpath('n:alternative_protein', namespaces = NS)
				for prot in altprots:
					altprot = prot.get('protein')
					try:
						altdesc = prot.get('protein_descr')
					except:
						altdesc = ''
					hit.add_alt_prot(altprot, altdesc)
			except:
				pass
			hits.append(hit)
	print('Found '+str(len(hits))+' significant hits from '+str(len(sqs))+' spectrum queries')
	return hits
	
def RunHits(pepclass, pepxml, fdr):
	''' A function to return a dict of hits '''
	tree = ReadXML(pepxml)
	ip = ProbCutoff(tree, fdr)
	msms_runs = SepAnalyses(tree)
	outdict = {}
	for run, sqs in msms_runs.iteritems():
		print('Parsing search hits from '+run)
		hits = HitParse(pepclass, sqs, ip)
		outdict[run] = hits
	return outdict

if __name__ == '__main__':
	RunHits(PepHit, PEPXML, FDR)
