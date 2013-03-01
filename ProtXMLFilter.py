#-- ProtXMLFilter
#
# Filters ProtXML based on FDR
# Maps resulting hits to class
# Author: Jason Winget
# Version: 0.1
#
# For speed, most XML parsing uses XPath expressions.
# These are likely to fail if the XML structure changes.
#--

from numpy import interp
from lxml import etree
from utils import Timer

# Just for testing
PROTXML = '1a_comet.interact.iproph.prot.xml'
FDR = 0.01

# ProtXML namespace, required by parser
NS = { 'n': 'http://regis-web.systemsbiology.net/protXML' }

class ProtHit(object):
	''' A class to hold hits from ProtXML '''
	def __init__(self, name, desc, probability, unique_peps, total_peps, coverage, perc_spec_ids):
		self.name = name
		self.desc = desc
		self.probability = probability
		self.unique_peps = unique_peps
		self.total_peps = total_peps
		self.coverage = coverage
		self.perc_spec_ids = perc_spec_ids

@Timer
def ReadXML(xmlfile):
	''' Read an XML file into memory '''
	with open(xmlfile, 'rb') as f:
		tree = etree.parse(xmlfile)
	return tree

def ProbCutoff(tree, fdr):
	''' Convert FDR to iprobability cutoff '''
	root = tree.getroot()
	exp = 'n:protein_summary_header/n:program_details/n:proteinprophet_details/n:protein_summary_data_filter'
	nodes = root.xpath(exp, namespaces = NS)
	prob = []
	error = []
	for node in nodes:
		prob.append(float(node.get('min_probability')))
		error.append(float(node.get('false_positive_error_rate')))
	error.reverse() # For numpy interp, x values must be increasing
	prob.reverse() # And therefore we flip these values as well
	pc = interp(fdr, error, prob) # This is the whole reason for requiring numpy
	print('Probability of '+str(pc)+' gives FDR of '+str(fdr))
	return pc

@Timer
def HitParse(tree, pc):
	''' Find hits above the probability cutoff
	Parse to ProtHit class '''
	# Get all spectrum queries
	proteins = tree.xpath('n:protein_group/n:protein', namespaces = NS)
	hits = []
	for p in proteins:
		prob = float(p.get('probability'))
		if prob > pc: # Parse this protein and add it to the output
			name = p.get('protein_name')
			probability = float(p.get('probability'))
			coverage = float(p.get('percent_coverage'))
			up = p.get('unique_stripped_peptides')
			ups = up.split('+')
			unique_peps = len([u for u in ups if not u == '+'])
			total_peps = int(p.get('total_number_peptides'))
			perc_spec_ids = float(p.get('pct_spectrum_ids'))
			try:
				desc = p.xpath('n:annotation', namespaces = NS)[0].get('protein_description')
			except:
				desc = ''
			hit = ProtHit(name, desc, probability, unique_peps, total_peps, coverage, perc_spec_ids)
			hits.append(hit)
	print('Found '+str(len(hits))+' significant hits from '+str(len(proteins))+' total proteins')
	return hits
	
def RunHits(pepxml, fdr):
	''' A function to return a list of hits '''
	tree = ReadXML(pepxml)
	p = ProbCutoff(tree, fdr)
	print('Parsing protein hits from '+pepxml)
	hits = HitParse(tree, p)
	return hits

if __name__ == '__main__':
	tree = ReadXML(PROTXML)
	prob = ProbCutoff(tree, FDR)
	hits = HitParse(tree, prob)
	print hits[0].name, hits[0].unique_peps, hits[0].total_peps, hits[0].desc
