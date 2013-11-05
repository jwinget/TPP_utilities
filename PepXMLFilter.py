#-- PepXMLFilter
#
# Filters PepXML based on FDR
# Author: Jason Winget
# Version: 1.1
#
#--

from numpy import interp
from lxml import etree
from utils import Timer

# Just for testing
PEPXML = '../2013-08-29_LatinSquare/LSC1.interact.iproph.pep.xml'
FDR = 0.01

# PepXML namespace, required by parser
NS = { 'n': 'http://regis-web.systemsbiology.net/pepXML' }

def IterParse(xmlfile, xmltag, fn, *args):
	''' Iteratively parse an XML file based on a given tag '''
	context = etree.iterparse(xmlfile, tag='{'+NS['n']+'}'+xmltag)
	for event, elem in context:
		result = fn(elem, *args)
		elem.clear()
	return result

@Timer
def ProbCutoff(pepxml, fdr):
	''' Convert FDR to iprobability cutoff '''
	error = []
	prob = []
	def genProbList(elem):
		''' Function to populate list of error points based on models '''
		if elem.get('charge') == 'all':
			roc_points = elem.xpath('n:roc_data_point', namespaces=NS)
			for rdp in roc_points:
				error.append(float(rdp.get('error')))
				prob.append(float(rdp.get('min_prob')))
		return error, prob
	IterParse(pepxml, 'roc_error_data', genProbList)
	ip = interp(fdr, error, prob)
	print('Probability of '+str(ip)+' give FDR of '+str(fdr))
	return ip

def GetAnalysisType(pepxml):
	''' Find if the data was processed with iProphet '''
	ats = []
	def getAnalysis(elem):
		ats.append(elem.get('analysis'))
		return ats
	IterParse(pepxml, 'analysis_summary', getAnalysis)
	return ats[0]

@Timer
def HitParse(pepxml, ip, at):
	''' Find hits above the iprobability cutoff
	Parse to PepHit class '''
	hits = {}
	print('Parsing significant peptides')
	def hitParser(*args):
		''' Function to populate dictionary of hits '''
		elem = args[0]
		basename = elem.get('base_name').split('/')[-1]
		if basename not in hits.keys():
			hits[basename] = []
		spectrum_queries = elem.xpath('n:spectrum_query', namespaces = NS)
		for sq in spectrum_queries:
			kw = {}
			expr = 'n:search_result/n:search_hit/n:analysis_result[@analysis="'+at+'"]/n:'+at+'_result'
			prob = float(sq.xpath(expr, namespaces = NS)[0].get('probability'))
			if prob > ip:
				bn = sq.get('spectrum').split('.')[0]
				kw['probability'] = prob
				scan = int(sq.get('start_scan'))
				kw['scan'] = scan
				mass = float(sq.get('precursor_neutral_mass'))
				kw['mass'] = mass
				charge = int(sq.get('assumed_charge'))
				kw['charge'] = charge
				rt = float(sq.get('retention_time_sec')) / 60.0 # Convert to minutes
				kw['rt'] = rt
				bh = sq.xpath('n:search_result/n:search_hit', namespaces = NS)[0]
				peptide = bh.get('peptide')
				kw['peptide'] = peptide
				protein = bh.get('protein')
				kw['protein'] = protein
				desc = bh.get('protein_descr')
				kw['desc'] = desc
				ntt = int(bh.get('num_tol_term'))
				kw['ntt'] = ntt
				mc = int(bh.get('num_missed_cleavages'))
				kw['missed_cleavages'] = mc
				try: # Get modified peptide if present
					mh = bh.xpath('n:modification_info', namespaces = NS[0])
					modpep = mh.get('modified_peptide')
					kw['modpep'] = modpep
				except:
					pass
				try: # Get compensation voltage if present (LOLFAIMS)
					cv = int(sq.get('compensation_voltage'))
					kw['compensation_voltage'] = cv
				except:
					pass
				try: # Get precursor intensity if present
					pi = int(sq.get('precursor_intensity'))
					kw['prec_intensity'] = pi
				except:
					pass
				if bn not in hits.keys():
					print('Error, base name mismatch')
				else:
					hits[bn].append(kw)
		return hits
	IterParse(pepxml, 'msms_run_summary', hitParser)
	return hits
	
if __name__ == '__main__':
	ip = ProbCutoff(PEPXML, FDR)
	at = GetAnalysisType(PEPXML)
	print(at)
	#hits = HitParse(PEPXML, ip, at)
	#i = 0
	#for k, v in hits.iteritems():
	#	for p in v:
	#		i += 1
	#print('Found '+str(i)+' (non-unique) peptides above cutoff')
