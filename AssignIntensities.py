#-- AssignIntensities
#
# Maps intensities calculated with Hardklor/Kronik
# to peptide hits filtered by FDR
# Author: Jason Winget
# Version: 0.1
#
#--

import argparse
import glob
import csv
import difflib
import operator
from PepXMLFilter import *

def ParseArgs():
	''' Handle command-line arguments '''
	args_dict = {}
	parser = argparse.ArgumentParser(
		description='Map Kronik feature intensities to peptide search hits'
		)
	parser.add_argument('pepXML',
						help = 'An iprophet pepXML file')
	parser.add_argument('-f', '--fdr', type = float, default = 0.01,
						help = 'An FDR cutoff value, e.g. 0.01.  Default is 0.01')
	parser.add_argument('-k', '--kro', default = '*.kro',
						help = 'An expression to match certain kronik files.  Default is *.kro')
	parser.add_argument('-m', '--masstol', type = int, default = 10,
						help = 'A mass tolerance in ppm.  Default is 10')
	parser.add_argument('-r', '--rttol', type = int, default = 2,
						help = 'A retention time tolerance in minutes.  Default is 2')
	args = parser.parse_args()
	args_dict['pepxml'] = args.pepXML
	args_dict['fdr'] = args.fdr
	args_dict['kexp'] = args.kro
	args_dict['masstol'] = args.masstol
	args_dict['rttol'] = args.rttol
	return args_dict

class QPep(PepHit):
	''' Subclass of PepHit with attr for kronik intensity '''
	def __init__(self, *args, **kw):
		super(QPep, self).__init__(*args, **kw)
		self.intensity = []
	def add_intensity(self, mean_int, stdev):
		self.intensity.append((mean_int, stdev))
	def __repr__(self):
		return '%s - %s' % (self.spectrum, self.scan)

class Kro(object):
	''' Kronik entry '''
	def __init__(self, fscan, lscan, numscans, charge, monomass, basepeak, bint, sint, frt, lrt, brt, bcorr, mods):
		self.scanrange = (fscan, lscan)
		self.numscans = numscans
		self.charge = charge
		self.monomass = monomass
		self.basepeak = basepeak
		self.bint = bint
		self.sint = sint
		self.rtrange = (frt, lrt)
		self.brt = brt
		self.bcorr = bcorr
		self.mods = mods
	def __repr__(self):
		return '%s - %s' % (self.monomass, self.charge)

@Timer
def KroParser(f):
	''' Parse kronik results and map to class '''
	khits = []
	with open(f, 'rb') as kf:
		reader = csv.reader(kf, delimiter='\t')
		reader.next() # Skip headers
		for row in reader:
			fscan = int(row[0])
			lscan = int(row[1])
			numscans = int(row[2])
			charge = int(row[3])
			monomass = float(row[4])
			basepeak = float(row[5])
			bint = float(row[6])
			sint = float(row[7])
			frt = float(row[8])
			lrt = float(row[9])
			brt = float(row[10])
			bcorr = float(row[11])
			mods = row[12]
			khit = Kro(fscan, lscan, numscans, charge, monomass, basepeak, bint, sint, frt, lrt, brt, bcorr, mods)
			khits.append(khit)
	return khits
	
def Matcher(iter1, iter2):
	''' Try to match filenames to one another
	This should _always_ find a match, which is not necessarily good'''
	print('-'*10)
	print('Matching MS/MS runs to Kronik files')
	print('-'*10)
	matches = []
	for i in iter1:
		c = 0.6
		closest_match = ''
		while (c > 0) and not closest_match:
			closest_match = difflib.get_close_matches(i, iter2, cutoff=c, n=1)
			if not closest_match:
				c -= 0.1
			else:
				try:
					matches.append((i, closest_match[0]))
				except:
					print 'No matches found for '+i
	for match in matches:
		print match[0]+' ==> '+match[1]
	print('-'*10)
	return matches

if __name__ == '__main__':
	args_dict = ParseArgs()
	rh = RunHits(QPep, args_dict['pepxml'], args_dict['fdr'])
	unique = {}
	print('Parsing peptides to remove those with indistinguishable RT or mass')
	for run, hits in rh.iteritems(): # Further winnow the hits to remove 'indistinguishable' peptides
		seen = []
		unique[run] = []
		for hit in hits:
			rt_lb = hit.rt - args_dict['rttol']
			rt_ub = hit.rt + args_dict['rttol']
			da_error = hit.mass * (args_dict['masstol'] / 1E6)
			mass_lb = hit.mass - da_error
			mass_ub = hit.mass + da_error
			found = 0
			while found == 0:
				for pep in seen:
					if hit.peptide == pep[0] and rt_lb <= pep[1] <= rt_ub and mass_lb <= pep[2] <= mass_ub:
						unique[run].append(hit)
						found = 1
				seen.append((hit.peptide, hit.rt, hit.mass))
				found = 1
	print('Parsing Kronik files')
	krofiles = glob.glob(args_dict['kexp'])
	kdict = {}
	for kf in krofiles:
		kp = KroParser(kf)
		kdict[kf] = kp
		print('Parsed '+str(len(kp))+' lines from '+kf)
	fm = Matcher(unique.keys(), kdict.keys())
	for m in fm: # Let's get down to brass tacks
		print('Assigning intensities to '+m[0])
		peps = unique[m[0]]
		kros = kdict[m[1]]
		for pep in peps:
			kromatches = []
			mass_lb = pep.mass - (pep.mass * (args_dict['masstol'] / 1E6))
			mass_ub = pep.mass + (pep.mass * (args_dict['masstol'] / 1E6))
			rt_lb = pep.rt - args_dict['rttol']
			rt_ub = pep.rt + args_dict['rttol']
			for k in kros:
				if k.rtrange[0] <= pep.rt <= k.rtrange[1] and mass_lb <= k.monomass <= mass_ub and k.charge == pep.charge:
					kromatches.append(k)
					kros.remove(k)
			if kromatches:
				intensities = [i.sint for i in kromatches]
				if len(intensities) > 1:
					int_mean = numpy.mean(intensities)
					int_stdev = numpy.std(intensities)
				else:
					int_mean = intensities[0]
					int_stdev = 0
				pep.add_intensity(int_mean, int_stdev)
		outfile = m[0]+'.csv'
		with open(outfile, 'wb') as f:
			fieldnames = ['rt', 'mass', 'charge', 'peptide',
							'iprob', 'protein', 'desc', 'alt_prots', 'intensity']
			writer = csv.DictWriter(f, delimiter=',', fieldnames=fieldnames)
			writer.writerow(dict((fn, fn) for fn in fieldnames))
			getter = operator.attrgetter('intensity')
			for pep in peps:
				if getter(pep):
					line = {'rt': pep.rt,
							'mass': pep.mass,
							'charge': pep.charge,
							'peptide': pep.peptide,
							'iprob': pep.iprob,
							'protein': pep.protein,
							'desc': pep.desc,
							'alt_prots': pep.alt_prots,
							'intensity': pep.intensity[0][0]}
					writer.writerow(line)
		print('Intensities written to '+outfile)
