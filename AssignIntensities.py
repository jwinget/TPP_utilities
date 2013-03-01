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
from PepXMLFilter import *

class QPep(PepHit):
	''' Subclass of PepHit with attr for kronik intensity '''
	def __init__(self, *args, **kw):
		super(QPep, self).__init__(*args, **kw)
		self.intensities = []
	def add_intensity(self, intensity):
		self.intensities.append(intensity)
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
	
def GetIntensities(PepHit, kros):
	''' Find the appropriate intensities for a peptide hit from available kronik entries '''

	return intensities

if __name__ == '__main__':
	# Handle arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('pepXML',
						help = 'An iprophet pepXML file')
	parser.add_argument('-f', '--fdr', type = float, default = 0.01,
						help = 'An FDR cutoff value, e.g. 0.01.  Default is 0.01')
	parser.add_argument('-k', '--kro', default = '*.kro',
						help = 'An expression to match certain kronik files.  Default is *.kro')
	args = parser.parse_args()
	pepxml = args.pepXML
	fdr = args.fdr
	kexp = args.kro

	# Do things
	krofiles = glob.glob(kexp)
	kdict = {}
	for kf in krofiles:
		kp = KroParser(kf)
		kdict[kf] = kp
		print('Parsed '+str(len(kp))+' lines from '+kf)
	rh = RunHits(QPep, pepxml, fdr)
	fm = Matcher(rh.keys(), kdict.keys())
