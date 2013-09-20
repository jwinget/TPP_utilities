import csv

infile = 'SiteGrinder_cleaned.csv'
outfile = 'TissueSpecProts.csv'

cr = csv.reader(open(infile, 'rb'))
cr.next()

prot_dict = {}
for row in cr:
	acc = row[1]
	tissue = row[0]
	if acc not in prot_dict.keys():
		prot_dict[acc] = [tissue]
	else:
		prot_dict[acc].append(tissue)

cout = csv.writer(open(outfile, 'wb'))
cout.writerow(['Accession','Tissue'])
for k,v in prot_dict.iteritems():
	if len(v) == 1:
		cout.writerow([k, v[0]])

print 'Complete'
