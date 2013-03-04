TPP Utils
=========

Various utilities for handling TPP output, written for python 2.7

Installation/Usage
------------------
`virtualenv venv`

`source venv/bin/activate`

`pip install -r REQUIREMENTS`

Scripts
-------

### utils.py

For storing generally useful functions.

### [Pep|Prot]XMLFilter.py

Filter pepXML or protXML for a given false discovery rate (FDR) cutoff

### [Pep|Prot]XML2Csv.py

Write filtered pepXML or protXML to a csv file

### AssignIntensities.py

Assign MS1 intensities (calculated by Kronik/Hardklor) to Peptide search results
