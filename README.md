TPP Utils
=========

Various utilities for handling TPP output, written in python

Installation/Usage
------------------
`virtualenv venv`

`source venv/bin/activate`

`pip install -r REQUIREMENTS`

Scripts
-------

### utils.py

For storing generally useful functions.  Currently only holds a timer decorator

### [Pep|Prot]XMLFilter.py

Filter pepXML or protXML for a given false discovery rate (FDR) cutoff

### [Pep|Prot]XML2Csv.py

Write filtered pepXML or protXML to a csv file
