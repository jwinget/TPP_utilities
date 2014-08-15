Tools for exporting TPP results to csv (Excel-compatible) format

## Installation (Requires python 2.7):
    virtualenv env
    source env/bin/activate
	pip install -r REQUIREMENTS.txt

## Usage example:
    python ProtXML2Csv.py -f 0.01 myfile.interact.iproph.prot.xml

This will generate myfile.interact.iproph.prot.csv, filtered at a 1% (model-based) FDR.
