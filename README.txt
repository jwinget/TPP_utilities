Dependencies:
=============
python2.7 with virtualenv
libxml2
libxslt

Installation:
=============
> virtualenv env
> source env/bin/activate
> pip install -r REQUIREMENTS.txt

Note on requirements and dependencies:
======================================
Currently requires numpy and lxml modules.
- numpy is only used for FDR interpolation,
  and can be easily replaced with a single
  function by someone with more math skill than I have
- lxml is the xml parser for Python.

Usage example:
==============
./ProtXML2Csv.py -f 0.01 myfile.interact.iproph.prot.xml # Proteins filtered at 1% FDR