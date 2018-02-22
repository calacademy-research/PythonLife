#!/usr/bin/env python

# usage: esummary.py "quoted_query_string"
#
# return the summary records that correspond to the query
# number of matching records written to stderr
# first line is column names of Accn TaxID Length and Title
# rest of lines are summaries with these fields for the records
# fields separated by | since that is not a likely character in Title

import sys, errno
from Bio import Entrez
Entrez.email = "ccg@calacademy.org"

def inJupyter(): # see if we are running in an iPython or Jupyter notebook
    try:
        get_ipython
    except:
        return False
    else:
        return True

def eSummary(qryStr):
    handle = Entrez.esearch(db="nucleotide", term=qryStr, retmax=10000)
    record = Entrez.read(handle)
    numIDs = int(record["Count"])
    sys.stderr.write(str(numIDs) + " results for \""+qryStr+"\"\n")
    if numIDs < 1:
        sys.exit(0)

    if inNoteBook and numIDs > maxToShow: # while testing limit to maxToShow
        numIDs = maxToShow

    id_list = record["IdList"]
    id_list.sort()

    print "Accn | TaxonID | Length | Title"
    startIx = 0; sliceSize = min(numIDs,200) # can't send too large a request at a time, so cut it up into this size slices
    while numIDs > 0:
        ids = ",".join( id_list[startIx : startIx+sliceSize] )

        handle = Entrez.esummary(db="nucleotide", id=ids)
        records = Entrez.read(handle)

        for rec in records:
            print rec['Caption'], '|', rec['TaxId'], '|', rec['Length'], '|', rec['Title']

        numIDs  -= sliceSize
        startIx += sliceSize
        
def usage():
    sys.stderr.write("Usage: esummary.py \"query_str\"\n   Ex: esummary.py \"Strix[Title] AND mitochondrial[Title]\"\n");
    sys.exit(0)

# set vars that allow us to run in Jupyter notebook or from command line file
inNoteBook = inJupyter()
if inNoteBook:
    maxToShow = 25
    qryStr = "Strix[Title] AND mitochondrion[Title]"
    #qryStr = "1246431282 1246431296 523582230"
    #qryStr = "KC953095 MF431746 MF431745 KU365899 MF001440"
elif __name__ == "__main__": # execute only if run as a script
    if len(sys.argv) < 2:
        usage()
    qryStr = sys.argv[1] # ex: "Strix[Title] AND mitochondrion[Title]"

    try:
        # get our summaries with the defined query string
        eSummary(qryStr)
    except (KeyboardInterrupt, SystemExit): # allow Ctrl+C to end output without traceback
        sys.exit(0)
    except IOError as e: # ditto for closing pipe on output; e.g. piping to head
        if e.errno == errno.EPIPE:
            sys.exit(0)
