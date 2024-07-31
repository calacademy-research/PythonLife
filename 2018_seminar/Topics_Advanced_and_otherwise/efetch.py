#!/usr/bin/env python

# usage: efetch.py IDs [ [-gb] | [-f <feature_name>|all [-q <qualifier_name>|all]] [-s max_seqlen] ]
#
# return the summary records that correspond to the Accession IDs 
# number of matching records written to stderr
# first line is column names of gID TaxID Length and Title
# rest of lines are summaries with these fields for the records
# fields separated by | since that is not a likely character in Title

import sys
from Bio import Entrez
Entrez.email = "ccg@calacademy.org"

def inJupyter(): # see if we are running in an iPython or Jupyter notebook
    try:
        get_ipython
    except:
        return False
    else:
        return True
    
def get_feature_qualifier_list(rec, feature_name, qualifier_name):
    result = []
    if not 'GBSeq_feature-table' in rec or feature_name == "":
        return result
    
    # build list of all qualifier values for feature_name[qualifier_name]
    # special case of feature_name of "all" for ALL feature_names
    # also special case for empty qualifier_name puts name of feature and name of qualifier in the list
    #   and if qualifier_name is "all" then also include value with feature and qualifier name
    special_case_qual = (qualifier_name == "" or qualifier_name.lower() == "all")
    
    for feat in rec['GBSeq_feature-table']:
        if (feature_name.lower() == "all" or feat['GBFeature_key'] == feature_name) and 'GBFeature_quals' in feat:
            for qual in feat['GBFeature_quals']:
                if qual['GBQualifier_name'] == qualifier_name:
                    result.append( qual['GBQualifier_value'] )
                elif special_case_qual: # special case: put feature name and qualifier name and maybe val in list
                    inf = feat['GBFeature_key'] + " \t" + qual['GBQualifier_name']
                    if qualifier_name.lower() == "all":
                        inf += " \t" + qual['GBQualifier_value']
                    result.append(inf)
    return result
    
def efetch(IDs, mode="xml", seq_end=-1, feat_name="", qual_name=""): #IDs is a string with comma separated values, each an ID
    global rec # only so we can inspect value (eg keys) later in a Jupyter notebook
    
    showing_feature = (feat_name != "")
    if mode == "xml" and feat_name == "":
        print "Accn | Length | UpdateDate | Title | Taxonomy"

    # if use a short seq_stop this speeds up retrieval quite a bit. but you won't get actual length values
    if seq_end < 0:
        handle = Entrez.efetch(db="nucleotide", id=IDs, rettype="gb", retmode=mode)
    else:
        handle = Entrez.efetch(db="nucleotide", id=IDs, rettype="gb", retmode=mode, strand=1,seq_start=1,seq_stop=seq_end)

    if mode == "text":
        print(handle.read())
    else:
        records = Entrez.read(handle)
        for rec in records:
            # set_trace() #enter ipdb debugger if we want to look at records interactively in Jupyter notebook
            if not showing_feature:
                print rec['GBSeq_primary-accession'], '|', rec['GBSeq_length'], '|', rec['GBSeq_update-date'], '|', \
                      rec['GBSeq_definition'],        '|', rec['GBSeq_taxonomy']
            else:
                features = get_feature_qualifier_list(rec, feat_name, qual_name)
                for qual in features:
                    print rec['GBSeq_primary-accession'] + " \t" + qual
        
def usage():
    sys.stderr.write('''
    Usage: efetch.py <ID1> [<ID2> ...] [ [-gb] | [-f <feature_name>|all [-q <qualifier_name>|all]] [-s max_seqlen] ]
    
    Examples:
        efetch.py KX773372 KX773303     # show Accn Length UpdateDate Title Taxonomy for each ID
        efetch.py KX773372 -gb          # show full textual GenBank output format
        efetch.py KX773372 -f source    # show name of a specific feature and each of its qualifier names
        efetch.py KX773372 -f source -q organelle # show specific feature and a specific qualifier value
        efetch.py MF431746 -f source -q all       # show specific feature and each of its qualifier names and value
        efetch.py MF431746 -f all                 # show names of each feature and each of its qualifier names
        efetch.py KX773372 -s 200       # limit the size of the sequence retrieved. this also limits which features are retrieved.
        
''');
    sys.exit(0)
 
def get_options(arg_list, min_args=2):
    class options:
        argv = arg_list
        id_list        = []
        mode           = "xml"
        sequence_max   = -1
        feature_name   = ""
        qualifier_name = ""
        
    num_args = len(options.argv)
    if num_args < min_args:
        usage()

    lst = options.id_list
    ixArg = 0; ixLast = num_args-1;
    while ixArg < ixLast:
        ixArg += 1
        idStr = arg_list[ixArg]
        if idStr[0] != "-":
            lst.append(idStr)
        elif idStr == "-gb" or idStr == "-text":
            options.mode = "text"
        else: # get options that have a value after the option flag
            optionType = idStr[1].lower()
            if optionType in "sfq" and ixArg < ixLast: # -s, -f or -q has a value after
                ixArg += 1
                optVal = arg_list[ixArg]
                if optionType == "f":
                    options.feature_name = optVal
                elif optionType == "q":
                    options.qualifier_name = optVal
                elif optionType == "s":
                    options.sequence_max = optVal
    return options

# main code entrypoint
if inJupyter():
    efetch("KC953095,KX773380,MF431746")
else:
    ops = get_options(sys.argv)
    ids = ",".join( ops.id_list )
    
    # call efetch with comma delimited set of IDs and optionally seq max, feature and qualifier name
    try:
        efetch(ids, ops.mode, seq_end=ops.sequence_max, feat_name=ops.feature_name, qual_name=ops.qualifier_name)
    except (KeyboardInterrupt, SystemExit): # allow Ctrl+C to end output without traceback
        sys.exit(0)
    except IOError as e: # ditto for closing pipe on output; e.g. piping to head
        if e.errno == errno.EPIPE:
            sys.exit(0)
