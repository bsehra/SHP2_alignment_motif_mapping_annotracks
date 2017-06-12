#!usr/bin/env python

"""Program to create a fasta input file for pattern searching using PatMan"""
#Written by Bhupinder Sehra; updated 8/15/16

import sys
import os
import pandas as pd

motifin = sys.argv[1]
patmaninput = sys.argv[2]
coldata1 = sys.argv[3]
coldata2 = sys.argv[4]

def makeIUPAC(motif):
    """Take in a degenerate consensus sequence and create motif with IUPAC codes"""
    codemotif = motif.upper().replace("[","(").replace("]",")").replace("(A/G)","R").replace("(G/A)","R").replace("(C/T)","Y").replace("(T/C)","Y").replace("(C/G)","S").replace("(G/C)","S").replace("(A/T)","W").replace("(T/A)","W").replace("(G/T)","K").replace("(T/G)","K").replace("(A/C)","M").replace("(C/A)","M").replace("(C/G/T)","B").replace("(A/G/T)","D").replace("(A/C/T)","H").replace("(A/C/G)","V").replace("(A)","A").replace("(T)","T").replace("(C)","C").replace("(G)","G")
    return codemotif

"""Main body. d is a dataframe and motifdict a dict to store columns from dataframe produced by reading in .xlsx file"""
temp = []
motifs = []
cmotif = ""
desc = ""
d = {}
motifdict = {}

"""read in tab delimited file with AGRIS binding sites from two columns. Column called 'name' or 'Name' contains fasta descriptor for motif
column 2 contains consensus motif."""
if  os.path.splitext(motifin)[-1].lower() == ".txt" or os.path.splitext(motifin)[-1].lower() == ".csv":
    with open(patmaninput, "w") as pin:
        with open(motifin, "r") as filein:
            for line in filein:
                line.strip()
                if not line.strip():
                    continue
                if line.startswith("Name") or line.startswith("name"):
                    continue
                else:
                    temp = line.split()
                    motifs = temp[:2]
                    descr = str(temp[1]).replace("[","(").replace("]",")")
                    cmotif = makeIUPAC(descr)
                    motifs.append(cmotif)
                    print str(motifs).replace("'","").replace(',',"").strip('[]')
                    if len(cmotif) > 0:
                        pin.write(">" + cmotif + "_" + str(temp[0]).replace("'","").replace(',',"").strip('[]') + "\n")
                        pin.write(cmotif + "\n")
                    else:
                        print "This motif is empty.", cmotifs.replace("'","").replace(',',"").strip('[]')
        filein.close()
    pin.close()

"""Read in column A and B from .xls or .xlsx file containing fasta descriptor and motif, repsectively"""
if  os.path.splitext(motifin)[-1].lower() == ".xls" or os.path.splitext(motifin)[-1].lower() == ".xlsx":
    d = pd.read_excel(motifin, parse_cols=[int(coldata1), int(coldata2)]) #leaves dataframe
    print d
    headers = []
    headers = d.columns.values
    print headers
    motifdict = d.set_index(headers[0]).to_dict()[headers[1]]
    #make a fasta file from this:
    #iterate over dict
    with open(patmaninput, "w") as pin:
        for key, val in motifdict.iteritems():
            cmotif = makeIUPAC(str(val))
            print cmotif, "This is the current motif"
            pin.write(">" + cmotif + "_" + key + "\n")
            pin.write(cmotif + "\n")
    pin.close()
