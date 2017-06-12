#!usr/bin/env python

"""Program to create a Patman input file for pattern searching from PLACE motifs"""
#Bhupinder Sehra 6/27/13

import sys
import os

motifin = sys.argv[1]
patmaninput = sys.argv[2]

"""Main body"""
temp = []
motifs = []
i = 0
#read in tab delimited file from PLACE database
numlines = sum(l for line in open(motifin))
print "the number of lines in motif input file is:", numlines
with open (patmaninput, "w") as pfile:
    with open (motifin, "r") as filein:
        for line in filein:
            line.strip()
        if not line.strip():
            continue
        if line.startswith("Name"):
            continue
        else:
            temp = line.split()
            motifs = temp[1:2]
            #write out to stdout and file
            print str(motifs).replace("'","").replace(',',"").strip('[]')
            pfile.write(">" + str(temp[1]).replace("'","").replace(',',"").strip('[]') + "_" + str(temp[2]).replace("'","").replace(',',"").strip('[]') + "\n")
            pfile.write(str(temp[2]).replace("'","").replace(',',"").strip('[]') + "\n")
    filein.close()
pfile.close()
