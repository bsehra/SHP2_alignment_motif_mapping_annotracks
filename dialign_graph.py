#!usr/bin/env python

#script to convert Dialign Chaos alignment to annotation track for viewing in GBrowse (TAIR10)
#for alignments which align multiple sequences to an Arabidopsis thaliana sequence which must be first in the fasta file input for Dialign Chaos and the alignment
#Written by Bhupinder Sehra; updated 8/14/16

import string
import sys
import os
import re


#dfilein = dialign alignment; ffilein = fasta file input; fout = file to write annotation track output; atlinenum = line of the alignment which belongs to A thaliana (this should be the 1st line)
#spnum = number of fasta sequences being aligned; genename = full AT locus ID in lower or uppercase; aligntype: "us" for upstream; "ds" for downstream; "long" for large genomic region, "introns" for coding sequence of gene to find alignment across introns
#gene dir = gene direction: enter 'r' for gene on negative strand, 'f' on positive strand.
dfilein = sys.argv[1]
ffilein = sys.argv[2]
fout = sys.argv[3]
atlinenum = int(sys.argv[4])
spnum = int(sys.argv[5])
genename = sys.argv[6]
aligntype = str(sys.argv[7])
gdir = str(sys.argv[8])

def genedir(dir, atype):
    """Search genefile for direction of gene ie 'f' or 'r' - returns string - works and tested. ETA: gene dir parsed changes as alignment type changes
    For 'f' gene, upstream regions are processed from R to L when mapping dialign-chaos scores to fasta files; 'long', 'intronic' and downstream regions
    are processed from L to R. For a reverse gene upstream regions aligned are processed from L to R while the reverse is true for other types of regions aligned."""
    process_dir = ""
    if dir == 'f':
        if atype == "us":
            process_dir = "f"
        else:
            process_dir = "r"
    if dir == 'r':
        if atype == "ds":
            process_dir = "f"
        else:
            process_dir = "r"
    return process_dir

def chrnum(atname):
    """take in AT ref of gene and return chromosome on which gene is based, according to A thaliana, TAIR10 version. - tested and works"""
    cnum = atname[2]
    return cnum

gfilein = "./Genealign/infiles/AIIgenedirection.txt"
atdirect = genedir(gdir, aligntype)

def getchrstartpos(fin, sdir):
    """method to get chromosome position from Fasta file input into Dialign-Chaos. The first sequence should belong to Athaliana. The fasta descriptor should be in the format as output from TAIR: >ChrX:123456789..234567890 etc with a space after the second chromosome number. Upstream regions of 'reverse' genes (on the negative strand), intronic sequences and longer genomic regions are not reverse complemented while upstream regions of 'forward' genes (on the positive strand) and d/s of reverse genes, sequences are reverse complemented. This method ensures that 'chrs' has the lowest value and 'chre' the highest - tested and works"""
    descriptor = ""
    temp = []
    i = 0
    fastain = open(fin, "r")
    for line in fastain:
        if not line.strip():
            continue
        i += 1
        if i == 1:
            if line.startswith(">"):
                line.strip()
                descriptor = `line`.replace(":"," ").replace(".."," ").replace("_"," ").replace('"',"").strip("\n")
                fulldescrip = line
                temp = descriptor.split()
        if i > 1:
            fastain.seek(-1, 2)
    fastain.close()
    chrs = 0
    chre = 0
    if sdir == "f":
        if int(temp[1]) > int(temp[2]):
            #reverse complemnted u/s forward seq or d/s reverse seq and order of chr positions goes from high to low.
            chrs = int(temp[2])
            chre = int(temp[1])
        if int(temp[1]) < int(temp[2]):
            chrs = int(temp[1])
            chre = int(temp[2])
    if sdir == "r":
        chrs = int(temp[1])
        chre = int(temp[2])
    chrlist = [chrs, chre, fulldescrip]
    return chrlist
    
def isgap(s):
    """take in string, 's' and check using isspace(), which returns True if s is only whitespace. 
    This method is used to filter() our lists of spaces if function returns False."""
    if s.isspace() == True:
        return False
    else:
        return True

def ischar(s):
    """take in string, 's', if '-' in s return False, else return True. 
    This method is used to filter our A thaliana seq list"""
    if '-' in `s`:
        return False
    else:
        return True

def findgap(arr):
    """Take in list (dialign score weightings), loop over and find position of '-' in athlist. Returns list"""
    poslist = []
    for y, z in enumerate(arr):
        if ischar(`z`) == False:
            poslist.append(y)
        else:
            continue
    return poslist

def checklen(list1, list2):
    if len(list1) == len(list2):
        return True
    else:
        return False

#y = holder of 'x' element as list is iterated over - shows what score remains constant over a stretch
#z = counter to add to chrs as iteration occurs - shows how many of 'y' there are
def calcsw(scorelist, gdir, fastain, dcpath):
    """take in scorearray, gene direction and the fasta file containing chromosome positions. Calculate chromosomal position as we iterate over array and the frequency of each score.
    Outputs a check file written to the same path as input Dialign alignment"""
    f = "Dcheck" + `genename`.replace("'","") + aligntype + ".txt"
    dcheckfile = os.path.join(dcpath, f)
    dfile = open(dcheckfile, "w")
    y = scorelist[0]
    #freq score counter
    z = 0
    chrposk = getchrstartpos(fastain, gdir)
    print "this is start and end chromosome number: ", chrposk
    cwlist = []
    if gdir == 'r':
        chrpos = chrposk[0]
        for i, x in enumerate(scorelist):
            if y == x:
                dfile.write("scorearray element = " + `y` + "\n")
                z += 1
            else:
                dfile.write("freq of this element = " + `z` + "\n")
                temptuple = (y, chrpos, chrpos+z-1)
                #print temptuple
                dfile.write("chr positions = " + `temptuple` + "\n")
                cwlist.append(temptuple)
                chrpos = chrpos + z
                #print "This is z: ", z
                y = x
                z = 1
        temptuple = (y, chrpos, chrpos+z-1)
        cwlist.append(temptuple)
        dfile.write("final scorearray element = " + `y` + "\n")
        dfile.write("freq of final element " + `z` + "\n")
        dfile.write("chr positions = " + `temptuple` + "\n")
        #print y
        #print z
        #print chrpos
        #print chrpos+z-1
    if gdir == 'f':
        chrpos = chrposk[1]
        for i, x in enumerate(scorelist):
            if y == x:
                dfile.write("scorearray element = " + `y` + "\n")
                z += 1
            else:
                dfile.write("freq of this element = " + `z` + "\n")
                temptuple = (y, chrpos, chrpos-z+1)
                dfile.write("chr positions = " + `temptuple` + "\n")
                cwlist.append(temptuple)
                chrpos = chrpos - z
                y = x
                z = 1
        temptuple = (y, chrpos, chrpos-z+1)
        cwlist.append(temptuple)
        dfile.write("final scorearray element = " + `y` + "\n")
        dfile.write("freq of final element " + `z` + "\n")
        dfile.write("chr positions = " + `temptuple` + "\n")
        #print y
        #print z
        #print chrpos
        #print chrpos-z+1
    dfile.close()
    return cwlist

#Read in Dialign alignment: ignores lines at the start and end of alignment
athstr = ""
scorestr = ""
atdescrip = ""
atharray = []
scorearray = []
chrinfo = []
chrinfo = getchrstartpos(ffilein, gdir)
print chrinfo
atdescrip = chrinfo[2]
print atdescrip, "this is atdescrip"
slinenum = atlinenum + spnum
tlist  = []
with open (dfilein, "r") as filein:
    for line in filein:
        if line.startswith("Sequence tree") or line.startswith("Tree constructed"):
            break
        if line.startswith(atdescrip[1:10]) and re.match("^[acgtACGT-]*$", "".join(line.split()[2:])):
                print line,  "".join(line.strip().split()[2:])
                athstr += "".join(line.strip().split()[2:])
        if "".join(line.split()).isdigit():
            print line 
            scorestr += line.replace("\n","").replace("\r","").replace(" ","").replace("'","").replace(",","")
        else:
            continue
filein.close()

print "this is the athstr: ", athstr, " and its len: ", len(athstr)
print "this is the score string: ",scorestr, "and its len: ", len(scorestr)
atharray = list(athstr)
scorearray = list(scorestr)
sarray = []
g = open(fout, "a")
print len(atharray)
print len(scorearray)
if checklen(atharray, scorearray) == True:
    print "both lists before removing gaps are equal in length"
else:
    print "Lists are unequal in length; an error has occurred in parsing the alignment"
gappos = []
gappos = findgap(atharray)
atharray = filter(ischar, atharray)
print "This is the ath string after gaps removed:" + "\n"
print "".join(`atharray`).strip('[]').replace("'","").replace(",","").replace(" ","")
for a, x in enumerate(scorearray):
    if a not in gappos:
        sarray.append(x)
    else:
        continue
print sarray, "\n"
print "this is the scorestring with gaps removed: " + "\n"
print "".join(`sarray`).strip('[]').replace("'","").replace(",","").replace(" ","")
if checklen(atharray, sarray) == True:
    print "both lists before after removing gaps are equal in length"
else:
    print "Lists are unequal; a problem has occured in parsing the alignment!!"
print len(atharray)
print len(sarray)
xsome = []
xsome = getchrstartpos(ffilein, atdirect)
dcheckpath = os.path.dirname(dfilein)
resultlist = calcsw(sarray, atdirect, ffilein, dcheckpath)
cref = chrnum(genename)
print "this is the chromosome number: ", cref

#write out annotation file = fout; the "key must match the line descriptor
dfout = open(fout, "w")
dfout.write("#" + genename + "_dialign" + aligntype + "\n")
dfout.write("[score]" + "\n")
dfout.write("glyph = xyplot" + "\n")
dfout.write("graph_type=boxes" + "\n")
dfout.write("fgcolor=black" + "\n")
if aligntype == "long":
    dfout.write("bgcolor=black" + "\n")
if aligntype == "introns":
    dfout.write("bgcolor=darkblue" + "\n")
if aligntype == "us":
    dfout.write("bgcolor=mediumblue" + "\n")
if aligntype == "ds":
    dfout.write("bgcolor=lightblue" + "\n")
dfout.write("height=100" + "\n")
dfout.write("min_score=0" + "\n")
dfout.write("max_score=9" + "\n")
dfout.write("label=1" + "\n")
dfout.write("key=Dialign_" + aligntype + "\n")
dfout.write("\n")
dfout.write("\n")
dfout.write("reference=Chr" + cref + "\n")
 
for i, x in enumerate(resultlist):
                dfout.write("score" + "    " + "Dialign_" + aligntype  + "   " + `x[1]`.strip('[]').replace(",","").replace("'","") + ".." + `x[2]`.strip('[]').replace(",","").replace("'","") + "    " + "score=" + `x[0]`.strip('[]').replace(",","").replace("'","") + "\n")  

dfout.write("\n")
dfout.write("\n")

dfout.close()
