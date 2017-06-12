#!usr/bin bash

#Created by Bhupinder Sehra
#Updated 1/10/17

<<EOF
This program will take the text file output (copied and pasted from the HTML page output) of a Dialign-Chaos multiple sequence alignment and produce an annotation file compatible wth TAIR Gbrowse
Command line arguments: 
'gene' = AT locus ID of gene in Arabidopsis thaliana
'dfile' = path of text file containing Dialign-Chaos output. The first line of this output should be the 1st line of the alignment
'fastafile' = path of fastafile containing Arabidopsis thaliana reference sequence
'outpath' = path of annotation track to be created
'annofile_name' = name of annotation track to be written (.txt)
'num_species_aligned' = number of species or sequences aligned 
'align_type' = 'us' for reference sequence that is on positive strand; 'ds' for reference sequence on the negative strand. 
'direction' = 'f' forward or 'r' reverse
EOF


gene=$1
dfile=$2
#Enter 
fastafile1=
#Enter path (dpath) where Dialign-Chaos text file containing the multiple sequence alignment is located
outpath=
annofile_name=
num_species_aligned= [enter number]
align_type = 
direction = 

mkdir ${outpath}/Dialigncheck/

echo "creating annotation files for ${gene}"
python ${dfile} ${fastafile1} ${outpath}${annofile_name} 1 5 ${gene} us r
echo "File created"
