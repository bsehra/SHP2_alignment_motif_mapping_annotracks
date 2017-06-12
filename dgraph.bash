#!usr/bin bash

gene=$1
dfile=$2
inpath1=/home/bsehra/scripts/Genealign/${gene}/downstream/
inpath2=/home/bsehra/scripts/Genealign/${gene}/longseq/
inpath3=/home/bsehra/scripts/Genealign/${gene}/intronic_analysis/
inpath4=/home/bsehra/scripts/Genealign/${gene}/upstream/
fastafile1=${inpath1}uAt_Al_Cr_Br_Th_${gene}_3kdown.fasta
fastafile2=${inpath2}uAt_Al_Cr_Br_Th_${gene}_long.fasta
fastafile3=${inpath3}uAt_Al_Cr_Br_Th_${gene}_introns.fasta
fastafile4=${inpath4}uAt_Al_Cr_Br_Th_${gene}_3kup.fasta

dpath=/home/bsehra/scripts/Genealign/${gene}/
#dpath=/home/bsehra/scripts/Genealign/${gene}/upstream/
#dfile=${dpath}Dialign_${gene}_allspecies_up.txt
#udfile=${dpath}uDialign_${gene}_allspecies_up.txt
#msamplerfile=${gene}msampler_w8_5motifs_100runs_M2_Ov1_BGM4.txt
#mrankingfile=${gene}mr_w8_m2_ov1_100runs_bgm4.txt
#memefoldersuffix=meme_bgm4_zoops
#weederfile=uAt_Al_Cr_Br_Th_${gene}_3kup_dustRM.fasta.wee
outpath=/home/bsehra/scripts/Genealign/${gene}/
script_path=/home/bsehra/scripts/Genealign/annofile_scripts/

mkdir ${outpath}/Dialigncheck/

<<EOF

cd ./Genealign/${gene}

rm -R dannotationfiles
mkdir dannotationfiles

rm uDialign*
EOF

#create unix files for set of Dialign alignments per gene
#bash ${script_path}filetounixformat.bash ${dfile} ${inpath4}

echo "creating anno files for ${gene}"

<<EOF

echo "downstream file for ${gene}"
python dialigngraph.py ${dpath}uDialign_${gene}_allspecies_down.txt ${fastafile1} ${outpath}${gene}dgraphanno_down.txt 1 7 ${gene} ds
echo "long file for ${gene}"
python dialigngraph.py ${dpath}uDialign_${gene}_allspecies_long.txt ${fastafile2} ${outpath}${gene}dgraphanno_long.txt 1 7 ${gene} long
echo "introns file for ${gene}"
python dialigngraph.py ${dpath}uDialign_${gene}_allspecies_introns.txt ${fastafile3} ${outpath}${gene}dgraphanno_introns.txt 1 7 ${gene} introns
EOF

echo "upstream file for ${gene}"
python ${script_path}dialigngraph_test2.py ${dfile} ${fastafile4} ${outpath}${gene}dgraphanno_up.txt 1 5 ${gene} us r

<<EOF
#to make one large annotation file
cd ${outpath}
cat *up* *introns* *down* *long* > ${gene}dgraphanno_merge.txt

for f in 
do
  echo "Processing file and converting to win compatibility..."
  awk 'sub("$", "\r")' $f > w${f}
done

mkdir winDannofiles_${gene}
mv w* winDannofiles_${gene}
EOF

<<EOF

python avid.py ${inpath}${gene}Avid_pwalign_long.txt ${fastafile} ${outpath}${gene}avidanno_long.txt ${gene}

python weederphylo.py ./Genealign/${gene}/${weederfile} ${fastafile} ${outpath}${gene}weederanno.txt ${gene}

python msamplerphylo.py ./Genealign/${gene}/${mrankingfile} ./Genealign/${gene}/${msamplerfile} ${fastafile} ${outpath}${gene}msampleranno.txt ${gene} w8M2Ov1BGM4_100runs

python memephylo.py ./Genealign/${gene}/${gene}${memefoldersuffix}/meme.txt ${fastafile} ${outpath}${gene}memeanno.txt ${gene} zoopsbgm4

cd ${outpath}

echo "Processing dgraph files: converting to win compatibility..."
awk 'sub("$", "\r")' ${gene}dgraphanno_long.txt > w${gene}dgraphanno_long.txt 
done
echo "Processing dgraph file: converting to win compatibility..."
awk 'sub("$", "\r")' ${gene}avidanno_long.txt > w$${gene}avidanno_long.txt

mkdir winannofiles${gene}_long
mv w* winannofiles${gene}_long/


echo "annotation files for Meme runs (31 gene input)"

python memeoverrep.py ./Genealign/Aiigenes/31align/memeout_1000bpBGM_500bpupstream/memeout_zoops_order4_dslrc500/meme.txt ./Genealign/Aiigenes/AIIgenes_500b_upstream_N_srlc_dust.fasta 500 zoops_BGM4.txt
python memeoverrep.py ./Genealign/Aiigenes/31align/memeout_1000bpBGM_1000bpupstream/memeout_zoops_order4_dslrc1000/meme.txt ./Genealign/Aiigenes/AIIgenes_1000b_upstream_N_srlc_dust.fasta 1000 zoops_BGM4.txt
python memeoverrep.py ./Genealign/Aiigenes/31align/memeout_1000bpBGM_3000bpupstream/memeout_zoops_order4_dslrc3000/meme.txt ./Genealign/Aiigenes/AIIgenes_3000b_upstream_N_srlc_dust.fasta 3000 zoops_BGM4.txt


cd ~/scripts/Genealign/Aiigenes/31align/memeannofiles

cd ~/scripts

echo "annotation files for Msampler runs (31gene input) : 10 runs"

python msampleroverrep.py ./Genealign/Aiigenes/31align/msampler_1000bpBGmodel_500bpup_100runs/BGM4/msampler_w8_5motifs_100runs_M2_Ov1_BGM4_500.txt ./Genealign/Aiigenes/31align/msampler_1000bpBGmodel_500bpup_100runs/MRanking_out/31genes_100runs_500bp_w8_m2_ov1_bgm4.txt ./Genealign/Aiigenes/AIIgenes_500b_upstream_N_srlc_dust.fasta 500 w8m2bgm4_100runs

python msampleroverrep.py ./Genealign/Aiigenes/31align/msampler_1000bpBGmodel_1000bpup_100runs/BGM4/msampler_w8_5motifs_100runs_M2_Ov1_BGM4_1000.txt ./Genealign/Aiigenes/31align/msampler_1000bpBGmodel_1000bpup_100runs/MRanking_out/31genes_100runs_1000bp_w8_m2_ov1_bgm4.txt ./Genealign/Aiigenes/AIIgenes_1000b_upstream_N_srlc_dust.fasta 1000 w8m2bgm4_100runs

python msampleroverrep.py ./Genealign/Aiigenes/31align/msampler_1000bpBGmodel_3000bpup_100runs/BGM4/msampler_w8_5motifs_100runs_M2_Ov1_BGM4_3000.txt ./Genealign/Aiigenes/31align/msampler_1000bpBGmodel_3000bpup_100runs/MRanking_out/31genes_100runs_3000bp_w8_m2_ov1_bgm4.txt ./Genealign/Aiigenes/AIIgenes_3000b_upstream_N_srlc_dust.fasta 3000 w8m2bgm4_100runs

cd ~/scripts/Genealign/Aiigenes/31align/msamplerannofiles

for f in `ls`
do
  echo "Processing file and converting to win compatibility..."
  awk 'sub("$", "\r")' $f > w${f}
done

mkdir msampleranno_win
mv w* msampleranno_win

mkdir msampleranno_unix
mv *.txt msampleranno_unix/

cd ~/scripts

echo "annotation files for weeder runs (31 gene input)"

python weederoverrep.py ./Genealign/Aiigenes/31align/weeder_500upstream_N_srlc/AIIgenes_500b_upstream_N_srlc_dust.fasta.wee ./Genealign/Aiigenes/AIIgenes_500b_upstream_N_srlc_dust.fasta 500 w
python weederoverrep.py ./Genealign/Aiigenes/31align/weeder_1000upstream_N_srlc/AIIgenes_1000b_upstream_N_srlc_dust.fasta.wee ./Genealign/Aiigenes/AIIgenes_1000b_upstream_N_srlc_dust.fasta 1000 w
python weederoverrep.py ./Genealign/Aiigenes/31align/weeder_3000upstream_N_srlc/approachIIgenes_3Kupstream_N_srlc_dust.fasta.wee ./Genealign/Aiigenes/AIIgenes_3000b_upstream_N_srlc_dust.fasta 3000 w

EOF
