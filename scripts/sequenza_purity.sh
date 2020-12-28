#!/bin/bash
printf "id_rep\tid_sample\tpurity\tploidy\tSLPP\n"
ls data/*/*/sequenza/R*_alternative_solutions.txt | \
while read f; do 
  awk -v OFS="\t" 'NR==2{split(FILENAME,a,"/");split(a[5],b,"_");print a[3],b[1],$0}' $f
done
