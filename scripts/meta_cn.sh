#!/bin/bash
for f in data/*/*/sim/bed/clone_*.bed; do 
  awk -v OFS="\t" 'BEGIN{split(ARGV[1],a,"/");r=a[3];c=gensub(/\.cn\.bed/,"","g",a[6])}{print r,c,$0}' $f
done
