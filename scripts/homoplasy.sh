#!/bin/bash
tool=$1
for f in data/???/*/${tool}/inf.snvs.csv; do 
  awk -F, -v OFS="," 'NR>1{a[$1]++}END{split(FILENAME,p,"/");r=p[3];h=0;for(x in a){c++;if(a[x]>1)h++}; print "'$(dirname $f)'",c,h}' $f; 
done
