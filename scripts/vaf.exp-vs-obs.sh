# vim: syntax=sh tabstop=4 expandtab
# coding: utf-8

if [[ $# -lt 1 ]]; then
  echo
  echo "usage: $0 <sample-id>"
  echo
  exit 1
fi

SAMPLE=$1
export LC_ALL=C

join \
  <(awk -v FS="," -v OFS="\t" '/^[0-9]/{print "s"$1,$2,$3,$4}' bed/${SAMPLE}.vaf.bed | sort) \
  <(awk -v OFS="\t" '/^s/{print $1,$3/$2}' bam/${SAMPLE}.vars.csv | sort) \
| sort -k3n | tr " " "\t" \
| awk -v OFS="\t" -v id=$SAMPLE '{print id,$0,$4-$5}'
