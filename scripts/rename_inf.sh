#!/bin/bash

BASEDIR=$1
TOOL=$2
PFX=${3:-"$2"}

#echo "$BASEDIR $TOOL $PFX"
for f in ${BASEDIR}/*/${TOOL}/${PFX}.*.csv; do
  fn=$(basename $f)
  mv $f $(dirname $f)/inf${fn#${PFX}}
done

<<'COMMENT'
for f in data/100/*/clonefinder*/clonefinder.*.csv; do
  fn=$(basename $f)
  mv $f $(dirname $f)/inf${fn#clonefinder}
done
COMMENT
