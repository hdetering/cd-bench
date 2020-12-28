#!/bin/bash
# Extract job id, job name, cohort, rep from SLURM log files
grep -o "data/.../[^/ ]*" log/*.*.out | uniq | awk -v OFS="," '{split($1,a,":");split(a[1],j,".");gsub("log/slurm-","",j[1]);split(a[2],r,"/");print j[1],j[2],r[2],r[3]}'
