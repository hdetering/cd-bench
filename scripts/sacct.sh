# NOTE: CPUTime is not reliable. It's just the product of assigned cpus and run time.
#       For "thinnodes", assigned cpus is 24 always, since entire node gets assigned to job.
sacct -P -j $(tail -n+2 slurm_jobs.csv | cut -d"," -f1 | tr "\n" ",") --format="JobId,JobName,Start,End,TotalCPU,MaxRSS,NodeList" \
| grep "batch" \
| tr '|' ',' \
| sed 's/\.batch//' \
| tee slurm_perf.$(date +"%FT%H%M%S")
