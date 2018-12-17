#$ -S /bin/sh
# First, activate all filesystems under /ihme
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_orig_umask=$(umask)
# ls /ihme/* 1>/dev/null
# ls /home/j 1>/dev/null
run_file="$1"; shift
singularity exec /ihme/singularity-images/hiv/hiv_11.img /usr/local/bin/R <$run_file  --no-save --args $@