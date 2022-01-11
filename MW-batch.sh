#!/bin/bash -l

# Run using the following command:
# qsub -t 1-[number of jobs] MW-batch.sh

# Set SCC project
#$ -P arithgeo

# Specify hard time limit for the job. 
#   The job will be aborted if it runs longer than this time.
#   The default time is 12 hours
#$ -l h_rt=2:00:00

# Require avx instruction set
#$ -l avx

# Exclude Haswell architecture CPUs, which are having problems with Magma
#$ -l cpu_arch=!haswell

# Join stdout and stderr streams
#$ -j y

# Specify location of stdout stream
#$ -o /projectnb/arithgeo/drhast/twocover-results/stdout

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m ea

module load sagemath/9.3
module load magma_usyd/2.26-9

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value

index=$(( $SGE_TASK_ID - 1 ))
sage --python ./twocover-processor.sage.py --index $index --output_directory /projectnb/arithgeo/drhast/twocover-results --stages search,locsolv,ainv,mw,reduce

