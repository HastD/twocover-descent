#!/bin/bash -l

# Set SCC project
#$ -P twocover_descent

# Submit an array job with 3 tasks 
#$ -t 1-20

# Specify hard time limit for the job. 
#   The job will be aborted if it runs longer than this time.
#   The default time is 12 hours
#$ -l h_rt=12:00:00

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m ea

# Use the SGE_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from SGE_TASK_ID value

module load sagemath/9.3
module load magma_usyd/2.26-9

sage --python ./twocover-processor.sage.py --index $(($SGE_TASK_ID-1)) --output_directory /projectnb/arithgeo/drhast/twocover-results --stages setup

