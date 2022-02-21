#!/bin/bash -l

# Run using the following command:
# qsub paper-examples.sh

# Set SCC project
#$ -P arithgeo

# Specify hard time limit for the job. 
#   The job will be aborted if it runs longer than this time.
#   The default time is 12 hours
#$ -l h_rt=12:00:00

# Require avx instruction set
#$ -l avx

# Run only on nodes with Magma license
#$ -l magma

# Join stdout and stderr streams
#$ -j y

# Specify location of stdout stream
#$ -o /projectnb/arithgeo/drhast/twocover-results/stdout

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m ea

module load sagemath/9.3
module load magma_usyd/2.26-10

outdir=/projectnb/arithgeo/drhast/twocover-results/examples

sage --python ./twocover-processor.sage.py --output_directory $outdir --stages all --label 141991.b.141991.1
sage --python ./twocover-processor.sage.py --output_directory $outdir --stages all --label 10681.a.117491.1
sage --python ./twocover-processor.sage.py --output_directory $outdir --stages all --label 7403.a.7403.1
sage --python ./twocover-processor.sage.py --output_directory $outdir --stages all --label 7211.a.7211.1

sage --python ./twocover-processor.sage.py --output_directory $outdir --stages all --label 6443.a.6443.1 --unconditional
magma "$outdir/curve-6443.a.6443.1-class-groups.m" > "$outdir/curve-6443.a.6443.1-class-groups.txt"

