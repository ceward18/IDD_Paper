#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 3 batch_knownE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 11-12 batch_knownE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 15 batch_knownE.job