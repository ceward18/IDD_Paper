#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 6-16 batch_knownE.job
qsub -pe smp 6 -cwd -q UI -V -e ~/err -o ~/out -t 1-5 batch_knownE.job