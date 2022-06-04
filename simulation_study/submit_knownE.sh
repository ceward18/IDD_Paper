#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 1-16 batch_knownE.job