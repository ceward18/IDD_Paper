#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 268 batch_estimatedE.job