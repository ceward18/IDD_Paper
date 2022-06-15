#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -q all.q -V -e ~/err -o ~/out -t 1-4800 batch_estimatedE.job