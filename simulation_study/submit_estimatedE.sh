#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -q UI -V -e ~/err -o ~/out -t 1-48 batch_estimatedE.job