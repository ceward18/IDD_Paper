#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -l mf=256G -q all.q -V -e ~/err -o ~/out -t 1-4400 batch_estimatedE.job
qsub -pe smp 6 -cwd -l mf=256G -q UI -V -e ~/err -o ~/out -t 1-4400 batch_estimatedE.job