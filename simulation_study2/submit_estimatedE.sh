#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -l mf=700G -q all.q -V -e ~/err -o ~/out -t 391-2645 batch_estimatedE.job
qsub -pe smp 6 -cwd -l mf=700G -q all.q -V -e ~/err -o ~/out -t 4400 batch_estimatedE.job