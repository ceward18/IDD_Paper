#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -l mf=700G -q all.q -V -e ~/err -o ~/out -t 201-400 batch_estimatedE.job
qsub -pe smp 6 -cwd -l mf=700G -q all.q -V -e ~/err -o ~/out -t 1001-1200 batch_estimatedE.job
qsub -pe smp 6 -cwd -l mf=700G -q all.q -V -e ~/err -o ~/out -t 2601-2800 batch_estimatedE.job