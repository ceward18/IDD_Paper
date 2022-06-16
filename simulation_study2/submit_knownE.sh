#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -l mf=256G -q all.q -V -e ~/err -o ~/out -t 3 batch_knownE.job
qsub -pe smp 6 -cwd -l mf=256G -q all.q -V -e ~/err -o ~/out -t 11-16 batch_knownE.job