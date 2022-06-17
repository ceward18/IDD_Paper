#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -l mf=700G -q all.q -V -e ~/err -o ~/out -t 1-3700 batch_estimatedE.job