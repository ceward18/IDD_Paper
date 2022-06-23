#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -l mf=700G -q all.q -V -e ~/err -o ~/out -t 390-4387 batch_estimatedE.job