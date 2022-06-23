#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -q UI -V -e ~/err -o ~/out -t 4164 batch_estimatedE.job
qsub -pe smp 6 -cwd -q UI -V -e ~/err -o ~/out -t 4161 batch_estimatedE.job
qsub -pe smp 6 -cwd -q UI -V -e ~/err -o ~/out -t 4154 batch_estimatedE.job
qsub -pe smp 6 -cwd -q UI -V -e ~/err -o ~/out -t 4150 batch_estimatedE.job