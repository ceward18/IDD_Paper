#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -q UI -V -e ~/err -o ~/out -t 4039 batch_estimatedE.job
qsub -pe smp 6 -cwd -q UI -V -e ~/err -o ~/out -t 4043 batch_estimatedE.job
qsub -pe smp 6 -cwd -q UI -V -e ~/err -o ~/out -t 4050 batch_estimatedE.job
qsub -pe smp 6 -cwd -q UI -V -e ~/err -o ~/out -t 4064 batch_estimatedE.job
qsub -pe smp 6 -cwd -q UI -V -e ~/err -o ~/out -t 4069 batch_estimatedE.job