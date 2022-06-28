#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 916 batch_estimatedE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 924 batch_estimatedE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 927 batch_estimatedE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 1700 batch_estimatedE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 1790 batch_estimatedE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 2596 batch_estimatedE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 3708 batch_estimatedE.job
qsub -pe smp 6 -cwd -l mf=700G -q UI -V -e ~/err -o ~/out -t 4143 batch_estimatedE.job