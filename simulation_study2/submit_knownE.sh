#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 4090 batch_estimatedE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 4092 batch_estimatedE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 4121 batch_estimatedE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 4143 batch_estimatedE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 4161 batch_estimatedE.job
qsub -pe smp 6 -cwd -q BIOSTAT -V -e ~/err -o ~/out -t 4164 batch_estimatedE.job