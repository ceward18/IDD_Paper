#!/bin/bash
# "-cwd": Run this job in the current working directory
qsub -pe smp 6 -cwd -q UI-HM -V -e ~/err -o ~/out -t 3 batch_knownE.job