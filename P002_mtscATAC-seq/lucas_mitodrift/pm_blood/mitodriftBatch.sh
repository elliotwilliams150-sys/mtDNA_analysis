#!/bin/bash
#SBATCH --account=rockhpc_spnmmd
#SBATCH --partition=default_free
#SBATCH --cpus-per-task=12 
#SBATCH --mem=90G
#SBATCH --time=10:00:00

./mitodrift.sh
