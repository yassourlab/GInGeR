#!/bin/tcsh
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
source /cs/usr/netta.barak/.tcshrc
conda activate evaluating_ginger
python3 -u compare_ginger_to_baseline.py