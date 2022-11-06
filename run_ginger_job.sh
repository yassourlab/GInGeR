#!/bin/tcsh
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
## conda env create -f /sci/labs/morani/morani/icore-data/lab/Projects/hybrid_assembler/GInGeR/ginger_conda_env.yml
conda activate evaluating_ginger
python3 -u /sci/labs/morani/morani/icore-data/lab/Projects/hybrid_assembler/GInGeR/run_ginger_from_python.py