#!/bin/csh
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --output=/cs/labs/morani/haimasree/GInGeR_test_run/O-%x.%j.out
#SBATCH --error=/cs/labs/morani/haimasree/GInGeR_test_run/E-%x.%j.err
conda activate GInGeR-minimal
python3 -u /cs/labs/morani/haimasree/GInGeR/run_ginger_from_python.py
