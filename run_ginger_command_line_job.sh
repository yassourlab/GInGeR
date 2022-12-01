#!/bin/tcsh
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
conda activate evaluating_ginger
python3 -u ginger.py tests/test_files/ecoli_1K_1.fq.gz tests/test_files/ecoli_1K_2.fq.gz tests/test_files/test_gene.fasta -t 16 --reads-ratio-th 0.1