#!/bin/tcsh
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=00:05:00
conda activate evaluating_ginger
python3 -u ginger.py tests/test_files/ecoli_1K_1.fq.gz tests/test_files/ecoli_1K_2.fq.gz tests/test_files/test_gene.fasta e2e_test_output --reads-ratio-th 0.1 --max-species-representatives 10