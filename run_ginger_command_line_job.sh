#!/bin/tcsh
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
conda activate evaluating_ginger
python3 -u ginger.py tests/test_files/short_reads_1.fastq.gz tests/test_files/short_reads_2.fastq.gz tests/test_files/genes.fasta -t 16