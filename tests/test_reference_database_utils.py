import unittest
import tempfile
import os
from unittest.mock import patch
import pandas as pd
from ginger import reference_database_utils as rdu
from tests import helper

TEST_FILES = helper.get_filedir()

class MyTestCase(unittest.TestCase):
    def test_get_paired_reads_seqkit_stats_parsing(self):
        reads_1 = f'{TEST_FILES}/ecoli_1K_1.fq.gz'
        reads_2 = f'{TEST_FILES}/ecoli_1K_2.fq.gz'

        class Dummy:
            def __init__(self, stdout, stderr='', returncode=0):
                self.stdout = stdout
                self.stderr = stderr
                self.returncode = returncode

        fake_out = (
            'file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n'
            f'{reads_1}\tFASTQ\tDNA\t1,000\t100000\t100\t100.0\t151\n'
            f'{reads_2}\tFASTQ\tDNA\t1,000\t100000\t100\t100.0\t150\n'
        )

        with patch('ginger.reference_database_utils.run', return_value=Dummy(fake_out)):
            a1, m1, a2, m2 = rdu.get_paired_reads_seqkit_stats(reads_1, reads_2)
        self.assertEqual(a1, 100.0)
        self.assertEqual(m1, 151)
        self.assertEqual(a2, 100.0)
        self.assertEqual(m2, 150)
    def test_get_kmer_length_options(self):
        kraken_db = f'{TEST_FILES}/fake_kraken_dir'
        kmer_length_options = rdu.get_kmer_length_options(kraken_db)
        self.assertEqual(set(kmer_length_options), set([100, 150]))

    def test_get_species_passing_coverage_threshold(self):
        bracken_output_path = f'{TEST_FILES}/bracken_out_coverage_test.txt'
        metadata_path = f'{TEST_FILES}/uhgg_metadata_coverage_test.tsv'
        avg_sum = 200.0
        passing = rdu.get_species_passing_coverage_threshold(
            bracken_output_path, avg_sum, metadata_path,
            max_refs_per_species=2,
            species_coverage_threshold=10,
        )

        # Coverage calculations:
        # EC: 250000*200/4000000 = 12.5 (pass)
        # SE: 200000*200/5000000 = 8.0 (fail)
        # YE: 200000*200/3000000 = 13.33.. (pass)
        self.assertEqual(set(passing), set(['Enterobacter cloacae', 'Yersinia enterocolitica']))

    def test_filter_kraken_report_by_distinct_kmer_count(self):
        kraken_report_path = f'{TEST_FILES}/kraken_report_filter_test.txt'
        with tempfile.TemporaryDirectory() as tmpdir:
            filtered_path = os.path.join(tmpdir, 'filtered_report.txt')
            rdu.filter_kraken_report_by_distinct_kmer_count(kraken_report_path, filtered_path, threshold=200000)
            filtered = pd.read_csv(filtered_path, sep='\t', header=None, names=rdu.KRAKEN_REPORT_COLS)

        # non-species rows (U, R, G) are kept regardless of distinct_kmer_count
        self.assertEqual(set(filtered.loc[filtered['rank'] != 'S', 'taxid']), {0, 1, 100})
        # only the species row above the threshold survives
        self.assertEqual(set(filtered.loc[filtered['rank'] == 'S', 'taxid']), {1001})


if __name__ == '__main__':
    unittest.main()
