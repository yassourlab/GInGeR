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

    def test_get_species_coverage_stats(self):
        bracken_output_path = f'{TEST_FILES}/bracken_out_coverage_test.txt'
        metadata_path = f'{TEST_FILES}/uhgg_metadata_coverage_test.tsv'
        avg_sum = 200.0
        stats = rdu.get_species_coverage_stats(
            bracken_output_path, avg_sum, metadata_path, max_refs_per_species=2,
        )
        stats_by_name = stats.set_index('name')
        # EC: median(4000000, 4000000) = 4000000; coverage = 250000*200/4000000 = 12.5
        self.assertEqual(stats_by_name.loc['Enterobacter cloacae', 'estimated_genome_length'], 4000000)
        self.assertAlmostEqual(stats_by_name.loc['Enterobacter cloacae', 'estimated_coverage'], 12.5)
        # SE: median(5000000, 5000000) = 5000000; coverage = 200000*200/5000000 = 8.0
        self.assertEqual(stats_by_name.loc['Salmonella enterica', 'estimated_genome_length'], 5000000)
        self.assertAlmostEqual(stats_by_name.loc['Salmonella enterica', 'estimated_coverage'], 8.0)

    def test_get_distinct_minimizers_by_species(self):
        kraken_report_path = f'{TEST_FILES}/kraken_report_filter_test.txt'
        distinct_minimizers = rdu.get_distinct_minimizers_by_species(kraken_report_path)
        self.assertEqual(distinct_minimizers, {'SpeciesPass': 300000, 'SpeciesFail': 50000})

    def test_get_species_included_in_analysis_df(self):
        bracken_output_path = f'{TEST_FILES}/bracken_out_coverage_test.txt'
        kraken_report_path = f'{TEST_FILES}/kraken_report_minimizers_test.txt'
        metadata_path = f'{TEST_FILES}/uhgg_metadata_coverage_test.tsv'
        avg_sum = 200.0
        top_species = ['Enterobacter cloacae', 'Yersinia enterocolitica']
        included = rdu.get_species_included_in_analysis_df(
            bracken_output_path, kraken_report_path, avg_sum, metadata_path,
            max_refs_per_species=2, top_species=top_species,
        )
        self.assertEqual(set(included['name']), set(top_species))
        included_by_name = included.set_index('name')
        self.assertEqual(included_by_name.loc['Enterobacter cloacae', 'distinct_minimizers'], 400000)
        self.assertEqual(included_by_name.loc['Yersinia enterocolitica', 'distinct_minimizers'], 300000)
        self.assertAlmostEqual(included_by_name.loc['Enterobacter cloacae', 'estimated_coverage'], 12.5)

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
