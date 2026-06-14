import unittest
from unittest.mock import patch
from shutil import rmtree, copytree
from click.testing import CliRunner
from ginger.ginger_runner import ginger_e2e_func, run_ginger_e2e
import os

from tests import helper

TEST_FILES = helper.get_filedir()

SPADES_OUTPUT = f'{TEST_FILES}/SPAdes'

"""
This test does not pass locally, but should pass on githubs CI.
When running it locally you will get the following error:
"FileNotFoundError: [Errno 2] No such file or directory: 'ginger/UHGG-metadata.tsv'"
If you do wish for it to pass locally (is a reasonable request), replace 'ginger/UHGG-metadata.tsv'
with '../ginger/UHGG-metadata.tsv' 
"""


def run_meta_or_hybrid_spades_mock(short_reads_1, short_reads_2, long_reads, output_folder, threads):
    copytree(f'{TEST_FILES}/SPAdes', output_folder, dirs_exist_ok=True)
    return output_folder


class GingerRunnerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.out_dir = 'e2e_test_output_skip_kraken'
        cls.short_reads_1 = f'{TEST_FILES}/ecoli_1K_1.fq.gz'
        cls.short_reads_2 = f'{TEST_FILES}/ecoli_1K_2.fq.gz'
        cls.threads = 1
        cls.genes_path = f'{TEST_FILES}/test_gene.faa'
        cls.merged_filtered_fasta = f'{TEST_FILES}/merged_filtered_ref_db.fasta.gz'
        cls.metadata_path = 'ginger/UHGG-metadata.tsv'  # for running on github CI
        # cls.metadata_path = '../ginger/UHGG-metadata.tsv'  # for running locally
        cls.max_species_representatives = 1
        cls.coverage_th = 10
        if os.path.exists(cls.out_dir):
            rmtree(cls.out_dir)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.out_dir):
            rmtree(cls.out_dir)

    @patch('ginger.assembly_utils.run_meta_or_hybrid_spades', run_meta_or_hybrid_spades_mock)
    def test_ginger_e2e_func(self):
        ginger_e2e_func(None, self.short_reads_1, self.short_reads_2, self.out_dir, None, self.threads, None,
                        None, self.coverage_th, self.metadata_path, 'references_dir', self.merged_filtered_fasta,
                        self.genes_path, 12, 1.5, 0, 2500, 0.9, 0.9, ['all'], False, self.max_species_representatives, False, 0.8,
                        add_plasmid_score=False)

        self.assertTrue(os.path.exists(f'{self.out_dir}/context_level_matches.csv'))
        self.assertTrue(os.path.exists(f'{self.out_dir}/species_level_matches.csv'))

        with open(f'{self.out_dir}/context_level_matches.csv', 'r') as test_out, open(
                f'{TEST_FILES}/context_level_matches.csv', 'r') as gt:
            lines_out = test_out.readlines()
            print(lines_out)
            lines_gt = gt.readlines()
            self.assertEqual(len(lines_out), len(lines_gt))
            self.assertEqual(len(lines_out[1].split('/')), len(lines_gt[1].split('/')))
            for field, field_gt in zip(lines_out[1].split(','), lines_gt[1].split(',')):
                self.assertEqual(field, field_gt)

            self.assertListEqual(lines_out, lines_gt)

        with open(f'{self.out_dir}/species_level_matches.csv', 'r') as test_out, open(
                f'{TEST_FILES}/species_level_matches.csv', 'r') as gt:
            lines_out = test_out.readlines()
            lines_gt = gt.readlines()
        self.assertListEqual(lines_out, lines_gt)

        no_match_csv = f'{self.out_dir}/genes_detected_in_graph_with_no_species_match.csv'
        if os.path.exists(no_match_csv):
            import pandas as pd
            df = pd.read_csv(no_match_csv)
            self.assertListEqual(list(df.columns), ['gene', 'contig', 'gene_match_score'])

    # @patch('ginger.assembly_utils.run_meta_or_hybrid_spades', run_meta_or_hybrid_spades_mock)
    # def test_ginger_e2e_func_with_subspecies_output(self):
    #     ginger_e2e_func(None, self.short_reads_1, self.short_reads_2, self.out_dir, None, 2,
    #                     None, self.read_ratio_th, f'{TEST_FILES}/e_coli_metadata.tsv', 'references_dir',
    #                     None, self.genes_path, 12, 1.5, 0,
    #                     2500, 0.9, 0.9, False, 3)
    
    #     self.assertTrue(os.path.exists(f'{self.out_dir}/context_level_matches.csv'))
    #     self.assertTrue(os.path.exists(f'{self.out_dir}/species_level_matches.csv'))

    @patch('ginger.assembly_utils.run_meta_or_hybrid_spades', run_meta_or_hybrid_spades_mock)
    def test_ginger_e2e_command_skip_kraken(self):
        runner = CliRunner()

        result = runner.invoke(run_ginger_e2e,
                               f'{self.short_reads_1} {self.short_reads_2} {self.genes_path} {self.out_dir} --sample-specific-references {self.merged_filtered_fasta} --species-coverage-threshold {self.coverage_th} --reference-genomes-metadata {self.metadata_path} --max-species-representatives 1 --no-add-plasmid-score'.split(
                                   ' '))
        self.assertEqual(result.exit_code, 0, str(result.exception))
        self.assertTrue(os.path.exists(f'{self.out_dir}/context_level_matches.csv'))
        self.assertTrue(os.path.exists(f'{self.out_dir}/species_level_matches.csv'))

        with open(f'{self.out_dir}/context_level_matches.csv', 'r') as test_out, open(
                f'{TEST_FILES}/context_level_matches.csv', 'r') as gt:
            lines_out = test_out.readlines()
            lines_gt = gt.readlines()
        self.assertListEqual(lines_out, lines_gt)

        with open(f'{self.out_dir}/species_level_matches.csv', 'r') as test_out, open(
                f'{TEST_FILES}/species_level_matches.csv', 'r') as gt:
            lines_out = test_out.readlines()
            lines_gt = gt.readlines()
        self.assertListEqual(lines_out, lines_gt)

        no_match_csv = f'{self.out_dir}/genes_detected_in_graph_with_no_species_match.csv'
        if os.path.exists(no_match_csv):
            import pandas as pd
            df = pd.read_csv(no_match_csv)
            self.assertListEqual(list(df.columns), ['gene', 'contig', 'gene_match_score'])


if __name__ == '__main__':
    unittest.main()
