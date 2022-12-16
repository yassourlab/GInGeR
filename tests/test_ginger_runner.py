import unittest
from unittest.mock import patch
from distutils.dir_util import copy_tree
from shutil import rmtree
from click.testing import CliRunner
from ginger.ginger_runner import run_ginger_e2e
from tests.helper import get_test_data_dir, get_test_output_dir
import os
from pathlib import Path


def run_meta_or_hybrid_spades_mock(short_reads_1, short_reads_2, long_reads, output_folder, threads):
    copy_tree(f'{get_test_data_dir()}/SPAdes', output_folder)
    return output_folder


class GingerRunnerTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.data_dir = get_test_data_dir()
        cls.ginger_out_dir = f'{get_test_output_dir()}/e2e_test_output_skip_kraken'

    @classmethod
    def tearDownClass(cls):
        rmtree(cls.ginger_out_dir)
        os.remove(f'{cls.data_dir}/merged_filtered_ref_db.mmi')

    @patch('ginger.assembly_utils.run_meta_or_hybrid_spades', run_meta_or_hybrid_spades_mock)
    def test_ginger_e2e_skip_kraken(self):
        runner = CliRunner()
        result = runner.invoke(run_ginger_e2e,
                               f'{self.data_dir}/ecoli_1K_1.fq.gz data/ecoli_1K_2.fq.gz {self.data_dir}/test_gene.fasta {self.ginger_out_dir} --merged-filtered-fasta {self.data_dir}/merged_filtered_ref_db.fasta.gz --metadata-path {Path(__file__).resolve().parent.parent}/ginger/UHGG-metadata.tsv --max-species-representatives 1'.split(
                                   ' '))
        if result.exit_code != 0:
            print(str(result.exception))
        self.assertEqual(result.exit_code, 0,
                         f'{self.data_dir}/ecoli_1K_1.fq.gz data/ecoli_1K_2.fq.gz {self.data_dir}/test_gene.fasta {self.ginger_out_dir} --merged-filtered-fasta {self.data_dir}/merged_filtered_ref_db.fasta.gz --metadata-path ginger/UHGG-metadata.tsv --max-species-representatives 1')
        self.assertTrue(os.path.exists(f'{self.ginger_out_dir}/context_level_matches.csv'))
        self.assertTrue(os.path.exists(f'{self.ginger_out_dir}/species_level_matches.csv'))

        with open(f'{self.ginger_out_dir}/context_level_matches.csv', 'r') as f1, open(
                f'{self.data_dir}/context_level_matches.csv', 'r') as f2:
            lines_1 = f1.readlines()
            lines_2 = f2.readlines()
        self.assertListEqual(lines_1, lines_2)

        with open(f'{self.ginger_out_dir}/species_level_matches.csv', 'r') as f1, open(
                f'{self.data_dir}/species_level_matches.csv', 'r') as f2:
            lines_1 = f1.readlines()
            lines_2 = f2.readlines()
        self.assertListEqual(lines_1, lines_2)


if __name__ == '__main__':
    unittest.main()
