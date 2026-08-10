import os
import tempfile
import unittest
from shutil import rmtree, copytree
from unittest.mock import patch

import pandas as pd
from click.testing import CliRunner

from ginger.ginger_runner import ginger_e2e_func, run_ginger_e2e
from tests import helper

TEST_FILES = helper.get_filedir()

# contexts are always exactly --context-len long, and the gene sits 337bp from one end of the 1000bp
# test contig, so anything longer than that leaves the test with no contexts at all
CONTEXT_LEN = 300
NO_SPECIES_MATCH_COLUMNS = ['gene', 'contig', 'gene_match_score']

def run_meta_or_hybrid_spades_mock(short_reads_1, short_reads_2, long_reads, output_folder, threads):
    copytree(f'{TEST_FILES}/SPAdes', output_folder, dirs_exist_ok=True)
    return output_folder


class GingerRunnerTest(unittest.TestCase):
    # the output dir is removed in tearDownClass, so a mismatch against the ground truth csv can only
    # be read off the failure message - print it in full rather than truncating it
    maxDiff = None

    @classmethod
    def setUpClass(cls) -> None:
        # a fresh temp dir, so a run leaves nothing behind in whatever directory it was started from
        cls.tmp_dir = tempfile.mkdtemp()
        cls.out_dir = f'{cls.tmp_dir}/e2e_test_output_skip_kraken'
        cls.references_dir = f'{cls.tmp_dir}/references_dir'
        cls.short_reads_1 = f'{TEST_FILES}/ecoli_1K_1.fq.gz'
        cls.short_reads_2 = f'{TEST_FILES}/ecoli_1K_2.fq.gz'
        cls.threads = 1
        cls.genes_path = f'{TEST_FILES}/test_gene.faa'
        cls.merged_filtered_fasta = f'{TEST_FILES}/merged_filtered_ref_db.fasta.gz'
        cls.metadata_path = helper.get_metadata_path()
        cls.max_species_representatives = 1
        cls.coverage_th = 10

    @classmethod
    def tearDownClass(cls):
        rmtree(cls.tmp_dir, ignore_errors=True)

    def setUp(self):
        # each test runs against an empty out_dir, the way a real run starts
        if os.path.exists(self.out_dir):
            rmtree(self.out_dir)

    GROUND_TRUTH_OUTPUTS = ['context_level_matches.csv', 'species_level_matches.csv']

    def _assert_outputs_match_ground_truth(self):
        """Compares every output csv against its ground truth, reporting all the mismatching files
        in one go and quoting what was produced.

        Stopping at the first failing assertion would hide the rest, and the output dir is deleted
        in tearDownClass - so when a deliberate change moves these numbers, the failure message is
        the only place the new ground truth can be read off.
        """
        mismatches = []
        for name in self.GROUND_TRUTH_OUTPUTS:
            with open(f'{self.out_dir}/{name}') as test_out, open(f'{TEST_FILES}/{name}') as gt:
                lines_out, lines_gt = test_out.readlines(), gt.readlines()
            if lines_out != lines_gt:
                mismatches.append(f'\n=== {name} - produced:\n{"".join(lines_out)}'
                                  f'=== {name} - ground truth:\n{"".join(lines_gt)}')
        self.assertEqual([], mismatches, ''.join(mismatches))

    def _assert_run_produced_the_expected_csvs(self):
        for name in self.GROUND_TRUTH_OUTPUTS:
            self.assertTrue(os.path.exists(f'{self.out_dir}/{name}'), f'{name} was not written')
        self._assert_outputs_match_ground_truth()

        # written only when some gene had no species level match, so its absence is not a failure
        no_match_csv = f'{self.out_dir}/genes_detected_in_graph_with_no_species_match.csv'
        if os.path.exists(no_match_csv):
            self.assertListEqual(list(pd.read_csv(no_match_csv).columns), NO_SPECIES_MATCH_COLUMNS)

    @patch('ginger.assembly_utils.run_meta_or_hybrid_spades', run_meta_or_hybrid_spades_mock)
    def test_ginger_e2e_func(self):
        ginger_e2e_func(None, self.short_reads_1, self.short_reads_2, self.out_dir, None, self.threads, None,
                        None, self.coverage_th, self.metadata_path, self.references_dir, self.merged_filtered_fasta,
                        self.genes_path, 12, 1.5, CONTEXT_LEN, 0.9, 0.9, ['all'], False,
                        self.max_species_representatives, False, 0.8, add_plasmid_score=False)

        self._assert_run_produced_the_expected_csvs()

    @patch('ginger.assembly_utils.run_meta_or_hybrid_spades', run_meta_or_hybrid_spades_mock)
    def test_ginger_e2e_command_skip_kraken(self):
        result = CliRunner().invoke(run_ginger_e2e, [
            self.short_reads_1, self.short_reads_2, self.genes_path, self.out_dir,
            '--sample-specific-references', self.merged_filtered_fasta,
            '--species-coverage-threshold', str(self.coverage_th),
            '--reference-genomes-metadata', self.metadata_path,
            '--max-species-representatives', '1',
            '--no-add-plasmid-score',
            '--context-len', str(CONTEXT_LEN)])

        self.assertEqual(result.exit_code, 0, str(result.exception))
        self._assert_run_produced_the_expected_csvs()


if __name__ == '__main__':
    unittest.main()
