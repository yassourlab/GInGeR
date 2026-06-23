import math
import unittest
import os
from types import SimpleNamespace
import pandas as pd
import pyfastg
from ginger import pipeline_utils as pu
from ginger.matches_classes import InOutPathsMatch

from tests import helper

TEST_FILES = helper.get_filedir()

"""
This test does not pass locally, but should pass on githubs CI.
When running it locally you will get the following error:
"FileNotFoundError: [Errno 2] No such file or directory: 'ginger/UHGG-metadata.tsv'"
If you do wish for it to pass locally (is a reasonable request), replace "cls.metadata_path = 'ginger/UHGG-metadata.tsv'"
with cls.metadata_path = '../ginger/UHGG-metadata.tsv' 
"""


class PipelineUtilsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.context_level_output_path = f'{TEST_FILES}/context_level_matches.csv'
        cls.context_level_output_path_with_dups = f'{TEST_FILES}/context_level_matches_with_dups.csv'  # dups means that the same gene is matched to the same genome twice
        cls.context_level_output_path_with_plasmid_score = f'{TEST_FILES}/context_level_matches_with_plasmid_score.csv'
        cls.context_level_output_path_with_context_species_diversity = f'{TEST_FILES}/context_level_matches_context_species_diversity.csv'
        cls.context_level_output_path_with_context_species_confidence_score = f'{TEST_FILES}/context_level_matches_context_species_confidence_score.csv'
        cls.context_species_diversity_metadata_path = f'{TEST_FILES}/context_species_diversity_metadata.tsv'
        cls.metadata_path = 'ginger/UHGG-metadata.tsv' # use this for running tests on github CI
        # cls.metadata_path = '../ginger/UHGG-metadata.tsv' # use this for running the test locally
        cls.species_level_output_path = 'test_species_level_matches.csv'
        cls.species_level_output_path_with_dups = f'test_species_level_matches_with_dups.csv'  # dups means that the same gene is matched to the same genome twice
        cls.species_level_output_path_with_plasmid_score = 'test_species_level_matches_with_plasmid_score.csv'
        cls.species_level_output_path_with_context_species_diversity = 'test_species_level_matches_context_species_diversity.csv'
        cls.species_level_output_path_with_context_species_confidence_score = 'test_species_level_matches_context_species_confidence_score.csv'
        # GT files
        cls.species_level_output_path_gt = f'{TEST_FILES}/species_level_matches.csv'
        cls.species_level_output_path_with_dups_gt = f'{TEST_FILES}/species_level_matches_with_dups.csv'

    @classmethod
    def tearDownClass(cls) -> None:
        if os.path.exists(cls.species_level_output_path):
            os.remove(cls.species_level_output_path)
        if os.path.exists(cls.species_level_output_path_with_plasmid_score):
            os.remove(cls.species_level_output_path_with_plasmid_score)
        if os.path.exists(cls.species_level_output_path_with_context_species_diversity):
            os.remove(cls.species_level_output_path_with_context_species_diversity)
        if os.path.exists(cls.species_level_output_path_with_context_species_confidence_score):
            os.remove(cls.species_level_output_path_with_context_species_confidence_score)

    def test_aggregate_context_level_output_to_species_level_output_and_write_csv(self):
        pu.aggregate_context_level_output_to_species_level_output_and_write_csv(self.context_level_output_path,
                                                                                self.metadata_path,
                                                                                self.species_level_output_path, 1)
        self.assertTrue(os.path.exists(self.species_level_output_path))

        with open(self.species_level_output_path, 'r') as out_f, open(self.species_level_output_path_gt, 'r') as gt_f:
            out_lines = out_f.readlines()
            gt_lines = gt_f.readlines()
        self.assertListEqual(out_lines, gt_lines)

        pu.aggregate_context_level_output_to_species_level_output_and_write_csv(
            self.context_level_output_path_with_dups,
            self.metadata_path,
            self.species_level_output_path_with_dups, 1)

        with open(self.species_level_output_path_with_dups, 'r') as out_f, open(
                self.species_level_output_path_with_dups_gt, 'r') as gt_f:
            out_lines = out_f.readlines()
            gt_lines = gt_f.readlines()
        self.assertListEqual(out_lines, gt_lines)

    def test_aggregate_context_level_output_to_species_level_output_and_write_csv_with_plasmid_score(self):
        species_level_output = pu.aggregate_context_level_output_to_species_level_output_and_write_csv(
            self.context_level_output_path_with_plasmid_score,
            self.metadata_path,
            self.species_level_output_path_with_plasmid_score, 1)

        self.assertTrue(os.path.exists(self.species_level_output_path_with_plasmid_score))
        self.assertIn('plasmid_score_mean', species_level_output.columns)
        self.assertIn('plasmid_score_most_common_context', species_level_output.columns)

        row = species_level_output.loc[('test_gene', 'Escherichia coli_D')]
        # plasmid_score_mean averages over the gene's two unique contexts (0.9 and 0.3)
        self.assertAlmostEqual(row['plasmid_score_mean'], 0.6)
        # plasmid_score_most_common_context is the score of the context matched to 2 genomes (0.9), not 1 (0.3)
        self.assertAlmostEqual(row['plasmid_score_most_common_context'], 0.9)

        # 2 distinct (in_context, out_context) trios: ctx1 and ctx2
        self.assertEqual(row['n_contexts'], 2)

    def test_aggregate_context_level_output_to_species_level_output_and_write_csv_with_context_species_diversity(self):
        species_level_output = pu.aggregate_context_level_output_to_species_level_output_and_write_csv(
            self.context_level_output_path_with_context_species_diversity,
            self.metadata_path,
            self.species_level_output_path_with_context_species_diversity, 1)

        self.assertTrue(os.path.exists(self.species_level_output_path_with_context_species_diversity))

        row = species_level_output.loc[('gene2', 'Escherichia coli_D')]
        # 5 distinct (in_context, out_context) trios: T1, T2, T3, T4, T6
        self.assertEqual(row['n_contexts'], 5)

    def test_aggregate_context_level_output_to_species_level_output_and_write_csv_with_context_species_confidence_score(self):
        species_level_output = pu.aggregate_context_level_output_to_species_level_output_and_write_csv(
            self.context_level_output_path_with_context_species_confidence_score,
            self.metadata_path,
            self.species_level_output_path_with_context_species_confidence_score, 1)

        self.assertTrue(os.path.exists(self.species_level_output_path_with_context_species_confidence_score))

        row = species_level_output.loc[('gene2', 'Escherichia coli_D')]
        # 5 distinct (in_context, out_context) trios (T1=0.1, T2=0.2, T3=0.3, T4=0.4, T6=0.5);
        # R1/R2 share T1 and R4/R5 share T4, so those duplicate rows must not be double-counted
        self.assertAlmostEqual(row['species_confidence_score'], 0.3)

    def test_compute_context_species_diversity(self):
        results_df = pd.DataFrame({
            'gene': ['g1', 'g1', 'g1', 'g1'],
            'in_context': ['in1', 'in1', 'in1', 'in2'],
            'out_context': ['out1', 'out1', 'out1', 'out2'],
            'Genome': ['G1', 'G2', 'G4', 'G3'],
            'species': ['species_A', 'species_A', 'species_B', 'species_A'],
        })
        species_reference_counts = pd.Series({'species_A': 2, 'species_B': 1})

        result = pu.compute_context_species_diversity(results_df, species_reference_counts)

        # trio (in1, out1) matches 2 genomes of species_A and 1 of species_B. Corrected counts
        # (2/2=1, 1/1=1) are equal, so each species gets probability 0.5 -> H' = ln(2)
        for _, row in result[result['in_context'] == 'in1'].iterrows():
            self.assertAlmostEqual(row['context_species_diversity'], math.log(2))

        # trio (in2, out2) matches only species_A -> single species -> H' = 0
        for _, row in result[result['in_context'] == 'in2'].iterrows():
            self.assertAlmostEqual(row['context_species_diversity'], 0.0)

    def test_compute_context_species_confidence_score(self):
        results_df = pd.DataFrame({
            'gene': ['g1', 'g1', 'g1', 'g1'],
            'in_context': ['in1', 'in1', 'in1', 'in2'],
            'out_context': ['out1', 'out1', 'out1', 'out2'],
            'Genome': ['G1', 'G2', 'G4', 'G3'],
            'species': ['species_A', 'species_A', 'species_B', 'species_A'],
        })
        species_reference_counts = pd.Series({'species_A': 2, 'species_B': 1})

        result = pu.compute_context_species_confidence_score(results_df, species_reference_counts)

        # trio (in1, out1) matches 2 genomes of species_A and 1 of species_B. Corrected counts
        # (2/2=1, 1/1=1) are equal, so each species gets confidence 0.5
        for _, row in result[result['in_context'] == 'in1'].iterrows():
            self.assertAlmostEqual(row['context_species_confidence_score'], 0.5)

        # trio (in2, out2) matches only species_A -> single species -> confidence 1.0
        for _, row in result[result['in_context'] == 'in2'].iterrows():
            self.assertAlmostEqual(row['context_species_confidence_score'], 1.0)

    def test_write_context_level_output_to_csv_context_species_diversity_and_confidence_score(self):
        def make_match(gene, in_query_name, out_query_name):
            in_path = SimpleNamespace(query_name=in_query_name, strand='+', ref_genome_start=100, ref_genome_end=200)
            out_path = SimpleNamespace(query_name=out_query_name, strand='+', ref_genome_start=300, ref_genome_end=400)
            return InOutPathsMatch(in_path, out_path, start=10, end=20, gap_ratio=0.0, score=0.9, gene_length=100,
                                   gene=gene, ref_genome='ref', gene_match_score=0.9, in_context_score=0.9,
                                   out_context_score=0.9)

        # trio (g1_in1, g1_out1) matches G1 and G2 (species_A) and G4 (species_B)
        # trio (g1_in2, g1_out2) matches only G3 (species_A)
        output = {
            ('g1', 'G1_1'): [make_match('g1', 'g1_in1', 'g1_out1')],
            ('g1', 'G2_1'): [make_match('g1', 'g1_in1', 'g1_out1')],
            ('g1', 'G4_1'): [make_match('g1', 'g1_in1', 'g1_out1')],
            ('g1', 'G3_1'): [make_match('g1', 'g1_in2', 'g1_out2')],
        }

        csv_path = 'test_context_level_matches_context_species_diversity_out.csv'
        try:
            # species_A has 3 genomes in the metadata, capped to 2 by max_species_representatives
            pu.write_context_level_output_to_csv(output, csv_path, self.context_species_diversity_metadata_path, 2)

            results_df = pd.read_csv(csv_path)
            diversity_by_contig = results_df.set_index('reference_contig')['context_species_diversity']
            confidence_by_contig = results_df.set_index('reference_contig')['context_species_confidence_score']

            for contig in ['G1_1', 'G2_1', 'G4_1']:
                self.assertAlmostEqual(diversity_by_contig[contig], math.log(2))
            self.assertAlmostEqual(diversity_by_contig['G3_1'], 0.0)

            # trio (g1_in1, g1_out1): species_A corrected count 2/2=1.0, species_B 1/1=1.0,
            # sum=2.0 -> confidence 0.5 for both species
            for contig in ['G1_1', 'G2_1', 'G4_1']:
                self.assertAlmostEqual(confidence_by_contig[contig], 0.5)
            # trio (g1_in2, g1_out2): species_A only -> confidence 1.0
            self.assertAlmostEqual(confidence_by_contig['G3_1'], 1.0)
        finally:
            if os.path.exists(csv_path):
                os.remove(csv_path)

    def test_parse_paths_file(self):
        paths_file = f'{TEST_FILES}/SPAdes/contigs.paths'
        paths_w_gaps_file = f'{TEST_FILES}/SPAdes/contigs_w_added_gaps.paths'
        assembly_graph =  pyfastg.parse_fastg(f'{TEST_FILES}/SPAdes/assembly_graph.fastg')

        parsed_paths_without_gaps, contigs_with_gaps = pu.parse_paths_file(paths_file, assembly_graph.nodes)
        self.assertEqual(len(parsed_paths_without_gaps), 2)
        self.assertEqual(len(contigs_with_gaps), 0)

        parsed_paths_without_gaps, contigs_with_gaps = pu.parse_paths_file(paths_w_gaps_file, assembly_graph.nodes)
        self.assertEqual(len(parsed_paths_without_gaps), 1)
        self.assertEqual(len(contigs_with_gaps), 1)


if __name__ == '__main__':
    unittest.main()
