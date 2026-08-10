import math
import unittest
import os
import shutil
import tempfile
from types import SimpleNamespace
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from ginger import extract_contexts_candidates as ecc
from ginger import pipeline_utils as pu
from ginger.matches_classes import GeneLocus, InOutPathsMatch

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

        # in1 is only ever paired with out1 (and in2 with out2), so the in_context-only and
        # out_context-only scores agree and their average equals either one.
        # trio (in1, out1) matches 2 genomes of species_A and 1 of species_B. Corrected counts
        # (2/2=1, 1/1=1) are equal, so each species gets confidence 0.5
        for _, row in result[result['in_context'] == 'in1'].iterrows():
            self.assertAlmostEqual(row['context_species_confidence_score'], 0.5)

        # trio (in2, out2) matches only species_A -> single species -> confidence 1.0
        for _, row in result[result['in_context'] == 'in2'].iterrows():
            self.assertAlmostEqual(row['context_species_confidence_score'], 1.0)

    def test_compute_context_species_confidence_score_in_out_context_disagree(self):
        # in1 is shared by two out_contexts that each match a single, different species, so the
        # in_context-only and out_context-only scores disagree and must be averaged.
        results_df = pd.DataFrame({
            'gene': ['g1', 'g1'],
            'in_context': ['in1', 'in1'],
            'out_context': ['out1', 'out2'],
            'Genome': ['G1', 'G2'],
            'species': ['species_A', 'species_B'],
        })
        species_reference_counts = pd.Series({'species_A': 1, 'species_B': 1})

        result = pu.compute_context_species_confidence_score(results_df, species_reference_counts)

        # in_context score (in1 matches both species equally) = 0.5 for both rows
        # out_context score (out1/out2 each match a single species) = 1.0 for both rows
        # average = 0.75
        for _, row in result.iterrows():
            self.assertAlmostEqual(row['context_species_confidence_score'], 0.75)

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
        parsed_paths, contigs_with_gaps = pu.parse_paths_file(paths_file)
        self.assertEqual(len(parsed_paths), 2)
        self.assertEqual(len(contigs_with_gaps), 0)
        # a contig assembled from a single graph path is a single segment
        self.assertEqual(parsed_paths['NODE_1_length_1000_cov_140.620106'], [['5+']])
        self.assertEqual(parsed_paths["NODE_1_length_1000_cov_140.620106'"], [['5-']])

        parsed_paths, contigs_with_gaps = pu.parse_paths_file(paths_w_gaps_file)
        self.assertEqual(len(parsed_paths), 2)
        self.assertEqual(contigs_with_gaps, {'NODE_1_length_1000_cov_140.620106'})
        # the path of a gap-containing contig is kept, split into one segment per part
        self.assertEqual(parsed_paths['NODE_1_length_1000_cov_140.620106'], [['5+'], ['6+']])

    def test_parse_paths_file_with_multiple_segments(self):
        paths_file = f'{TEST_FILES}/test_contigs.paths'
        parsed_paths, contigs_with_gaps = pu.parse_paths_file(paths_file)
        self.assertEqual(contigs_with_gaps, {'NODE_7_length_41181_cov_4.618952',
                                             "NODE_7_length_41181_cov_4.618952'",
                                             'NODE_9_length_39416_cov_5.216737',
                                             "NODE_9_length_39416_cov_5.216737'"})
        self.assertEqual(parsed_paths['NODE_7_length_41181_cov_4.618952'],
                         [['597920+'], ['6335-'],
                          ['99310+', '66769-', '66771+', '316649+', '153255+', '424761-', '424763+', '47+']])
        self.assertEqual(parsed_paths['NODE_9_length_39416_cov_5.216737'], [['74101+'], ['69851+']])
        # the last line of the file has no trailing newline
        self.assertEqual(parsed_paths["NODE_10_length_36078_cov_6.312495'"], [['76999-']])


class WriteInGeneOutContextsFastaTest(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _write_fasta(self, name, records):
        path = os.path.join(self.tmp_dir, name)
        with open(path, 'w') as f:
            for header, seq in records.items():
                f.write(f'>{header}\n{seq}\n')
        return path

    def test_writes_context_trios_and_unmatched_contigs(self):
        in_paths_fasta = self._write_fasta('in.fasta', {'in_ctx1': 'IIIIIIII'})
        out_paths_fasta = self._write_fasta('out.fasta', {'out_ctx1': 'OOOOOOOO'})
        contigs_fasta = self._write_fasta('contigs.fasta', {
            'contig1': 'ACGTACGTAC',
            'contig2': 'GGGGCCCCAA',
        })

        genes_with_location_in_graph = [
            helper.FakeGeneMatch('geneA', 'contig1', 1.0),
            helper.FakeGeneMatch('geneB', 'contig2', 1.0),
        ]
        context_level_results = {
            ('geneA', 'ref_genome1'): [helper.FakeInOutMatch('geneA', 'in_ctx1', 'out_ctx1',
                                                             GeneLocus('contig1', 2, 6))],
        }
        matched_genes = {'geneA'}

        output_fasta_path = os.path.join(self.tmp_dir, 'in_gene_out_contexts.fasta')
        result_path, records_by_seq_id = pu.write_in_gene_out_contexts_fasta(
            context_level_results, genes_with_location_in_graph, matched_genes,
            in_paths_fasta, out_paths_fasta, contigs_fasta, output_fasta_path)

        self.assertEqual(result_path, output_fasta_path)
        records = pu._fasta_to_dict(output_fasta_path)
        # geneA's context: in-path + gene sequence (contig1[2:6] = "GTAC") + out-path
        self.assertEqual(records_by_seq_id,
                         {'ctx0000000': pu.ContextSeqRecord('ctx0000000', 'geneA', 'in_ctx1', 'out_ctx1', 8, 12)})
        self.assertEqual(records['ctx0000000'], 'IIIIIIIIGTACOOOOOOOO')
        # geneB has no context match, so its full (unmatched) contig is included as-is
        self.assertEqual(records['contig2'], 'GGGGCCCCAA')
        self.assertNotIn('contig1', records)

    def test_gene_offsets_cut_the_gene_back_out_of_the_record(self):
        # the offsets are what a user slices the record with to highlight the gene, so they have to
        # land on the gene and not on either flank
        in_paths_fasta = self._write_fasta('in.fasta', {'in_ctx1': 'IIIIIIII'})
        out_paths_fasta = self._write_fasta('out.fasta', {'out_ctx1': 'OOOO'})
        contigs_fasta = self._write_fasta('contigs.fasta', {'contig1': 'ACGTACGTAC'})
        context_level_results = {
            ('geneA', 'ref_genome1'): [helper.FakeInOutMatch('geneA', 'in_ctx1', 'out_ctx1',
                                                             GeneLocus('contig1', 2, 6))],
        }

        output_fasta_path = os.path.join(self.tmp_dir, 'in_gene_out_contexts.fasta')
        _, records_by_seq_id = pu.write_in_gene_out_contexts_fasta(
            context_level_results, [], set(), in_paths_fasta, out_paths_fasta, contigs_fasta,
            output_fasta_path)

        record = records_by_seq_id['ctx0000000']
        seq = pu._fasta_to_dict(output_fasta_path)['ctx0000000']
        self.assertEqual(seq[record.gene_start:record.gene_end], 'GTAC')
        self.assertEqual(seq[:record.gene_start], 'IIIIIIII')
        self.assertEqual(seq[record.gene_end:], 'OOOO')

    def test_gene_segment_is_taken_forward_out_of_the_contig(self):
        # in/out path sequences are always extracted in the contig's forward orientation, so the
        # spliced gene segment stays forward too, however the gene itself is oriented
        in_paths_fasta = self._write_fasta('in.fasta', {'in_ctx1': 'IIIIIIII'})
        out_paths_fasta = self._write_fasta('out.fasta', {'out_ctx1': 'OOOOOOOO'})
        contigs_fasta = self._write_fasta('contigs.fasta', {'contig1': 'AAAAACCCCC'})

        genes_with_location_in_graph = [helper.FakeGeneMatch('geneA', 'contig1', 1.0)]
        context_level_results = {
            ('geneA', 'ref_genome1'): [helper.FakeInOutMatch('geneA', 'in_ctx1', 'out_ctx1',
                                                             GeneLocus('contig1', 2, 6))],
        }

        output_fasta_path = os.path.join(self.tmp_dir, 'in_gene_out_contexts.fasta')
        pu.write_in_gene_out_contexts_fasta(
            context_level_results, genes_with_location_in_graph, set(),
            in_paths_fasta, out_paths_fasta, contigs_fasta, output_fasta_path)

        records = pu._fasta_to_dict(output_fasta_path)
        # contig1[2:6] = "AAAC", kept forward - NOT its reverse complement "GTTT", since
        # in/out paths are already extracted in the contig's forward orientation.
        self.assertEqual(records['ctx0000000'], 'IIIIIIIIAAACOOOOOOOO')

    def test_the_gene_segment_comes_from_the_copy_the_contexts_flank(self):
        # geneA is on two contigs. the trio's contexts were cut from the copy on contig2, so that is
        # the copy whose sequence goes between them - not the copy that happens to match best
        in_paths_fasta = self._write_fasta('in.fasta', {'in_ctx1': 'IIIIIIII'})
        out_paths_fasta = self._write_fasta('out.fasta', {'out_ctx1': 'OOOOOOOO'})
        contigs_fasta = self._write_fasta('contigs.fasta', {'contig1': 'AAAAAAAAAA',
                                                            'contig2': 'CCCCGGGGTT'})

        genes_with_location_in_graph = [
            helper.FakeGeneMatch('geneA', 'contig1', 1.0),
            helper.FakeGeneMatch('geneA', 'contig2', 0.5),
        ]
        context_level_results = {
            ('geneA', 'ref_genome1'): [helper.FakeInOutMatch('geneA', 'in_ctx1', 'out_ctx1',
                                                             GeneLocus('contig2', 4, 8))],
        }

        output_fasta_path = os.path.join(self.tmp_dir, 'in_gene_out_contexts.fasta')
        pu.write_in_gene_out_contexts_fasta(
            context_level_results, genes_with_location_in_graph, {'geneA'},
            in_paths_fasta, out_paths_fasta, contigs_fasta, output_fasta_path)

        records = pu._fasta_to_dict(output_fasta_path)
        # contig2[4:8], not contig1[2:6] - splicing the other copy's sequence in would produce a
        # sequence that is on neither contig
        self.assertEqual(records['ctx0000000'], 'IIIIIIIIGGGGOOOOOOOO')

    def test_returns_none_and_removes_file_when_nothing_to_write(self):
        in_paths_fasta = self._write_fasta('in.fasta', {})
        out_paths_fasta = self._write_fasta('out.fasta', {})
        contigs_fasta = self._write_fasta('contigs.fasta', {'contig1': 'ACGTACGTAC'})

        genes_with_location_in_graph = [helper.FakeGeneMatch('geneA', 'contig1', 1.0)]
        matched_genes = {'geneA'}  # geneA is matched, so its contig is not included

        output_fasta_path = os.path.join(self.tmp_dir, 'in_gene_out_contexts.fasta')
        result_path, records_by_seq_id = pu.write_in_gene_out_contexts_fasta(
            {}, genes_with_location_in_graph, matched_genes,
            in_paths_fasta, out_paths_fasta, contigs_fasta, output_fasta_path)

        self.assertIsNone(result_path)
        self.assertEqual(records_by_seq_id, {})
        self.assertFalse(os.path.exists(output_fasta_path))

    def test_sequences_are_numbered_the_same_way_on_a_rerun(self):
        # ids are not stable across samples - they follow the sorted trios - but the same input has
        # to number them the same way twice, or a rerun's fasta and csv would not agree
        in_paths_fasta = self._write_fasta('in.fasta', {'in_ctx1': 'IIIIIIII', 'in_ctx2': 'IIII'})
        out_paths_fasta = self._write_fasta('out.fasta', {'out_ctx1': 'OOOOOOOO', 'out_ctx2': 'OOOO'})
        contigs_fasta = self._write_fasta('contigs.fasta', {'contig1': 'ACGTACGTAC'})
        locus = GeneLocus('contig1', 2, 6)
        context_level_results = {
            ('geneA', 'ref_genome1'): [helper.FakeInOutMatch('geneA', 'in_ctx1', 'out_ctx1', locus),
                                       helper.FakeInOutMatch('geneA', 'in_ctx2', 'out_ctx2', locus)],
        }

        numbering = []
        for run in range(2):
            _, records_by_seq_id = pu.write_in_gene_out_contexts_fasta(
                context_level_results, [], set(), in_paths_fasta, out_paths_fasta, contigs_fasta,
                os.path.join(self.tmp_dir, f'in_gene_out_contexts_{run}.fasta'))
            numbering.append(records_by_seq_id)

        self.assertEqual(numbering[0], numbering[1])
        self.assertEqual(sorted(numbering[0]), ['ctx0000000', 'ctx0000001'])


class StitchedContextGeneContextTest(unittest.TestCase):
    """The invariant the whole context-gene-context fasta rests on: the two contexts abut the gene
    with no overlap and no gap.

    Any drift at either junction - a base of the gene duplicated into the incoming context, or a
    base of the contig skipped before the outgoing one - breaks the reading frame of whatever sits
    across that junction, which for a full-length gene is its own stop codon.

    test_gene has one path through the graph on each side, and those paths are the contig's own, so
    here the invariant takes its sharpest form: the record has to be a byte-exact substring of the
    contig. A gene with a branch on either side also gets records read off the alternative paths,
    and those are not on any single contig - the junctions still have to be exact.
    """

    CONTEXT_LEN = 300
    CONTIG_NAME = 'NODE_1_length_1000_cov_140.620106'

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.in_paths_fasta = os.path.join(self.tmp_dir, 'in.fasta')
        self.out_paths_fasta = os.path.join(self.tmp_dir, 'out.fasta')

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _stitch(self):
        contigs_path = f'{TEST_FILES}/SPAdes/contigs.fasta'
        genes_with_location_in_graph = helper.get_genes_with_location_in_graph()
        ecc.extract_all_in_out_paths_and_write_them_to_fastas(helper.get_assembly_graph(), helper.get_geometry(),
                                                              genes_with_location_in_graph, 12, self.CONTEXT_LEN,
                                                              self.in_paths_fasta, self.out_paths_fasta,
                                                              contigs_path)
        in_contexts = pu._fasta_to_dict(self.in_paths_fasta)
        out_contexts = pu._fasta_to_dict(self.out_paths_fasta)
        self.assertEqual(len(in_contexts), 1, 'expected exactly one incoming context for test_gene')
        self.assertEqual(len(out_contexts), 1, 'expected exactly one outgoing context for test_gene')
        in_context, out_context = next(iter(in_contexts)), next(iter(out_contexts))

        gene_match = genes_with_location_in_graph[0]
        locus = GeneLocus(gene_match.contig, gene_match.start, gene_match.end)
        output_fasta_path = os.path.join(self.tmp_dir, 'in_gene_out_contexts.fasta')
        _, records_by_seq_id = pu.write_in_gene_out_contexts_fasta(
            {('test_gene', 'ref_genome1'): [helper.FakeInOutMatch('test_gene', in_context, out_context, locus)]},
            genes_with_location_in_graph, {'test_gene'},
            self.in_paths_fasta, self.out_paths_fasta, contigs_path, output_fasta_path)

        self.record = records_by_seq_id['ctx0000000']
        stitched = pu._fasta_to_dict(output_fasta_path)['ctx0000000']
        return stitched, in_contexts[in_context], out_contexts[out_context]

    def test_the_gene_offsets_land_on_the_gene(self):
        # the same boundaries the two tests below check on the contig, but expressed the way a user
        # gets them - as the columns of context_level_matches.csv that slice this record
        stitched, in_context_seq, _ = self._stitch()

        self.assertEqual((self.record.gene_start, self.record.gene_end),
                         (self.CONTEXT_LEN, self.CONTEXT_LEN + 279))
        self.assertEqual(stitched[self.record.gene_start:self.record.gene_end],
                         stitched[len(in_context_seq):len(in_context_seq) + 279])

    def test_stitched_record_is_an_exact_substring_of_the_contig(self):
        stitched, _, _ = self._stitch()

        contig_seq = helper.get_contig_seq(self.CONTIG_NAME)
        self.assertIn(stitched, contig_seq,
                      'the stitched context-gene-context is not a contiguous stretch of the contig')
        self.assertEqual(len(stitched), 2 * self.CONTEXT_LEN + 279)

    def test_the_gene_segment_sits_between_the_contexts_in_frame(self):
        stitched, in_context_seq, out_context_seq = self._stitch()

        gene_seq = stitched[len(in_context_seq):len(stitched) - len(out_context_seq)]
        with open(f'{TEST_FILES}/test_gene.faa') as f:
            protein = str(next(SeqIO.parse(f, 'fasta')).seq)
        self.assertEqual(str(Seq(gene_seq).translate()), protein,
                         'the segment between the two contexts is not the gene mmseqs2 aligned there')

    def test_the_contexts_abut_the_gene_with_no_overlap_and_no_gap(self):
        # the sharpest form of the invariant: cutting the two contexts back off the contig has to
        # leave exactly the gene, at the offset mmseqs2 reported for it
        stitched, in_context_seq, out_context_seq = self._stitch()

        contig_seq = helper.get_contig_seq(self.CONTIG_NAME)
        gene_start = contig_seq.index(stitched) + len(in_context_seq)
        self.assertEqual((gene_start, gene_start + 279), (336, 615),
                         'the gene does not sit where mmseqs2 aligned it (1-based 337-615)')


class AddContextSeqIdsToContextLevelCsvTest(unittest.TestCase):
    """The join from a context level row to the sequence it was concluded from."""

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.csv_path = os.path.join(self.tmp_dir, 'context_level_matches.csv')

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _write_csv(self, rows):
        pd.DataFrame(rows).to_csv(self.csv_path, index=False)

    def test_adds_the_seq_id_and_the_gene_offsets(self):
        # geneA matched two reference genomes off the same pair of contexts - one sequence was
        # written for the pair, so both rows point at it
        self._write_csv({'gene': ['geneA', 'geneA'],
                          'in_context': ['in_ctx1', 'in_ctx1'],
                          'out_context': ['out_ctx1', 'out_ctx1'],
                          'reference_contig': ['ref1', 'ref2']})

        pu.add_context_seq_ids_to_context_level_csv(
            self.csv_path,
            {'ctx0000000': pu.ContextSeqRecord('ctx0000000', 'geneA', 'in_ctx1', 'out_ctx1', 300, 579)})

        result = pd.read_csv(self.csv_path)
        self.assertEqual(len(result), 2, 'a row was duplicated or dropped by the merge')
        self.assertListEqual(list(result[pu.CONTEXT_SEQ_ID_COLUMN]), ['ctx0000000', 'ctx0000000'])
        self.assertListEqual(list(result['gene_start_in_context_seq']), [300, 300])
        self.assertListEqual(list(result['gene_end_in_context_seq']), [579, 579])

    def test_offsets_are_written_as_integers(self):
        # a float column would write every offset as "300.0", which is not a slice index
        self._write_csv({'gene': ['geneA'], 'in_context': ['in_ctx1'], 'out_context': ['out_ctx1']})

        pu.add_context_seq_ids_to_context_level_csv(
            self.csv_path,
            {'ctx0000000': pu.ContextSeqRecord('ctx0000000', 'geneA', 'in_ctx1', 'out_ctx1', 300, 579)})

        with open(self.csv_path) as f:
            self.assertIn('ctx0000000,300,579', f.read())

    def test_a_row_with_no_record_is_left_empty(self):
        self._write_csv({'gene': ['geneA', 'geneB'],
                          'in_context': ['in_ctx1', 'in_ctx2'],
                          'out_context': ['out_ctx1', 'out_ctx2']})

        pu.add_context_seq_ids_to_context_level_csv(
            self.csv_path,
            {'ctx0000000': pu.ContextSeqRecord('ctx0000000', 'geneA', 'in_ctx1', 'out_ctx1', 300, 579)})

        result = pd.read_csv(self.csv_path)
        geneb = result[result['gene'] == 'geneB'].iloc[0]
        self.assertTrue(pd.isna(geneb[pu.CONTEXT_SEQ_ID_COLUMN]))
        self.assertTrue(pd.isna(geneb['gene_start_in_context_seq']))


if __name__ == '__main__':
    unittest.main()
