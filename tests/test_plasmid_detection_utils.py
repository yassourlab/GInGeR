import os
import tempfile
import shutil
import unittest
from unittest.mock import patch

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from ginger import extract_contexts_candidates as ecc
from ginger import plasmid_detection_utils as pdu
from tests import helper

TEST_FILES = helper.get_filedir()
CONTIGS_PATH = f'{TEST_FILES}/SPAdes/contigs.fasta'
CONTIG_NAME = 'NODE_1_length_1000_cov_140.620106'


class FakeGeneMatch:
    """start/end are 0-based half-open, as they are on matches_classes.GeneContigMatch."""

    def __init__(self, gene, contig, score, start, end, strand):
        self.gene = gene
        self.contig = contig
        self.score = score
        self.start = start
        self.end = end
        self.strand = strand


class FakePathMatch:
    def __init__(self, query_name):
        self.query_name = query_name


class FakeInOutMatch:
    def __init__(self, gene, in_context, out_context):
        self.gene = gene
        self.in_path = FakePathMatch(in_context)
        self.out_path = FakePathMatch(out_context)


class WritePlasmidDetectionInputFastaTest(unittest.TestCase):
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
            FakeGeneMatch('geneA', 'contig1', 1.0, 2, 6, '+'),
            FakeGeneMatch('geneB', 'contig2', 1.0, 0, 4, '-'),
        ]
        context_level_results = {
            ('geneA', 'ref_genome1'): [FakeInOutMatch('geneA', 'in_ctx1', 'out_ctx1')],
        }
        matched_genes = {'geneA'}

        output_fasta_path = os.path.join(self.tmp_dir, 'plasmid_input.fasta')
        result_path = pdu.write_plasmid_detection_input_fasta(
            context_level_results, genes_with_location_in_graph, matched_genes,
            in_paths_fasta, out_paths_fasta, contigs_fasta, output_fasta_path)

        self.assertEqual(result_path, output_fasta_path)
        records = pdu._fasta_to_dict(output_fasta_path)
        # geneA's context: in-path + gene sequence (contig1[2:6] = "GTAC") + out-path
        self.assertEqual(records['geneA|in_ctx1|out_ctx1'], 'IIIIIIIIGTACOOOOOOOO')
        # geneB has no context match, so its full (unmatched) contig is included as-is
        self.assertEqual(records['contig2'], 'GGGGCCCCAA')
        self.assertNotIn('contig1', records)

    def test_minus_strand_gene_sequence_is_not_reverse_complemented(self):
        # in/out path sequences are always extracted in the contig's forward orientation, so
        # the spliced gene segment must stay forward too, even when gene_match.strand is '-'.
        in_paths_fasta = self._write_fasta('in.fasta', {'in_ctx1': 'IIIIIIII'})
        out_paths_fasta = self._write_fasta('out.fasta', {'out_ctx1': 'OOOOOOOO'})
        contigs_fasta = self._write_fasta('contigs.fasta', {'contig1': 'AAAAACCCCC'})

        genes_with_location_in_graph = [
            FakeGeneMatch('geneA', 'contig1', 1.0, 2, 6, '-'),
        ]
        context_level_results = {
            ('geneA', 'ref_genome1'): [FakeInOutMatch('geneA', 'in_ctx1', 'out_ctx1')],
        }

        output_fasta_path = os.path.join(self.tmp_dir, 'plasmid_input.fasta')
        pdu.write_plasmid_detection_input_fasta(
            context_level_results, genes_with_location_in_graph, set(),
            in_paths_fasta, out_paths_fasta, contigs_fasta, output_fasta_path)

        records = pdu._fasta_to_dict(output_fasta_path)
        # contig1[2:6] = "AAAC", kept forward - NOT its reverse complement "GTTT", since
        # in/out paths are already extracted in the contig's forward orientation.
        self.assertEqual(records['geneA|in_ctx1|out_ctx1'], 'IIIIIIIIAAACOOOOOOOO')

    def test_returns_none_and_removes_file_when_nothing_to_write(self):
        in_paths_fasta = self._write_fasta('in.fasta', {})
        out_paths_fasta = self._write_fasta('out.fasta', {})
        contigs_fasta = self._write_fasta('contigs.fasta', {'contig1': 'ACGTACGTAC'})

        genes_with_location_in_graph = [
            FakeGeneMatch('geneA', 'contig1', 1.0, 2, 6, '+'),
        ]
        matched_genes = {'geneA'}  # geneA is matched, so its contig is not included

        output_fasta_path = os.path.join(self.tmp_dir, 'plasmid_input.fasta')
        result_path = pdu.write_plasmid_detection_input_fasta(
            {}, genes_with_location_in_graph, matched_genes,
            in_paths_fasta, out_paths_fasta, contigs_fasta, output_fasta_path)

        self.assertIsNone(result_path)
        self.assertFalse(os.path.exists(output_fasta_path))


class StitchedContextGeneContextTest(unittest.TestCase):
    """The invariant the whole context-gene-context fasta rests on: a stitched record is a
    byte-exact substring of the contig the gene was found on.

    Any drift at either junction - a base of the gene duplicated into the incoming context, or a
    base of the contig skipped before the outgoing one - breaks the reading frame of whatever sits
    across that junction, which for a full-length gene is its own stop codon.
    """

    CONTEXT_LEN = 300

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.in_paths_fasta = os.path.join(self.tmp_dir, 'in.fasta')
        self.out_paths_fasta = os.path.join(self.tmp_dir, 'out.fasta')

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _stitch(self):
        genes_with_location_in_graph = helper.get_genes_with_location_in_graph()
        ecc.extract_all_in_out_paths_and_write_them_to_fastas(helper.get_assembly_graph(),
                                                              helper.get_assembly_graph_nodes(),
                                                              genes_with_location_in_graph, 12, self.CONTEXT_LEN,
                                                              self.in_paths_fasta, self.out_paths_fasta,
                                                              CONTIGS_PATH)
        in_contexts = pdu._fasta_to_dict(self.in_paths_fasta)
        out_contexts = pdu._fasta_to_dict(self.out_paths_fasta)
        self.assertEqual(len(in_contexts), 1, 'expected exactly one incoming context for test_gene')
        self.assertEqual(len(out_contexts), 1, 'expected exactly one outgoing context for test_gene')
        in_context, out_context = next(iter(in_contexts)), next(iter(out_contexts))

        output_fasta_path = os.path.join(self.tmp_dir, 'plasmid_input.fasta')
        pdu.write_plasmid_detection_input_fasta(
            {('test_gene', 'ref_genome1'): [FakeInOutMatch('test_gene', in_context, out_context)]},
            genes_with_location_in_graph, {'test_gene'},
            self.in_paths_fasta, self.out_paths_fasta, CONTIGS_PATH, output_fasta_path)

        stitched = pdu._fasta_to_dict(output_fasta_path)[f'test_gene|{in_context}|{out_context}']
        return stitched, in_contexts[in_context], out_contexts[out_context]

    def test_stitched_record_is_an_exact_substring_of_the_contig(self):
        stitched, _, _ = self._stitch()

        contig_seq = helper.get_contig_seq(CONTIG_NAME)
        self.assertIn(stitched, contig_seq,
                      'the stitched context-gene-context is not a contiguous stretch of the contig')
        self.assertEqual(len(stitched), 2 * self.CONTEXT_LEN + 279)

    def test_the_gene_segment_sits_between_the_contexts_in_frame(self):
        stitched, in_context_seq, out_context_seq = self._stitch()

        gene_seq = stitched[len(in_context_seq):len(stitched) - len(out_context_seq)]
        protein = str(next(SeqIO.parse(f'{TEST_FILES}/test_gene.faa', 'fasta')).seq)
        self.assertEqual(str(Seq(gene_seq).translate()), protein,
                         'the segment between the two contexts is not the gene mmseqs2 aligned there')

    def test_the_contexts_abut_the_gene_with_no_overlap_and_no_gap(self):
        # the sharpest form of the invariant: cutting the two contexts back off the contig has to
        # leave exactly the gene, at the offset mmseqs2 reported for it
        stitched, in_context_seq, out_context_seq = self._stitch()

        contig_seq = helper.get_contig_seq(CONTIG_NAME)
        gene_start = contig_seq.index(stitched) + len(in_context_seq)
        self.assertEqual((gene_start, gene_start + 279), (336, 615),
                         'the gene does not sit where mmseqs2 aligned it (1-based 337-615)')


class ReadPlasmidScoresTest(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_no_context_rows_returns_empty_context_df(self):
        """GeNomad output with only contig rows (no '|' in seq_name) must not crash."""
        summary_path = os.path.join(self.tmp_dir, 'plasmid_summary.tsv')
        pd.DataFrame({
            'seq_name': ['contig1', 'contig2'],
            'length': [1000, 2000],
            'plasmid_score': [0.9, 0.3],
        }).to_csv(summary_path, sep='\t', index=False)

        context_plasmid_scores, contig_plasmid_scores = pdu.read_plasmid_scores(summary_path)

        self.assertListEqual(list(context_plasmid_scores.columns), ['gene', 'in_context', 'out_context', 'plasmid_score'])
        self.assertEqual(len(context_plasmid_scores), 0)
        self.assertListEqual(list(contig_plasmid_scores.columns), ['contig', 'plasmid_score'])
        self.assertEqual(len(contig_plasmid_scores), 2)

    def test_no_contig_rows_returns_empty_contig_df(self):
        """GeNomad output with only context rows must not crash."""
        summary_path = os.path.join(self.tmp_dir, 'plasmid_summary.tsv')
        pd.DataFrame({
            'seq_name': ['geneA|in_ctx1|out_ctx1', 'geneB|in_ctx2|out_ctx2'],
            'length': [500, 600],
            'plasmid_score': [0.8, 0.5],
        }).to_csv(summary_path, sep='\t', index=False)

        context_plasmid_scores, contig_plasmid_scores = pdu.read_plasmid_scores(summary_path)

        self.assertListEqual(list(context_plasmid_scores.columns), ['gene', 'in_context', 'out_context', 'plasmid_score'])
        self.assertEqual(len(context_plasmid_scores), 2)
        self.assertListEqual(list(contig_plasmid_scores.columns), ['contig', 'plasmid_score'])
        self.assertEqual(len(contig_plasmid_scores), 0)

    def test_splits_context_and_contig_rows(self):
        summary_path = os.path.join(self.tmp_dir, 'plasmid_summary.tsv')
        summary_df = pd.DataFrame({
            'seq_name': ['geneA|in_ctx1|out_ctx1', 'contig2'],
            'length': [1234, 5678],
            'plasmid_score': [0.9, 0.1],
        })
        summary_df.to_csv(summary_path, sep='\t', index=False)

        context_plasmid_scores, contig_plasmid_scores = pdu.read_plasmid_scores(summary_path)

        self.assertListEqual(list(context_plasmid_scores.columns), ['gene', 'in_context', 'out_context', 'plasmid_score'])
        self.assertEqual(len(context_plasmid_scores), 1)
        row = context_plasmid_scores.iloc[0]
        self.assertEqual(row['gene'], 'geneA')
        self.assertEqual(row['in_context'], 'in_ctx1')
        self.assertEqual(row['out_context'], 'out_ctx1')
        self.assertEqual(row['plasmid_score'], 0.9)

        self.assertListEqual(list(contig_plasmid_scores.columns), ['contig', 'plasmid_score'])
        self.assertEqual(len(contig_plasmid_scores), 1)
        row = contig_plasmid_scores.iloc[0]
        self.assertEqual(row['contig'], 'contig2')
        self.assertEqual(row['plasmid_score'], 0.1)


class AddPlasmidScoresToCsvTest(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_add_plasmid_scores_to_context_level_csv(self):
        csv_path = os.path.join(self.tmp_dir, 'context_level_matches.csv')
        pd.DataFrame({
            'gene': ['geneA', 'geneB'],
            'in_context': ['in_ctx1', 'in_ctx2'],
            'out_context': ['out_ctx1', 'out_ctx2'],
            'score': [0.9, 0.8],
        }).to_csv(csv_path, index=False)

        context_plasmid_scores = pd.DataFrame({
            'gene': ['geneA'],
            'in_context': ['in_ctx1'],
            'out_context': ['out_ctx1'],
            'plasmid_score': [0.7],
        })

        pdu.add_plasmid_scores_to_context_level_csv(csv_path, context_plasmid_scores)

        result = pd.read_csv(csv_path)
        self.assertEqual(result.loc[result['gene'] == 'geneA', 'plasmid_score'].iloc[0], 0.7)
        # genes with no plasmid score reported by GeNomad default to 0
        self.assertEqual(result.loc[result['gene'] == 'geneB', 'plasmid_score'].iloc[0], 0.0)

    def test_add_plasmid_scores_to_genes_no_species_match_csv(self):
        csv_path = os.path.join(self.tmp_dir, 'genes_detected_in_graph_with_no_species_match.csv')
        pd.DataFrame({
            'gene': ['geneA', 'geneB'],
            'contig': ['contig1', 'contig2'],
            'gene_match_score': [1.0, 0.95],
        }).to_csv(csv_path, index=False)

        contig_plasmid_scores = pd.DataFrame({
            'contig': ['contig1'],
            'plasmid_score': [0.4],
        })

        pdu.add_plasmid_scores_to_genes_no_species_match_csv(csv_path, contig_plasmid_scores)

        result = pd.read_csv(csv_path)
        self.assertEqual(result.loc[result['gene'] == 'geneA', 'plasmid_score'].iloc[0], 0.4)
        self.assertEqual(result.loc[result['gene'] == 'geneB', 'plasmid_score'].iloc[0], 0.0)


class RunGenomadTest(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    class Dummy:
        def __init__(self, stdout='', stderr='', returncode=0):
            self.stdout = stdout
            self.stderr = stderr
            self.returncode = returncode

    def test_raises_if_genomad_db_missing(self):
        fasta_path = os.path.join(self.tmp_dir, 'plasmid_input.fasta')
        with open(fasta_path, 'w') as f:
            f.write('>seq1\nACGT\n')

        with self.assertRaises(Exception):
            pdu.run_genomad(fasta_path, os.path.join(self.tmp_dir, 'out'),
                            os.path.join(self.tmp_dir, 'missing_genomad_db'), threads=1)

    def test_runs_genomad_and_returns_summary_path(self):
        fasta_path = os.path.join(self.tmp_dir, 'plasmid_input.fasta')
        with open(fasta_path, 'w') as f:
            f.write('>seq1\nACGT\n')
        genomad_db = os.path.join(self.tmp_dir, 'genomad_db')
        os.makedirs(genomad_db)
        output_dir = os.path.join(self.tmp_dir, 'genomad_output')

        with patch('ginger.plasmid_detection_utils.run', return_value=self.Dummy()) as mock_run:
            summary_path = pdu.run_genomad(fasta_path, output_dir, genomad_db, threads=4)

        self.assertEqual(summary_path,
                        os.path.join(output_dir, 'plasmid_input_summary', 'plasmid_input_plasmid_summary.tsv'))
        command = mock_run.call_args[0][0]
        self.assertIn('genomad end-to-end', command)
        self.assertIn('--threads 4', command)
        self.assertIn(fasta_path, command)
        self.assertIn(output_dir, command)
        self.assertIn(genomad_db, command)

    def test_raises_if_genomad_fails(self):
        fasta_path = os.path.join(self.tmp_dir, 'plasmid_input.fasta')
        with open(fasta_path, 'w') as f:
            f.write('>seq1\nACGT\n')
        genomad_db = os.path.join(self.tmp_dir, 'genomad_db')
        os.makedirs(genomad_db)
        output_dir = os.path.join(self.tmp_dir, 'genomad_output')

        with patch('ginger.plasmid_detection_utils.run', return_value=self.Dummy(stderr='boom', returncode=1)):
            with self.assertRaises(Exception):
                pdu.run_genomad(fasta_path, output_dir, genomad_db, threads=1)


if __name__ == '__main__':
    unittest.main()
