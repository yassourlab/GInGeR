import os
import unittest
from unittest.mock import patch

import pandas as pd

from ginger import matches_classes as mc
from ginger import pipeline_utils as pu
from ginger import plasmid_detection_utils as pdu
from tests import helper

# the two frames read_plasmid_scores splits GeNomad's summary into
CONTEXT_SCORE_COLUMNS = ['gene', 'in_context', 'out_context', 'plasmid_score']
CONTIG_SCORE_COLUMNS = ['contig', 'plasmid_score']


class ReadPlasmidScoresTest(helper.TempDirTestCase):
    def test_no_context_rows_returns_empty_context_df(self):
        """GeNomad output with only contig rows must not crash."""
        summary_path = os.path.join(self.tmp_dir, 'plasmid_summary.tsv')
        pd.DataFrame({
            'seq_name': ['contig1', 'contig2'],
            'length': [1000, 2000],
            'plasmid_score': [0.9, 0.3],
        }).to_csv(summary_path, sep='\t', index=False)

        context_plasmid_scores, contig_plasmid_scores = pdu.read_plasmid_scores(summary_path, {})

        self.assertListEqual(list(context_plasmid_scores.columns), CONTEXT_SCORE_COLUMNS)
        self.assertEqual(len(context_plasmid_scores), 0)
        self.assertListEqual(list(contig_plasmid_scores.columns), CONTIG_SCORE_COLUMNS)
        self.assertEqual(len(contig_plasmid_scores), 2)

    def test_no_contig_rows_returns_empty_contig_df(self):
        """GeNomad output with only context rows must not crash."""
        summary_path = os.path.join(self.tmp_dir, 'plasmid_summary.tsv')
        pd.DataFrame({
            'seq_name': ['ctx0000000', 'ctx0000001'],
            'length': [500, 600],
            'plasmid_score': [0.8, 0.5],
        }).to_csv(summary_path, sep='\t', index=False)

        context_plasmid_scores, contig_plasmid_scores = pdu.read_plasmid_scores(
            summary_path, {'ctx0000000': pu.ContextSeqRecord('ctx0000000', 'geneA', 'in_ctx1', 'out_ctx1', 8, 12),
                           'ctx0000001': pu.ContextSeqRecord('ctx0000001', 'geneB', 'in_ctx2', 'out_ctx2', 8, 12)})

        self.assertListEqual(list(context_plasmid_scores.columns), CONTEXT_SCORE_COLUMNS)
        self.assertEqual(len(context_plasmid_scores), 2)
        self.assertListEqual(list(contig_plasmid_scores.columns), CONTIG_SCORE_COLUMNS)
        self.assertEqual(len(contig_plasmid_scores), 0)

    def test_splits_context_and_contig_rows(self):
        summary_path = os.path.join(self.tmp_dir, 'plasmid_summary.tsv')
        summary_df = pd.DataFrame({
            'seq_name': ['ctx0000000', 'contig2'],
            'length': [1234, 5678],
            'plasmid_score': [0.9, 0.1],
        })
        summary_df.to_csv(summary_path, sep='\t', index=False)

        context_plasmid_scores, contig_plasmid_scores = pdu.read_plasmid_scores(
            summary_path, {'ctx0000000': pu.ContextSeqRecord('ctx0000000', 'geneA', 'in_ctx1', 'out_ctx1', 8, 12)})

        self.assertListEqual(list(context_plasmid_scores.columns), CONTEXT_SCORE_COLUMNS)
        self.assertEqual(len(context_plasmid_scores), 1)
        row = context_plasmid_scores.iloc[0]
        self.assertEqual(row['gene'], 'geneA')
        self.assertEqual(row['in_context'], 'in_ctx1')
        self.assertEqual(row['out_context'], 'out_ctx1')
        self.assertEqual(row['plasmid_score'], 0.9)

        self.assertListEqual(list(contig_plasmid_scores.columns), CONTIG_SCORE_COLUMNS)
        self.assertEqual(len(contig_plasmid_scores), 1)
        row = contig_plasmid_scores.iloc[0]
        self.assertEqual(row['contig'], 'contig2')
        self.assertEqual(row['plasmid_score'], 0.1)


class GeneNamesWithSeparatorsTest(helper.TempDirTestCase):
    """Gene names come from a fasta the caller supplies, so no character is safe to build a
    composite sequence name out of. SARG's are pipe delimited and CARD's contain both '|' and ':'.

    The round trip that has to survive one: writing the contexts fasta, then reading GeNomad's scores
    for it back onto the gene the sequence was built for.
    """

    SARG_GENE = 'SARG|multidrug@MFS|emrB|WP_145513356.1'

    def test_a_gene_name_full_of_separators_survives_the_round_trip(self):
        in_paths_fasta = self.write_fasta('in.fasta', {'in_ctx1': 'IIIIIIII'})
        out_paths_fasta = self.write_fasta('out.fasta', {'out_ctx1': 'OOOOOOOO'})
        contigs_fasta = self.write_fasta('contigs.fasta', {'contig1': 'ACGTACGTAC'})
        context_level_results = {
            (self.SARG_GENE, 'ref_genome1'): [helper.FakeInOutMatch(self.SARG_GENE, 'in_ctx1', 'out_ctx1',
                                                                    mc.GeneLocus('contig1', 2, 6))],
        }
        output_fasta_path = os.path.join(self.tmp_dir, 'in_gene_out_contexts.fasta')

        _, records_by_seq_id = pu.write_in_gene_out_contexts_fasta(
            context_level_results, [helper.FakeGeneMatch(self.SARG_GENE, 'contig1', 1.0)], {self.SARG_GENE},
            in_paths_fasta, out_paths_fasta, contigs_fasta, output_fasta_path)

        # the record's own name holds none of this, so nothing about the gene reaches GeNomad's input
        self.assertEqual(list(pu._fasta_to_dict(output_fasta_path)), ['ctx0000000'])

        summary_path = os.path.join(self.tmp_dir, 'plasmid_summary.tsv')
        pd.DataFrame({'seq_name': ['ctx0000000'], 'plasmid_score': [0.7]}).to_csv(summary_path, sep='\t',
                                                                                  index=False)
        context_plasmid_scores, _ = pdu.read_plasmid_scores(summary_path, records_by_seq_id)

        # the gene comes back whole. splitting a "{gene}|{in}|{out}" name on '|' would have made
        # this gene 'SARG', its incoming context 'multidrug@MFS' and its outgoing one 'emrB', so
        # the merge onto the context level csv would miss every row and leave the score at 0
        row = context_plasmid_scores.iloc[0]
        self.assertEqual(row['gene'], self.SARG_GENE)
        self.assertEqual(row['in_context'], 'in_ctx1')
        self.assertEqual(row['out_context'], 'out_ctx1')
        self.assertEqual(row['plasmid_score'], 0.7)


class KeepOnlyPlasmidSummaryTest(helper.TempDirTestCase):
    def setUp(self):
        super().setUp()
        self.genomad_out_dir = os.path.join(self.tmp_dir, 'genomad_output')
        self.summary_path = os.path.join(self.genomad_out_dir, 'input_summary', 'input_plasmid_summary.tsv')
        os.makedirs(os.path.dirname(self.summary_path))
        with open(self.summary_path, 'w') as f:
            f.write('seq_name\tplasmid_score\nctx0000000\t0.9\n')
        # the ~80MB per sample that nothing downstream reads
        os.makedirs(os.path.join(self.genomad_out_dir, 'input_annotate'))
        with open(os.path.join(self.genomad_out_dir, 'input_annotate', 'input_proteins.faa'), 'w') as f:
            f.write('>ctx0000000_1\nMKV\n')

    def test_keeps_the_summary_and_removes_the_tree(self):
        kept_path = os.path.join(self.tmp_dir, 'plasmid_summary.tsv')

        self.assertEqual(pdu.keep_only_plasmid_summary(self.genomad_out_dir, self.summary_path, kept_path),
                         kept_path)

        self.assertFalse(os.path.exists(self.genomad_out_dir))
        with open(kept_path) as f:
            self.assertEqual(f.read(), 'seq_name\tplasmid_score\nctx0000000\t0.9\n')

    def test_refuses_to_keep_the_summary_inside_the_tree_it_removes(self):
        # the move would succeed and the rmtree would then delete the only file worth keeping
        with self.assertRaises(ValueError):
            pdu.keep_only_plasmid_summary(self.genomad_out_dir, self.summary_path,
                                          os.path.join(self.genomad_out_dir, 'plasmid_summary.tsv'))
        self.assertTrue(os.path.exists(self.summary_path))


class AddPlasmidScoresToCsvTest(helper.TempDirTestCase):
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


class RunGenomadTest(helper.TempDirTestCase):
    class Dummy:
        def __init__(self, stdout='', stderr='', returncode=0):
            self.stdout = stdout
            self.stderr = stderr
            self.returncode = returncode

    def test_raises_if_genomad_db_missing(self):
        fasta_path = os.path.join(self.tmp_dir, 'in_gene_out_contexts.fasta')
        with open(fasta_path, 'w') as f:
            f.write('>seq1\nACGT\n')

        with self.assertRaises(Exception):
            pdu.run_genomad(fasta_path, os.path.join(self.tmp_dir, 'out'),
                            os.path.join(self.tmp_dir, 'missing_genomad_db'), threads=1)

    def test_runs_genomad_and_returns_summary_path(self):
        fasta_path = os.path.join(self.tmp_dir, 'in_gene_out_contexts.fasta')
        with open(fasta_path, 'w') as f:
            f.write('>seq1\nACGT\n')
        genomad_db = os.path.join(self.tmp_dir, 'genomad_db')
        os.makedirs(genomad_db)
        output_dir = os.path.join(self.tmp_dir, 'genomad_output')

        with patch('ginger.pipeline_utils.run', return_value=self.Dummy()) as mock_run:
            summary_path = pdu.run_genomad(fasta_path, output_dir, genomad_db, threads=4)

        self.assertEqual(summary_path,
                        os.path.join(output_dir, 'in_gene_out_contexts_summary',
                                     'in_gene_out_contexts_plasmid_summary.tsv'))
        command = mock_run.call_args[0][0]
        self.assertIn('genomad end-to-end', command)
        self.assertIn('--threads 4', command)
        self.assertIn(fasta_path, command)
        self.assertIn(output_dir, command)
        self.assertIn(genomad_db, command)

    def test_raises_if_genomad_fails(self):
        fasta_path = os.path.join(self.tmp_dir, 'in_gene_out_contexts.fasta')
        with open(fasta_path, 'w') as f:
            f.write('>seq1\nACGT\n')
        genomad_db = os.path.join(self.tmp_dir, 'genomad_db')
        os.makedirs(genomad_db)
        output_dir = os.path.join(self.tmp_dir, 'genomad_output')

        with patch('ginger.pipeline_utils.run', return_value=self.Dummy(stderr='boom', returncode=1)):
            with self.assertRaises(Exception):
                pdu.run_genomad(fasta_path, output_dir, genomad_db, threads=1)


if __name__ == '__main__':
    unittest.main()
