import unittest
from shutil import rmtree
import pickle
import os
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pafpy import PafRecord
import ginger.extract_contexts_candidates as ecc
from ginger import matches_classes as mc
import pyfastg
from tests import helper

TEST_FILES = helper.get_filedir()
CONTIGS_PATH = f'{TEST_FILES}/SPAdes/contigs.fasta'
CONTIG_NAME = 'NODE_1_length_1000_cov_140.620106'
# the node the gene sits on in build_short_dead_end_locus is shorter than this, so the graph can't
# supply a context, while the contig has 400bp of flanking sequence on either side of the gene
FALLBACK_MIN_CONTEXT_LEN = 50
FALLBACK_MAX_CONTEXT_LEN = 100


class TestExtractContextsCandidates(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.test_outputs_dir = 'extraction_candidates_tests_output'
        cls.in_paths_fasta = f'{cls.test_outputs_dir}/test_in_paths_fasta.fasta'
        cls.out_paths_fasta = f'{cls.test_outputs_dir}/test_out_paths_fasta.fasta'
        cls.fallback_in_paths_fasta = f'{cls.test_outputs_dir}/test_fallback_in_paths_fasta.fasta'
        cls.fallback_out_paths_fasta = f'{cls.test_outputs_dir}/test_fallback_out_paths_fasta.fasta'
        if os.path.exists(cls.test_outputs_dir):
            rmtree(cls.test_outputs_dir)
        os.mkdir(cls.test_outputs_dir)
        cls.min_context_len = 10
        cls.max_context_len = 100

    @classmethod
    def tearDownClass(cls):
        pass
        rmtree(cls.test_outputs_dir)

    def test_extract_all_in_out_paths_and_write_them_to_fastas(self):
        # with open(f'{TEST_FILES}/assembly_graph.pkl', 'rb') as f:
            # assembly_graph = pickle.load(f)
        assembly_graph= pyfastg.parse_fastg(f'{TEST_FILES}/SPAdes/assembly_graph.fastg')
        with open(f'{TEST_FILES}/nodes_with_edges_and_sequences.pkl', 'rb') as f:
            nodes_with_edges_and_sequences = pickle.load(f)
        with open(f'{TEST_FILES}/genes_with_location_in_graph.pkl', 'rb') as f:
            genes_with_location_in_graph = pickle.load(f)

        gene_lengths = ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph,
                                                                             nodes_with_edges_and_sequences,
                                                                             genes_with_location_in_graph,
                                                                             12, self.min_context_len,
                                                                             self.max_context_len,
                                                                             self.in_paths_fasta, self.out_paths_fasta,
                                                                             CONTIGS_PATH)
        self.assertTrue(os.path.exists(self.in_paths_fasta))
        self.assertTrue(os.path.exists(self.out_paths_fasta))

        with open(self.in_paths_fasta) as f:
            lines = f.readlines()
            # the -1 is beacuse of the \n in the end of the line
            self.assertEqual(len(lines[1]) - 1, self.max_context_len, 'Incoming context are not of correct length')
            self.assertListEqual(lines,
                                 ['>test_gene_nodes_5+_match_1.0000_path_5+\n',
                                  'GGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCA\n'],
                                 'Wrong incoming contexts extracted')
            

        with open(self.out_paths_fasta) as f:
            lines = f.readlines()
            self.assertEqual(len(lines[1]) - 1, self.max_context_len, 'Outgoing context are not of correct length')
            self.assertListEqual(lines,
                                 ['>test_gene_nodes_5+_match_1.0000_path_5+\n',
                                  'TCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGA\n'],
                                 'Wrong outgoing contexts extracted')
            

        self.assertDictEqual(gene_lengths, {'test_gene': 279}, 'Genes length dictionary is incorrect')


    @staticmethod
    def build_short_dead_end_locus():
        """A gene located on a graph node that is a dead end and is too short to supply a context of
        min_context_len on either side, on a 1000bp contig that has plenty of flanking sequence -
        the situation the contig fallback exists for."""
        gene_contig_match = mc.GeneContigMatch(f'{CONTIG_NAME}\tfallback_gene\t401\t700\t100\t100')
        gene_contig_match.nodes_list = ['1+']
        gene_contig_match.start_in_first_node = 10

        assembly_graph = nx.DiGraph()
        assembly_graph.add_node('1+', length=40)
        nodes_with_edges_and_sequences = {'1+': SeqRecord(Seq('A' * 40), id='1+')}
        return assembly_graph, nodes_with_edges_and_sequences, [gene_contig_match]

    def test_contig_context_fallback(self):
        assembly_graph, nodes_with_edges_and_sequences, genes_to_contigs = self.build_short_dead_end_locus()

        gene_lengths = ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph,
                                                                             nodes_with_edges_and_sequences,
                                                                             genes_to_contigs, 12,
                                                                             FALLBACK_MIN_CONTEXT_LEN,
                                                                             FALLBACK_MAX_CONTEXT_LEN,
                                                                             self.fallback_in_paths_fasta,
                                                                             self.fallback_out_paths_fasta,
                                                                             CONTIGS_PATH)

        contig_seq = str(SeqIO.index(CONTIGS_PATH, 'fasta')[CONTIG_NAME].seq)
        expected_header = f'>fallback_gene_nodes_1+_match_1.0000_path_contigfallback_{CONTIG_NAME}_401_700\n'
        with open(self.fallback_in_paths_fasta) as f:
            self.assertListEqual(f.readlines(), [expected_header, f'{contig_seq[301:401]}\n'],
                                 'Incoming context was not sliced out of the contig')
        with open(self.fallback_out_paths_fasta) as f:
            self.assertListEqual(f.readlines(), [expected_header, f'{contig_seq[700:800]}\n'],
                                 'Outgoing context was not sliced out of the contig')
        self.assertDictEqual(gene_lengths, {'fallback_gene': 300}, 'Genes length dictionary is incorrect')

    def test_contig_context_fallback_disabled(self):
        assembly_graph, nodes_with_edges_and_sequences, genes_to_contigs = self.build_short_dead_end_locus()

        ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph, nodes_with_edges_and_sequences,
                                                              genes_to_contigs, 12, FALLBACK_MIN_CONTEXT_LEN,
                                                              FALLBACK_MAX_CONTEXT_LEN, self.fallback_in_paths_fasta,
                                                              self.fallback_out_paths_fasta, CONTIGS_PATH,
                                                              contig_context_fallback=False)

        # the fastas are still created (downstream minimap2 needs them to exist), but the graph
        # supplies no context here, so with the fallback off they must stay empty
        with open(self.fallback_in_paths_fasta) as f:
            self.assertEqual(f.read(), '', 'Incoming context should not have been written')
        with open(self.fallback_out_paths_fasta) as f:
            self.assertEqual(f.read(), '', 'Outgoing context should not have been written')

    def test_contig_fallback_context_name_is_parseable(self):
        query_name = f'fallback_gene_nodes_1+_match_1.0000_path_contigfallback_{CONTIG_NAME}_401_700'
        paf_line = f'{query_name}\t100\t0\t100\t+\tref_contig\t5000\t1000\t1100\t100\t100\t60'

        path_match = mc.PathRefGenomeMatch(PafRecord.from_str(paf_line), {})

        self.assertEqual(path_match.query_name, query_name)
        self.assertEqual(path_match.gene, 'fallback_gene')
        self.assertEqual(path_match.nodes_list, '1+')
        self.assertEqual(path_match.gene_match_score, 1.0)
        self.assertEqual(path_match.path, f'contigfallback_{CONTIG_NAME}_401_700')


class FakeSeqRecord:
    def __init__(self, seq):
        self.seq = seq


class TestSavePathsToFastaIoPathsApproach(unittest.TestCase):
    """Regression tests for the min_context_len/max_context_len bug: the length written to the
    fasta (after covered_by_gene is trimmed off) must itself satisfy
    min_context_len < len(seq) <= max_context_len, not the untrimmed node-path length."""

    def setUp(self):
        self.test_outputs_dir = 'save_paths_tests_output'
        os.mkdir(self.test_outputs_dir)
        self.fasta_path = f'{self.test_outputs_dir}/test.fasta'

    def tearDown(self):
        rmtree(self.test_outputs_dir)

    def _read_written_seqs(self):
        if not os.path.exists(self.fasta_path):
            return []
        with open(self.fasta_path) as f:
            lines = f.readlines()
        return [lines[i + 1].strip() for i in range(0, len(lines), 2)]

    def test_path_dropped_when_trimmed_length_is_below_min_context_len(self):
        # node is 30bp, but 25bp of it is covered by the gene -> only 5bp of real context remain,
        # which is below min_context_len=10. Before the fix this path was incorrectly written
        # because the min_context_len check ran on the untrimmed 30bp length.
        node = '1+'
        records_dict = {node: FakeSeqRecord('A' * 30)}
        in_paths_lengths, node_locations = ecc.save_paths_to_fasta_io_paths_approach(
            [[node]], self.fasta_path, records_dict,
            max_context_len=100, min_context_len=10, in_or_out='in', covered_by_gene=25)

        self.assertEqual(self._read_written_seqs(), [], 'A too-short (post-trim) context should not be written')
        self.assertEqual(in_paths_lengths, {})
        self.assertEqual(node_locations, {})

    def test_path_kept_when_trimmed_length_is_within_bounds(self):
        node = '1+'
        records_dict = {node: FakeSeqRecord('A' * 50)}
        in_paths_lengths, _ = ecc.save_paths_to_fasta_io_paths_approach(
            [[node]], self.fasta_path, records_dict,
            max_context_len=100, min_context_len=10, in_or_out='in', covered_by_gene=10)

        written = self._read_written_seqs()
        self.assertEqual(len(written), 1)
        self.assertEqual(len(written[0]), 40)
        self.assertEqual(in_paths_lengths['1+'], 40)

    def test_trimmed_length_is_capped_at_max_context_len(self):
        node = '1+'
        records_dict = {node: FakeSeqRecord('A' * 200)}
        in_paths_lengths, _ = ecc.save_paths_to_fasta_io_paths_approach(
            [[node]], self.fasta_path, records_dict,
            max_context_len=100, min_context_len=10, in_or_out='in', covered_by_gene=10)

        written = self._read_written_seqs()
        self.assertEqual(len(written), 1)
        self.assertEqual(len(written[0]), 100)
        self.assertEqual(in_paths_lengths['1+'], 100)

    def test_out_path_dropped_when_trimmed_length_is_below_min_context_len(self):
        node = '1+'
        records_dict = {node: FakeSeqRecord('A' * 30)}
        in_paths_lengths, _ = ecc.save_paths_to_fasta_io_paths_approach(
            [[node]], self.fasta_path, records_dict,
            max_context_len=100, min_context_len=10, in_or_out='out', covered_by_gene=25)

        self.assertEqual(self._read_written_seqs(), [])
        self.assertEqual(in_paths_lengths, {})

    def test_covered_by_gene_zero_does_not_produce_empty_sequence(self):
        # regression for seq[:-0] evaluating to '' when covered_by_gene == 0
        node = '1+'
        records_dict = {node: FakeSeqRecord('A' * 50)}
        in_paths_lengths, _ = ecc.save_paths_to_fasta_io_paths_approach(
            [[node]], self.fasta_path, records_dict,
            max_context_len=100, min_context_len=10, in_or_out='in', covered_by_gene=0)

        written = self._read_written_seqs()
        self.assertEqual(len(written), 1)
        self.assertEqual(len(written[0]), 50)
        self.assertEqual(in_paths_lengths['1+'], 50)


if __name__ == '__main__':
    unittest.main()
