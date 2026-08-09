import unittest
from shutil import rmtree
import os
import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pafpy import PafRecord
import ginger.extract_contexts_candidates as ecc
from ginger import matches_classes as mc
from tests import helper

TEST_FILES = helper.get_filedir()
CONTIGS_PATH = f'{TEST_FILES}/SPAdes/contigs.fasta'
CONTIG_NAME = 'NODE_1_length_1000_cov_140.620106'
# where test_gene sits on that contig, 0-based half-open. mmseqs reports it as 1-based 337-615 in
# genes_to_contigs.m8, and 615 - 337 + 1 == 279 == 93aa * 3
GENE_START, GENE_END = 336, 615
# the node the gene sits on in build_gene_on_gappy_contig is shorter than this, so the graph can't
# supply a context, while the contig has 400bp of flanking sequence on either side of the gene
FALLBACK_CONTEXT_LEN = 100


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
        cls.context_len = 100

    @classmethod
    def tearDownClass(cls):
        pass
        rmtree(cls.test_outputs_dir)

    def test_extract_all_in_out_paths_and_write_them_to_fastas(self):
        assembly_graph = helper.get_assembly_graph()
        nodes_with_edges_and_sequences = helper.get_assembly_graph_nodes()
        genes_with_location_in_graph = helper.get_genes_with_location_in_graph()

        gene_lengths, _ = ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph,
                                                                             nodes_with_edges_and_sequences,
                                                                             genes_with_location_in_graph,
                                                                             12, self.context_len,
                                                                             self.in_paths_fasta, self.out_paths_fasta,
                                                                             CONTIGS_PATH)
        self.assertTrue(os.path.exists(self.in_paths_fasta))
        self.assertTrue(os.path.exists(self.out_paths_fasta))

        # test_gene is aligned to 1-based 337-615 of the contig, which is GENE_START:GENE_END
        # 0-based. the contexts are the context_len bases that flank exactly that interval - no base
        # of the gene may appear in them, and no base of the contig may be skipped between them.
        contig_seq = helper.get_contig_seq(CONTIG_NAME)
        with open(self.in_paths_fasta) as f:
            lines = f.readlines()
            # the -1 is beacuse of the \n in the end of the line
            self.assertEqual(len(lines[1]) - 1, self.context_len, 'Incoming context are not of correct length')
            self.assertListEqual(lines,
                                 ['>test_gene_nodes_5+_match_1.0000_path_5+\n',
                                  f'{contig_seq[GENE_START - self.context_len:GENE_START]}\n'],
                                 'Wrong incoming contexts extracted')

        with open(self.out_paths_fasta) as f:
            lines = f.readlines()
            self.assertEqual(len(lines[1]) - 1, self.context_len, 'Outgoing context are not of correct length')
            self.assertListEqual(lines,
                                 ['>test_gene_nodes_5+_match_1.0000_path_5+\n',
                                  f'{contig_seq[GENE_END:GENE_END + self.context_len]}\n'],
                                 'Wrong outgoing contexts extracted')

        self.assertDictEqual(gene_lengths, {'test_gene': 279}, 'Genes length dictionary is incorrect')

    def test_the_outgoing_context_starts_at_the_end_of_the_alignment_not_of_the_reference_gene(self):
        # test_gene is a 93aa fragment whose alignment covers the contig exactly, so pad gene_length
        # out to a reference protein that is 60bp longer than the aligned span. the outgoing context
        # must not move: it follows the alignment, and gene_length says nothing about the contig.
        assembly_graph = helper.get_assembly_graph()
        nodes_with_edges_and_sequences = helper.get_assembly_graph_nodes()
        genes_with_location_in_graph = helper.get_genes_with_location_in_graph()
        for gene_contigs_match in genes_with_location_in_graph:
            gene_contigs_match.gene_length += 60

        ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph, nodes_with_edges_and_sequences,
                                                              genes_with_location_in_graph, 12, self.context_len,
                                                              self.in_paths_fasta, self.out_paths_fasta, CONTIGS_PATH)

        contig_seq = helper.get_contig_seq(CONTIG_NAME)
        with open(self.out_paths_fasta) as f:
            self.assertEqual(f.readlines()[1].strip(), contig_seq[GENE_END:GENE_END + self.context_len],
                             'The outgoing context was placed using gene_length instead of aligned_length')


    def test_contexts_are_mapped_to_the_locus_they_were_cut_from(self):
        contexts_to_loci_path = f'{self.test_outputs_dir}/contexts_to_loci.tsv'

        _, contexts_to_loci = ecc.extract_all_in_out_paths_and_write_them_to_fastas(
            helper.get_assembly_graph(), helper.get_assembly_graph_nodes(),
            helper.get_genes_with_location_in_graph(), 12, self.context_len,
            self.in_paths_fasta, self.out_paths_fasta, CONTIGS_PATH,
            contexts_to_loci_path=contexts_to_loci_path)

        context_name = 'test_gene_nodes_5+_match_1.0000_path_5+'
        self.assertEqual(contexts_to_loci, {context_name: mc.GeneLocus(CONTIG_NAME, GENE_START, GENE_END)})

        with open(contexts_to_loci_path) as f:
            header, *rows = [line.rstrip('\n').split('\t') for line in f]
        self.assertEqual(header, ecc.CONTEXTS_TO_LOCI_COLUMNS)
        # the gene sits on a single node, so both of its contexts are named after the same one-node
        # path - one name, but a row per side
        self.assertEqual(rows, [
            [context_name, 'in', CONTIG_NAME, str(GENE_START), str(GENE_END), '5+', '1.0', 'graph'],
            [context_name, 'out', CONTIG_NAME, str(GENE_START), str(GENE_END), '5+', '1.0', 'graph'],
        ])

    def test_contig_contexts_are_mapped_to_their_locus_without_any_nodes(self):
        # the gene could not be located in the graph, so there are no nodes in its context names at
        # all - nothing about the graph identifies which copy of the gene they flank
        _, contexts_to_loci = self._extract_for_gappy_contig(
            *self.build_gene_on_gappy_contig(located_in_graph=False))

        name = f'fallback_gene_nodes__match_1.0000_path_contigfallback_{CONTIG_NAME}_400_700'
        self.assertEqual(contexts_to_loci, {name: mc.GeneLocus(CONTIG_NAME, 400, 700)})

    def test_each_copy_of_a_gene_gets_its_own_locus(self):
        # two copies of one gene on one contig. their context names share everything up to _path_,
        # which is the part any key built out of the gene and its nodes would keep
        first = mc.GeneContigMatch(f'{CONTIG_NAME}\tfallback_gene\t401\t700\t100\t100')
        second = mc.GeneContigMatch(f'{CONTIG_NAME}\tfallback_gene\t701\t1000\t100\t100')

        _, contexts_to_loci = self._extract_for_gappy_contig(nx.DiGraph(), {}, [first, second])

        self.assertEqual(sorted(contexts_to_loci.values()),
                         [mc.GeneLocus(CONTIG_NAME, 400, 700), mc.GeneLocus(CONTIG_NAME, 700, 1000)])

    @staticmethod
    def build_gene_on_gappy_contig(located_in_graph=True):
        """A gene on a 1000bp gap-containing contig that has plenty of flanking sequence, located on
        a graph node that is a dead end and is too short to supply a context of context_len on
        either side - the situation the contig fallback exists for. With located_in_graph=False the
        gene could not be located in the graph at all."""
        gene_contig_match = mc.GeneContigMatch(f'{CONTIG_NAME}\tfallback_gene\t401\t700\t100\t100')
        if located_in_graph:
            gene_contig_match.nodes_list = ['1+']
            gene_contig_match.start_in_first_node = 10

        assembly_graph = nx.DiGraph()
        assembly_graph.add_node('1+', length=40)
        nodes_with_edges_and_sequences = {'1+': SeqRecord(Seq('A' * 40), id='1+')}
        return assembly_graph, nodes_with_edges_and_sequences, [gene_contig_match]

    def _extract_for_gappy_contig(self, assembly_graph, nodes_with_edges_and_sequences, genes_to_contigs,
                                  contigs_with_gaps=frozenset({CONTIG_NAME})):
        """Returns (gene_lengths, contexts_to_loci)."""
        return ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph,
                                                                     nodes_with_edges_and_sequences,
                                                                     genes_to_contigs, 12, FALLBACK_CONTEXT_LEN,
                                                                     self.fallback_in_paths_fasta,
                                                                     self.fallback_out_paths_fasta,
                                                                     CONTIGS_PATH, contigs_with_gaps)

    def _assert_contexts_were_sliced_out_of_the_contig(self, expected_header):
        # the gene is at 1-based 401-700 of the contig, so 400:700 0-based half-open. the flanks
        # abut it exactly on both sides.
        contig_seq = helper.get_contig_seq(CONTIG_NAME)
        with open(self.fallback_in_paths_fasta) as f:
            self.assertListEqual(f.readlines(), [expected_header, f'{contig_seq[300:400]}\n'],
                                 'Incoming context was not sliced out of the contig')
        with open(self.fallback_out_paths_fasta) as f:
            self.assertListEqual(f.readlines(), [expected_header, f'{contig_seq[700:800]}\n'],
                                 'Outgoing context was not sliced out of the contig')

    def test_contig_context_for_gene_on_gappy_contig(self):
        gene_lengths, contexts_to_loci = self._extract_for_gappy_contig(*self.build_gene_on_gappy_contig())

        self._assert_contexts_were_sliced_out_of_the_contig(
            f'>fallback_gene_nodes_1+_match_1.0000_path_contigfallback_{CONTIG_NAME}_400_700\n')
        self.assertDictEqual(gene_lengths, {'fallback_gene': 300}, 'Genes length dictionary is incorrect')

    def test_contig_context_for_gene_not_located_in_the_graph(self):
        # the gene has no nodes_list, so the graph is not consulted at all - being on a gap-containing
        # contig is enough to get contexts, and the nodes part of their name is left empty
        gene_lengths, contexts_to_loci = self._extract_for_gappy_contig(*self.build_gene_on_gappy_contig(located_in_graph=False))

        self._assert_contexts_were_sliced_out_of_the_contig(
            f'>fallback_gene_nodes__match_1.0000_path_contigfallback_{CONTIG_NAME}_400_700\n')
        self.assertDictEqual(gene_lengths, {'fallback_gene': 300}, 'Genes length dictionary is incorrect')

    def test_no_contig_context_when_the_contig_has_no_gaps(self):
        self._extract_for_gappy_contig(*self.build_gene_on_gappy_contig(), contigs_with_gaps=frozenset())

        # the fastas are still created (downstream minimap2 needs them to exist), but the graph
        # supplies no context here and the contig is not consulted, so they must stay empty
        with open(self.fallback_in_paths_fasta) as f:
            self.assertEqual(f.read(), '', 'Incoming context should not have been written')
        with open(self.fallback_out_paths_fasta) as f:
            self.assertEqual(f.read(), '', 'Outgoing context should not have been written')

    def test_contig_context_is_added_next_to_the_graph_contexts(self):
        # a gene on a gap-containing contig gets its contexts from the contig even when the graph
        # supplied contexts of its own
        assembly_graph = helper.get_assembly_graph()
        nodes_with_edges_and_sequences = helper.get_assembly_graph_nodes()
        genes_with_location_in_graph = helper.get_genes_with_location_in_graph()

        ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph, nodes_with_edges_and_sequences,
                                                              genes_with_location_in_graph, 12, self.context_len,
                                                              self.fallback_in_paths_fasta,
                                                              self.fallback_out_paths_fasta, CONTIGS_PATH,
                                                              frozenset({CONTIG_NAME}))

        contig_seq = helper.get_contig_seq(CONTIG_NAME)
        with open(self.fallback_in_paths_fasta) as f:
            records = dict(zip(*[iter(line.strip() for line in f)] * 2))
        self.assertIn('>test_gene_nodes_5+_match_1.0000_path_5+', records, 'The graph context is missing')
        self.assertEqual(
            records.get(f'>test_gene_nodes_5+_match_1.0000_path_contigfallback_{CONTIG_NAME}'
                        f'_{GENE_START}_{GENE_END}'),
            contig_seq[GENE_START - self.context_len:GENE_START], 'The contig context is missing')

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
    """The length written to the fasta is measured after covered_by_gene is trimmed off, and a
    context is written only if it is exactly context_len long."""

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

    def test_path_dropped_when_trimmed_length_is_below_context_len(self):
        # node is 30bp, but 25bp of it is covered by the gene -> only 5bp of real context remain,
        # which is short of context_len=10. The length check must not run on the untrimmed 30bp.
        node = '1+'
        records_dict = {node: FakeSeqRecord('A' * 30)}
        in_paths_lengths, node_locations = ecc.save_paths_to_fasta_io_paths_approach(
            [[node]], self.fasta_path, records_dict, 10, in_or_out='in', covered_by_gene=25)

        self.assertEqual(self._read_written_seqs(), [], 'A too-short (post-trim) context should not be written')
        self.assertEqual(in_paths_lengths, {})
        self.assertEqual(node_locations, {})

    def test_path_kept_when_trimmed_length_is_exactly_context_len(self):
        node = '1+'
        records_dict = {node: FakeSeqRecord('A' * 50)}
        in_paths_lengths, _ = ecc.save_paths_to_fasta_io_paths_approach(
            [[node]], self.fasta_path, records_dict, 40, in_or_out='in', covered_by_gene=10)

        written = self._read_written_seqs()
        self.assertEqual(len(written), 1)
        self.assertEqual(len(written[0]), 40)
        self.assertEqual(in_paths_lengths['1+'], 40)

    def test_trimmed_length_is_cut_down_to_context_len(self):
        node = '1+'
        records_dict = {node: FakeSeqRecord('A' * 200)}
        in_paths_lengths, _ = ecc.save_paths_to_fasta_io_paths_approach(
            [[node]], self.fasta_path, records_dict, 100, in_or_out='in', covered_by_gene=10)

        written = self._read_written_seqs()
        self.assertEqual(len(written), 1)
        self.assertEqual(len(written[0]), 100)
        self.assertEqual(in_paths_lengths['1+'], 100)

    def test_out_path_dropped_when_trimmed_length_is_below_context_len(self):
        node = '1+'
        records_dict = {node: FakeSeqRecord('A' * 30)}
        in_paths_lengths, _ = ecc.save_paths_to_fasta_io_paths_approach(
            [[node]], self.fasta_path, records_dict, 10, in_or_out='out', covered_by_gene=25)

        self.assertEqual(self._read_written_seqs(), [])
        self.assertEqual(in_paths_lengths, {})

    def test_covered_by_gene_zero_does_not_produce_empty_sequence(self):
        # regression for seq[:-0] evaluating to '' when covered_by_gene == 0
        node = '1+'
        records_dict = {node: FakeSeqRecord('A' * 50)}
        in_paths_lengths, _ = ecc.save_paths_to_fasta_io_paths_approach(
            [[node]], self.fasta_path, records_dict, 50, in_or_out='in', covered_by_gene=0)

        written = self._read_written_seqs()
        self.assertEqual(len(written), 1)
        self.assertEqual(len(written[0]), 50)
        self.assertEqual(in_paths_lengths['1+'], 50)

    def test_the_flanking_end_of_the_path_is_the_one_that_is_kept(self):
        # 'in' contexts end at the gene, so the last context_len bases are kept; 'out' contexts
        # start at the gene, so the first context_len bases are kept
        node = '1+'
        records_dict = {node: FakeSeqRecord('C' * 20 + 'G' * 20)}
        ecc.save_paths_to_fasta_io_paths_approach([[node]], self.fasta_path, records_dict, 10, in_or_out='in')
        ecc.save_paths_to_fasta_io_paths_approach([[node]], self.fasta_path, records_dict, 10, in_or_out='out')

        self.assertEqual(self._read_written_seqs(), ['G' * 10, 'C' * 10])


if __name__ == '__main__':
    unittest.main()
