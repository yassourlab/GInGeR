import unittest
from shutil import rmtree
import os
import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pafpy import PafRecord
import ginger.extract_contexts_candidates as ecc
from ginger import matches_classes as mc
from ginger import pipeline_utils as pu
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

GRAPH_CONTEXT_NAME = f'test_gene|{CONTIG_NAME}|{GENE_START}|{GENE_END}|1.0000|5+|5+'


def geometry_for(node_sequences):
    return pu.PathGeometry(node_sequences)


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
        rmtree(cls.test_outputs_dir)

    def test_extract_all_in_out_paths_and_write_them_to_fastas(self):
        gene_lengths = ecc.extract_all_in_out_paths_and_write_them_to_fastas(
            helper.get_assembly_graph(), helper.get_geometry(), helper.get_genes_with_location_in_graph(),
            12, self.context_len, self.in_paths_fasta, self.out_paths_fasta, CONTIGS_PATH)
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
                                 [f'>{GRAPH_CONTEXT_NAME}|in\n',
                                  f'{contig_seq[GENE_START - self.context_len:GENE_START]}\n'],
                                 'Wrong incoming contexts extracted')

        with open(self.out_paths_fasta) as f:
            lines = f.readlines()
            self.assertEqual(len(lines[1]) - 1, self.context_len, 'Outgoing context are not of correct length')
            self.assertListEqual(lines,
                                 [f'>{GRAPH_CONTEXT_NAME}|out\n',
                                  f'{contig_seq[GENE_END:GENE_END + self.context_len]}\n'],
                                 'Wrong outgoing contexts extracted')

        self.assertDictEqual(gene_lengths, {'test_gene': 279}, 'Genes length dictionary is incorrect')

    def test_the_outgoing_context_starts_at_the_end_of_the_alignment_not_of_the_reference_gene(self):
        # test_gene is a 93aa fragment whose alignment covers the contig exactly, so pad gene_length
        # out to a reference protein that is 60bp longer than the aligned span. the outgoing context
        # must not move: it follows the alignment, and gene_length says nothing about the contig.
        genes_with_location_in_graph = helper.get_genes_with_location_in_graph()
        for gene_contigs_match in genes_with_location_in_graph:
            gene_contigs_match.gene_length += 60

        ecc.extract_all_in_out_paths_and_write_them_to_fastas(
            helper.get_assembly_graph(), helper.get_geometry(), genes_with_location_in_graph, 12, self.context_len,
            self.in_paths_fasta, self.out_paths_fasta, CONTIGS_PATH)

        contig_seq = helper.get_contig_seq(CONTIG_NAME)
        with open(self.out_paths_fasta) as f:
            self.assertEqual(f.readlines()[1].strip(), contig_seq[GENE_END:GENE_END + self.context_len],
                             'The outgoing context was placed using gene_length instead of aligned_length')

    # the locus a context was cut from, which is carried in its name

    def test_a_contexts_name_carries_the_locus_it_was_cut_from(self):
        ecc.extract_all_in_out_paths_and_write_them_to_fastas(
            helper.get_assembly_graph(), helper.get_geometry(), helper.get_genes_with_location_in_graph(),
            12, self.context_len, self.in_paths_fasta, self.out_paths_fasta, CONTIGS_PATH)

        with open(self.in_paths_fasta) as f:
            header = f.readline().strip().lstrip('>')
        match = mc.PathRefGenomeMatch(PafRecord.from_str(f'{header}\t100\t0\t100\t+\tref\t5000\t10\t110\t100\t100\t60'),
                                      {})
        self.assertEqual(match.gene, 'test_gene')
        self.assertEqual(match.locus, mc.GeneLocus(CONTIG_NAME, GENE_START, GENE_END))
        self.assertEqual(match.side, 'in')

    def test_contig_contexts_carry_their_locus_without_any_nodes(self):
        # the gene could not be located in the graph, so there are no nodes in its context names at
        # all - nothing about the graph identifies which copy of the gene they flank
        self._extract_for_gappy_contig(*self.build_gene_on_gappy_contig(located_in_graph=False))

        self._assert_contexts_were_sliced_out_of_the_contig(
            f'>fallback_gene|{CONTIG_NAME}|400|700|1.0000||contigfallback')

    def test_each_copy_of_a_gene_gets_its_own_contexts(self):
        # two copies of one gene on one contig, located on the same node. everything a name used to
        # be built from - the gene and its nodes - is identical between them, so only the locus tells
        # them apart. their contexts flank different places and both have to be written
        assembly_graph, node_sequences, [first] = self.build_gene_on_gappy_contig()
        second = mc.GeneContigMatch(f'{CONTIG_NAME}\tfallback_gene\t601\t900\t100\t100')
        second.nodes_list, second.start_in_first_node = first.nodes_list, first.start_in_first_node

        self._extract_for_gappy_contig(assembly_graph, node_sequences, [first, second])

        with open(self.fallback_in_paths_fasta) as f:
            headers = [line.strip().lstrip('>') for line in f if line.startswith('>')]
        self.assertEqual(headers, [f'fallback_gene|{CONTIG_NAME}|400|700|1.0000|1+|contigfallback|in',
                                   f'fallback_gene|{CONTIG_NAME}|600|900|1.0000|1+|contigfallback|in'])

    def test_a_gene_name_containing_the_separator_survives(self):
        # SARG's and CARD's gene names have '|' in them, and the name is taken apart from the right
        gene_contig_match = mc.GeneContigMatch(f'{CONTIG_NAME}\tgb|AAA|blaTEM-1\t401\t700\t100\t100')
        name = ecc.context_name(gene_contig_match, '5+', 'in')

        match = mc.PathRefGenomeMatch(PafRecord.from_str(f'{name}\t100\t0\t100\t+\tref\t5000\t10\t110\t100\t100\t60'),
                                      {})
        self.assertEqual(match.gene, 'gb|AAA|blaTEM-1')
        self.assertEqual(match.locus, mc.GeneLocus(CONTIG_NAME, 400, 700))

    # contexts sliced straight out of a gap-containing contig

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
        node_sequences = {'1+': SeqRecord(Seq('A' * 40), id='1+')}
        return assembly_graph, node_sequences, [gene_contig_match]

    def _extract_for_gappy_contig(self, assembly_graph, node_sequences, genes_to_contigs,
                                  contigs_with_gaps=frozenset({CONTIG_NAME})):
        return ecc.extract_all_in_out_paths_and_write_them_to_fastas(
            assembly_graph, geometry_for(node_sequences), genes_to_contigs, 12, FALLBACK_CONTEXT_LEN,
            self.fallback_in_paths_fasta, self.fallback_out_paths_fasta, CONTIGS_PATH, contigs_with_gaps)

    def _assert_contexts_were_sliced_out_of_the_contig(self, expected_name):
        # the gene is at 1-based 401-700 of the contig, so 400:700 0-based half-open. the flanks
        # abut it exactly on both sides.
        contig_seq = helper.get_contig_seq(CONTIG_NAME)
        with open(self.fallback_in_paths_fasta) as f:
            self.assertListEqual(f.readlines(), [f'{expected_name}|in\n', f'{contig_seq[300:400]}\n'],
                                 'Incoming context was not sliced out of the contig')
        with open(self.fallback_out_paths_fasta) as f:
            self.assertListEqual(f.readlines(), [f'{expected_name}|out\n', f'{contig_seq[700:800]}\n'],
                                 'Outgoing context was not sliced out of the contig')

    def test_contig_context_for_gene_on_gappy_contig(self):
        gene_lengths = self._extract_for_gappy_contig(*self.build_gene_on_gappy_contig())

        self._assert_contexts_were_sliced_out_of_the_contig(
            f'>fallback_gene|{CONTIG_NAME}|400|700|1.0000|1+|contigfallback')
        self.assertDictEqual(gene_lengths, {'fallback_gene': 300}, 'Genes length dictionary is incorrect')

    def test_contig_context_for_gene_not_located_in_the_graph(self):
        # the gene has no nodes_list, so the graph is not consulted at all - being on a gap-containing
        # contig is enough to get contexts, and the nodes part of their name is left empty
        gene_lengths = self._extract_for_gappy_contig(*self.build_gene_on_gappy_contig(located_in_graph=False))

        self._assert_contexts_were_sliced_out_of_the_contig(
            f'>fallback_gene|{CONTIG_NAME}|400|700|1.0000||contigfallback')
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
        ecc.extract_all_in_out_paths_and_write_them_to_fastas(
            helper.get_assembly_graph(), helper.get_geometry(), helper.get_genes_with_location_in_graph(),
            12, self.context_len, self.fallback_in_paths_fasta, self.fallback_out_paths_fasta, CONTIGS_PATH,
            frozenset({CONTIG_NAME}))

        contig_seq = helper.get_contig_seq(CONTIG_NAME)
        with open(self.fallback_in_paths_fasta) as f:
            records = dict(zip(*[iter(line.strip() for line in f)] * 2))
        self.assertIn(f'>{GRAPH_CONTEXT_NAME}|in', records, 'The graph context is missing')
        self.assertEqual(
            records.get(f'>test_gene|{CONTIG_NAME}|{GENE_START}|{GENE_END}|1.0000|5+|contigfallback|in'),
            contig_seq[GENE_START - self.context_len:GENE_START], 'The contig context is missing')


class TestSavePathsToFastaIoPathsApproach(unittest.TestCase):
    """The length written to the fasta is measured after covered_by_gene is trimmed off, and a
    context is written only if it is exactly context_len long."""

    def setUp(self):
        self.test_outputs_dir = 'save_paths_tests_output'
        os.mkdir(self.test_outputs_dir)
        self.fasta_path = f'{self.test_outputs_dir}/test.fasta'
        self.gene_contigs_match = mc.GeneContigMatch(f'{CONTIG_NAME}\tgene\t401\t700\t100\t100')
        self.gene_contigs_match.nodes_list = ['1+']

    def tearDown(self):
        rmtree(self.test_outputs_dir)

    def _write(self, node_seq, context_len, in_or_out, covered_by_gene):
        geometry = geometry_for({'1+': SeqRecord(Seq(node_seq), id='1+')})
        return ecc.save_paths_to_fasta_io_paths_approach([['1+']], self.fasta_path, geometry, context_len,
                                                         in_or_out, covered_by_gene, self.gene_contigs_match)

    def _read_written_seqs(self):
        if not os.path.exists(self.fasta_path):
            return []
        with open(self.fasta_path) as f:
            lines = f.readlines()
        return [lines[i + 1].strip() for i in range(0, len(lines), 2)]

    def test_path_dropped_when_trimmed_length_is_below_context_len(self):
        # node is 30bp, but 25bp of it is covered by the gene -> only 5bp of real context remain,
        # which is short of context_len=10. The length check must not run on the untrimmed 30bp.
        written = self._write('A' * 30, 10, 'in', 25)

        self.assertEqual(self._read_written_seqs(), [], 'A too-short (post-trim) context should not be written')
        self.assertEqual(written, [])

    def test_path_kept_when_trimmed_length_is_exactly_context_len(self):
        written = self._write('A' * 50, 40, 'in', 10)

        self.assertEqual([len(seq) for seq in self._read_written_seqs()], [40])
        self.assertEqual(len(written), 1)

    def test_trimmed_length_is_cut_down_to_context_len(self):
        self._write('A' * 200, 100, 'in', 10)

        self.assertEqual([len(seq) for seq in self._read_written_seqs()], [100])

    def test_out_path_dropped_when_trimmed_length_is_below_context_len(self):
        written = self._write('A' * 30, 10, 'out', 25)

        self.assertEqual(self._read_written_seqs(), [])
        self.assertEqual(written, [])

    def test_covered_by_gene_zero_does_not_produce_empty_sequence(self):
        # regression for seq[:-0] evaluating to '' when covered_by_gene == 0
        self._write('A' * 50, 50, 'in', 0)

        self.assertEqual([len(seq) for seq in self._read_written_seqs()], [50])

    def test_the_flanking_end_of_the_path_is_the_one_that_is_kept(self):
        # 'in' contexts end at the gene, so the last context_len bases are kept; 'out' contexts
        # start at the gene, so the first context_len bases are kept
        self._write('C' * 20 + 'G' * 20, 10, 'in', 0)
        self._write('C' * 20 + 'G' * 20, 10, 'out', 0)

        self.assertEqual(self._read_written_seqs(), ['G' * 10, 'C' * 10])


if __name__ == '__main__':
    unittest.main()
