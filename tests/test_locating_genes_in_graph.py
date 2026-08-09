import unittest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ginger import locating_genes_in_graph as lg
from ginger import matches_classes as mc

from tests import helper

TEST_FILES = helper.get_filedir()

CONTIG_NAME = 'NODE_1_length_1000_cov_1'
# node 2 and node 3 are consecutive in the graph and overlap by 55bp
NODE_SEQUENCES = {'1+': SeqRecord(Seq('T' * 100), id='1+'),
                  '2+': SeqRecord(Seq('A' * 145 + 'G' * 55), id='2+'),
                  '3+': SeqRecord(Seq('G' * 55 + 'C' * 245), id='3+')}
# the alignments of the nodes to the contig, as computed by map_nodes_to_contigs_w_gaps. node 1 sits
# at the start of the contig and node 2 (the first node of the second segment) 300bp into it - after
# a gap of unknown length
NODES_TO_CONTIGS_DF = pd.DataFrame([
    {'contig': CONTIG_NAME, 'contig_start': 0, 'contig_end': 100, 'node': 'EDGE_1_length_100_cov_1', 'score': 1.0,
     'strand': '+'},
    {'contig': CONTIG_NAME, 'contig_start': 300, 'contig_end': 500, 'node': 'EDGE_2_length_200_cov_1', 'score': 1.0,
     'strand': '+'},
    {'contig': CONTIG_NAME, 'contig_start': 445, 'contig_end': 745, 'node': 'EDGE_3_length_300_cov_1', 'score': 1.0,
     'strand': '+'}])
GAPPY_CONTIG_SEGMENTS = [['1+'], ['2+', '3+']]


def gene_contig_match(start, end):
    """start and end are 1-based inclusive, the way mmseqs2 reports them. GeneContigMatch turns them
    into 0-based half-open, so start_in_first_node below is one less than the number passed here
    minus the offset of the node the gene starts on."""
    return mc.GeneContigMatch(f'{CONTIG_NAME}\tgene\t{start}\t{end}\t100\t100')


class LocatingGenesInGraphTest(unittest.TestCase):
    # find_genes_in_contigs is covered by tests in sequence_alignment_utils.py
    # locate_genes_in_graph is covered by testing add_node_list_to_genes_to_contigs which is the main function here

    def test_add_node_list_to_genes_to_contigs_single_segment(self):
        matches = lg.add_node_list_to_genes_to_contigs([gene_contig_match(50, 250)],
                                                       {CONTIG_NAME: [['2+', '3+']]}, NODE_SEQUENCES,
                                                       NODES_TO_CONTIGS_DF)
        self.assertEqual([match.nodes_list for match in matches], [['2+', '3+']])
        self.assertEqual([match.start_in_first_node for match in matches], [49])

    def test_add_node_list_to_genes_to_contigs_gappy_contig(self):
        # the gene is at 350-500 in the contig, which is inside the second segment (anchored at 300,
        # 445bp long) and covered by both of its nodes
        matches = lg.add_node_list_to_genes_to_contigs([gene_contig_match(350, 500)],
                                                       {CONTIG_NAME: GAPPY_CONTIG_SEGMENTS}, NODE_SEQUENCES,
                                                       NODES_TO_CONTIGS_DF)
        self.assertEqual([match.nodes_list for match in matches], [['2+', '3+']])
        self.assertEqual([match.start_in_first_node for match in matches], [49])

    def test_add_node_list_to_genes_to_contigs_gappy_contig_first_segment(self):
        matches = lg.add_node_list_to_genes_to_contigs([gene_contig_match(20, 70)],
                                                       {CONTIG_NAME: GAPPY_CONTIG_SEGMENTS}, NODE_SEQUENCES,
                                                       NODES_TO_CONTIGS_DF)
        self.assertEqual([match.nodes_list for match in matches], [['1+']])
        self.assertEqual([match.start_in_first_node for match in matches], [19])

    def test_add_node_list_to_genes_to_contigs_gene_over_a_gap(self):
        # a gene spanning the join between the segments is covered by no segment, so it falls back to
        # the single node that aligned to its location in the contig
        matches = lg.add_node_list_to_genes_to_contigs([gene_contig_match(90, 350)],
                                                       {CONTIG_NAME: GAPPY_CONTIG_SEGMENTS}, NODE_SEQUENCES,
                                                       NODES_TO_CONTIGS_DF)
        self.assertEqual([match.nodes_list for match in matches], [['1+']])
        self.assertEqual([match.start_in_first_node for match in matches], [89])

    def test_add_node_list_to_genes_to_contigs_unanchored_segments(self):
        # none of the segments' first nodes aligned to the contig, so the single graph node that
        # aligned to the gene's location is used, as it was before segments were kept
        matches = lg.add_node_list_to_genes_to_contigs([gene_contig_match(350, 500)],
                                                       {CONTIG_NAME: [['4+'], ['5+']]}, NODE_SEQUENCES,
                                                       NODES_TO_CONTIGS_DF)
        self.assertEqual([match.nodes_list for match in matches], [['2+']])
        self.assertEqual([match.start_in_first_node for match in matches], [49])

    def test_add_node_list_to_genes_to_contigs_keeps_unlocated_gene_on_a_gappy_contig(self):
        # the gene is at 800-900 in the contig, where no segment could be anchored and no node
        # aligned - it is kept anyway, because its context will be taken from the contig
        matches = lg.add_node_list_to_genes_to_contigs([gene_contig_match(800, 900)],
                                                       {CONTIG_NAME: [['4+'], ['5+']]}, NODE_SEQUENCES,
                                                       NODES_TO_CONTIGS_DF, {CONTIG_NAME})
        self.assertEqual([match.nodes_list for match in matches], [None])
        self.assertEqual([match.start_in_first_node for match in matches], [None])

    def test_add_node_list_to_genes_to_contigs_drops_unlocated_gene_on_a_gap_free_contig(self):
        matches = lg.add_node_list_to_genes_to_contigs([gene_contig_match(800, 900)],
                                                       {CONTIG_NAME: [['4+'], ['5+']]}, NODE_SEQUENCES,
                                                       NODES_TO_CONTIGS_DF)
        self.assertEqual(matches, [])

    # def test_get_nodes_dict_from_fastg_file(self):
    #     nodes_with_edges_and_sequences = lg.get_nodes_dict_from_fastg_file(f'{TEST_FILES}/SPAdes/assembly_graph.fastg')
    #     self.assertEqual(len(nodes_with_edges_and_sequences), 2)


if __name__ == '__main__':
    unittest.main()
