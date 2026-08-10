import unittest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ginger import locating_genes_in_graph as lg
from ginger import matches_classes as mc
from ginger import pipeline_utils as pu

from tests import helper

TEST_FILES = helper.get_filedir()

CONTIG_NAME = 'NODE_1_length_1000_cov_1'
# nodes 2 and 3 are consecutive in the graph and overlap by 55bp. node 4 is 56bp - too short for
# minimap2 to report, which is why it never appears in the placements below
NODE_SEQUENCES = {'1+': SeqRecord(Seq('T' * 100), id='1+'),
                  '2+': SeqRecord(Seq('A' * 145 + 'G' * 55), id='2+'),
                  '3+': SeqRecord(Seq('G' * 55 + 'C' * 245), id='3+'),
                  '4+': SeqRecord(Seq('C' * 245 + 'T' * 55), id='4+'),
                  '5+': SeqRecord(Seq('T' * 55 + 'A' * 145), id='5+')}
GAPPY_CONTIG_SEGMENTS = [['1+'], ['2+', '3+']]


def geometry():
    return pu.PathGeometry(NODE_SEQUENCES)


def placements(*nodes_and_origins):
    """{node: [NodePlacement]}, the shape map_nodes_to_contigs_w_gaps hands the locator."""
    by_node = {}
    for node, origin in nodes_and_origins:
        by_node.setdefault(node, []).append(mc.NodePlacement(node=node, origin=origin, score=1.0))
    return by_node


# node 1 sits at the start of the contig and node 2 (the first node of the second segment) 300bp
# into it - after a gap of unknown length. node 3 follows node 2, overlapping it by 55bp
CONTIG_PLACEMENTS = placements(('1+', 0), ('2+', 300), ('3+', 445))


def gene_contig_match(start, end):
    """start and end are 1-based inclusive, the way mmseqs2 reports them. GeneContigMatch turns them
    into 0-based half-open, so start_in_first_node below is one less than the number passed here
    minus the offset of the node the gene starts on."""
    return mc.GeneContigMatch(f'{CONTIG_NAME}\tgene\t{start}\t{end}\t100\t100')


def locate(match, segments, placements_by_node=None, contigs_with_gaps=frozenset()):
    return lg.add_node_list_to_genes_to_contigs(
        [match], {CONTIG_NAME: segments}, geometry(),
        {CONTIG_NAME: CONTIG_PLACEMENTS if placements_by_node is None else placements_by_node},
        contigs_with_gaps)


class LocatingGenesInGraphTest(unittest.TestCase):
    # find_genes_in_contigs is covered by tests in sequence_alignment_utils.py
    # locate_genes_in_graph is covered by testing add_node_list_to_genes_to_contigs which is the main function here

    # case 1 - a gene on a contig SPAdes assembled from a single graph path

    def test_gene_on_a_gap_free_contig_is_located_on_the_contigs_path(self):
        matches = locate(gene_contig_match(50, 250), [['2+', '3+']])
        self.assertEqual([match.nodes_list for match in matches], [['2+', '3+']])
        self.assertEqual([match.start_in_first_node for match in matches], [49])

    def test_gene_running_off_the_end_of_the_contigs_path_is_dropped(self):
        # the path spells out 445bp and the gene reaches 600, so it is not on these nodes and the
        # contexts cut from them would not flank it
        matches = locate(gene_contig_match(300, 600), [['2+', '3+']])
        self.assertEqual(matches, [])

    def test_gene_on_a_contig_with_no_path_is_dropped(self):
        matches = locate(gene_contig_match(50, 250), [])
        self.assertEqual(matches, [])

    # case 2 - a gene on a contig assembled from several graph paths joined by paired-end evidence

    def test_gene_on_a_gappy_contig_is_located_in_the_segment_that_contains_it(self):
        # the gene is at 350-500 in the contig, which is inside the second segment (anchored at 300,
        # 445bp long) and covered by both of its nodes
        matches = locate(gene_contig_match(350, 500), GAPPY_CONTIG_SEGMENTS, contigs_with_gaps={CONTIG_NAME})
        self.assertEqual([match.nodes_list for match in matches], [['2+', '3+']])
        self.assertEqual([match.start_in_first_node for match in matches], [49])

    def test_gene_on_a_gappy_contig_is_located_in_the_first_segment(self):
        matches = locate(gene_contig_match(20, 70), GAPPY_CONTIG_SEGMENTS, contigs_with_gaps={CONTIG_NAME})
        self.assertEqual([match.nodes_list for match in matches], [['1+']])
        self.assertEqual([match.start_in_first_node for match in matches], [19])

    def test_gene_spanning_a_gap_is_kept_but_not_located(self):
        # a gene crossing the join between two segments lies on nodes that are not adjacent in the
        # graph, so no segment holds its whole span. it is kept for the contig fallback to context
        matches = locate(gene_contig_match(90, 350), GAPPY_CONTIG_SEGMENTS, contigs_with_gaps={CONTIG_NAME})
        self.assertEqual([match.nodes_list for match in matches], [None])
        self.assertEqual([match.start_in_first_node for match in matches], [None])

    def test_unanchored_segments_leave_the_gene_unlocated(self):
        # none of the segments' nodes aligned to the contig, so there is nothing to place them by
        matches = locate(gene_contig_match(350, 500), [['4+'], ['5+']], placements_by_node={},
                         contigs_with_gaps={CONTIG_NAME})
        self.assertEqual([match.nodes_list for match in matches], [None])
        self.assertEqual([match.start_in_first_node for match in matches], [None])

    def test_gene_on_a_gappy_contig_outside_every_segment_is_kept_unlocated(self):
        # the gene is at 800-900 in the contig, past everything the segments cover - it is kept
        # anyway, because its context will be taken from the contig
        matches = locate(gene_contig_match(800, 900), GAPPY_CONTIG_SEGMENTS, contigs_with_gaps={CONTIG_NAME})
        self.assertEqual([match.nodes_list for match in matches], [None])
        self.assertEqual([match.start_in_first_node for match in matches], [None])

    # anchoring a segment in contig coordinates

    def test_segment_is_anchored_on_a_later_node_when_the_first_is_too_short_to_align(self):
        # node 4 is 56bp, below what minimap2 reports, so only node 5 has a placement. it starts
        # 245bp into the segment, so the segment starts 245bp before where node 5 aligned
        origin = lg.anchor_segment_in_contig(geometry(), ['4+', '5+'], placements(('5+', 545)))
        self.assertEqual(origin, 300)

    def test_segment_with_no_aligned_node_cannot_be_anchored(self):
        self.assertIsNone(lg.anchor_segment_in_contig(geometry(), ['4+', '5+'], placements(('1+', 0))))

    def test_agreeing_nodes_anchor_the_segment_together(self):
        # both nodes of the segment place it at 300, one directly and one through its 145bp offset
        origin = lg.anchor_segment_in_contig(geometry(), ['2+', '3+'], placements(('2+', 300), ('3+', 445)))
        self.assertEqual(origin, 300)

    def test_a_node_aligned_to_a_repeat_elsewhere_does_not_move_the_anchor(self):
        # node 3 also aligned 4000bp away; the two nodes that agree outvote it
        origin = lg.anchor_segment_in_contig(geometry(), ['2+', '3+'],
                                             placements(('2+', 300), ('3+', 445), ('3+', 4445)))
        self.assertEqual(origin, 300)

    def test_segment_whose_only_node_aligned_twice_cannot_be_anchored(self):
        # nothing says which of the two copies of the repeat this segment is
        self.assertIsNone(lg.anchor_segment_in_contig(geometry(), ['2+', '3+'],
                                                      placements(('3+', 445), ('3+', 4445))))

    # turning minimap2's output into placements

    def test_placements_resolve_orientation_and_undo_the_alignment_clipping(self):
        nodes_to_contigs_df = pd.DataFrame([
            # the forward record aligned on the minus strand, so the reverse complement is the node
            # that runs with the contig. its first 12bp were clipped, so it starts 12bp before the
            # alignment does
            {'qname': "EDGE_4_length_100_cov_1", 'qlen': 100, 'qstart': 12, 'strand': '-',
             'tname': CONTIG_NAME, 'tstart': 600, 'mlen': 100}])
        by_contig = lg.placements_from_minimap_results(nodes_to_contigs_df)
        self.assertEqual(list(by_contig[CONTIG_NAME]), ['4-'])
        self.assertEqual(by_contig[CONTIG_NAME]['4-'][0].origin, 588)

    def test_partially_aligned_nodes_are_not_placed(self):
        nodes_to_contigs_df = pd.DataFrame([
            {'qname': "EDGE_4_length_100_cov_1", 'qlen': 100, 'qstart': 0, 'strand': '+',
             'tname': CONTIG_NAME, 'tstart': 600, 'mlen': 50}])
        self.assertEqual(lg.placements_from_minimap_results(nodes_to_contigs_df), {})

    def test_node_oriented_with_contig(self):
        # every edge is in the fastg twice, as itself and as its reverse complement, and either can
        # be the record that aligned to the contig. what decides the orientation to walk the graph
        # in is the combination of which record it was and which strand it aligned on
        forward, reverse_complement = 'EDGE_4_length_100_cov_1', "EDGE_4_length_100_cov_1'"
        self.assertEqual(lg.node_oriented_with_contig(forward, '+'), '4+')
        self.assertEqual(lg.node_oriented_with_contig(forward, '-'), '4-')
        self.assertEqual(lg.node_oriented_with_contig(reverse_complement, '+'), '4-')
        self.assertEqual(lg.node_oriented_with_contig(reverse_complement, '-'), '4+')


if __name__ == '__main__':
    unittest.main()
