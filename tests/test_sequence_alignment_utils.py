import unittest
from ginger import sequence_alignment_utils as sau
from ginger import matches_classes as mc


class SequenceAlignmentUtilsTest(unittest.TestCase):
    def test_read_and_filter_minimap_matches(self):
        genes_to_contigs_path = 'data/genes_to_contigs.paf'
        pident_filtering_th = 0.9
        genes_to_contigs = sau.read_and_filter_minimap_matches(mc.GeneContigMatch, genes_to_contigs_path,
                                                               pident_filtering_th)
        self.assertEqual(len(genes_to_contigs), 1)
        match = genes_to_contigs[0]
        self.assertEqual(match.contig, 'NODE_1_length_1000_cov_140.620106')
        self.assertEqual(match.contig_start, 400)
        self.assertEqual(match.contig_end, 600)
        self.assertEqual(match.gene, 'test_gene')
        self.assertEqual(match.gene_start, 0)
        self.assertEqual(match.gene_end, 200)
        self.assertEqual(match.match_score, 1)


if __name__ == '__main__':
    unittest.main()
