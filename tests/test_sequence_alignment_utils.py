import unittest
from ginger import sequence_alignment_utils as sau
from ginger import matches_classes as mc

from tests import helper

TEST_FILES = helper.get_filedir()


class SequenceAlignmentUtilsTest(unittest.TestCase):

    def test_read_and_filter_mmseq2_matches(self):
        genes_to_contigs_path = f'{TEST_FILES}/genes_to_contigs.m8'
        pident_filtering_th = 0.9
        genes_to_contigs = sau.read_and_filter_mmseq2_matches(mc.GeneContigMatch, genes_to_contigs_path,
                                                               pident_filtering_th)
        self.assertEqual(len(genes_to_contigs), 1)
        match = genes_to_contigs[0]
        self.assertEqual(match.contig, 'NODE_1_length_1000_cov_140.620106')
        self.assertEqual(match.start, 337)
        self.assertEqual(match.end, 615)
        self.assertEqual(match.gene, 'test_gene')
        self.assertEqual(match.score, 1)

    def test_read_and_filter_minimap_matches(self):
        genes_to_contigs_path = f'{TEST_FILES}/in_paths_to_reference.paf'
        pident_filtering_th = 0.9
        genes_to_contigs = sau.read_and_filter_minimap_matches(mc.PathRefGenomeMatch, genes_to_contigs_path,
                                                               pident_filtering_th, {})
        self.assertEqual(len(genes_to_contigs), 2)




if __name__ == '__main__':
    unittest.main()
