import unittest
import sequence_alignment_utils as sau
import matches_classes as mc


class MyTestCase(unittest.TestCase):
    def test_read_and_filter_minimap_matches(self):
        genes_to_contigs_path = '/Volumes/netta.barak/lab/Projects/hybrid_assembler/GInGeR/tests/test_files/temp_dir/genes_to_contigs.paf'
        pident_filtering_th = 0.9
        genes_to_contigs = sau.read_and_filter_minimap_matches(mc.GeneContigMatch, genes_to_contigs_path,
                                                               pident_filtering_th,
                                                               verbose=True)
        genes_to_contigs_list = list(genes_to_contigs)
        self.assertEqual(len(genes_to_contigs_list), 4)


if __name__ == '__main__':
    unittest.main()
