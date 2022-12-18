import unittest
from ginger import locating_genes_in_graph as lg

from tests import helper

TEST_FILES = helper.get_filedir()

class LocatingGenesInGraphTest(unittest.TestCase):

    # find_genes_in_contigs is covered by tests in sequence_alignment_utils.py
    # locate_genes_in_graph is covered by testing get_genes_to_contigs_with_nodes_list which is the main function here
    def test_get_nodes_dict_from_fastg_file(self):
        nodes_with_edges_and_sequences = lg.get_nodes_dict_from_fastg_file(f'{TEST_FILES}/SPAdes/assembly_graph.fastg')
        self.assertEqual(len(nodes_with_edges_and_sequences), 2)

if __name__ == '__main__':
    unittest.main()
