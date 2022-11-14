import unittest
import locating_genes_in_graph as lg

ASSEMBLY_DIR = 'test_files/assembly_dir'
GENE_PIDENT_FILTERING_TH = 0.9
GENES_PATH = 'test_files/genes.fasta'
TEMP_FOLDER = 'temp_folder'


class MyTestCase(unittest.TestCase):

    # def test_locate_genes_in_graph(self):
    #     lg.locate_genes_in_graph(ASSEMBLY_DIR, GENE_PIDENT_FILTERING_TH, GENES_PATH, 1, TEMP_FOLDER)
    def test_get_records_dict_from_assembly_graph(self):
        nodes_with_edges_and_sequences = lg.get_nodes_dict_from_fastg_file(ASSEMBLY_DIR + '/assembly_graph.fastg')
        self.assertEqual(len(nodes_with_edges_and_sequences), 55897)

    def test_find_genes_in_contigs(self):
        pass
    # TODO you are here - continue this test


if __name__ == '__main__':
    unittest.main()
