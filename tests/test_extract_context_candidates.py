import unittest
from shutil import rmtree
import pickle
import os
from tests import helper

import ginger.extract_contexts_candidates as ecc


class TestExtractContextsCandidates(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data_dir = helper.get_test_data_dir()
        cls.test_outputs_dir = 'extraction_candidates_tests_output'
        cls.in_paths_fasta = f'{cls.test_outputs_dir}/test_in_paths_fasta.fasta'
        cls.out_paths_fasta = f'{cls.test_outputs_dir}/test_out_paths_fasta.fasta'
        os.mkdir(cls.test_outputs_dir)
        cls.context_len = 100

    @classmethod
    def tearDownClass(cls):
        pass
        rmtree(cls.test_outputs_dir)

    def test_extract_all_in_out_paths_and_write_them_to_fastas(self):
        print('os.getcwd()', os.getcwd())
        with open(f'{self.data_dir}/assembly_graph.pkl', 'rb') as f:
            assembly_graph = pickle.load(f)
        with open(f'{self.data_dir}/nodes_with_edges_and_sequences.pkl', 'rb') as f:
            nodes_with_edges_and_sequences = pickle.load(f)
        with open(f'{self.data_dir}/genes_with_location_in_graph.pkl', 'rb') as f:
            genes_with_location_in_graph = pickle.load(f)

        gene_lengths = ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph,
                                                                             nodes_with_edges_and_sequences,
                                                                             genes_with_location_in_graph,
                                                                             12,
                                                                             self.context_len,
                                                                             self.in_paths_fasta, self.out_paths_fasta)
        self.assertTrue(os.path.exists(self.in_paths_fasta))
        self.assertTrue(os.path.exists(self.out_paths_fasta))

        with open(self.in_paths_fasta) as f:
            lines = f.readlines()
            self.assertListEqual(lines,
                                 ['>test_gene_nodes_5+_path_5+\n',
                                  'CTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCG\n'],
                                 'Wrong incoming contexts extracted')
            # the -1 is beacuse of the \n in the end of the line
            self.assertEqual(len(lines[1]) - 1, self.context_len, 'Incoming context are not of correct length')

        with open(self.out_paths_fasta) as f:
            lines = f.readlines()

            self.assertListEqual(lines,
                                 ['>test_gene_nodes_5+_path_5+\n',
                                  'CAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGC\n'],
                                 'Wrong outgoing contexts extracted')
            self.assertEqual(len(lines[1]) - 1, self.context_len, 'Outgoing context are not of correct length')

        self.assertDictEqual(gene_lengths, {'test_gene': 200}, 'Genes length dictionary is incorrect')


if __name__ == '__main__':
    unittest.main()
