import unittest
from ginger import locating_genes_in_graph as lg
import pyfastg

from tests import helper

TEST_FILES = helper.get_filedir()

class LocatingGenesInGraphTest(unittest.TestCase):
    pass
    # find_genes_in_contigs is covered by tests in sequence_alignment_utils.py
    # locate_genes_in_graph is covered by testing get_genes_to_contigs_with_nodes_list which is the main function here
    # def test_get_nodes_dict_from_fastg_file(self):
    #     nodes_with_edges_and_sequences = lg.get_nodes_dict_from_fastg_file(f'{TEST_FILES}/SPAdes/assembly_graph.fastg')
    #     self.assertEqual(len(nodes_with_edges_and_sequences), 2)

    # def test_locate_genes_in_graph(self):
    #     assembly_dir = f'/sci/backup/morani/lab/Projects/GInGeR/analysis_data/SPAdes_short_3M_long_0/'
    #     assembly_graph_path = '/sci/backup/morani/lab/Projects/GInGeR/analysis_data/SPAdes_short_3M_long_0/assembly_graph.fastg'
    #     genes_to_contigs_path = f'/sci/backup/morani/lab/Projects/GInGeR/mlflow/mlruns/611452708858167350/a86cebb0c7184c05826386ecf50f2207/artifacts/genes_to_contigs.m8'
    #     assembly_graph_nodes = lg.get_nodes_dict_from_fastg_file(assembly_graph_path)
    #     genes_to_contigs = lg.find_genes_in_contigs(None, None, None, 0.9, genes_to_contigs_path)
    #     print(f'found {len(set([m.gene for m in genes_to_contigs]))} genes in the assembly graph')
    #
    #     if genes_to_contigs is None:
    #         return None, None, None
    #
    #     assembly_graph = pyfastg.parse_fastg(assembly_graph_path)
    #     problematic_genes = {'411481.BIFADO_00245', 'V1.UC39-4_GL0000256','MH0012_GL0211130'}
    #     problematic_genes_to_contigs = [m for m in genes_to_contigs if m.gene in problematic_genes]
    #     genes_with_location_in_graph = lg.get_genes_to_contigs_with_nodes_list(problematic_genes_to_contigs, assembly_graph,
    #                                                                         assembly_graph_nodes, assembly_dir)
    #     print(len(genes_with_location_in_graph))
    #

if __name__ == '__main__':
    unittest.main()
