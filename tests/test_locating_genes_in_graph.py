import unittest
from ginger import locating_genes_in_graph as lg

ASSEMBLY_DIR = 'test_files/assembly_dir'
GENE_PIDENT_FILTERING_TH = 0.9
GENES_PATH = 'test_files/genes.fasta'
TEMP_FOLDER = 'temp_folder'


class MyTestCase(unittest.TestCase):

    # find_genes_in_contigs is covered by tests in sequence_alignment_utils.py
    # locate_genes_in_graph is covered by testing get_genes_to_contigs_with_nodes_list which is the main function here
    def test_get_nodes_dict_from_fastg_file(self):
        nodes_with_edges_and_sequences = lg.get_nodes_dict_from_fastg_file(ASSEMBLY_DIR + '/assembly_graph.fastg')
        self.assertEqual(len(nodes_with_edges_and_sequences), 55897)

    # def test_get_genes_to_contigs_with_nodes_list(self): # TODO you are here
    #     # # Subject_10042055__Scaffold_1__Start_21742__End_22848	1107	0	1107	+	NODE_2634_length_3645_cov_4.273259	3645	850	1957	1094	1107	0	NM:i:13	ms:i:1042	AS:i:1042	nn:i:0	tp:A:S	cm:i:151	s1:i:962	de:f:0.0117	rl:i:0	cg:Z:1107M
    #     genes_to_contigs = [
    #         mc.GeneContigMatch(
    #             mc.PafLine("Subject_10042055__Scaffold_1__Start_21742__End_22848", "1107", "0", "1107", "+",
    #                        "NODE_2634_length_3645_cov_4", "3645", "850", "1957", "1094", "1107", "0", ""))
    #     ]

    #     genes_to_contigs
    #     lg.get_genes_to_contigs_with_nodes_list(genes_to_contigs: Iterator[
    #         mc.GeneContigMatch], assembly_graph: nx.DiGraph,
    #     nodes_sequences_dict: Dict[str, SeqIO.SeqRecord], assembly_dir: str)

    if __name__ == '__main__':
        unittest.main()
