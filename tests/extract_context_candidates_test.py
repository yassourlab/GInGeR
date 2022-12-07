import unittest
import extract_contexts_candidates as ecc
import locating_genes_in_graph as lg
import sequence_alignment_utils as sau
import constants as c
import matches_classes as mc
import pyfastg


class TestExtractContextsCandidates(unittest.TestCase):
    def test_extract_all_in_out_paths_and_write_them_to_fastas(self):
        # TODO - you are here - find this bug
        # out_dir = 'test_outputs'
        # assembly_graph_path = 'assembly_dir/assembly_graph.fastg'
        out_dir = '/ginger_output_'
        assembly_dir = '/ginger_output_/SPAdes'
        assembly_graph_path = f'{assembly_dir}/assembly_graph.fastg'
        genes_to_contigs_path = '/ginger_output_/genes_to_contigs.paf'
        assembly_graph = pyfastg.parse_fastg(assembly_graph_path)
        nodes_with_edges_and_sequences = lg.get_nodes_dict_from_fastg_file(assembly_graph_path)
        genes_to_contigs = sau.read_and_filter_minimap_matches(mc.GeneContigMatch, genes_to_contigs_path,
                                                               0.9)
        genes_with_location_in_graph = lg.get_genes_to_contigs_with_nodes_list(genes_to_contigs, assembly_graph,
                                                                               nodes_with_edges_and_sequences,
                                                                               assembly_dir)
        in_paths_fasta = c.IN_PATHS_FASTA_TEMPLATE.format(temp_folder=out_dir)
        out_paths_fasta = c.OUT_PATHS_FASTA_TEMPLATE.format(temp_folder=out_dir)
        gene_lengths = ecc.extract_all_in_out_paths_and_write_them_to_fastas(assembly_graph,
                                                                             nodes_with_edges_and_sequences,
                                                                             genes_with_location_in_graph,
                                                                             12,
                                                                             2500,
                                                                             in_paths_fasta, out_paths_fasta)
        print(gene_lengths)


if __name__ == '__main__':
    unittest.main()
