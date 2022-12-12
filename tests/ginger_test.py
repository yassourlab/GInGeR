import unittest
from ginger import pipeline_utils as pu
import pickle

LAB = '/Volumes/netta.barak/lab'
METADATA_PATH = f'{LAB}/Tools/UnifiedHumanGastrointestinalGenome/genomes-all_metadata.tsv'


class GingerTest(unittest.TestCase):
    def test_write_context_level_output_to_csv(self):
        with open(f'{LAB}/Projects/hybrid_assembler/test_run/output_pickled.pkl', 'rb') as f:
            output = pickle.load(f)
        test_dict = {}
        counter = 0
        for gene_species_tuple, matches_list in output.items():
            test_dict[gene_species_tuple] = matches_list
            counter += 1
            if counter > 20:
                break
        with open('context_level_output.pkl', 'wb') as f:
            pickle.dump(test_dict, f)

        pu.write_context_level_output_to_csv(test_dict,
                                                 f'{LAB}/Projects/hybrid_assembler/GInGeR/test_outputs/context_level_output_test.csv')

    def test_aggregate_context_level_output_to_species_level_output_and_write_csv(self):
        pu.aggregate_context_level_output_to_species_level_output_and_write_csv(
            'test_files/context_level_output_test.csv', METADATA_PATH,
            f'{LAB}/Projects/hybrid_assembler/GInGeR/test_outputs/species_level_output_test.csv')

    # def test_get_command_parser(self):
    #     parser = ginger.get_command_parser()
    #     parsed_args = vars(
    #         parser.parse_args('short_reads_1 short_reads_2 genes_path output_folder --threads=32'.split(' ')))
    #     self.assertEqual(parsed_args['short_reads_1'], 'short_reads_1')
    #     self.assertEqual(parsed_args['short_reads_2'], 'short_reads_2')
    #     self.assertEqual(parsed_args['genes_path'], 'genes_path')
    #     self.assertEqual(parsed_args['output_folder'], 'output_folder')
    #     self.assertEqual(parsed_args['threads'], 32)


if __name__ == '__main__':
    unittest.main()
