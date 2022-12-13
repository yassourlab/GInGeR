import unittest
import pipeline_utils as pu


class MyTestCase(unittest.TestCase):
    def test_aggregate_context_level_output_to_species_level_output_and_write_csv(self):
        context_level_output_path = '../e2e_test_output/context_level_matches.csv'
        metadata_path = '../UHGG-metadata.tsv'
        species_level_output_path = '../e2e_test_output/context_level_matches.csv'
        pu.aggregate_context_level_output_to_species_level_output_and_write_csv(context_level_output_path,
                                                                                metadata_path,
                                                                                species_level_output_path, 10)


if __name__ == '__main__':
    unittest.main()
