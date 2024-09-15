import unittest
from ginger import reference_database_utils as rdu
from tests import helper

TEST_FILES = helper.get_filedir()

class MyTestCase(unittest.TestCase):
    def test_get_max_read_len(self):
        gz_fastq_file = f'{TEST_FILES}/ecoli_1K_1.fq.gz'
        fastq_file = f'{TEST_FILES}/ecoli_1K_1_3_reads.fq'

        # test non gzipped file
        num_reads, max_read_len = rdu.get_max_read_len(fastq_file)
        self.assertEqual(max_read_len, 100)
        self.assertEqual(num_reads, 3)

        # test gzipped file
        num_reads, max_read_len = rdu.get_max_read_len(gz_fastq_file)
        self.assertEqual(max_read_len, 100)
    def test_get_kmer_length_options(self):
        kraken_db = f'{TEST_FILES}/fake_kraken_dir'
        kmer_length_options = rdu.get_kmer_length_options(kraken_db)
        self.assertEqual(set(kmer_length_options), set([100, 150]))

    def test_get_list_of_top_species_by_bracken(self):
        bracken_output_path = f'{TEST_FILES}/bracken_out_test.txt'
        reads_fraction = 0.01
        top_species = rdu.get_list_of_top_species_by_bracken(bracken_output_path, reads_fraction)
        self.assertEqual(set(top_species), set(['Enterobacter cloacae','Salmonella enterica','Yersinia enterocolitica']))

if __name__ == '__main__':
    unittest.main()
