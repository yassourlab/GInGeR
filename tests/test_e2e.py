import unittest
import ginger.pipeline_utils as pu

class MyTestCase(unittest.TestCase):
    def test_ginger_e2e_skip_kraken(self):
        pu.is_contig_name_func('fafdasfds')
        # ginger.run_ginger_e2e('test_files/ecoli_1K_1.fq.gz test_files/ecoli_1K_2.fq.gz test_files/test_gene.fasta e2e_temp_output --merged-filtered-fasta tests/test_files/merged_filtered_ref_db.fasta.gz'.split(' '))


if __name__ == '__main__':
    unittest.main()
