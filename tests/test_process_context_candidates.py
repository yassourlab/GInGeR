from Bio import SeqIO
import unittest
import pickle

from ginger import verify_context_candidates as pcc


class ProcessContextCandidatesTest(unittest.TestCase):

    def test_process_in_and_out_paths_to_results(self):
        in_path_mapping_to_bugs = f'test_files/in_paths_to_reference.paf'
        out_path_mapping_to_bugs = f'test_files/out_paths_to_reference.paf'

        genes_lengths = {'test_gene': 200}
        paths_pident_filtering_th = 0.9
        minimal_gap_ratio = 0
        maximal_gap_ratio = 1.5
        matches_per_gene_no_overlaps = pcc.process_in_and_out_paths_to_results(in_path_mapping_to_bugs,
                                                                               out_path_mapping_to_bugs, genes_lengths,
                                                                               paths_pident_filtering_th,
                                                                               minimal_gap_ratio, maximal_gap_ratio)

        with open('test_files/matches_per_gene_no_overlaps.pkl', 'rb') as f:
            expected_matches_per_gene_no_overlaps = pickle.load(f)
        self.assertEqual(len(matches_per_gene_no_overlaps), len(expected_matches_per_gene_no_overlaps))
        self.assertListEqual(list(matches_per_gene_no_overlaps.keys()),
                             list(expected_matches_per_gene_no_overlaps.keys()))
        match = matches_per_gene_no_overlaps['test_gene', 'MGYG000260594_1'][0]
        self.assertEqual(match.bug, 'MGYG000260594_1')
        self.assertEqual(match.gene, 'test_gene')
        self.assertEqual(match.gene_length, 200)
        self.assertEqual(match.match_score, 0.9775)
        self.assertEqual(match.start, 5077972)
        self.assertEqual(match.end, 5078172)


if __name__ == '__main__':
    unittest.main()