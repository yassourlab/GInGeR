import unittest
import ginger


class GingerTest(unittest.TestCase):
    def test_get_command_parser(self):
        parser = ginger.get_command_parser()
        parsed_args = vars(
            parser.parse_args('short_reads_1 short_reads_2 genes_path output_folder --threads=32'.split(' ')))
        self.assertEqual(parsed_args['short_reads_1'], 'short_reads_1')
        self.assertEqual(parsed_args['short_reads_2'], 'short_reads_2')
        self.assertEqual(parsed_args['genes_path'], 'genes_path')
        self.assertEqual(parsed_args['output_folder'], 'output_folder')
        self.assertEqual(parsed_args['threads'], 32)


if __name__ == '__main__':
    unittest.main()
