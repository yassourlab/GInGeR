import unittest
import tempfile
import os
import zipfile
from types import SimpleNamespace
from unittest.mock import patch
import pandas as pd
from ginger import reference_database_utils as rdu
from tests import helper

TEST_FILES = helper.get_filedir()

class MyTestCase(unittest.TestCase):
    def test_get_paired_reads_seqkit_stats_parsing(self):
        reads_1 = f'{TEST_FILES}/ecoli_1K_1.fq.gz'
        reads_2 = f'{TEST_FILES}/ecoli_1K_2.fq.gz'

        class Dummy:
            def __init__(self, stdout, stderr='', returncode=0):
                self.stdout = stdout
                self.stderr = stderr
                self.returncode = returncode

        fake_out = (
            'file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n'
            f'{reads_1}\tFASTQ\tDNA\t1,000\t100000\t100\t100.0\t151\n'
            f'{reads_2}\tFASTQ\tDNA\t1,000\t100000\t100\t100.0\t150\n'
        )

        with patch('ginger.reference_database_utils.run', return_value=Dummy(fake_out)):
            a1, m1, a2, m2 = rdu.get_paired_reads_seqkit_stats(reads_1, reads_2)
        self.assertEqual(a1, 100.0)
        self.assertEqual(m1, 151)
        self.assertEqual(a2, 100.0)
        self.assertEqual(m2, 150)
    def test_get_kmer_length_options(self):
        kraken_db = f'{TEST_FILES}/fake_kraken_dir'
        kmer_length_options = rdu.get_kmer_length_options(kraken_db)
        self.assertEqual(set(kmer_length_options), set([100, 150]))

    def test_get_species_passing_coverage_threshold(self):
        bracken_output_path = f'{TEST_FILES}/bracken_out_coverage_test.txt'
        metadata_path = f'{TEST_FILES}/uhgg_metadata_coverage_test.tsv'
        avg_sum = 200.0
        stats = rdu.get_species_coverage_stats(bracken_output_path, avg_sum, metadata_path, max_refs_per_species=2)
        passing = rdu.get_species_passing_coverage_threshold(stats, species_coverage_threshold=10)

        # Coverage calculations:
        # EC: 250000*200/4000000 = 12.5 (pass)
        # SE: 200000*200/5000000 = 8.0 (fail)
        # YE: 200000*200/3000000 = 13.33.. (pass)
        self.assertEqual(set(passing), set(['Enterobacter cloacae', 'Yersinia enterocolitica']))

    def test_get_species_coverage_stats(self):
        bracken_output_path = f'{TEST_FILES}/bracken_out_coverage_test.txt'
        metadata_path = f'{TEST_FILES}/uhgg_metadata_coverage_test.tsv'
        avg_sum = 200.0
        stats = rdu.get_species_coverage_stats(
            bracken_output_path, avg_sum, metadata_path, max_refs_per_species=2,
        )
        stats_by_name = stats.set_index('name')
        # EC: median(4000000, 4000000) = 4000000; coverage = 250000*200/4000000 = 12.5
        self.assertEqual(stats_by_name.loc['Enterobacter cloacae', 'estimated_genome_length'], 4000000)
        self.assertAlmostEqual(stats_by_name.loc['Enterobacter cloacae', 'estimated_coverage'], 12.5)
        # SE: median(5000000, 5000000) = 5000000; coverage = 200000*200/5000000 = 8.0
        self.assertEqual(stats_by_name.loc['Salmonella enterica', 'estimated_genome_length'], 5000000)
        self.assertAlmostEqual(stats_by_name.loc['Salmonella enterica', 'estimated_coverage'], 8.0)

    def test_get_distinct_minimizers_by_species(self):
        kraken_report_path = f'{TEST_FILES}/kraken_report_filter_test.txt'
        distinct_minimizers = rdu.get_distinct_minimizers_by_species(kraken_report_path)
        self.assertEqual(distinct_minimizers, {'SpeciesPass': 300000, 'SpeciesFail': 50000})

    def test_get_species_included_in_analysis_df(self):
        bracken_output_path = f'{TEST_FILES}/bracken_out_coverage_test.txt'
        kraken_report_path = f'{TEST_FILES}/kraken_report_minimizers_test.txt'
        metadata_path = f'{TEST_FILES}/uhgg_metadata_coverage_test.tsv'
        avg_sum = 200.0
        top_species = ['Enterobacter cloacae', 'Yersinia enterocolitica']
        stats = rdu.get_species_coverage_stats(bracken_output_path, avg_sum, metadata_path, max_refs_per_species=2)
        included = rdu.get_species_included_in_analysis_df(stats, kraken_report_path, top_species)
        self.assertEqual(set(included['name']), set(top_species))
        included_by_name = included.set_index('name')
        self.assertEqual(included_by_name.loc['Enterobacter cloacae', 'distinct_minimizers'], 400000)
        self.assertEqual(included_by_name.loc['Yersinia enterocolitica', 'distinct_minimizers'], 300000)
        self.assertAlmostEqual(included_by_name.loc['Enterobacter cloacae', 'estimated_coverage'], 12.5)

    def test_filter_kraken_report_by_distinct_kmer_count(self):
        kraken_report_path = f'{TEST_FILES}/kraken_report_filter_test.txt'
        metadata_path = f'{TEST_FILES}/kraken_report_filter_test_metadata.tsv'
        with tempfile.TemporaryDirectory() as tmpdir:
            filtered_path = os.path.join(tmpdir, 'filtered_report.txt')
            rdu.filter_kraken_report_by_distinct_kmer_count(kraken_report_path, filtered_path,
                                                             metadata_path, max_refs_per_species=1)
            filtered = pd.read_csv(filtered_path, sep='\t', header=None, names=rdu.KRAKEN_REPORT_COLS)

        # SpeciesPass: 300000 / 3000000 = 0.1 (pass); SpeciesFail: 50000 / 5000000 = 0.01 (fail)

        # non-species rows (U, R, G) are kept regardless of distinct_kmer_count
        self.assertEqual(set(filtered.loc[filtered['rank'] != 'S', 'taxid']), {0, 1, 100})
        # only the species row above the threshold survives
        self.assertEqual(set(filtered.loc[filtered['rank'] == 'S', 'taxid']), {1001})


class DownloadRetryTest(unittest.TestCase):
    """The reference download retries a failing `datasets` call N_ATTEMPTS times.

    An earlier version called time.sleep without importing time, so the first failure raised
    NameError and no retry ever happened - which on a multi-hour run threw away the whole reference
    database step over one hiccup. Retries are per chunk of accessions, because `datasets` fetches a
    chunk in a single request.
    """

    def _download_with_failing_datasets(self, n_failures, already_downloaded=(), genomes=('GCF_000001.1',)):
        """Run download_missing_references against a `datasets` that fails its first n_failures calls.

        Returns the accessions files each call was given, and the error that escaped (if any).
        """
        calls = []

        def failing_run(command, **kwargs):
            # the accessions file is written before the call and deleted after, so its contents have
            # to be captured here rather than inspected afterwards
            accessions_file = command.split('--inputfile ')[1].split(' ')[0]
            with open(accessions_file) as f:
                calls.append([line.strip() for line in f if line.strip()])
            if len(calls) <= n_failures:
                return SimpleNamespace(returncode=1, stdout='', stderr='NCBI is down')
            zip_path = command.split('--filename ')[1].strip()
            _write_datasets_zip(zip_path, calls[-1])
            return SimpleNamespace(returncode=0, stdout='', stderr='')

        with tempfile.TemporaryDirectory() as tmpdir:
            for genome in already_downloaded:
                with open(rdu.reference_fasta_path(tmpdir, genome), 'w') as f:
                    f.write(f'>{genome}_contig_1\nACGT\n')
            with patch.object(rdu, 'N_ATTEMPTS', 3), patch.object(rdu, 'SLEEP_SECS', 0), \
                    patch.object(rdu, 'run', failing_run):
                error = None
                try:
                    rdu.download_missing_references(list(genomes), tmpdir)
                except Exception as e:
                    error = e
        return calls, error

    def test_retries_until_the_download_succeeds(self):
        calls, error = self._download_with_failing_datasets(n_failures=2)
        self.assertIsNone(error)
        self.assertEqual(len(calls), 3)  # two failures, then the one that worked

    def test_raises_the_download_error_after_the_last_attempt(self):
        calls, error = self._download_with_failing_datasets(n_failures=99)
        self.assertIsInstance(error, RuntimeError)
        self.assertEqual(len(calls), 3)  # N_ATTEMPTS, patched down from 10

    def test_does_not_download_a_reference_that_is_already_there(self):
        calls, error = self._download_with_failing_datasets(
            n_failures=99, already_downloaded=['GCF_000001.1'])
        self.assertIsNone(error)
        self.assertEqual(calls, [])

    def test_only_the_missing_accessions_are_requested(self):
        calls, error = self._download_with_failing_datasets(
            n_failures=0, already_downloaded=['GCF_000001.1'],
            genomes=['GCF_000001.1', 'GCF_000002.1'])
        self.assertIsNone(error)
        self.assertEqual(calls, [['GCF_000002.1']])

    def test_accessions_are_requested_in_chunks(self):
        genomes = [f'GCF_{i:06d}.1' for i in range(5)]
        with patch.object(rdu, 'DOWNLOAD_CHUNK_SIZE', 2):
            calls, error = self._download_with_failing_datasets(n_failures=0, genomes=genomes)
        self.assertIsNone(error)
        # 5 accessions in chunks of 2, and every accession asked for exactly once
        self.assertEqual([len(c) for c in calls], [2, 2, 1])
        self.assertEqual(sorted(a for call in calls for a in call), genomes)


def _write_datasets_zip(zip_path, accessions):
    """A stand-in for what `datasets download genome` produces: one .fna per accession under
    ncbi_dataset/data/, alongside metadata members that must be ignored."""
    with zipfile.ZipFile(zip_path, 'w') as archive:
        archive.writestr('README.md', 'not a genome')
        archive.writestr('ncbi_dataset/data/dataset_catalog.json', '{}')
        for accession in accessions:
            archive.writestr(f'ncbi_dataset/data/{accession}/{accession}_genomic.fna',
                             f'>{accession}_contig_1 some description\nACGT\n')


class DatasetsArchiveTest(unittest.TestCase):
    """Unpacking a datasets archive into the one-file-per-accession layout the rest of the pipeline
    expects."""

    def test_extracts_one_fasta_per_accession(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, 'chunk.zip')
            _write_datasets_zip(zip_path, ['GCF_000001.1', 'GCA_000002.1'])
            extracted = rdu.extract_genomes_from_datasets_zip(zip_path, tmpdir)

        self.assertEqual(extracted, {'GCF_000001.1', 'GCA_000002.1'})

    def test_concatenates_an_assembly_split_over_several_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, 'chunk.zip')
            with zipfile.ZipFile(zip_path, 'w') as archive:
                archive.writestr('ncbi_dataset/data/GCF_000001.1/chr1.fna', '>c1\nAAAA\n')
                archive.writestr('ncbi_dataset/data/GCF_000001.1/chr2.fna', '>c2\nCCCC\n')
            rdu.extract_genomes_from_datasets_zip(zip_path, tmpdir)
            with open(rdu.reference_fasta_path(tmpdir, 'GCF_000001.1')) as f:
                content = f.read()

        self.assertEqual(content, '>c1\nAAAA\n>c2\nCCCC\n')

    def test_an_accession_ncbi_does_not_return_is_skipped_rather_than_failing(self):
        """GTDB's metadata outlives NCBI's suppressions, so an accession that comes back empty is
        dropped with a warning instead of killing a multi-hour run."""
        calls = []

        def run_returning_one_of_two(command, **kwargs):
            calls.append(command)
            zip_path = command.split('--filename ')[1].strip()
            _write_datasets_zip(zip_path, ['GCF_000001.1'])  # GCF_000002.1 is not in the archive
            return SimpleNamespace(returncode=0, stdout='', stderr='')

        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.object(rdu, 'run', run_returning_one_of_two):
                available = rdu.download_missing_references(['GCF_000001.1', 'GCF_000002.1'], tmpdir)

        self.assertEqual(available, {'GCF_000001.1'})


class ContigToGenomeMapTest(unittest.TestCase):
    """The merged reference fasta and the contig->genome map are written together, so that every
    contig in the database can be traced back to the genome it came from."""

    def test_records_every_contig_of_every_genome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(rdu.reference_fasta_path(tmpdir, 'GCF_000001.1'), 'w') as f:
                f.write('>NZ_CP007265.1 Escherichia coli chromosome\nACGT\n>NZ_CP007266.1 plasmid\nTTTT\n')
            merged = os.path.join(tmpdir, 'merged.fasta')
            map_path = os.path.join(tmpdir, 'map.tsv')
            with open(merged, 'w') as merged_f, open(map_path, 'w') as map_f:
                rdu.write_genome_to_merged_fasta('GCF_000001.1', tmpdir, merged_f, map_f)
            with open(map_path) as f:
                rows = [line.strip().split('\t') for line in f if line.strip()]
            with open(merged) as f:
                merged_content = f.read()

        # keyed by the first whitespace-delimited token, which is what minimap2 reports as tname
        self.assertEqual(rows, [['NZ_CP007265.1', 'GCF_000001.1'], ['NZ_CP007266.1', 'GCF_000001.1']])
        self.assertIn('>NZ_CP007265.1 Escherichia coli chromosome', merged_content)


if __name__ == '__main__':
    unittest.main()
