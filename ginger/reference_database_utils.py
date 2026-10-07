import gzip
import shutil
import tempfile
import urllib.request
import zipfile
from subprocess import run

import numpy as np
import pandas as pd
import logging
import csv
from ginger import pipeline_utils as pu
import os
import re
import time

log = logging.getLogger(__name__)
KRAKEN_COMMAND = 'kraken2 --db {kraken_db} {extra_args} --paired {reads_1} {reads_2} --threads {threads} --output {kraken_output} --report {kraken_report} --confidence 0.1 --use-names --report-minimizer-data'  # --report {report}
BRACKEN_COMMAND = 'bracken -d {kraken_db} -i {kraken_report} -o {bracken_output} -w {bracken_report} -r {read_len} -l S -t {min_reads_for_bracken}'
# NCBI's datasets CLI, which fetches reference genomes by assembly accession. An accessions
# file rather than one accession per invocation: the top references of every selected species
# come to a few thousand genomes, and that many round trips is both slow and far more likely to
# be throttled.
DATASETS_COMMAND = ('datasets download genome accession --inputfile {accessions_file} '
                    '--include genome --no-progressbar --filename {zip_path}')
# how many accessions go into one `datasets` call. Small enough that a failure re-fetches
# little, large enough that a few thousand genomes take tens of calls rather than thousands
DOWNLOAD_CHUNK_SIZE = 200
N_ATTEMPTS = 10
SLEEP_SECS = 60
URLOPEN_TIMEOUT = 60

# The reference catalogs GInGeR knows about. They differ in one thing that matters - where genome
# sequence comes from - plus the database and defaults that go with it. Everything downstream of the
# download is column-driven and shared. Paths are relative; ginger_runner resolves them.
REFERENCE_SOURCES = {
    'gtdb': dict(metadata='GTDB-metadata.tsv', kraken_db='kraken2_db_gtdb_r226',
                 references_dir='references_dir_gtdb',
                 kraken_extra_args='--memory-mapping', # as the database is very large
                 distinct_kmer_ratio_threshold=0.01), # GTDB has more species so the number of k-mers mapped uniqely to the specific clade is smaller in average
    'uhgg': dict(metadata='UHGG-metadata.tsv', kraken_db='kraken2_db_uhgg_v2.0.2',
                 references_dir='references_dir_uhgg',
                 kraken_extra_args='', # UHGG's is 15.5GB, so loading it into RAM is both possible and much faster than memory-mapping it would be
                 distinct_kmer_ratio_threshold=0.05),
}
DEFAULT_REFERENCE_SOURCE = 'gtdb'
BRACKEN_MIN_READS_RELAXATION_FACTOR = 0.5
DISTINCT_KMER_RATIO_THRESHOLD = 0.01
KRAKEN_REPORT_COLS = ['pct', 'reads_clade', 'reads_direct', 'kmer_count', 'distinct_kmer_count', 'rank', 'taxid', 'name']

def get_paired_reads_seqkit_stats(reads_1: str, reads_2: str):
    """Return (avg_len_r1, max_len_r1, avg_len_r2, max_len_r2) from `seqkit stats -T`."""
    out = run(['seqkit', 'stats', '-T', reads_1, reads_2], capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"seqkit stats failed. stderr: {out.stderr.strip()}")
    lines = [ln for ln in out.stdout.splitlines() if ln.strip()]
    if not lines:
        raise RuntimeError('seqkit stats returned empty output')

    reader = csv.DictReader(lines, delimiter='\t')
    if reader.fieldnames is None:
        raise RuntimeError('seqkit stats output missing header')
    for col in ['file', 'avg_len', 'max_len']:
        if col not in reader.fieldnames:
            raise RuntimeError(f"seqkit stats output missing column '{col}'. Header: {reader.fieldnames}")

    by_file = {row['file']: row for row in reader if row.get('file')}
    if reads_1 not in by_file or reads_2 not in by_file:
        raise RuntimeError(f"seqkit stats output missing one of the input files. Returned: {list(by_file.keys())}")

    avg1 = float(by_file[reads_1]['avg_len'])
    max1 = int(by_file[reads_1]['max_len'].replace(',', ''))
    avg2 = float(by_file[reads_2]['avg_len'])
    max2 = int(by_file[reads_2]['max_len'].replace(',', ''))
    return avg1, max1, avg2, max2


def run_kraken(reads_1, reads_2, threads, output_path, report_path, kraken_db, extra_args=''):
    # if kraken db does not exist, raise an error
    if not os.path.exists(kraken_db):
        raise Exception(f'Kraken database does not exist in {kraken_db}')

    pu.stream_tool('Kraken2', KRAKEN_COMMAND.format(kraken_db=kraken_db, extra_args=extra_args, reads_1=reads_1,
                                                    reads_2=reads_2, threads=threads, kraken_output=output_path,
                                                    kraken_report=report_path))


def filter_kraken_report_by_distinct_kmer_count(kraken_report_path, filtered_kraken_report_path,
                                               metadata_path, max_refs_per_species,
                                               threshold=DISTINCT_KMER_RATIO_THRESHOLD):
    """Drop low-confidence species rows (low distinct_kmer_count / genome_length ratio) from a Kraken2 report.

    Species with no estimable genome length (e.g. missing from the reference metadata) are dropped,
    since the ratio can't be computed. Non-species rows are kept untouched, since Bracken needs them
    for its tree walk.
    """
    report = pd.read_csv(kraken_report_path, sep='\t', header=None, names=KRAKEN_REPORT_COLS)
    species_mask = report['rank'] == 'S'
    species_names = report.loc[species_mask, 'name'].str.strip().tolist()

    metadata = pd.read_csv(metadata_path, sep='\t')
    genome_length_by_species = get_species_median_genome_length_by_quality(metadata, species_names, max_refs_per_species)

    genome_lengths = report['name'].str.strip().map(genome_length_by_species)
    kmer_ratio = report['distinct_kmer_count'] / genome_lengths
    keep_species_mask = species_mask & (kmer_ratio > threshold)
    filtered_report = report[~species_mask | keep_species_mask]
    filtered_report.to_csv(filtered_kraken_report_path, sep='\t', header=False, index=False, float_format='%.2f')


def get_kmer_length_options(kraken_db):
    """The read lengths a Kraken2 database has Bracken distributions for.

    Keyed on the .kmer_distrib files, which are what Bracken reads at run time. The
    database{N}mers.kraken files are bracken-build intermediates a database need not ship - the
    prebuilt GTDB one does not, and matching those would leave run_bracken with nothing to pick from.
    """
    read_lengths = []
    for filename in os.listdir(kraken_db):
        match = re.fullmatch(r'database(\d+)mers\.kmer_distrib', filename)
        if match:
            read_lengths.append(int(match.group(1)))
    return read_lengths


def get_min_reads_for_bracken(metadata_path: str, species_coverage_threshold: float, avg_sum: float, bracken_relaxation_factor:float =BRACKEN_MIN_READS_RELAXATION_FACTOR) -> int:
    """Minimum reads Bracken requires before re-estimating a taxon's abundance.

    Derived from `species_coverage_threshold` (the same absolute-coverage threshold used in
    `get_species_passing_coverage_threshold`) using the smallest genome length in the reference
    metadata, so this pre-filter never excludes a species that could pass the coverage threshold
    downstream.
    The relaxation factor allows for an even more permissive threshold, so that in case that many reads are mapped only to the genus level, we still give the species some chance to be considered by bracken and potentially by GInGeR.
    """
    metadata = pd.read_csv(metadata_path, sep='\t')
    lengths = pd.to_numeric(metadata['Length'], errors='coerce')
    lengths = lengths[lengths > 0]
    if len(lengths) == 0 or avg_sum <= 0:
        return 0
    min_genome_length = lengths.min()
    return max(0, int(np.ceil(species_coverage_threshold * min_genome_length / avg_sum) * bracken_relaxation_factor)) 


def run_bracken(kraken_report, bracken_output, bracken_report, kraken_db, min_reads_for_bracken, max_read_len: int):
    kmer_length_options = get_kmer_length_options(kraken_db)
    # get the kmer length that is closest to the read length
    read_len = min(kmer_length_options, key=lambda x: abs(int(x) - max_read_len))
    pu.stream_tool('Bracken', BRACKEN_COMMAND.format(kraken_db=kraken_db, kraken_report=kraken_report,
                                                     bracken_output=bracken_output, bracken_report=bracken_report,
                                                     read_len=read_len,
                                                     min_reads_for_bracken=min_reads_for_bracken))


def get_list_of_top_species_by_bracken(bracken_output_path, fraction_of_reads):
    bracken_out = pd.read_csv(bracken_output_path, sep='\t')
    top_species = bracken_out[bracken_out['fraction_total_reads'] > fraction_of_reads]['name'].tolist()
    log.info(f'Top species detected: {top_species}')
    return top_species


def compute_quality(metadata: pd.DataFrame) -> pd.Series:
    """A reference genome's quality score: Completeness - 5 * Contamination + ln(N50).

    Both the references that get downloaded and the genome lengths the coverage estimate is based on
    are the top scorers by this, so they have to score them the same way.
    """
    return metadata['Completeness'] - 5 * metadata['Contamination'] + \
        metadata['N50'].apply(lambda n50: 0 if n50 <= 0 else np.log(n50))


def get_species_median_genome_length_by_quality(metadata: pd.DataFrame, species_list, max_refs_per_species: int):
    """Compute median genome length per species based on the top references by Quality.

    For each species, take the top `max_refs_per_species` references by Quality (ties broken by Genome)
    and return the median of the `Length` column for these selected references.
    """
    if metadata is None or len(metadata) == 0:
        return {}

    for col in ['Genome', 'Completeness', 'Contamination', 'N50', 'species', 'Length']:
        if col not in metadata.columns:
            raise ValueError(f"Metadata is missing required column '{col}'")

    df = metadata[metadata['species'].isin(species_list)].copy()
    # UHGG lists references it has no download URL for; they cannot be used, so they should not
    # steer the median either. GTDB's table has no such column - every row is fetchable by accession
    if 'FTP_download' in df.columns:
        df = df[df['FTP_download'].astype(str).str.startswith('ftp')]
    if len(df) == 0:
        return {}

    df['Completeness'] = pd.to_numeric(df['Completeness'], errors='coerce')
    df['Contamination'] = pd.to_numeric(df['Contamination'], errors='coerce')
    df['N50'] = pd.to_numeric(df['N50'], errors='coerce')
    df['Length'] = pd.to_numeric(df['Length'], errors='coerce')

    df = df.dropna(subset=['Completeness', 'Contamination', 'N50', 'Length', 'Genome', 'species'])
    df = df[df['Length'] > 0]
    if len(df) == 0:
        return {}

    df['Quality'] = compute_quality(df)

    # Select the top references per species by Quality, breaking ties by Genome
    df = df.sort_values(['species', 'Quality', 'Genome'])
    top_refs = df.groupby('species', group_keys=False).tail(max_refs_per_species)

    medians = top_refs.groupby('species')['Length'].median().to_dict()
    return medians


def get_species_coverage_stats(bracken_output_path: str,
                              avg_sum: float,
                              metadata_path: str,
                              max_refs_per_species: int) -> pd.DataFrame:
    """Bracken output enriched with per-species estimated genome length and estimated coverage.

    estimated_coverage = new_est_reads * (avg_len_r1 + avg_len_r2) / estimated_genome_length
    """
    bracken_out = pd.read_csv(bracken_output_path, sep='\t')
    if 'name' not in bracken_out.columns or 'new_est_reads' not in bracken_out.columns:
        raise ValueError('Bracken output must include columns: name, new_est_reads')

    bracken_out['name'] = bracken_out['name'].astype(str)
    bracken_out['new_est_reads'] = pd.to_numeric(bracken_out['new_est_reads'], errors='coerce').fillna(0)
    detected_species = bracken_out['name'].dropna().tolist()

    metadata = pd.read_csv(metadata_path, sep='\t')
    median_len_by_species = get_species_median_genome_length_by_quality(metadata, detected_species, max_refs_per_species)

    bracken_out['estimated_genome_length'] = bracken_out['name'].map(median_len_by_species)
    bracken_out['estimated_coverage'] = (bracken_out['new_est_reads'] * avg_sum) / bracken_out['estimated_genome_length']
    return bracken_out


def get_species_passing_coverage_threshold(stats: pd.DataFrame, species_coverage_threshold: float):
    """The species of a `get_species_coverage_stats` table whose estimated coverage clears the
    threshold - the ones GInGeR will download references for and look for genes in."""
    passing_mask = stats['estimated_genome_length'].notna() & (stats['estimated_genome_length'] > 0) & \
        (stats['estimated_coverage'] > species_coverage_threshold)
    passing = stats.loc[passing_mask, 'name'].tolist()

    log.info(f"Detected {len(stats)} species by Bracken, keeping {len(passing)} with estimated_coverage > {species_coverage_threshold}")
    return passing


def get_distinct_minimizers_by_species(kraken_report_path: str) -> dict:
    """Map species name -> distinct minimizer count, from a Kraken2 report's species-rank rows."""
    kraken_report = pd.read_csv(kraken_report_path, sep='\t', header=None, names=KRAKEN_REPORT_COLS)
    species_rows = kraken_report[kraken_report['rank'] == 'S'].copy()
    species_rows['name'] = species_rows['name'].str.strip()
    return species_rows.set_index('name')['distinct_kmer_count'].to_dict()


def get_species_included_in_analysis_df(stats: pd.DataFrame, kraken_report_path: str, top_species) -> pd.DataFrame:
    """A `get_species_coverage_stats` table restricted to the species included in the analysis, with
    each one's distinct minimizer count added."""
    distinct_minimizers_by_species = get_distinct_minimizers_by_species(kraken_report_path)

    included = stats[stats['name'].isin(top_species)].copy()
    included['distinct_minimizers'] = included['name'].map(distinct_minimizers_by_species)
    return included


def reference_fasta_path(references_folder: str, genome: str) -> str:
    """Where a downloaded reference genome lives - one file per assembly accession, so that
    "do we have this genome already?" is an existence check."""
    return os.path.join(references_folder, f'{genome}.fna')


def run_datasets_download(accessions, zip_path: str):
    """Fetch a chunk of assembly accessions with NCBI's datasets CLI, retrying the whole chunk.

    Retried per chunk rather than per genome, because `datasets` fetches a chunk in one request.
    """
    with tempfile.NamedTemporaryFile('w', suffix='.txt') as accessions_file:
        accessions_file.write('\n'.join(accessions) + '\n')
        accessions_file.flush()
        command = DATASETS_COMMAND.format(accessions_file=accessions_file.name, zip_path=zip_path)
        for attempt in range(N_ATTEMPTS):
            out = run(command, shell=True, capture_output=True, text=True)
            # the zip check catches a datasets that reports success without writing anything, which
            # would otherwise surface as a confusing failure to unpack the archive
            if out.returncode == 0 and os.path.exists(zip_path):
                return
            error = out.stderr.strip() or out.stdout.strip()
            # the last attempt raises rather than sleeping through a retry it will not make
            if attempt == N_ATTEMPTS - 1:
                raise RuntimeError(f'datasets download failed for {len(accessions)} accessions in '
                                   f'{N_ATTEMPTS} attempts. stderr: {error}')
            log.error(f'datasets download failed: {error}. Retrying in {SLEEP_SECS} seconds')
            time.sleep(SLEEP_SECS)


def extract_genomes_from_datasets_zip(zip_path: str, references_folder: str, sample_tag: str) -> set:
    """Unpack one datasets archive into {accession}.fna files, returning the accessions found.

    A datasets archive lays genomes out as ncbi_dataset/data/{accession}/{something}.fna, with an
    assembly's sequence possibly split over several files. Each accession's files are concatenated
    into a `.part` file tagged with sample_tag and moved onto the real name so an interrupted
    extraction cannot leave a partial file that the next run treats as already downloaded, and two
    samples fetching the same accession into a shared references_folder don't interleave their
    writes into the same `.part` file (references_folder is deliberately shared across concurrent
    runs of the same cohort, since they tend to need the same species).
    """
    downloaded = set()
    with zipfile.ZipFile(zip_path) as archive:
        by_accession = {}
        for name in archive.namelist():
            parts = name.split('/')
            # ncbi_dataset/data/<accession>/<file>.fna - anything else is the archive's own
            # metadata (dataset_catalog.json, README.md, ...)
            if len(parts) >= 4 and parts[0] == 'ncbi_dataset' and parts[1] == 'data' \
                    and name.endswith('.fna'):
                by_accession.setdefault(parts[2], []).append(name)

        for accession, members in sorted(by_accession.items()):
            target = reference_fasta_path(references_folder, accession)
            part_path = f'{target}.{sample_tag}.part'
            with open(part_path, 'wb') as out_f:
                for member in sorted(members):
                    with archive.open(member) as in_f:
                        shutil.copyfileobj(in_f, out_f)
            os.replace(part_path, target)
            downloaded.add(accession)
    return downloaded


def download_missing_references(genomes, references_folder: str, sample_tag: str):
    """Make sure every genome in `genomes` has a .fna in references_folder, fetching what is absent.

    Returns the genomes that are available afterwards. NCBI suppresses assemblies over time, and
    GTDB's metadata outlives those removals, so an accession that cannot be fetched is dropped with
    a warning rather than failing the run - the references that did make it are recorded in
    references_used.csv.
    """
    missing = {genome for genome in genomes
               if not os.path.exists(reference_fasta_path(references_folder, genome))}
    available = set(genomes) - missing
    log.info(f'{len(available)} reference genomes already downloaded, fetching {len(missing)}')

    missing = sorted(missing)
    for chunk_start in range(0, len(missing), DOWNLOAD_CHUNK_SIZE):
        chunk = missing[chunk_start:chunk_start + DOWNLOAD_CHUNK_SIZE]
        # tagged with sample_tag so two samples sharing references_folder (deliberate, for cohorts
        # with overlapping species) don't collide on the same zip name - one run's cleanup
        # `os.remove` would otherwise delete the other run's in-flight download
        zip_path = os.path.join(references_folder, f'datasets_chunk_{chunk_start}_{sample_tag}.zip')
        try:
            run_datasets_download(chunk, zip_path)
            available |= extract_genomes_from_datasets_zip(zip_path, references_folder, sample_tag)
        finally:
            if os.path.exists(zip_path):
                os.remove(zip_path)
        log.info(f'downloaded {min(chunk_start + DOWNLOAD_CHUNK_SIZE, len(missing))}/{len(missing)}')

    unavailable = set(missing) - available
    if unavailable:
        log.warning(f'{len(unavailable)} reference genomes could not be downloaded from NCBI and '
                    f'are excluded from the analysis: {sorted(unavailable)}')
    return available


def gffgz_to_fasta(gff_gz_path, fasta_f):
    """Write the FASTA half of a UHGG .gff.gz - everything after its ##FASTA marker - to fasta_f."""
    with gzip.open(gff_gz_path, 'rt') as gzip_fin:
        fasta_part = False
        for line in gzip_fin:
            if fasta_part:
                fasta_f.write(line)
            elif line.startswith('##FASTA'):
                fasta_part = True


def download_uhgg_references(selected_df, references_folder: str, sample_tag: str):
    """The UHGG counterpart of download_missing_references: fetch what is missing, leave {Genome}.fna.

    UHGG serves one .gff.gz per genome over FTP rather than assemblies by accession, so this needs
    the FTP_download column rather than just the ids. The .gff.gz is kept as the cache, since shared
    reference directories are already full of them, and converted to the same {Genome}.fna layout the
    NCBI path produces - which is what lets everything downstream stay catalog-agnostic.
    """
    if 'FTP_download' not in selected_df.columns:
        raise ValueError("--reference-source uhgg needs a metadata table with an FTP_download column; "
                         "this one has none. A GTDB table carrying UHGG species names is still a GTDB "
                         "table - run it with --reference-source gtdb.")

    available, missing = set(), []
    for genome, url in zip(selected_df['Genome'], selected_df['FTP_download']):
        if os.path.exists(reference_fasta_path(references_folder, genome)):
            available.add(genome)
        else:
            missing.append((genome, url))
    log.info(f'{len(available)} reference genomes already downloaded, fetching {len(missing)}')

    unavailable = []
    for n, (genome, url) in enumerate(missing, start=1):
        gff_gz_path = os.path.join(references_folder, os.path.basename(str(url)))
        try:
            if not os.path.exists(gff_gz_path):
                download_with_retries(str(url), gff_gz_path, sample_tag)
            target = reference_fasta_path(references_folder, genome)
            part_path = f'{target}.{sample_tag}.part'
            with open(part_path, 'w') as fasta_f:
                gffgz_to_fasta(gff_gz_path, fasta_f)
            os.replace(part_path, target)
            available.add(genome)
        except Exception as e:
            # UHGG's FTP outlives individual files going missing; one unreachable genome should not
            # end a multi-hour run, and references_used.csv records what actually went in
            log.error(f'Failed to fetch {genome} from {url}: {e}')
            unavailable.append(genome)
        if n % 100 == 0:
            log.info(f'downloaded {n}/{len(missing)}')

    if unavailable:
        log.warning(f'{len(unavailable)} reference genomes could not be downloaded from UHGG and are '
                    f'excluded from the analysis: {sorted(unavailable)}')
    return available


def download_with_retries(url: str, target_path: str, sample_tag: str):
    """Fetch one URL to target_path, retrying N_ATTEMPTS times.

    Written through a sample-tagged .part file and moved onto the real name, so an interrupted
    download leaves nothing rather than a truncated file the next run treats as cached.
    """
    part_path = f'{target_path}.{sample_tag}.part'
    for attempt in range(N_ATTEMPTS):
        try:
            data = urllib.request.urlopen(url, timeout=URLOPEN_TIMEOUT).read()
            with open(part_path, 'wb') as f:
                f.write(data)
            os.replace(part_path, target_path)
            return
        except Exception as e:
            # the last attempt raises rather than sleeping through a retry it will not make
            if attempt == N_ATTEMPTS - 1:
                raise
            log.error(f'Failed to download {url}: {e}. Retrying in {SLEEP_SECS} seconds')
            time.sleep(SLEEP_SECS)


def write_genome_to_merged_fasta(genome: str, references_folder: str, merged_filtered_fasta_f,
                                 contig_to_genome_f):
    """Append one reference genome to the merged fasta, recording which genome each contig came from.

    The map is what attributes a match on a reference contig to a genome and from there to a species,
    which the contig's name cannot do: NCBI contigs are nucleotide accessions (NZ_CP007265.1).
    """
    with open(reference_fasta_path(references_folder, genome)) as genome_f:
        for line in genome_f:
            if line.startswith('>'):
                # minimap2 reports the first whitespace-delimited token as the target name, so that
                # is the key the PAF will have to be looked up by
                contig_to_genome_f.write(f'{line[1:].split()[0]}\t{genome}\n')
            merged_filtered_fasta_f.write(line)


def generate_filtered_minimap_db_according_to_selected_species(top_species, metadata_path, references_folder,
                                                               merged_filtered_fasta, max_refs_per_species,
                                                               contig_to_genome_path, sample_tag,
                                                               source=DEFAULT_REFERENCE_SOURCE):
    """Build the sample-specific reference database out of the top references of every selected species.

    In two phases - pick the references and fetch what is missing in batches, then stream what is on
    disk into the merged fasta - because `datasets` is fed many accessions per call.
    """
    metadata = pd.read_csv(metadata_path, sep='\t')
    metadata['Quality'] = compute_quality(metadata)
    # max_refs_per_species applies per subspecies when the metadata names them, and to the species as
    # a whole when it does not
    has_subspecies = 'subspecies' in metadata.columns
    selected_samples_dfs_list = []
    has_ftp = 'FTP_download' in metadata.columns
    for species in top_species:
        single_species_table = metadata[metadata.species == species]
        if has_ftp:
            single_species_table = single_species_table[
                single_species_table.FTP_download.astype(str).str.startswith('ftp')]
        tables_to_take_from = ([table for _, table in single_species_table.groupby('subspecies')]
                               if has_subspecies else [single_species_table])
        for table in tables_to_take_from:
            selected_samples_dfs_list.append(take_top_references_per_species(max_refs_per_species, table))

    if not selected_samples_dfs_list:
        raise ValueError(f'None of the {len(top_species)} selected species has references in '
                         f'{metadata_path}')
    # a genome listed twice would go into the merged fasta twice, leaving minimap2 with duplicate
    # sequence names
    selected_samples_df = pd.concat(selected_samples_dfs_list).drop_duplicates(subset='Genome')

    # the one place the catalogs genuinely differ
    if source == 'uhgg':
        available = download_uhgg_references(selected_samples_df, references_folder, sample_tag)
    else:
        available = download_missing_references(selected_samples_df['Genome'].tolist(), references_folder, sample_tag)
    selected_samples_df = selected_samples_df[selected_samples_df['Genome'].isin(available)]

    with open(merged_filtered_fasta, 'w') as merged_filtered_fasta_f, \
            open(contig_to_genome_path, 'w') as contig_to_genome_f:
        contig_to_genome_f.write('contig\tGenome\n')
        for genome in selected_samples_df['Genome']:
            write_genome_to_merged_fasta(genome, references_folder, merged_filtered_fasta_f,
                                         contig_to_genome_f)
    log.info(f'built a reference database of {len(selected_samples_df)} genomes for '
             f'{len(top_species)} species')
    return selected_samples_df


def take_top_references_per_species(max_refs_per_species, single_species_table):
    # Take top X references according to Quality score, breaking ties alphabetically by Genome
    return single_species_table.sort_values(['Quality', 'Genome']).tail(max_refs_per_species)


@pu.step_timing
def get_filtered_references_database(reads_1, reads_2, threads, kraken_output_path, kraken_report_path,
                                     bracken_output,
                                     bracken_report, species_coverage_threshold, metadata_path, references_folder,
                                     merged_filtered_fasta, references_used_path, max_species_representatives, kraken_db,
                                     species_included_in_analysis_path, contig_to_genome_path,
                                     reuse_existing_kraken_output=False, source=DEFAULT_REFERENCE_SOURCE,
                                     distinct_kmer_ratio_threshold=DISTINCT_KMER_RATIO_THRESHOLD):
    pu.check_and_makedir(kraken_output_path)
    pu.check_and_make_dir_no_file_name(references_folder)
    # Kraken2 is by far the slowest step - when asked to, reuse a previous run's output instead of redoing it
    if not (reuse_existing_kraken_output and os.path.exists(kraken_output_path) and os.path.exists(kraken_report_path)):
        run_kraken(reads_1, reads_2, threads, kraken_output_path, kraken_report_path, kraken_db,
                   extra_args=REFERENCE_SOURCES[source]['kraken_extra_args'])
    filtered_kraken_report_path = f'{kraken_report_path}.distinct_kmer_filtered'
    filter_kraken_report_by_distinct_kmer_count(kraken_report_path, filtered_kraken_report_path,
                                                metadata_path, max_species_representatives,
                                                threshold=distinct_kmer_ratio_threshold)

    avg1, max1, avg2, max2 = get_paired_reads_seqkit_stats(reads_1, reads_2)
    avg_sum = avg1 + avg2
    max_read_len = max(max1, max2)

    min_reads_for_bracken = get_min_reads_for_bracken(metadata_path, species_coverage_threshold, avg_sum)
    run_bracken(filtered_kraken_report_path, bracken_output, bracken_report, kraken_db, min_reads_for_bracken, max_read_len)
    # both of the steps below read the same coverage table, which costs a pass over the reference
    # metadata to estimate every detected species' genome length - so it is built once here
    coverage_stats = get_species_coverage_stats(bracken_output, avg_sum, metadata_path, max_species_representatives)
    top_species = get_species_passing_coverage_threshold(coverage_stats, species_coverage_threshold)
    species_included_in_analysis_df = get_species_included_in_analysis_df(coverage_stats, filtered_kraken_report_path,
                                                                          top_species)
    species_included_in_analysis_df.to_csv(species_included_in_analysis_path, index=False, sep='\t')
    # identifies this sample's in-flight downloads within a references_folder shared across a
    # cohort's concurrent runs, so they don't collide on the same temp file
    sample_tag = os.path.basename(reads_1)
    selected_species_df = generate_filtered_minimap_db_according_to_selected_species(top_species, metadata_path,
                                                                                     references_folder,
                                                                                     merged_filtered_fasta,
                                                                                     max_refs_per_species=max_species_representatives,
                                                                                     contig_to_genome_path=contig_to_genome_path,
                                                                                     sample_tag=sample_tag,
                                                                                     source=source)
    selected_species_df.to_csv(references_used_path, index=False, sep='\t')
    return merged_filtered_fasta
