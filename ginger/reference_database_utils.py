import urllib.request
import time
from subprocess import run, Popen, PIPE
from collections import defaultdict

import numpy as np
import pandas as pd
from glob import glob
import urllib
import gzip
import logging
import csv
from ginger import pipeline_utils as pu
from tqdm import tqdm
import os
import re

log = logging.getLogger(__name__)
KRAKEN_COMMAND = 'kraken2 --db {kraken_db} --paired {reads_1} {reads_2} --threads {threads} --output {kraken_output} --report {kraken_report} --confidence 0.1 --use-names'  # --report {report}
BRACKEN_COMMAND = 'bracken -d {kraken_db} -i {kraken_report} -o {bracken_output} -w {bracken_report} -r {read_len} -l S -t {threads}'
URLOPEN_TIMEOUT = 60
N_ATTEMPTS = 10
SLEEP_SECS = 60


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


def run_kraken(reads_1, reads_2, threads, output_path, report_path, kraken_db):
    command = KRAKEN_COMMAND.format(kraken_db=kraken_db, reads_1=reads_1, reads_2=reads_2,
                                    threads=threads, kraken_output=output_path, kraken_report=report_path)
    # if kraken db does not exist, raise an error
    if not os.path.exists(kraken_db):
        raise Exception(f'Kraken database does not exist in {kraken_db}')

    log.info(f'Running Kraken2 - {command}')
    # command_output = run(command, shell=True, capture_output=True)
    with Popen(command.split(' '), stdout=PIPE) as kraken_process:
        output_lines = [output_line for output_line in tqdm(iter(lambda: kraken_process.stdout.readline(), b""))]
        if kraken_process.returncode:
            log.error(kraken_process.stderr)
            raise Exception('Kraken2 failed - GInGeR aborted')


def get_kmer_length_options(kraken_db):
    pattern = r'database(.*?)mers\.kraken'
    # List to store the extracted *** parts
    extracted_parts = []
    # Walk through the directory
    for filename in os.listdir(kraken_db):
        # Check if the filename matches the pattern
        match = re.match(pattern, filename)
        if match:
            # Extract the *** part and add it to the list
            extracted_parts.append(int(match.group(1)))
    return extracted_parts


def run_bracken(kraken_report, bracken_output, bracken_report, kraken_db, threads, max_read_len: int):
    kmer_length_options = get_kmer_length_options(kraken_db)
    # get the kmer length that is closest to the read length
    read_len = min(kmer_length_options, key=lambda x: abs(int(x) - max_read_len))
    command = BRACKEN_COMMAND.format(kraken_db=kraken_db, kraken_report=kraken_report, bracken_output=bracken_output,
                                     bracken_report=bracken_report, read_len=read_len,
                                     threads=threads)
    log.info(f'Running Bracken - {command}')
    # command_output = run(command, shell=True, capture_output=True)
    with Popen(command.split(' '), stdout=PIPE) as bracken_process:
        output_lines = [output_line for output_line in tqdm(iter(lambda: bracken_process.stdout.readline(), b""))]
        if bracken_process.returncode:
            log.error(bracken_process.stderr)
            raise Exception('Bracken failed - GInGeR aborted')
            log.info(bracken_process.stdout)


def get_list_of_top_species_by_bracken(bracken_output_path, fraction_of_reads):
    bracken_out = pd.read_csv(bracken_output_path, sep='\t')
    top_species = bracken_out[bracken_out['fraction_total_reads'] > fraction_of_reads]['name'].tolist()
    log.info(f'Top species detected: {top_species}')
    return top_species


def get_species_median_genome_length_by_quality(metadata: pd.DataFrame, species_list, max_refs_per_species: int):
    """Compute median genome length per species based on the top references by Quality.

    For each species, take the top `max_refs_per_species` references by Quality (ties broken by Genome)
    and return the median of the `Length` column for these selected references.
    """
    if metadata is None or len(metadata) == 0:
        return {}

    for col in ['Genome', 'Completeness', 'Contamination', 'N50', 'FTP_download', 'species', 'Length']:
        if col not in metadata.columns:
            raise ValueError(f"Metadata is missing required column '{col}'")

    df = metadata[metadata['species'].isin(species_list)].copy()
    df = df[df['FTP_download'].astype(str).str.startswith('ftp')]
    if len(df) == 0:
        return {}

    df['Completeness'] = pd.to_numeric(df['Completeness'], errors='coerce')
    df['Contamination'] = pd.to_numeric(df['Contamination'], errors='coerce')
    df['N50'] = pd.to_numeric(df['N50'], errors='coerce')
    df['Length'] = pd.to_numeric(df['Length'], errors='coerce')

    df = df.dropna(subset=['Completeness', 'Contamination', 'N50', 'Length', 'Genome', 'species', 'FTP_download'])
    df = df[df['Length'] > 0]
    if len(df) == 0:
        return {}

    df['Quality'] = df['Completeness'] - 5 * df['Contamination'] + df['N50'].apply(lambda x: 0 if x <= 0 else np.log(x))

    # Select the top references per species by Quality, breaking ties by Genome
    df = df.sort_values(['species', 'Quality', 'Genome'])
    top_refs = df.groupby('species', group_keys=False).tail(max_refs_per_species)

    medians = top_refs.groupby('species')['Length'].median().to_dict()
    return medians


def get_species_passing_coverage_threshold(bracken_output_path: str,
                                          avg_sum: float,
                                          metadata_path: str,
                                          max_refs_per_species: int,
                                          species_coverage_threshold: float):
    """Return species list filtered by estimated coverage.

    estimated_coverage = new_est_reads * (avg_len_r1 + avg_len_r2) / median_genome_length
    """
    bracken_out = pd.read_csv(bracken_output_path, sep='\t')
    if 'name' not in bracken_out.columns or 'new_est_reads' not in bracken_out.columns:
        raise ValueError('Bracken output must include columns: name, new_est_reads')

    bracken_out['new_est_reads'] = pd.to_numeric(bracken_out['new_est_reads'], errors='coerce').fillna(0)
    detected_species = bracken_out['name'].dropna().astype(str).tolist()

    metadata = pd.read_csv(metadata_path, sep='\t')
    median_len_by_species = get_species_median_genome_length_by_quality(metadata, detected_species, max_refs_per_species)

    passing = []
    for _, row in bracken_out.iterrows():
        species = str(row['name'])
        reads_est = float(row['new_est_reads'])
        med_len = median_len_by_species.get(species)
        if not med_len or med_len <= 0:
            continue
        coverage = (reads_est * avg_sum) / float(med_len)
        if coverage > species_coverage_threshold:
            passing.append(species)

    log.info(f"Detected {len(detected_species)} species by Bracken, keeping {len(passing)} with estimated_coverage > {species_coverage_threshold}")
    return passing


def download_and_write_content_to_file(references_folder, references_folder_content, ftp_download_str: str,
                                       merged_filtered_fasta_f):
    mgyg_file = ftp_download_str.split('/')[-1]
    local_tar_gz_path = f'{references_folder}/{mgyg_file}'
    # download file from FTP if needed

    if mgyg_file in references_folder_content:
        log.debug(f'{mgyg_file} already exists in {references_folder} - skipping download')
    else:
        log.debug(f'{mgyg_file} not found in {references_folder} - downloading file')
        for attempt in range(N_ATTEMPTS):
            try:
                data = urllib.request.urlopen(ftp_download_str, timeout=URLOPEN_TIMEOUT).read()
                with open(local_tar_gz_path, 'wb') as f:
                    f.write(data)
                break
            except Exception as e:
                log.error(f'Failed to download {ftp_download_str}: {e}. Retrying in {SLEEP_SECS} seconds')
                time.sleep(SLEEP_SECS)
                if attempt == N_ATTEMPTS - 1:
                    raise e

    # read tar.gt file and add it's content to the merged filtered fasta
    gffgz_to_fasta(local_tar_gz_path, merged_filtered_fasta_f)


# TODO sed -n '/>MGYG000005036.fa_1/,$p' MGYG000005036.gff > MGYG000005036.fasta works in command line. Can I use it here?
def gffgz_to_fasta(local_tar_gz_path, merged_filtered_fasta_f):
    with gzip.open(local_tar_gz_path, 'rt') as gzip_fin:
        fasta_part = False
        for line in gzip_fin:
            if fasta_part:
                merged_filtered_fasta_f.write(line)
            elif line.startswith('##FASTA'):
                fasta_part = True


def generate_filtered_minimap_db_according_to_selected_species(top_species, metadata_path, references_folder,
                                                               merged_filtered_fasta, max_refs_per_species):
    metadata = pd.read_csv(metadata_path, sep='\t')
    # Create column 'Quality' as Completeness - 5 * Contamination + ln(N50)
    metadata['Quality'] = metadata.Completeness - 5 * metadata.Contamination + \
                                      metadata.N50.apply(lambda x: 0 if x <= 0 else np.log(x))
    # metadata['species'] = metadata.Lineage.apply(lambda x: x.split('s__')[-1])
    references_folder_content = [x.split('/')[-1] for x in glob(references_folder + '/*')]
    selected_samples_dfs_list = []
    with open(merged_filtered_fasta, 'w') as merged_filtered_fasta_f:
        for species in top_species:
            single_species_table = metadata[
                (metadata.species == species) & (metadata.FTP_download.str.startswith('ftp'))]
            if 'subspecies' in single_species_table.columns:
                for subspecies, subspecies_table in single_species_table.groupby('subspecies'):
                    top_x_df = take_top_species_and_download_to_file(max_refs_per_species, subspecies_table,
                                                                     references_folder, references_folder_content,
                                                                     merged_filtered_fasta_f)
                    selected_samples_dfs_list.append(top_x_df)
            else:
                top_x_df = take_top_species_and_download_to_file(max_refs_per_species, single_species_table,
                                                                 references_folder,
                                                                 references_folder_content, merged_filtered_fasta_f)
                selected_samples_dfs_list.append(top_x_df)

    return pd.concat(selected_samples_dfs_list)


def take_top_species_and_download_to_file(max_refs_per_species, single_species_table, references_folder,
                                          references_folder_content, merged_filtered_fasta_f):
    # Take top X references according to Quality score, breaking ties alphabetically by Genome
    top_x_df = single_species_table.sort_values(['Quality', 'Genome']).tail(max_refs_per_species)
    top_x_df.apply(lambda x: download_and_write_content_to_file(references_folder,
                                                                references_folder_content,
                                                                x.FTP_download,
                                                                merged_filtered_fasta_f), axis=1)
    return top_x_df


@pu.step_timing
def get_filtered_references_database(reads_1, reads_2, threads, kraken_output_path, kraken_report_path,
                                     bracken_output,
                                     bracken_report, species_coverage_threshold, metadata_path, references_folder,
                                     merged_filtered_fasta, references_used_path, max_species_representatives, kraken_db):
    pu.check_and_makedir(kraken_output_path)
    pu.check_and_make_dir_no_file_name(references_folder)
    run_kraken(reads_1, reads_2, threads, kraken_output_path, kraken_report_path, kraken_db)

    avg1, max1, avg2, max2 = get_paired_reads_seqkit_stats(reads_1, reads_2)
    avg_sum = avg1 + avg2
    max_read_len = max(max1, max2)

    run_bracken(kraken_report_path, bracken_output, bracken_report, kraken_db, threads, max_read_len)
    top_species = get_species_passing_coverage_threshold(bracken_output, avg_sum, metadata_path,
                                                         max_species_representatives, species_coverage_threshold)
    selected_species_df = generate_filtered_minimap_db_according_to_selected_species(top_species, metadata_path,
                                                                                     references_folder,
                                                                                     merged_filtered_fasta,
                                                                                     max_refs_per_species=max_species_representatives)
    selected_species_df.to_csv(references_used_path, index=False, sep='\t')
    return merged_filtered_fasta
