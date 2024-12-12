import urllib.request
from subprocess import run, Popen, PIPE
from collections import defaultdict
import pandas as pd
from glob import glob
import urllib
import gzip
import logging
from ginger import pipeline_utils as pu
from tqdm import tqdm
from Bio import SeqIO
import os
import re

log = logging.getLogger(__name__)
KRAKEN_DB = f'/sci/labs/morani/morani/icore-data/lab/Tools/kraken2_db/UnifiedHumanGastrointestinalGenome'
KRAKEN_COMMAND = 'kraken2 --db {kraken_db} --paired {reads_1} {reads_2} --threads {threads} --output {kraken_output} --report {kraken_report} --confidence 0.1 --use-names'  # --report {report}
BRACKEN_COMMAND = 'bracken -d {kraken_db} -i {kraken_report} -o {bracken_output} -w {bracken_report} -r {read_len} -l S -t {min_reads_for_bracken}'
KRAKEN_OUTPUT_HEADER = ['classified', 'read', 'genome', 'reads_len', 'mapping_str']

READ_STATUS_INDEX = 0
SPECIES_NAME_INDEX = 2
NUM_READS_RATIO = 0.005
URLOPEN_TIMEOUT = 60
N_ATTEMPTS = 10
SLEEP_SECS = 60


def run_kraken(reads_1, reads_2, threads, output_path, report_path):
    command = KRAKEN_COMMAND.format(kraken_db=KRAKEN_DB, reads_1=reads_1, reads_2=reads_2,
                                    threads=threads, kraken_output=output_path, kraken_report=report_path)
    # if kraken db does not exist, raise an error
    if not os.path.exists(KRAKEN_DB):
        raise Exception(f'Kraken database does not exist in {KRAKEN_DB}')

    log.info(f'running Kraken2: {command}')
    # command_output = run(command, shell=True, capture_output=True)
    with Popen(command.split(' '), stdout=PIPE) as kraken_process:
        output_lines = [output_line for output_line in tqdm(iter(lambda: kraken_process.stdout.readline(), b""))]
        if kraken_process.returncode:
            log.error(kraken_process.stderr)
            raise Exception('GInGeR failed to run Kraken2. The pipeline will abort')


def get_max_read_len(fastq_file):
    max_length = 0
    num_reads = 0
    if fastq_file.endswith(".gz"):
        fastq_file = gzip.open(fastq_file, "rt")
    for record in SeqIO.parse(fastq_file, "fastq"):
        num_reads += 1
        read_length = len(record.seq)
        if read_length > max_length:
            max_length = read_length

    # Output the maximum read length
    return num_reads, max_length


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


def run_bracken(reads_1, kraken_report, bracken_output, bracken_report):
    num_reads, max_read_len = get_max_read_len(reads_1)
    kmer_length_options = get_kmer_length_options(KRAKEN_DB)
    # get the kmer length that is closest to the read length
    read_len = min(kmer_length_options, key=lambda x: abs(int(x) - max_read_len))
    command = BRACKEN_COMMAND.format(kraken_db=KRAKEN_DB, kraken_report=kraken_report, bracken_output=bracken_output,
                                     bracken_report=bracken_report, read_len=read_len,
                                     min_reads_for_bracken=int(num_reads * NUM_READS_RATIO))
    log.info(f'running Bracken: {command}')
    # command_output = run(command, shell=True, capture_output=True)
    with Popen(command.split(' '), stdout=PIPE) as bracken_process:
        output_lines = [output_line for output_line in tqdm(iter(lambda: bracken_process.stdout.readline(), b""))]
        if bracken_process.returncode:
            log.error(bracken_process.stderr)
            raise Exception('GInGeR failed to run bracken_process. The pipeline will abort')
            log.info(bracken_process.stdout)


def get_list_of_top_species_by_bracken(bracken_output_path, fraction_of_reads):
    bracken_out = pd.read_csv(bracken_output_path, sep='\t')
    top_species = bracken_out[bracken_out['fraction_total_reads'] > fraction_of_reads]['name'].tolist()
    log.info(top_species)
    return top_species


def download_and_write_content_to_file(references_folder, references_folder_content, ftp_download_str: str,
                                       merged_filtered_fasta_f):
    mgyg_file = ftp_download_str.split('/')[-1]
    local_tar_gz_path = f'{references_folder}/{mgyg_file}'
    # download file from FTP if needed

    if mgyg_file in references_folder_content:
        log.debug(f'{mgyg_file} already in {references_folder}. File will not be downloaded')
    else:
        log.debug(f'{mgyg_file} not in {references_folder}. File will be downloaded')
        for attempt in range(N_ATTEMPTS):
            try:
                data = urllib.request.urlopen(ftp_download_str, timeout=URLOPEN_TIMEOUT).read()
                with open(local_tar_gz_path, 'wb') as f:
                    f.write(data)
                break
            except Exception as e:
                log.error(f'Failed to download {ftp_download_str} with error: {e} trying again in {SLEEP_SECS} seconds')
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
    top_x_df = single_species_table.sort_values(['Completeness', 'Genome'], ascending=[False, True]).head(
        max_refs_per_species)
    top_x_df.apply(lambda x: download_and_write_content_to_file(references_folder,
                                                                references_folder_content,
                                                                x.FTP_download,
                                                                merged_filtered_fasta_f), axis=1)
    return top_x_df


def take_top_refs_and_download_to_file():
    @pu.step_timing
    def get_filtered_references_database(reads_1, reads_2, threads, kraken_output_path, kraken_report_path,
                                         bracken_output,
                                         bracken_report, reads_ratio_th, metadata_path, references_folder,
                                         merged_filtered_fasta, references_used_path, max_species_representatives):
        pu.check_and_makedir(kraken_output_path)
        pu.check_and_make_dir_no_file_name(references_folder)
        run_kraken(reads_1, reads_2, threads, kraken_output_path, kraken_report_path)
        run_bracken(reads_1, kraken_report_path, bracken_output, bracken_report)
        top_species = get_list_of_top_species_by_bracken(bracken_output, reads_ratio_th)
        selected_species_df = generate_filtered_minimap_db_according_to_selected_species(top_species, metadata_path,
                                                                                         references_folder,
                                                                                         merged_filtered_fasta,
                                                                                         max_refs_per_species=max_species_representatives)
        selected_species_df.to_csv(references_used_path, index=False, sep='\t')
        return merged_filtered_fasta
