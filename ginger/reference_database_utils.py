import urllib.request
from subprocess import run,Popen,PIPE
from collections import defaultdict
import pandas as pd
from glob import glob
import urllib
import gzip
import logging
from ginger import pipeline_utils as pu
from tqdm import tqdm

log = logging.getLogger(__name__)
TOOLS = '/sci/labs/morani/morani/icore-data/lab/Tools'
KRAKEN_PATH = f'kraken2'
KRAKEN_DB = f'{TOOLS}/kraken2_db/UnifiedHumanGastrointestinalGenome'
KRAKEN_COMMAND = '{kraken_path} --db {kraken_db} --paired {reads_1} {reads_2} --threads {threads} --output {output_path} --use-names'  # --report {report}
KRAKEN_OUTPUT_HEADER = ['classified', 'read', 'genome', 'reads_len', 'mapping_str']

READ_STATUS_INDEX = 0
SPECIES_NAME_INDEX = 2


# TODO do I need to constantly log minimap's output (see assembly_utils) or is it enough to just log it in the end
def run_kraken(reads_1, reads_2, threads, output_path):
    command = KRAKEN_COMMAND.format(kraken_path=KRAKEN_PATH, kraken_db=KRAKEN_DB, reads_1=reads_1, reads_2=reads_2,
                                    threads=threads, output_path=output_path)
    log.info(f'running kraken: {command}')
    # command_output = run(command, shell=True, capture_output=True)
    with Popen(command.split(' '), stdout=PIPE) as kraken_process:
        output_lines = [output_line for output_line in tqdm(iter(lambda: kraken_process.stdout.readline(), b""))]
        if kraken_process.returncode:
            log.error(kraken_process.stderr)
            raise Exception('GInGeR failed to run Kraken2. The pipeline will abort')
        else:
            log.info(kraken_process.stdout)


def get_list_of_top_species_by_kraken(kraken_output_path, reads_ratio_th):
    species_dict = defaultdict(int)
    reads_count = 0
    with open(kraken_output_path) as kraken_f:
        for line in kraken_f:
            line_split = line.split('\t')
            if line_split[READ_STATUS_INDEX] == 'C':
                reads_count += 1
                species_dict[line_split[SPECIES_NAME_INDEX]] += 1
    top_species = [species.split(' (')[0] for species, c in species_dict.items() if
                   c / reads_count > reads_ratio_th]
    log.info(top_species)
    return top_species


def download_and_write_content_to_file(references_folder, references_folder_content, ftp_download_str: str,
                                       merged_filtered_fasta_f):
    mgyg_file = ftp_download_str.split('/')[-1]
    local_tar_gz_path = f'{references_folder}/{mgyg_file}'
    # download file from FTP if needed

    if mgyg_file in references_folder_content:
        log.info(f'{mgyg_file} already in {references_folder}. File will not be downloaded')
    else:
        log.info(f'{mgyg_file} not in {references_folder}. File will be downloaded')

        data = urllib.request.urlopen(ftp_download_str).read()
        with open(local_tar_gz_path, 'wb') as f:
            f.write(data)

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
                                                               merged_filtered_fasta, max_refs_per_species=10):
    # TODO - move max refs per species to external parameters
    metadata = pd.read_csv(metadata_path, sep='\t')
    # metadata['species'] = metadata.Lineage.apply(lambda x: x.split('s__')[-1])
    references_folder_content = [x.split('/')[-1] for x in glob(references_folder + '/*')]
    with open(merged_filtered_fasta, 'w') as merged_filtered_fasta_f:
        for s in top_species:
            single_species_table = metadata[(metadata.species == s) & (metadata.FTP_download.str.startswith('ftp'))].sort_values('Completeness', ascending=False).head(
                max_refs_per_species)
            single_species_table.apply(
                lambda x: download_and_write_content_to_file(references_folder, references_folder_content,
                                                             x.FTP_download,
                                                             merged_filtered_fasta_f), axis=1)

    return merged_filtered_fasta


@pu.step_timing
def get_filtered_references_database(reads_1, reads_2, threads, kraken_output_path, reads_ratio_th, metadata_path,
                                     references_folder, merged_filtered_fasta, max_species_representatives):
    pu.check_and_makedir(kraken_output_path)
    pu.check_and_make_dir_no_file_name(references_folder)
    run_kraken(reads_1, reads_2, threads, kraken_output_path)
    top_species = get_list_of_top_species_by_kraken(kraken_output_path, reads_ratio_th)
    generate_filtered_minimap_db_according_to_selected_species(top_species, metadata_path, references_folder,
                                                               merged_filtered_fasta,
                                                               max_refs_per_species=max_species_representatives)
    return merged_filtered_fasta
