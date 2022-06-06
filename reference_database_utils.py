import urllib.request
from subprocess import run
import timeit
from collections import defaultdict
import pandas as pd
import sequence_alignment_utils as sau
from glob import glob
import urllib
import gzip

# /vol/sci/bio/data/moran.yassour/lab/Projects/hybrid_assembler/pipeline_tree/zymo_standard/data/short_reads/ERR2984773_Illumina_MiSeq_paired_end_sequencing_1.fastq.gz
# /vol/sci/bio/data/moran.yassour/lab/Projects/hybrid_assembler/pipeline_tree/zymo_standard/data/short_reads/ERR2984773_Illumina_MiSeq_paired_end_sequencing_2.fastq.gz
# output /vol/sci/bio/data/moran.yassour/lab/Projects/hybrid_assembler/arg-to-bug-pipeline/kraken_in_200k_genes.txt
# report /vol/sci/bio/data/moran.yassour/lab/Projects/hybrid_assembler/arg-to-bug-pipeline/kraken_in_200k_genes_report.txt


KRAKEN_PATH = '/vol/sci/bio/data/moran.yassour/lab/Tools/kraken2/kraken2'
KRAKEN_DB = '/vol/sci/bio/data/moran.yassour/lab/Tools/kraken2_db/UnifiedHumanGastrointestinalGenome'
KRAKEN_COMMAND = '{kraken_path} --db {kraken_db} --paired {reads_1} {reads_2} --threads {threads} --output {output} --use-names'  # --report {report}
KRAKEN_OUTPUT_HEADER = ['classified', 'read', 'genome', 'reads_len', 'mapping_str']

READ_STATUS_INDEX = 0
SPECIES_NAME_INDEX = 2


def run_kraken(reads_1, reads_2, threads, output):
    command = KRAKEN_COMMAND.format(kraken_path=KRAKEN_PATH, kraken_db=KRAKEN_DB, reads_1=reads_1, reads_2=reads_2,
                                    threads=threads, output=output)
    print(command)
    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        print(command_output.stderr)
    else:
        print(command_output.stdout)


def get_list_of_top_species_by_kraken(kraken_output_path, reads_ratio_th):
    start = timeit.default_timer()
    species_dict = defaultdict(int)
    reads_count = 0
    with open(kraken_output_path) as kraken_f:
        for l in kraken_f:
            line_split = l.split('\t')
            reads_count += 1
            if line_split[READ_STATUS_INDEX] == 'C':
                species_dict[line_split[SPECIES_NAME_INDEX]] += 1
    top_species = [species.split(' (')[0] for species, c in species_dict.items() if
                   c / reads_count > reads_ratio_th]
    print(top_species)
    stop = timeit.default_timer()
    print('Time get_list_of_top_species_by_kraken: ', stop - start)
    return top_species


def download_and_write_content_to_file(references_folder, references_folder_content, ftp_download_str: str,
                                       merged_filtered_fasta):
    mgyg_file = ftp_download_str.split('/')[-1]
    local_tar_gz_path = f'{references_folder}/{mgyg_file}'
    # download file from FTP if needed

    if mgyg_file in references_folder_content:
        print(f'{mgyg_file} in {references_folder}')
    else:
        print(f'{mgyg_file} not in {references_folder}')

        data = urllib.request.urlopen(ftp_download_str).read()
        with open(local_tar_gz_path, 'wb') as f:
            f.write(data)

    # read tar.gt file and add it's content to the merged filtered fasta
    with gzip.open(local_tar_gz_path, 'rt') as gzip_fin:
        fasta_part = False
        for line in gzip_fin:
            if fasta_part:
                merged_filtered_fasta.write(line)
            elif line.startswith('##FASTA'):
                fasta_part = True


def generate_filtered_minimap_db_according_to_selected_species(top_species, metadata_path, references_folder,
                                                               merged_filtered_fasta):
    metadata = pd.read_csv(metadata_path, sep='\t')
    metadata['species'] = metadata.Lineage.apply(lambda x: x.split('s__')[-1])
    references_folder_content = [x.split('/')[-1] for x in glob(references_folder + '/*')]
    with open(merged_filtered_fasta, 'w'):
        for s in top_species:
            single_species_table = metadata[metadata.species == s].sort_values('Length', ascending=False).head(100)
            single_species_table.apply(
                lambda x: download_and_write_content_to_file(references_folder, references_folder_content,
                                                             x.FTP_download,
                                                             merged_filtered_fasta), axis=1)

    return merged_filtered_fasta


def get_filtered_references_database(reads_1, reads_2, threads, kraken_output_path, reads_ratio_th, metadata_path,
                                     references_folder, merged_filtered_fasta):
    run_kraken(reads_1, reads_2, threads, kraken_output_path)
    top_species = get_list_of_top_species_by_kraken(kraken_output_path, reads_ratio_th)
    generate_filtered_minimap_db_according_to_selected_species(top_species, metadata_path, references_folder,
                                                               merged_filtered_fasta)


