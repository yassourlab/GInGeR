from subprocess import Popen, PIPE
import logging
import os
from tqdm import tqdm
from ginger import pipeline_utils as pu
from ginger import constants as c

log = logging.getLogger(__name__)
META_SPADES_COMMAND = "spades.py --meta -1 {short_reads_1} -2 {short_reads_2} {optional_long_reads}-o {output_folder} -t {threads}"
OPTIONAL_LONG_READS_ADDITION = "--nanopore {long_reads} "
SECS_TO_MIN = 60
TQDM_INTERVAL_MINS = 10


def validate_spades_output(output_folder):
    """
    Validates that SPAdes generated the required output files.
    
    :param output_folder: The SPAdes output directory
    :raises RuntimeError: If required files are missing
    """
    contigs_path = c.CONTIGS_PATH_TEMPLATE.format(assembly_dir=output_folder)
    assembly_graph_path = c.ASSEMBLY_GRAPH_PATH_TEMPLATE.format(assembly_dir=output_folder)
    spades_log_path = os.path.join(output_folder, 'spades.log')
    
    missing_files = []
    if not os.path.exists(contigs_path):
        missing_files.append('contigs.fasta')
    if not os.path.exists(assembly_graph_path):
        missing_files.append('assembly_graph.fastg')
    
    if missing_files:
        error_msg = f"SPAdes assembly failed. Missing required output files: {', '.join(missing_files)}. "
        if os.path.exists(spades_log_path):
            error_msg += f"Please check the SPAdes log file for details: {spades_log_path}"
        else:
            error_msg += f"SPAdes log file not found at: {spades_log_path}"
        log.error(error_msg)
        raise RuntimeError(error_msg)
    

@pu.step_timing
def run_meta_or_hybrid_spades(short_reads_1, short_reads_2, long_reads, output_folder, threads):
    long_reads_str = ''
    if long_reads:
        long_reads_str = OPTIONAL_LONG_READS_ADDITION.format(long_reads=long_reads)
    command = META_SPADES_COMMAND.format(short_reads_1=short_reads_1, short_reads_2=short_reads_2,
                                         optional_long_reads=long_reads_str, output_folder=output_folder,
                                         threads=threads)
    log.info(f'running MetaSPAdes - {command}')
    with Popen(command.split(' '), stdout=PIPE) as p:
        tqdm_boject = tqdm(iter(lambda: p.stdout.readline(), b""), mininterval=TQDM_INTERVAL_MINS * SECS_TO_MIN)
        output_lines = [output_line for output_line in tqdm_boject]
    
    validate_spades_output(output_folder)
    return output_folder
