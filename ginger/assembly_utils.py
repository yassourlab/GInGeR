from subprocess import run, Popen, PIPE
import logging
from tqdm import tqdm
from ginger import pipeline_utils as pu

log = logging.getLogger(__name__)
META_SPADES_COMMAND = "spades.py --meta -1 {short_reads_1} -2 {short_reads_2} {optional_long_reads}-o {output_folder} -t {threads}"
OPTIONAL_LONG_READS_ADDITION = "--nanopore {long_reads} "

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
        output_lines = [output_line for output_line in tqdm(iter(lambda: p.stdout.readline(), b""))]
    return output_folder
