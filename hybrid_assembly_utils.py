from subprocess import run
import logging

log = logging.getLogger(__name__)

# TODO use logging instead of prints
# TODO change path to be in the same directory + add links to my directory
SPADES = 'spades.py'
META_SPADES_COMMAND = "{spades} --meta -1 {short_reads_1} -2 {short_reads_2} {optional_long_reads}-o {output_folder} -t {threads}"
OPTIONAL_LONG_READS_ADDITION = "--nanopore {long_reads} "


def run_meta_or_hybrid_spades(short_reads_1, short_reads_2, long_reads, output_folder, threads):
    long_reads_str = ''
    if long_reads:
        long_reads_str = OPTIONAL_LONG_READS_ADDITION.format(long_reads=long_reads)
    command = META_SPADES_COMMAND.format(spades=SPADES, short_reads_1=short_reads_1, short_reads_2=short_reads_2,
                                         optional_long_reads=long_reads_str, output_folder=output_folder,
                                         threads=threads)
    log.info(f'running MetaSPAdes - {command}')
    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        raise Exception(f'Hybrid spades failed to run - {command_output}')
    return output_folder
