from subprocess import run

# TODO use logging instead of prints
# TODO change path to be in the same directory + add links to my directory
SPADES = '/sci/labs/morani/morani/icore-data/lab/Tools/SPAdes-3.15.2-Linux/bin/spades.py'
HYBRID_SPADES_COMMAND = "{spades} --meta -1 {short_reads_1} -2 {short_reads_2} --nanopore {long_reads} -o {output_folder} -t {threads}"
META_SPADES_COMMAND = "{spades} --meta -1 {short_reads_1} -2 {short_reads_2} -o {output_folder} -t {threads}"


def run_meta_or_hybrid_spades(short_reads_1, short_reads_2, long_reads, output_folder, threads):
    if long_reads:
        command = HYBRID_SPADES_COMMAND.format(spades=SPADES, short_reads_1=short_reads_1, short_reads_2=short_reads_2,
                                               long_reads=long_reads, output_folder=output_folder, threads=threads)
    else:
        command = META_SPADES_COMMAND.format(spades=SPADES, short_reads_1=short_reads_1, short_reads_2=short_reads_2,
                                             output_folder=output_folder, threads=threads)
    print(command)
    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        raise Exception(f'Hybrid spades failed to run - {command_output}')
    return output_folder
