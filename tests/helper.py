from pathlib import Path

import pandas as pd
import pyfastg
from Bio import SeqIO

from ginger import locating_genes_in_graph as lg
from ginger import matches_classes as mc
from ginger import pipeline_utils as pu


def get_filedir() -> str:
    currentdir = Path(__file__).resolve().parent
    return f"{currentdir}/test_files"


def get_contig_seq(contig_name, contigs_path=None):
    """The sequence of a single contig, read without leaving the fasta handle open the way
    SeqIO.index does."""
    contigs_path = contigs_path or f'{get_filedir()}/SPAdes/contigs.fasta'
    with open(contigs_path) as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id == contig_name:
                return str(record.seq)
    raise KeyError(f'{contig_name} is not in {contigs_path}')


# the columns map_nodes_to_contigs_w_gaps produces, needed so that filtering works on an empty frame
NODE_TO_CONTIG_COLUMNS = ['contig', 'contig_start', 'contig_end', 'node', 'score', 'strand']


def get_assembly_graph():
    return pyfastg.parse_fastg(f'{get_filedir()}/SPAdes/assembly_graph.fastg')


def get_assembly_graph_nodes():
    return lg.get_nodes_dict_from_fastg_file(f'{get_filedir()}/SPAdes/assembly_graph.fastg')


def get_genes_with_location_in_graph():
    """The test_gene match from genes_to_contigs.m8, located in the assembly graph exactly the way
    locate_genes_in_graph does it.

    Built from the fixtures on every call rather than unpickled, so that it can never go stale
    against the coordinate convention in matches_classes.GeneContigMatch.
    """
    files = get_filedir()
    assembly_graph = get_assembly_graph()
    parsed_paths, _ = pu.parse_paths_file(f'{files}/SPAdes/contigs.paths', assembly_graph.nodes)
    with open(f'{files}/genes_to_contigs.m8') as f:
        next(f)  # skip header
        genes_to_contigs = [mc.GeneContigMatch(line) for line in f]
    return lg.add_node_list_to_genes_to_contigs(genes_to_contigs, parsed_paths, get_assembly_graph_nodes(),
                                                pd.DataFrame(columns=NODE_TO_CONTIG_COLUMNS))
