import os
import shutil
import tempfile
import unittest
from pathlib import Path

import pyfastg
from Bio import SeqIO

from ginger import locating_genes_in_graph as lg
from ginger import matches_classes as mc
from ginger import pipeline_utils as pu


class TempDirTestCase(unittest.TestCase):
    """A test case with a temp dir of its own in self.tmp_dir, removed when the test ends.

    Everything a test writes belongs in there - a test that writes a relative path instead leaves its
    output in whatever directory the suite was started from, and one of them used to do exactly that.
    """

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        # addCleanup rather than tearDown, so it still runs when a subclass's own setUp raises
        self.addCleanup(shutil.rmtree, self.tmp_dir, ignore_errors=True)

    def tmp_path(self, name) -> str:
        return os.path.join(self.tmp_dir, name)

    def write_fasta(self, name, records) -> str:
        """A fasta of {header: sequence} written into the temp dir, returning its path."""
        path = self.tmp_path(name)
        with open(path, 'w') as f:
            for header, seq in records.items():
                f.write(f'>{header}\n{seq}\n')
        return path


def get_filedir() -> str:
    currentdir = Path(__file__).resolve().parent
    return f"{currentdir}/test_files"


def get_metadata_path() -> str:
    """The UHGG metadata table that ships inside the package.

    Resolved from where ginger is installed rather than from the working directory, so that the tests
    that need it pass whether they are run from the repo root the way CI does or from anywhere else.
    """
    return str(Path(pu.__file__).resolve().parent / 'UHGG-metadata.tsv')


class FakeGeneMatch:
    """A gene located on a contig. Only the gene and the contig matter to the contexts fasta - a trio
    takes its sequence from the locus its match carries, so this is left for the unmatched contig
    listing."""

    def __init__(self, gene, contig, score):
        self.gene = gene
        self.contig = contig
        self.score = score


class FakePathMatch:
    def __init__(self, query_name):
        self.query_name = query_name


class FakeInOutMatch:
    def __init__(self, gene, in_context, out_context, locus=None):
        self.gene = gene
        self.in_path = FakePathMatch(in_context)
        self.out_path = FakePathMatch(out_context)
        self.locus = locus


def get_contig_seq(contig_name, contigs_path=None):
    """The sequence of a single contig, read without leaving the fasta handle open the way
    SeqIO.index does."""
    contigs_path = contigs_path or f'{get_filedir()}/SPAdes/contigs.fasta'
    with open(contigs_path) as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id == contig_name:
                return str(record.seq)
    raise KeyError(f'{contig_name} is not in {contigs_path}')


def get_assembly_graph():
    return pyfastg.parse_fastg(f'{get_filedir()}/SPAdes/assembly_graph.fastg')


def get_assembly_graph_nodes():
    return lg.get_nodes_dict_from_fastg_file(f'{get_filedir()}/SPAdes/assembly_graph.fastg')


def get_geometry():
    return pu.PathGeometry(get_assembly_graph_nodes())


def get_genes_with_location_in_graph():
    """The test_gene match from genes_to_contigs.m8, located in the assembly graph exactly the way
    locate_genes_in_graph does it.

    Built from the fixtures on every call rather than unpickled, so that it can never go stale
    against the coordinate convention in matches_classes.GeneContigMatch.
    """
    files = get_filedir()
    parsed_paths, _ = pu.parse_paths_file(f'{files}/SPAdes/contigs.paths')
    with open(f'{files}/genes_to_contigs.m8') as f:
        next(f)  # skip header
        genes_to_contigs = [mc.GeneContigMatch(line) for line in f]
    return lg.add_node_list_to_genes_to_contigs(genes_to_contigs, parsed_paths, get_geometry(), {})
