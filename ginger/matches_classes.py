from pafpy import PafRecord
import re

class PathRefGenomeMatch:
    def __init__(self, paf_line: PafRecord, contigs_to_species: dict):
        query_name_splt = paf_line.qname.split('_path_')
        gene_and_nodes_list = query_name_splt[0].split('_nodes_')
        self.query_name = paf_line.qname
        self.path = query_name_splt[1] if len(query_name_splt) > 1 else None
        self.gene = gene_and_nodes_list[0].split('|')[0]
        self.nodes_list = gene_and_nodes_list[1] if len(gene_and_nodes_list) > 1 else None
        self.path_length = paf_line.qlen
        self.path_start = paf_line.qstart
        self.path_end = paf_line.qend

        self.strand = str(paf_line.strand)

        self.ref_genome = paf_line.tname
        self.species = contigs_to_species.get(re.split(self.ref_genome, '[\._]')[0], f'unknown_{self.ref_genome}')
        self.ref_genome_length = paf_line.tlen
        self.ref_genome_start = paf_line.tstart
        self.ref_genome_end = paf_line.tend

        self.score = paf_line.mlen / min(self.path_length, self.ref_genome_end, self.ref_genome_length - self.ref_genome_start)

    def __str__(self):
        return f'{self.ref_genome} {self.species} {self.gene} nodes: {self.nodes_list} path: {self.path} match score:{self.score} strand: {self.strand} {self.ref_genome_start} to {self.ref_genome_end}'


class GeneContigMatch:
    def __init__(self, mmseq_line: str):
        target, query, tstart, tend, nident, qlen = mmseq_line.split('\t')
        self.gene_length = int(qlen) * 3
        self.gene = query
        start_int = int(tstart)
        end_int = int(tend)

        self.strand = '+' if start_int < end_int else '-'

        self.contig = target
        self.start = min(int(tstart), int(tend))
        self.end = max(int(tstart), int(tend))

        self.score = int(nident) / int(qlen)
        self.nodes_list = None
        self.start_in_first_node = None

    def __str__(self):
        return f'{self.gene} {self.contig} {self.score} {self.nodes_list}'


class InOutPathsMatch:
    def __init__(self, in_path, out_path, start, end, gap_ratio, score, gene_length, gene=None, ref_genome=None):
        self.in_path = in_path
        self.out_path = out_path
        if gene is None:
            self.gene = self.in_path.gene
        else:
            self.gene = gene

        if ref_genome is None:
            self.ref_genome = self.in_path.ref_genome
        else:
            self.ref_genome = ref_genome

        self.start = start
        self.end = end
        self.gap_ratio = gap_ratio
        self.score = score
        self.gene_length = gene_length

    def to_dict(self):
        out_dict = dict(gene=self.gene, ref_genome=self.ref_genome, start=self.start, end=self.end, score=self.score)
        if self.in_path:
            out_dict.update(dict(in_path_name=self.in_path.query_name, out_path_name=self.out_path.query_name,
                                 in_path_score=self.in_path.score, out_path_score=self.out_path.score))

        return out_dict
