from pafpy import PafRecord


class PathBugMatch:
    def __init__(self, paf_line: PafRecord):
        query_name_splt = paf_line.qname.split('_path_')
        gene_and_nodes_list = query_name_splt[0].split('_nodes_')
        self.query_name = paf_line.qname
        self.path = query_name_splt[1]
        self.gene = gene_and_nodes_list[0]
        self.nodes_list = gene_and_nodes_list[1]
        self.path_length = paf_line.qlen
        self.path_start = paf_line.qstart
        self.path_end = paf_line.qend

        self.strand = paf_line.strand

        self.bug = paf_line.tname
        self.bug_length = paf_line.tlen
        self.bug_start = paf_line.tstart
        self.bug_end = paf_line.tend

        self.match_score = paf_line.mlen / self.path_length

    def __str__(self):
        return f'{self.gene} nodes: {self.nodes_list} path: {self.path} match score:{self.match_score} strand: {self.strand} {self.bug_start} to {self.bug_end}'


# class GeneRefGTMatch:
#     def __init__(self, paf_line: PafRecord):
#         self.gene_species_start_end = paf_line.qname
#         self.gene_with_context_length = paf_line.qlen
#         self.gene_with_context_start = paf_line.qstart
#         self.gene_with_context_end = paf_line.qend
#
#         self.gene = self.gene_species_start_end.split('|')[0]
#
#         self.strand = paf_line.strand
#
#         self.bug = paf_line.tname
#         self.bug_length = paf_line.tlen
#         self.start = paf_line.tstart
#         self.end = paf_line.tend
#
#         self.match_score = paf_line.mlen / self.gene_with_context_length
#
#     def __str__(self):
#         return f'{self.gene_species_start_end} {self.bug} {self.match_score}'
#

class GeneContigMatch:
    def __init__(self, mmseq_line: str):
        # TODO make the list of fields shared with the one in sequence_alignment_utils.py
        target, query, tstart, tend, nident, qlen = mmseq_line.split('\t')
        # self.gene_species_start_end = query
        self.gene_length = int(qlen)
        # self.gene_with_context_start = paf_line.qstart
        # self.gene_with_context_end = paf_line.qend

        self.gene = query

        self.strand = '+' if int(tstart) < int(tend) else '-'

        self.contig = target
        # self.bug_length = paf_line.tlen
        self.start = int(tstart)
        self.end = int(tend)

        self.match_score = int(nident) / int(qlen)
        self.nodes_list = None
        self.start_in_first_node = None

    def __str__(self):
        return f'{self.gene} {self.contig} {self.match_score} {self.nodes_list}'


class InOutPathsMatch:
    def __init__(self, in_path, out_path, start, end, gap_ratio, match_score, gene_length, gene=None, bug=None):
        self.in_path = in_path
        self.out_path = out_path
        if gene is None:
            self.gene = self.in_path.gene
        else:
            self.gene = gene

        if bug is None:
            self.bug = self.in_path.bug
        else:
            self.bug = bug

        self.start = start
        self.end = end
        self.gap_ratio = gap_ratio
        self.match_score = match_score
        self.gene_length = gene_length

    def to_dict(self):
        out_dict = dict(gene=self.gene, bug=self.bug, start=self.start, end=self.end, match_score=self.match_score)
        if self.in_path:
            out_dict.update(dict(in_path_name=self.in_path.query_name, out_path_name=self.out_path.query_name,
                                 in_path_score=self.in_path.match_score, out_path_score=self.out_path.match_score))

        return out_dict
