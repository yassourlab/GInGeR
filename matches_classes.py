class PathBugMatch:
    def __init__(self, paf_line):
        query_name_splt = paf_line.query_name.split('_path_')
        gene_and_nodes_list = query_name_splt[0].split('_nodes_')
        self.query_name = paf_line.query_name
        self.path = query_name_splt[1]
        self.gene = gene_and_nodes_list[0]
        self.nodes_list = gene_and_nodes_list[1]
        self.path_length = paf_line.query_length
        self.path_start = paf_line.query_start
        self.path_end = paf_line.query_end

        self.strand = paf_line.strand

        self.bug = paf_line.target_name
        self.bug_length = paf_line.target_length
        self.bug_start = paf_line.target_start
        self.bug_end = paf_line.target_end

        self.match_score = paf_line.residue_matches / self.path_length

    def __str__(self):
        return f'{self.gene} nodes: {self.nodes_list} path: {self.path} match score:{self.match_score} strand: {self.strand} {self.bug_start} to {self.bug_end}'


class GeneRefGTMatch:
    def __init__(self, paf_line):
        self.gene_species_start_end = paf_line.query_name
        self.gene_with_context_length = paf_line.query_length
        self.gene_with_context_start = paf_line.query_start
        self.gene_with_context_end = paf_line.query_end

        self.gene = self.gene_species_start_end.split('|')[0]

        self.strand = paf_line.strand

        self.bug = paf_line.target_name
        self.bug_length = paf_line.target_length
        self.bug_start = paf_line.target_start
        self.bug_end = paf_line.target_end

        self.match_score = paf_line.residue_matches / self.gene_with_context_length

    def __str__(self):
        return f'{self.gene_species_start_end} {self.bug} {self.match_score}'


class GeneContigMatch:
    def __init__(self, paf_line):
        self.gene = paf_line.query_name
        self.gene_length = paf_line.query_length
        self.gene_start = paf_line.query_start
        self.gene_end = paf_line.query_end

        self.strand = paf_line.strand

        self.contig = paf_line.target_name
        self.contig_length = paf_line.target_length
        self.contig_start = paf_line.target_start
        self.contig_end = paf_line.target_end

        self.match_score = paf_line.residue_matches / self.gene_length
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
