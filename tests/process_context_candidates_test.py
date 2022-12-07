import verify_context_candidates as pcc
from Bio import SeqIO


def extract_genes_lengths(genes_path):
    gene_lengths = {}
    with open(genes_path) as f:
        fasta_sequences = SeqIO.parse(f, 'fasta')
        for s in fasta_sequences:
            gene_length = len(s.seq)
            gene_name = s.id

            gene_lengths[gene_name] = gene_length
    return gene_lengths


if __name__ == '__main__':
    LAB = '/Volumes/netta.barak/lab'
    genes_path = f'{LAB}/Projects/hybrid_assembler/pipeline_tree/genes_from_amir_erez/processed/genes_without_dup_locations.fasta'
    in_path_mapping_to_bugs = f'{LAB}/Projects/hybrid_assembler/test_run/in_paths_to_reference.paf'
    out_path_mapping_to_bugs = f'{LAB}/Projects/hybrid_assembler/test_run/out_paths_to_reference.paf'

    genes_lengths = extract_genes_lengths(genes_path)
    paths_pident_filtering_th = 0.9
    minimal_gap_ratio = 0
    maximal_gap_ratio = 1.5
    pcc.process_in_and_out_paths_to_results(in_path_mapping_to_bugs, out_path_mapping_to_bugs, genes_lengths,
                                            paths_pident_filtering_th, minimal_gap_ratio, maximal_gap_ratio)
