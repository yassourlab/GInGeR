import os
import logging
from subprocess import run
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from ginger import pipeline_utils as pu

log = logging.getLogger(__name__)

GENOMAD_COMMAND = 'genomad end-to-end --threads {threads} --cleanup {fasta_path} {output_dir} {genomad_db} --min-score 0'


def _fasta_to_dict(fasta_path: str) -> dict:
    with open(fasta_path) as f:
        return {rec.id: str(rec.seq) for rec in SeqIO.parse(f, 'fasta')}


def _get_gene_sequence(contig_seq: str, gene_match) -> str:
    seq = contig_seq[gene_match.start - 1:gene_match.end]
    if gene_match.strand == '-':
        seq = str(Seq(seq).reverse_complement())
    return seq


def write_plasmid_detection_input_fasta(context_level_results, genes_with_location_in_graph, matched_genes,
                                        in_paths_fasta, out_paths_fasta, contigs_fasta, output_fasta_path):
    """Writes a FASTA file to be used as GeNomad's input, containing:
    - for every unique (gene, in_context, out_context) trio in context_level_results, the
      concatenation of the in-path, gene and out-path sequences, with header "{gene}|{in_context}|{out_context}"
    - for every gene in genes_with_location_in_graph that is not in matched_genes, the full
      sequence of the contig it was found on, with header "{contig}"

    Returns the path to the written fasta, or None if there was nothing to write.
    """
    best_match_by_gene = {}
    for gene_match in genes_with_location_in_graph:
        current_best = best_match_by_gene.get(gene_match.gene)
        if current_best is None or gene_match.score > current_best.score:
            best_match_by_gene[gene_match.gene] = gene_match

    contig_seq_by_id = _fasta_to_dict(contigs_fasta)

    wrote_any = False
    with open(output_fasta_path, 'w') as f:
        if context_level_results:
            in_seq_by_id = _fasta_to_dict(in_paths_fasta)
            out_seq_by_id = _fasta_to_dict(out_paths_fasta)

            trios = set()
            for matches_list in context_level_results.values():
                for match in matches_list:
                    trios.add((match.gene, match.in_path.query_name, match.out_path.query_name))

            for gene, in_context, out_context in trios:
                gene_match = best_match_by_gene.get(gene)
                if gene_match is None:
                    continue
                gene_seq = _get_gene_sequence(contig_seq_by_id[gene_match.contig], gene_match)
                full_seq = in_seq_by_id[in_context] + gene_seq + out_seq_by_id[out_context]
                f.write(f'>{gene}|{in_context}|{out_context}\n{full_seq}\n')
                wrote_any = True

        written_contigs = set()
        for gene_match in genes_with_location_in_graph:
            if gene_match.gene not in matched_genes and gene_match.contig not in written_contigs:
                f.write(f'>{gene_match.contig}\n{contig_seq_by_id[gene_match.contig]}\n')
                written_contigs.add(gene_match.contig)
                wrote_any = True

    if not wrote_any:
        os.remove(output_fasta_path)
        return None
    return output_fasta_path


@pu.step_timing
def run_genomad(fasta_path, output_dir, genomad_db, threads):
    if not os.path.exists(genomad_db):
        raise Exception(f'GeNomad database does not exist in {genomad_db}')

    pu.check_and_make_dir_no_file_name(output_dir)
    command = GENOMAD_COMMAND.format(threads=threads, fasta_path=fasta_path, output_dir=output_dir,
                                     genomad_db=genomad_db)
    log.info(f'Running GeNomad - {command}')
    command_output = run(command, shell=True, capture_output=True)
    if command_output.returncode:
        log.error(f'GeNomad failed: {command_output.stderr}')
        raise Exception('GeNomad failed - GInGeR aborted')
    log.info('GeNomad completed successfully')

    fasta_stem = os.path.splitext(os.path.basename(fasta_path))[0]
    return os.path.join(output_dir, f'{fasta_stem}_summary', f'{fasta_stem}_plasmid_summary.tsv')


def read_plasmid_scores(plasmid_summary_path):
    """Reads GeNomad's plasmid_summary.tsv and splits the results into context-level and
    contig-level plasmid scores based on whether seq_name is a "{gene}|{in_context}|{out_context}"
    composite header or a plain contig name.

    Returns a tuple (context_plasmid_scores, contig_plasmid_scores):
    - context_plasmid_scores has columns [gene, in_context, out_context, plasmid_score]
    - contig_plasmid_scores has columns [contig, plasmid_score]
    """
    summary_df = pd.read_csv(plasmid_summary_path, sep='\t')[['seq_name', 'plasmid_score']]
    is_context = summary_df['seq_name'].str.contains('|', regex=False)

    context_plasmid_scores = summary_df[is_context].copy()
    if context_plasmid_scores.empty:
        context_plasmid_scores = pd.DataFrame(columns=['gene', 'in_context', 'out_context', 'plasmid_score'])
    else:
        context_split = context_plasmid_scores['seq_name'].str.split('|', expand=True)
        context_plasmid_scores['gene'] = context_split[0]
        context_plasmid_scores['in_context'] = context_split[1]
        context_plasmid_scores['out_context'] = context_split[2]
        context_plasmid_scores = context_plasmid_scores[['gene', 'in_context', 'out_context', 'plasmid_score']]

    contig_plasmid_scores = summary_df[~is_context].rename(columns={'seq_name': 'contig'})[['contig', 'plasmid_score']]

    return context_plasmid_scores, contig_plasmid_scores


def add_plasmid_scores_to_context_level_csv(context_level_csv_path, context_plasmid_scores):
    context_level_df = pd.read_csv(context_level_csv_path)
    context_level_df = context_level_df.merge(context_plasmid_scores, on=['gene', 'in_context', 'out_context'],
                                              how='left')
    context_level_df['plasmid_score'] = context_level_df['plasmid_score'].fillna(0)
    context_level_df.to_csv(context_level_csv_path, index=False)


def add_plasmid_scores_to_genes_no_species_match_csv(csv_path, contig_plasmid_scores):
    genes_df = pd.read_csv(csv_path)
    genes_df = genes_df.merge(contig_plasmid_scores, on='contig', how='left')
    genes_df['plasmid_score'] = genes_df['plasmid_score'].fillna(0)
    genes_df.to_csv(csv_path, index=False)
