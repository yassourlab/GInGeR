import os
import logging
from subprocess import run
import pandas as pd
from Bio import SeqIO

from ginger import pipeline_utils as pu

log = logging.getLogger(__name__)

GENOMAD_COMMAND = 'genomad end-to-end --threads {threads} --cleanup {fasta_path} {output_dir} {genomad_db} --min-score 0'


def _fasta_to_dict(fasta_path: str) -> dict:
    with open(fasta_path) as f:
        return {rec.id: str(rec.seq) for rec in SeqIO.parse(f, 'fasta')}


def _get_gene_sequence(contig_seq: str, locus) -> str:
    """The gene segment to splice between a pair of contexts, in the contig's forward orientation.

    Every route that produces a context produces it in that orientation: the ones that read
    contigs.paths get their nodes in contig order, and the contig fallback slices the contig itself.
    So the segment stays forward - reverse complementing it would splice a flipped middle into
    forward-oriented flanks. A locus carries no strand, so there is nothing here to be tempted by.

    Its coordinates are 0-based half-open (see matches_classes.GeneLocus), so this slice is exactly
    the aligned part of the gene.
    """
    return contig_seq[locus.start:locus.end]


TRIOS_TABLE_COLUMNS = ['seq_id', 'gene', 'in_context', 'out_context']
CONTEXT_SEQ_ID_PREFIX = 'ctx'


def _write_trios_table(trios_by_seq_id, trios_table_path):
    with open(trios_table_path, 'w') as f:
        f.write('\t'.join(TRIOS_TABLE_COLUMNS) + '\n')
        for seq_id, trio in trios_by_seq_id.items():
            f.write('\t'.join((seq_id,) + trio) + '\n')


def write_plasmid_detection_input_fasta(context_level_results, genes_with_location_in_graph, matched_genes,
                                        in_paths_fasta, out_paths_fasta, contigs_fasta, output_fasta_path,
                                        trios_table_path=None):
    """Writes a FASTA file to be used as GeNomad's input, containing:
    - for every unique (gene, in_context, out_context) trio in context_level_results, the
      concatenation of the in-path, gene and out-path sequences, named "ctx0000001" and up
    - for every gene in genes_with_location_in_graph that is not in matched_genes, the full
      sequence of the contig it was found on, named after the contig

    The gene sequence comes from the copy of the gene the trio's two contexts were cut from, which
    the match carries: they are only ever paired within one copy, so splicing that copy's sequence
    between them reproduces a stretch of the contig exactly.

    A trio used to be named "{gene}|{in_context}|{out_context}", which cannot be taken apart again
    when the gene's own name contains a '|' - as SARG's and CARD's do. Gene names come from a fasta
    the caller supplies, so no separator is safe; the names here are opaque, and what they stand for
    is returned as a mapping and written to trios_table_path when one is given.

    Returns the path to the written fasta and that mapping, or (None, {}) if there was nothing to
    write.
    """
    contig_seq_by_id = _fasta_to_dict(contigs_fasta)
    trios_by_seq_id = {}

    wrote_any = False
    with open(output_fasta_path, 'w') as f:
        if context_level_results:
            in_seq_by_id = _fasta_to_dict(in_paths_fasta)
            out_seq_by_id = _fasta_to_dict(out_paths_fasta)

            trios = set()
            for matches_list in context_level_results.values():
                for match in matches_list:
                    trios.add((match.gene, match.in_path.query_name, match.out_path.query_name, match.locus))

            # sorted, so that a rerun on the same input numbers the sequences the same way - and so
            # that GeNomad's own per-sequence gene numbering stays comparable between runs
            for n, (gene, in_context, out_context, locus) in enumerate(sorted(trios)):
                gene_seq = _get_gene_sequence(contig_seq_by_id[locus.contig], locus)
                full_seq = in_seq_by_id[in_context] + gene_seq + out_seq_by_id[out_context]
                seq_id = f'{CONTEXT_SEQ_ID_PREFIX}{n:07d}'
                trios_by_seq_id[seq_id] = (gene, in_context, out_context)
                f.write(f'>{seq_id}\n{full_seq}\n')
                wrote_any = True

        written_contigs = set()
        for gene_match in genes_with_location_in_graph:
            if gene_match.gene not in matched_genes and gene_match.contig not in written_contigs:
                f.write(f'>{gene_match.contig}\n{contig_seq_by_id[gene_match.contig]}\n')
                written_contigs.add(gene_match.contig)
                wrote_any = True

    if not wrote_any:
        os.remove(output_fasta_path)
        return None, {}
    if trios_table_path is not None:
        _write_trios_table(trios_by_seq_id, trios_table_path)
    return output_fasta_path, trios_by_seq_id


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


def read_plasmid_scores(plasmid_summary_path, trios_by_seq_id):
    """Reads GeNomad's plasmid_summary.tsv and splits the results into context-level and
    contig-level plasmid scores.

    A sequence is a gene context if trios_by_seq_id - as returned by
    write_plasmid_detection_input_fasta - names it, and a contig otherwise. Membership is exact, so
    no contig can be mistaken for a context whatever it is called.

    Returns a tuple (context_plasmid_scores, contig_plasmid_scores):
    - context_plasmid_scores has columns [gene, in_context, out_context, plasmid_score]
    - contig_plasmid_scores has columns [contig, plasmid_score]
    """
    summary_df = pd.read_csv(plasmid_summary_path, sep='\t')[['seq_name', 'plasmid_score']]
    is_context = summary_df['seq_name'].isin(trios_by_seq_id)

    if trios_by_seq_id and len(summary_df) and not is_context.any():
        raise ValueError(f'none of the {len(summary_df)} sequences GeNomad reported on is one of the '
                         f'{len(trios_by_seq_id)} gene contexts written for it - the trios were not '
                         f'produced by this run, so every context would silently score 0')

    context_plasmid_scores = summary_df[is_context].copy()
    if context_plasmid_scores.empty:
        context_plasmid_scores = pd.DataFrame(columns=['gene', 'in_context', 'out_context', 'plasmid_score'])
    else:
        trios = context_plasmid_scores['seq_name'].map(trios_by_seq_id)
        context_plasmid_scores['gene'] = [trio[0] for trio in trios]
        context_plasmid_scores['in_context'] = [trio[1] for trio in trios]
        context_plasmid_scores['out_context'] = [trio[2] for trio in trios]
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
