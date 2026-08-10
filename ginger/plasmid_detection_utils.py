import os
import logging
import shutil
import pandas as pd

from ginger import pipeline_utils as pu

log = logging.getLogger(__name__)

GENOMAD_COMMAND = 'genomad end-to-end --threads {threads} --cleanup {fasta_path} {output_dir} {genomad_db} --min-score 0'


@pu.step_timing
def run_genomad(fasta_path, output_dir, genomad_db, threads):
    if not os.path.exists(genomad_db):
        raise Exception(f'GeNomad database does not exist in {genomad_db}')

    pu.check_and_make_dir_no_file_name(output_dir)
    pu.run_tool('GeNomad', GENOMAD_COMMAND.format(threads=threads, fasta_path=fasta_path, output_dir=output_dir,
                                                  genomad_db=genomad_db))

    fasta_stem = os.path.splitext(os.path.basename(fasta_path))[0]
    return os.path.join(output_dir, f'{fasta_stem}_summary', f'{fasta_stem}_plasmid_summary.tsv')


def keep_only_plasmid_summary(genomad_out_dir, plasmid_summary_path, kept_summary_path):
    """Moves the one file GInGeR reads out of GeNomad's output tree and removes the rest of it.

    An end-to-end run leaves ~80MB per sample of which the plasmid summary is ~100KB; everything else
    is per-module intermediates - protein fastas, mmseqs2 hits, feature matrices - that nothing
    downstream reads. `--cleanup`, which GENOMAD_COMMAND already passes, does not remove any of them.

    Call this only once the summary has been read, so that a run which failed anywhere inside GeNomad
    keeps its whole output tree to be debugged.
    """
    genomad_out_dir, kept_summary_path = os.path.abspath(genomad_out_dir), os.path.abspath(kept_summary_path)
    if os.path.commonpath([genomad_out_dir, kept_summary_path]) == genomad_out_dir:
        raise ValueError(f'{kept_summary_path} is inside the GeNomad output tree {genomad_out_dir} that is '
                         f'about to be removed - the summary has to be kept outside of it')

    shutil.move(plasmid_summary_path, kept_summary_path)
    shutil.rmtree(genomad_out_dir)
    log.info(f"kept GeNomad's plasmid summary as {kept_summary_path} and removed the rest of {genomad_out_dir}")
    return kept_summary_path


def read_plasmid_scores(plasmid_summary_path, records_by_seq_id):
    """Reads GeNomad's plasmid_summary.tsv and splits the results into context-level and
    contig-level plasmid scores.

    A sequence is a gene context if records_by_seq_id - as returned by
    pipeline_utils.write_in_gene_out_contexts_fasta - names it, and a contig otherwise. Membership is
    exact, so no contig can be mistaken for a context whatever it is called.

    Returns a tuple (context_plasmid_scores, contig_plasmid_scores):
    - context_plasmid_scores has columns [gene, in_context, out_context, plasmid_score]
    - contig_plasmid_scores has columns [contig, plasmid_score]
    """
    summary_df = pd.read_csv(plasmid_summary_path, sep='\t')[['seq_name', 'plasmid_score']]
    is_context = summary_df['seq_name'].isin(records_by_seq_id)

    if records_by_seq_id and len(summary_df) and not is_context.any():
        raise ValueError(f'none of the {len(summary_df)} sequences GeNomad reported on is one of the '
                         f'{len(records_by_seq_id)} gene contexts written for it - the contexts were '
                         f'not produced by this run, so every context would silently score 0')

    context_plasmid_scores = summary_df[is_context].copy()
    if context_plasmid_scores.empty:
        context_plasmid_scores = pd.DataFrame(columns=['gene', 'in_context', 'out_context', 'plasmid_score'])
    else:
        records = context_plasmid_scores['seq_name'].map(records_by_seq_id)
        context_plasmid_scores['gene'] = [record.gene for record in records]
        context_plasmid_scores['in_context'] = [record.in_context for record in records]
        context_plasmid_scores['out_context'] = [record.out_context for record in records]
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
