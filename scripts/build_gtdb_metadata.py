"""Build ginger/GTDB-metadata.tsv from GTDB's bac120 + ar53 metadata tables.

GInGeR's reference metadata table is a slim view of the reference catalog: one row per
reference genome, with the columns the pipeline actually reads. This turns GTDB's ~110-column
release tables into that view, so the table can be regenerated for a future GTDB release.

    python scripts/build_gtdb_metadata.py --out ginger/GTDB-metadata.tsv

Source tables default to GTDB's release site; pass --bac120/--ar53 to use local copies.

Note on quality columns: GTDB r226 carries both checkm and checkm2 estimates. checkm2 is the
one GTDB itself uses for its quality filtering from r207 on, so it is what feeds GInGeR's
Quality score (Completeness - 5 * Contamination + ln(N50)).
"""
import argparse
import csv
import gzip
import logging
import os
import sys
import urllib.request

log = logging.getLogger(__name__)

GTDB_RELEASE_BASE = 'https://data.gtdb.ecogenomic.org/releases/release226/226.0'
BAC120_URL = f'{GTDB_RELEASE_BASE}/bac120_metadata_r226.tsv.gz'
AR53_URL = f'{GTDB_RELEASE_BASE}/ar53_metadata_r226.tsv.gz'

# The slim table's header. Matches the contract documented in README for
# --reference-genomes-metadata, minus UHGG's FTP_download - GTDB genomes are fetched from NCBI
# by accession, so there is no per-genome URL to carry.
OUTPUT_COLUMNS = ['Genome', 'Completeness', 'Contamination', 'N50', 'Length', 'species']

# slim column -> GTDB column
SOURCE_COLUMNS = {
    'Completeness': 'checkm2_completeness',
    'Contamination': 'checkm2_contamination',
    'N50': 'n50_contigs',
    'Length': 'genome_size',
}
ACCESSION_COLUMN = 'accession'
TAXONOMY_COLUMN = 'gtdb_taxonomy'

# GTDB prefixes an accession by the archive it came from; NCBI's datasets CLI wants it without
ACCESSION_PREFIXES = ('GB_', 'RS_')


def strip_accession_prefix(accession: str) -> str:
    for prefix in ACCESSION_PREFIXES:
        if accession.startswith(prefix):
            return accession[len(prefix):]
    return accession


def species_from_taxonomy(taxonomy: str) -> str:
    """The species name of a GTDB lineage: the s__ rank of 'd__Bacteria;...;s__Escherichia coli'.

    Returned without the 's__' prefix, which is how Kraken2 reports species names for a
    GTDB-derived database - so these join straight onto Bracken's output.
    """
    for rank in reversed(taxonomy.split(';')):
        if rank.startswith('s__'):
            return rank[len('s__'):].strip()
    return ''


def open_maybe_download(path_or_url: str, download_dir: str):
    """A readable text handle on a local .tsv.gz, downloading it first if given a URL."""
    if not path_or_url.startswith(('http://', 'https://')):
        return gzip.open(path_or_url, 'rt')

    local_path = os.path.join(download_dir, os.path.basename(path_or_url))
    if os.path.exists(local_path):
        log.info(f'using already downloaded {local_path}')
    else:
        log.info(f'downloading {path_or_url} -> {local_path}')
        os.makedirs(download_dir, exist_ok=True)
        # via a .part file so an interrupted download cannot leave a truncated table that the
        # next run would treat as complete
        urllib.request.urlretrieve(path_or_url, f'{local_path}.part')
        os.replace(f'{local_path}.part', local_path)
    return gzip.open(local_path, 'rt')


def slim_rows(handle, source_name: str):
    """The slim rows of one GTDB metadata table, skipping rows missing anything required."""
    # QUOTE_NONE: a '"' inside a free-text field (isolation source, organism name) is data, not
    # a quote character, and letting csv treat it as one would swallow the rows that follow
    reader = csv.DictReader(handle, delimiter='\t', quoting=csv.QUOTE_NONE)
    missing = [c for c in [ACCESSION_COLUMN, TAXONOMY_COLUMN, *SOURCE_COLUMNS.values()]
               if c not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f'{source_name} is missing required columns: {missing}')

    kept = skipped = 0
    for row in reader:
        genome = strip_accession_prefix(row[ACCESSION_COLUMN].strip())
        species = species_from_taxonomy(row[TAXONOMY_COLUMN])
        values = {'Genome': genome, 'species': species}
        values.update({slim: row[src].strip() for slim, src in SOURCE_COLUMNS.items()})
        if not all(values[c] for c in OUTPUT_COLUMNS):
            skipped += 1
            continue
        kept += 1
        yield [values[c] for c in OUTPUT_COLUMNS]
    log.info(f'{source_name}: kept {kept} genomes, skipped {skipped} with missing fields')


def build(bac120: str, ar53: str, out_path: str, download_dir: str):
    species = set()
    total = 0
    tmp_path = f'{out_path}.part'
    with open(tmp_path, 'w', newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t', lineterminator='\n')
        writer.writerow(OUTPUT_COLUMNS)
        for source, name in [(bac120, 'bac120'), (ar53, 'ar53')]:
            with open_maybe_download(source, download_dir) as handle:
                for row in slim_rows(handle, name):
                    writer.writerow(row)
                    species.add(row[OUTPUT_COLUMNS.index('species')])
                    total += 1
    os.replace(tmp_path, out_path)
    log.info(f'wrote {out_path}: {total} genomes across {len(species)} species')
    return total, len(species)


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bac120', default=BAC120_URL, help='bac120 metadata .tsv.gz (path or URL)')
    parser.add_argument('--ar53', default=AR53_URL, help='ar53 metadata .tsv.gz (path or URL)')
    parser.add_argument('--out', default=os.path.join(os.path.dirname(__file__), '..', 'ginger',
                                                      'GTDB-metadata.tsv'),
                        help="where to write the slim table")
    parser.add_argument('--download-dir', default='.',
                        help='where to keep source tables downloaded from GTDB')
    args = parser.parse_args(argv)

    build(args.bac120, args.ar53, args.out, args.download_dir)
    return 0


if __name__ == '__main__':
    sys.exit(main())
