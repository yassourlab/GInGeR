#!/bin/bash
#SBATCH -c 4
#SBATCH -t 2-0
#SBATCH --mem=8G
#SBATCH -J gtdb_k2_download
#SBATCH -o gtdb_k2_download_%j.out
#
# Downloads the prebuilt GTDB r226 Kraken2 database (k2_gtdb_genome_reps_20250609),
# file by file rather than as the 538GB .tar.gz: same bytes on disk, resumable per file,
# and nothing to unpack and then delete.
#
#   sbatch scripts/download_gtdb_kraken_db.sbatch
#
# Safe to re-run: curl -C - resumes partial files, and files that already match their
# md5 are skipped entirely. ~645 GiB, so expect a few hours.

set -uo pipefail

BASE_URL='https://genome-idx.s3.amazonaws.com/kraken/gtdb_genome_reps_20250609'
DEST='/cs/usr/netta.barak/lab/Tools/kraken2_db/gtdb_genome_reps_20250609'

# library_report.tsv (5.66 GiB) is a provenance record listing the genomes that went into
# the build. Kraken2 and Bracken never read it. Set to 1 if you want it.
WANT_LIBRARY_REPORT=0

# hash/opts/taxo are what kraken2 classifies with; seqid2taxid maps sequences to taxa;
# the kmer_distrib files are what bracken reads; the rest is small and useful to keep.
FILES=(
  opts.k2d
  taxo.k2d
  inspect.txt
  ktaxonomy.tsv
  nodes.dmp
  names.dmp
  database50mers.kmer_distrib
  database75mers.kmer_distrib
  database100mers.kmer_distrib
  database150mers.kmer_distrib
  database200mers.kmer_distrib
  database250mers.kmer_distrib
  database300mers.kmer_distrib
  seqid2taxid.map
  hash.k2d            # 644 GiB - last, so everything small is in place if this needs a retry
)
(( WANT_LIBRARY_REPORT )) && FILES+=( library_report.tsv )

mkdir -p "$DEST" || exit 1
cd "$DEST" || exit 1

NEEDED_GIB=646
(( WANT_LIBRARY_REPORT )) && NEEDED_GIB=652
AVAIL_GIB=$(df -BG --output=avail "$DEST" | tail -1 | tr -dc '0-9')
echo "=== destination: $DEST"
echo "=== free space: ${AVAIL_GIB} GiB, need ~${NEEDED_GIB} GiB"
if (( AVAIL_GIB < NEEDED_GIB )); then
    echo "ERROR: not enough free space. Aborting before writing anything." >&2
    exit 1
fi

echo "=== fetching md5 manifest"
curl -fsS --retry 5 --retry-delay 10 -o gtdb_genome_reps.md5 "$BASE_URL/gtdb_genome_reps.md5" || exit 1

# One entry per file, so a file can be checked without running md5sum over all 645 GiB.
expected_md5() { awk -v f="$1" '$2==f {print $1}' gtdb_genome_reps.md5; }

verify() {
    local f="$1" want got
    want=$(expected_md5 "$f")
    [[ -z "$want" ]] && { echo "  no md5 listed for $f - skipping check"; return 0; }
    [[ -f "$f" ]] || return 1
    echo "  verifying $f (this reads the whole file)"
    got=$(md5sum "$f" | awk '{print $1}')
    [[ "$got" == "$want" ]]
}

FAILED=()
for f in "${FILES[@]}"; do
    echo
    echo "=== $f"
    if verify "$f"; then
        echo "  already present and md5 matches - skipping"
        continue
    fi

    ok=0
    for attempt in 1 2 3; do
        echo "  download attempt $attempt"
        # -C - resumes a partial file; on attempt 3 start clean, in case the partial
        # file is corrupt rather than merely truncated
        (( attempt == 3 )) && rm -f "$f"
        if curl -fsS -C - --retry 5 --retry-delay 30 --speed-limit 1024 --speed-time 300 \
                -o "$f" "$BASE_URL/$f"; then
            if verify "$f"; then ok=1; break; fi
            echo "  md5 mismatch after download"
        else
            echo "  curl failed"
        fi
        sleep 60
    done

    if (( ok )); then echo "  OK"; else echo "  FAILED"; FAILED+=("$f"); fi
done

echo
echo "================================================"
df -h "$DEST" | tail -1
du -sh "$DEST"
if (( ${#FAILED[@]} )); then
    echo "FAILED: ${FAILED[*]}"
    echo "Re-run this script - completed files are skipped."
    exit 1
fi
echo "All files downloaded and verified."
echo "Kraken DB ready at: $DEST"
