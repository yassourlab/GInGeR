#!/bin/bash
#SBATCH -c 16
#SBATCH -t 1-0
#SBATCH -p short
#SBATCH --mem=750G
#SBATCH -J kraken_mmap_bench
#SBATCH -o kraken_mmap_bench_%j.out
#
# How much does --memory-mapping cost on the GTDB r226 Kraken2 DB?
#
# GTDB's hash.k2d is 644 GiB against UHGG's 15.5 GiB. --memory-mapping lets kraken2 read it
# off disk so GInGeR keeps needing only ~16-32GB of RAM, but the DB sits on NFS
# (isi-02-sci), so classification does random reads over the network. This measures whether
# that is acceptable before committing the 36-task subsampling grid to it.
#
#   sbatch scripts/benchmark_kraken_memory_mapping.sbatch
#
# --mem=750G is needed only for the no-mmap arm, which loads the whole 644 GiB hash into
# RAM. It restricts this job to the ~1008 GiB nodes (firefoot-*, gringolet-02).
#
# Four arms, on the 1M-read pair from the subsampling grid:
#   1. gtdb_mmap_cold  --memory-mapping, first touch of the DB in this job
#   2. gtdb_mmap_warm  --memory-mapping again, now that the page cache is populated
#   3. gtdb_ram        no flag: the full hash read into RAM
#   4. uhgg_ram        the current production setup, for scale
#
# CACHE CAVEAT, which matters when reading the numbers: page cache cannot be dropped
# without root, so the arms are not mutually independent. Arm 1 is the only arm that sees a
# genuinely cold-ish cache, and it is the one that predicts a fresh job on a fresh node -
# which is what each task of the subsampling array actually is. Arms 2 and 3 run against a
# cache that arm 1 warmed, so arm 3's load time is an OPTIMISTIC lower bound on what a cold
# no-mmap run would cost. Do not read arm 3 as "the RAM route is only this much faster".

set -uo pipefail

GTDB_DB='/cs/usr/netta.barak/lab/Tools/kraken2_db/gtdb_genome_reps_20250609'
UHGG_DB='/cs/usr/netta.barak/lab/Tools/kraken2_db/UnifiedHumanGastrointestinalGenome'
READS_DIR='/sci/labs/morani/morani/icore-data/lab/Data/mock_communities/atcc_1006_2006/short_reads'
READS_1="$READS_DIR/MSA-1006_1_R1_1M.fastq.gz"
READS_2="$READS_DIR/MSA-1006_1_R2_1M.fastq.gz"

OUT_DIR="kraken_mmap_bench_${SLURM_JOB_ID:-local}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

source /cs/usr/netta.barak/.bashrc
mamba activate ginger_env_312 || { echo "ERROR: could not activate ginger_env_312" >&2; exit 1; }

echo "=== node:    $(hostname)"
echo "=== threads: $THREADS"
echo "=== kraken2: $(which kraken2)"
kraken2 --version | head -1
echo "=== memory:"; free -g | head -2
echo

for f in "$READS_1" "$READS_2"; do
    [[ -f "$f" ]] || { echo "ERROR: missing reads $f" >&2; exit 1; }
done
for d in "$GTDB_DB" "$UHGG_DB"; do
    [[ -f "$d/hash.k2d" ]] || { echo "ERROR: no hash.k2d in $d" >&2; exit 1; }
done
echo "=== GTDB hash.k2d: $(du -h "$GTDB_DB/hash.k2d" | cut -f1)"
echo "=== UHGG hash.k2d: $(du -h "$UHGG_DB/hash.k2d" | cut -f1)"
echo

mkdir -p "$OUT_DIR" || exit 1
RESULTS="$OUT_DIR/results.tsv"
printf 'arm\tdb\tmemory_mapping\twall_seconds\twall_hms\tmax_rss_gib\tclassified_pct\n' > "$RESULTS"

# GInGeR's own kraken invocation (reference_database_utils.KRAKEN_COMMAND), so the numbers
# transfer directly to a real run rather than measuring a different command.
run_arm() {
    local arm="$1" db="$2" mmap_flag="$3"
    local out="$OUT_DIR/${arm}.output" report="$OUT_DIR/${arm}.report"
    local timing="$OUT_DIR/${arm}.time" log="$OUT_DIR/${arm}.log"

    echo "############################################################"
    echo "### arm: $arm   (memory-mapping: ${mmap_flag:-no})"
    echo "### started: $(date '+%F %T')"

    /usr/bin/time -v -o "$timing" \
        kraken2 --db "$db" $mmap_flag --paired "$READS_1" "$READS_2" \
                --threads "$THREADS" --output "$out" --report "$report" \
                --confidence 0.1 --use-names --report-minimizer-data \
        > "$log" 2>&1
    local rc=$?

    local wall_s wall_hms rss_gib pct
    # "Elapsed (wall clock) time (h:mm:ss or m:ss)" - normalise both shapes to seconds
    wall_hms=$(awk -F': ' '/Elapsed \(wall clock\)/{print $NF}' "$timing")
    wall_s=$(awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1 }' <<<"$wall_hms")
    rss_gib=$(awk '/Maximum resident set size/{printf "%.1f", $NF/1048576}' "$timing")
    # kraken2 prints "  N sequences classified (P%)" to stderr, captured in $log
    pct=$(grep -oE '\(([0-9.]+)%\)' "$log" | head -1 | tr -d '()')

    if (( rc != 0 )); then
        echo "### FAILED (exit $rc). kraken2 output:"
        cat "$log"
        printf '%s\t%s\t%s\tFAILED\tFAILED\t%s\t-\n' \
            "$arm" "$(basename "$db")" "${mmap_flag:-no}" "${rss_gib:--}" >> "$RESULTS"
        return 1
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$arm" "$(basename "$db")" "${mmap_flag:-no}" "$wall_s" "$wall_hms" "$rss_gib" "${pct:--}" >> "$RESULTS"
    echo "### wall: $wall_hms (${wall_s}s)   max RSS: ${rss_gib} GiB   classified: ${pct:--}"
    echo "### finished: $(date '+%F %T')"
    echo
}

# Order matters: the cold-cache mmap arm must run before anything else has touched the DB.
run_arm gtdb_mmap_cold "$GTDB_DB" --memory-mapping
run_arm gtdb_mmap_warm "$GTDB_DB" --memory-mapping
run_arm gtdb_ram       "$GTDB_DB" ""
run_arm uhgg_ram       "$UHGG_DB" ""

echo "############################################################"
echo "### RESULTS"
column -t -s$'\t' "$RESULTS"
echo
echo "Reminder: gtdb_ram ran against a cache warmed by the two mmap arms, so its wall time"
echo "is a lower bound. gtdb_mmap_cold is the number that predicts a real grid task."
echo
echo "Full logs and reports: $OUT_DIR"
