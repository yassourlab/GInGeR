# GInGeR

GInGeR (Gene-in-Genomic-Region) is a Python pipeline that finds the genomic context of
genes of interest in metagenomic samples (Illumina paired-end reads, optionally + Oxford
Nanopore long reads), given the genes as protein (amino acid) FASTA sequences. It builds
a metagenomic assembly graph, locates the genes in the graph, extracts candidate flanking
("in"/"out") contexts, verifies them against a reference genome database, and assigns
matches to species.

CLI entry point: `run_ginger <SHORT_READS_1> <SHORT_READS_2> <GENES_PATH> <OUT_DIR>`
(registered in [pyproject.toml](pyproject.toml) as `ginger.ginger_runner:run_ginger_e2e`).
See [README.md](README.md) for full CLI options and output format documentation.

## Pipeline stages and key modules

1. **Reference filtering** — [ginger/reference_database_utils.py](ginger/reference_database_utils.py):
   runs Kraken2/Bracken to detect dominant species and estimate coverage, then downloads/builds a
   sample-specific UHGG reference subset (`merged_filtered_ref_db.fasta`).
2. **Assembly** — [ginger/assembly_utils.py](ginger/assembly_utils.py): runs MetaSPAdes (or
   HybridSPAdes if long reads are given) to produce `contigs.fasta`, `contigs.paths`, and
   `assembly_graph.fastg`.
3. **Gene locating** — [ginger/locating_genes_in_graph.py](ginger/locating_genes_in_graph.py):
   MMseqs2 maps genes -> contigs, then contigs are mapped to assembly-graph nodes (via
   `contigs.paths`, falling back to minimap2 for contigs with gaps) to locate genes within graph nodes.
4. **Context extraction** — [ginger/extract_contexts_candidates.py](ginger/extract_contexts_candidates.py):
   bidirectional DFS over the assembly graph (`pyfastg`/`networkx`) from each gene location to
   enumerate "in" (upstream) and "out" (downstream) path candidates, bounded by `--depth-limit`
   and `--max-context-len`. Writes `all_in_paths.fasta` / `all_out_paths.fasta`.
5. **Context verification** — [ginger/verify_context_candidates.py](ginger/verify_context_candidates.py):
   minimap2-aligns in/out paths to the reference DB (PAF output), then pairs in/out matches per
   gene+reference-genome subject to gap-ratio constraints (`--max-gap-ratio`).
6. **Aggregation** — [ginger/pipeline_utils.py](ginger/pipeline_utils.py): writes
   `context_level_matches.csv`, aggregates to `species_level_matches.csv` (and
   `subspecies_level_matches.csv` when subspecies metadata is available), and writes
   `genes_detected_in_graph_with_no_species_match.csv` for genes located in the graph but with
   no downstream match.

Orchestration of all stages happens in [ginger/ginger_runner.py](ginger/ginger_runner.py)
(`ginger_e2e_func`), which also handles intermediate-file cleanup (`--keep-intermediate`).

Other modules:
- [ginger/matches_classes.py](ginger/matches_classes.py) — core match classes:
  `GeneContigMatch` (gene->contig, from MMseqs2 `.m8`), `PathRefGenomeMatch` (path->reference
  genome, from minimap2 PAF), `InOutPathsMatch` (paired in/out context match -> CSV row).
- [ginger/sequence_alignment_utils.py](ginger/sequence_alignment_utils.py) — wrappers for
  running/parsing minimap2 and MMseqs2, plus NMS (non-max-suppression) for overlapping matches.
- [ginger/constants.py](ginger/constants.py) — centralized path templates for intermediate/output
  files. **Don't hardcode paths** — add/use templates here.
- [ginger/run_ginger_from_python.py](ginger/run_ginger_from_python.py) — example/profiling
  script for invoking the pipeline directly (with cProfile) instead of via the CLI.

## Conventions and gotchas

- Assembly graph nodes are identified with explicit orientation suffixes `+`/`-` (e.g. `123+`,
  `123-`); this must be preserved through all traversal and path-stitching logic.
- `.fastg` is the assembly graph format (parsed via `pyfastg`); MMseqs2 produces `.m8`; minimap2
  produces PAF — each has dedicated parsing in `sequence_alignment_utils.py` /
  `matches_classes.py`.
- Sequence overlaps between adjacent graph nodes are resolved via k-mer overlap detection
  (`pipeline_utils.get_sequence_overlap`, k in `{55, 33, 21, 43}`) — get this wrong and stitched
  context sequences will be corrupted.
- Empty results at any stage (no genes found, no in/out matches, etc.) are valid early-exit
  conditions, not crashes — `ginger_e2e_func` logs and returns early in these cases.
- External tools invoked via subprocess: `kraken2`, `bracken`, `seqkit`, `spades.py`, `mmseqs`,
  `minimap2`. Failures raise exceptions with the tool's stderr.

## Environment and testing

- Use the existing conda env `ginger_env_312` for running tests/tools; don't create new virtual
  environments (venv/poetry/pipenv) unless explicitly requested.
- **Don't run the test suite yourself** — the user runs tests (`python -m unittest discover -s
  tests`, matching CI in [.github/workflows/build-and-test.yml](.github/workflows/build-and-test.yml))
  and shares results for you to interpret/debug.
- Tests live in [tests/](tests/), one file per module roughly mirroring `ginger/`.
