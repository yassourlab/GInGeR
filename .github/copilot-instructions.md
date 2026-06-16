---
name: ginger-repository-overview
description: "WORKFLOW SKILL - Explain how the GInGeR analysis repository works, including pipeline flow, key modules, conventions, and common pitfalls. Use when onboarding, planning changes, or debugging experiment behavior."
---

# GInGeR Repository Overview Skill

## Purpose
Provide a consistent way to explain the GInGeR ecosystem end to end: what the pipeline does, where each stage lives, what artifacts are produced, and which conventions matter when changing code.

## When To Use
- Onboarding to GInGeR analysis and experiments code
- Tracing failures or low precision/recall across pipeline stages
- Planning edits in assembly, mapping, context extraction, and metrics paths
- Reviewing whether a change is consistent with repository conventions

## Inputs
- Target scenario or question (onboarding, debugging, architecture review, change planning)
- Optional focus area (assembly, ground truth, baselines, MLflow, metrics)
- Optional run artifacts or paths (for concrete debugging walkthroughs)

## Workflow
1. State the high-level objective and expected outputs.
2. Map the request to pipeline stage(s) and owning module(s).
3. Explain data flow, intermediate artifacts, and key parameters.
4. Highlight decision points and assumptions for the selected stage(s).
5. Run quality checks against conventions and known gotchas.
6. Produce a concise summary with next debugging or implementation steps.

## Pipeline Summary
1. Reference filtering: Kraken2/Bracken species filtering and UHGG subset creation.
2. Assembly: MetaSPAdes or HybridSPAdes graph generation.
3. Gene locating: MMseqs2 mapping from genes to contigs to graph nodes.
4. Context extraction: bidirectional DFS over graph for in/out candidate contexts.
5. Context verification: Minimap2 alignment and in/out pairing by gap constraints.
6. Aggregation: context-level and species-level outputs.

## Key Modules
- `ginger_runner.py`: top-level orchestration.
- `reference_database_utils.py`: sample-specific reference DB generation.
- `assembly_utils.py`: assembly execution and graph-related prep.
- `locating_genes_in_graph.py`: gene-to-node localization.
- `extract_contexts_candidates.py`: context candidate traversal logic.
- `verify_context_candidates.py`: alignment-based verification and pairing.
- `pipeline_utils.py`: aggregation and pipeline utilities.
- `matches_classes.py`: core match classes (`GeneContigMatch`, `PathRefGenomeMatch`, `InOutPathsMatch`).

## Data And Conventions
- Assembly graph uses `.fastg`, parsed to directed graph representation.
- Node orientation is explicit (`+` and `-`) and must be preserved through traversal.
- Paths are bounded by depth and context-length controls.
- Intermediate path templates are centralized in `constants.py`; do not hardcode paths.
- MMseqs2 outputs `.m8`; Minimap2 outputs PAF; parsing is stage-specific.

## Decision Points
- Whether to run full species filtering or use a pre-built filtered reference.
- Whether to run assembly or reuse existing SPAdes artifacts with skip flags.
- Which thresholds to tune first (NMS IoU, context lengths, gap ratio).
- Whether low recall is from missing genes in assembly, weak alignment, or pairing filters.

## Quality Checks
- Verify all expected SPAdes outputs exist before downstream steps.
- Validate node naming/orientation conversions before graph traversal conclusions.
- Confirm metadata/reference compatibility for UHGG IDs and headers.
- Check that output CSVs are generated and split/aggregated as expected.
- Separate expected empty-result behavior from real pipeline failures.

## Common Pitfalls
- Metadata path mismatches between CI and local runs.
- Insufficient memory for Kraken2-based filtering.
- Orientation errors when stitching graph paths.
- Incorrect sequence overlap handling when concatenating nodes.
- Misinterpreting early empty-result exits as crashes.

## Output Format
Return answers in this order:
1. High-level flow (1 short paragraph).
2. Stage-by-stage mapping (numbered list).
3. Decision points and failure modes (bullets).
4. Concrete next steps (numbered list).

## Notes
- Keep explanations grounded in actual repository modules and artifacts.
- Prefer short, practical guidance over theoretical bioinformatics detail.
- Environment: Use the user's existing conda environments; never create new virtual environments (venv/virtualenv/poetry/pipenv) unless explicitly requested.If you want to run a tests using a conda env, you can use ginger_env_312
- Tests: don't run the tests yourself. I'll run them and share the results with you. You can then help me interpret the results and debug any failures.