# Private release review

Version **0.2.0**, reviewed 2026-10-28. The repository remains private.
The eight findings from the 0.1.0 review have regression coverage and are
resolved in the reviewed implementation. See the [finding-by-finding closure
report](review-fixes.md). No blocking scientific correctness defect was found
in the follow-up review of the tested workflows. The limits below still apply.

Public visibility, hosted documentation and PyPI publication remain separate
owner-approved release steps. None was performed during this review.

## Verification

- **85 tests** pass locally, including independent exhaustive RDKit scores,
  780 fingerprint candidate-completeness checks, exact cutoff decisions
  immediately below/at/above 0.8 for both operators, stereo-aware nondefault
  fingerprints, aliases, invalid protected inputs, process preparation,
  task-local reads, compact indexes, pipeline incompatibility and interruptions.
- The [six-job CI matrix](https://github.com/Jozefov/MolecularSimilarityToolkit/actions/runs/40767367300) passed the
  suite, Python/CLI examples, integrity checks, strict documentation build,
  wheel/source builds and clean-wheel example on Python 3.11–3.13 for Linux
  x86-64 and macOS arm64. This is a private release workflow, not deployment.
- Ruff, strict MkDocs generation and package builds pass. A separate wheel-only
  environment ran the quickstart and process preparation from a plain Python
  command without a multiprocessing main guard.
- The generated package description renders with absolute links and PNG assets.
  The GitHub README was inspected at desktop 1,280 px and mobile 390 px, in
  light and dark themes. All responsive artwork loaded without page overflow.
- Full verification passed for the 1M index, search and filter outputs.
  A worst-case filter rejecting **all 1M real molecules** completed in
  8.81 s with 1639 MiB parent peak RSS;
  its output also passed full verification. Unchanged index resume took
  0.97 s, including input/artifact checksums.

## Real-molecule measurements

These are **nonrandom first-row subsets of a local MSnGym pretraining pool**,
snapshot `20260401_msngym_pretrain_pool_split.tsv`. They are not a whole-PubChem
benchmark. Each subset had the stated number of valid unique structures.
No input molecular collection is redistributed.

Machine: macOS 26.5.1 arm64, 10 logical CPUs, 32 GiB RAM; Python 3.12.9,
RDKit 2026.9.1, FPSim2 0.7.4. Default stereo-free Morgan radius 2 / 2,048 bits,
strict `> 0.8` exclusion, top-k 20 and 4 GiB configured memory budget.
The ten-query cases use ten protected molecules. The larger case uses
**1,000 queries and 1,000 protected molecules over four 250K-molecule shards**.

| Run | Maximum shard size | Build | Complete search | Complete filter | Build / search-filter peak tree RSS |
|---|---:|---:|---:|---:|---:|
| real-100k-1cpu | 1,000,000 | 20.02 s | 0.58 s | 0.63 s | 501 / 429 MiB |
| real-100k-4cpu | 1,000,000 | 8.56 s | 0.54 s | 0.63 s | 628 / 430 MiB |
| real-1m | 1,000,000 | 92.97 s | 1.49 s | 3.43 s | 1467 / 1712 MiB |
| real-1m-1000queries | 250,000 | 88.67 s | 7.59 s | 5.35 s | 1486 / 1055 MiB |

On the 100K subset, four requested preparation workers used four child PIDs.
Chemistry CPU time was 17.14 s over
5.39 s of staging wall time, demonstrating
actual parallel preparation. The one-worker staging time was
16.72 s. Whole-build speedup in this pair was
2.34×, not a promise of
linear CPU scaling.

Every real-data correctness flag passed. The ten-query cases compared full
query coverage, ranked molecule identities, integer intersection/union counts
and rejection identities with both packed NumPy and direct FPSim2 plus RDKit
rescoring. The 1,000-query case checked all 20,000 ranked matches and all 2,002
rejected molecule IDs against the independent direct path; the slower NumPy
baseline was explicitly omitted for that larger-query run.

Direct FPSim2 kernel timings omit preparation, boundary expansion, exact
rescoring, validation, checkpoints and output. Package timings include those
services. NumPy's filter mask reuses search similarities. These fields measure
different amounts of work and are labeled in the reports.

## Synthetic 10M fingerprint stress test

Ten 1M-row shards of uniform random packed fingerprints passed the actual
FPSim2 adapter versus a chunked NumPy score oracle. Fingerprint-file build time
was 293.42 s; loading/searching all ten shards and
rescoring candidates for one query took 0.87 s.
Parent peak RSS was 1235 MiB.

These are **synthetic fingerprints, not real molecules**. The measurement
excludes chemistry, provenance catalogs and complete 10M-row package scheduling.
It supports bounded shard-resident search, not a full 10M-molecule workflow claim.

## Interpretation and unresolved scope

The cache was uncontrolled and freshly built files were likely cached. The
measurements ran sequentially on a development workstation. Small local QA
commands overlapped the synthetic run; this was not an isolated performance
host. Process-tree RSS was sampled every 50 ms and includes preparation
subprocesses. It may count shared pages more than once and miss short peaks.
Parent OS high-water marks are retained separately. Memory budgets are planning
limits, not hard OS enforcement.

One preliminary direct-backend phase using four internal FPSim2 workers
segfaulted; an identical diagnostic rerun passed. The root cause is unknown.
The package uses serial native kernels with outer query concurrency. The final
benchmark explicitly uses the same native mode, has a regression guarding it,
and was rerun from fresh directories. The failed preliminary run is excluded
from the table above; see [performance notes](performance.md).

Whole-PubChem/ZINC throughput, real 10M+ molecule preparation, very large
reference collections and production PBS/Slurm network filesystems require
separate measurements. Local and CI cluster-task tests exercise the same work
units, but do not establish a scheduler deployment benchmark. No fixed full-corpus
runtime or linear CPU scaling is claimed.

**Migration:** 0.1 indexes and runs must be rebuilt in new directories. Scientific
pipeline versions prevent silent reuse. Portable compact indexes preserve
provenance and runtime receipts without staging files; the original build remains
available for repair.

Machine-readable evidence (input/code hashes, versions, settings, timing, memory
and complete correctness flags): [real-100k-1cpu](benchmarks/v020-real-100k-1cpu.json), [real-100k-4cpu](benchmarks/v020-real-100k-4cpu.json), [real-1m](benchmarks/v020-real-1m.json), [real-1m-1000queries](benchmarks/v020-real-1m-1000queries.json), [synthetic-10m](benchmarks/v020-synthetic-10m.json),
[full verification, worst-case filtering and resume](benchmarks/v020-extra-checks.json).
Reports from 0.1 remain historical artifacts; they are superseded by these 0.2
measurements. Wheel and source artifacts are built while private. No PyPI upload,
public documentation deployment or GitHub visibility change has occurred.
