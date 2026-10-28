# Private release review

Version **0.1.0**, reviewed 2026-10-28. The repository remains private.
Public visibility, hosted documentation and PyPI publication are separate
owner-approved release steps.

## Verification

- 49 automated tests cover independent exhaustive RDKit agreement, exact 0.8
  boundaries, cross-shard ties, repeated IDs, identity exclusion, chemistry
  policies, input formats, malformed plans, corrupt/missing receipts and
  interruptions during staging, shard creation, search and reduction.
- Python 3.11, 3.12 and 3.13 passed on Linux x86-64 and macOS arm64 in the
  [six-job CI matrix](https://github.com/Jozefov/MolecularSimilarityToolkit/actions/runs/40752464603).
- Ruff, the strict documentation build, source distribution and wheel build pass.
  The wheel installs in a separate environment and runs the complete example.
- Python and CLI quickstarts run in CI. The example retains both ethanol source
  records, deduplicates seven input rows to six structures, and keeps four
  structures when the two validation/test molecules are protected.
- README artwork was rendered at 1,280-pixel desktop and 390-pixel mobile widths
  in light and dark themes; no page-width overflow or missing images was found.
  The infographic uses actual molecule drawings and calculated example scores.

## Real-molecule measurements

Input: the first 100K and first 1M rows of the user's local MSnGym pretraining
pool snapshot `20260401_msngym_pretrain_pool_split.tsv`. These are **nonrandom
subsets of that file**, not a whole-PubChem benchmark. Each subset contained the
stated number of unique valid canonical structures and no invalid rows.

Machine: macOS 26.5.1, arm64, 10 logical CPUs, 32 GiB RAM. Python 3.12.9,
RDKit 2026.9.1, FPSim2 0.7.4. Four requested workers, 4 GiB configured budget,
up to 1M molecules per shard. Each test used **10 query/protected molecules**,
selected in stable-ID order, top-k 20, and a strict `> 0.8` cutoff. Queries can
differ between subsets. The OS cache was uncontrolled; this was a development
workstation, not an isolated performance host.

| Measurement | 100K molecules | 1M molecules |
|---|---:|---:|
| Index build | 22.64 s | 279.08 s |
| Complete package search, 10 queries | 0.65 s | 1.56 s |
| Complete package filter, 10 protected molecules | 0.82 s | 4.01 s |
| Package search/filter peak RSS | 426 MiB | 1,697 MiB |
| Build peak RSS | 496 MiB | 1,481 MiB |
| Direct FPSim2 load | 0.065 s | 0.553 s |
| Direct FPSim2 search kernel, 10 queries | 0.009 s | 0.044 s |
| Direct FPSim2 filter kernel, 10 queries | 0.005 s | 0.020 s |
| Historical NumPy scan and top-k, 10 queries | 0.84 s | 9.09 s |
| Removed molecules | 22 | 24 |

The package's extra time relative to direct FPSim2 is expected: query index
preparation, input validation, boundary handling, exact rescoring, deterministic
witnesses, checkpoint receipts, provenance joins and output writing are included.
The filter writes almost the entire kept pool and its aliases. These are different
services from a raw search kernel. Full package search was faster than the
historical NumPy baseline in both runs. Exact scores and rejected counts agreed
with NumPy; direct backend scores agreed within their float32 precision.

Machine-readable reports include input subset checksums, resolved settings,
versions, per-phase RSS and correctness flags:
[100K report](benchmarks/real-100k.json), [1M report](benchmarks/real-1m.json).
The input collections are not committed or redistributed.

## Synthetic 10M-row stress test

Ten shards of one million **synthetic random fingerprints** each passed the
adapter-versus-NumPy exact-score checks. Peak RSS was **1,203 MiB**. Total
fingerprint-file build time was 355.48 seconds; loading/searching all ten shards
for one query and rescoring candidates took 1.11 seconds. NumPy oracle scanning
and random-data generation are outside those timing fields.

These are uniform random bytes, approximately 50% on-bits, not real compounds.
The run validates bounded shard storage, the FPSim2 adapter and exact rescoring;
it does not measure molecular normalization, aliases or the complete 10M-row
package workflow. Local documentation generation and small checks overlapped
this stress run, so it is not an isolated speed measurement.
[Full synthetic report](benchmarks/synthetic-10m.json).

## Practical limits and release boundary

The evidence covers up to 1M real molecules through the full workflow and 10M
synthetic fingerprints through the search adapter. Whole-PubChem performance,
large protected sets and network-filesystem behavior need their own measurements.
No fixed runtime or linear CPU scaling is promised.

Standalone distributed build tasks validate the shared catalog independently,
so each can incur a catalog checksum read. Local execution reuses validated
context and resident shards safely. Avoid duplicate concurrent attempts of the
same task. Memory settings are planning limits, not a hard OS sandbox.

The repository includes the MIT license, citation metadata, contribution guide,
locked development environment, editable SVGs and PNG exports. Only small
hand-written molecule fixtures are bundled. The wheel and source archive are
built while private; there is no PyPI upload, public documentation deployment,
or visibility change in this release review.
