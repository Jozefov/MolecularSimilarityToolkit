# Changelog

## 0.1.0 — unreleased

Initial private release candidate: reusable molecular fingerprint indexes,
exact ranked candidate search, protected training pools, resumable shard jobs,
and documented Python and command-line interfaces.

Verified with 49 tests across Python 3.11–3.13 on Linux and macOS. Private
release evidence includes 100K/1M real-molecule subsets and a separate 10M
synthetic fingerprint stress test. README artwork includes real molecular
examples, mobile layouts, and light/dark variants.

## 0.2.0 — private review candidate

- Fail closed on malformed protected SDF records during fold selection; preserve
  original SDF identifier/property strings and reject missing explicit ID mappings.
- Remove ignored annotations before normalization; strip atom maps by default
  with an explicit retention option. Correct isotope/tautomer interactions.
- Version scientific pipelines and reject legacy caches. Reuse completed staging
  without input reparsing; expose preparation progress and compact index export.
- Use bounded process workers for chemistry. Isolate build partitions and search
  task inputs; share a complete dependency contract across merge and verify.
- Strengthen benchmark validation to include every query, exact ranked identities
  and rejection sets, and measure process-tree memory and larger multi-shard runs.
