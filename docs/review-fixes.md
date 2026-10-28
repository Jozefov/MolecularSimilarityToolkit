# Review findings and regression coverage

Version 0.2.0 addresses the eight findings from the private 0.1.0 review.
The repository remains private. Visibility, hosted documentation and PyPI
publication are separate release decisions.

## Eight findings

| Finding | Correction | Regression evidence |
|---|---|---|
| Malformed protected SDF rows disappeared during fold selection | Reader failures reach quarantine before fold selection; unknown fold membership is invalid; protection fails by default | `test_invalid_selected_sdf_is_never_silently_dropped`, including gzip and malformed first/last records; `test_unknown_fold_is_invalid_and_extra_receipts_are_rejected` |
| Isotope removal happened after tautomer selection | Remove ignored isotope, map and stereo annotations before canonical choices | `test_ignored_isotopes_do_not_influence_tautomer_choice`, including interactions with parent normalization and end-to-end filtering |
| Every search task read the global query file; build tasks read the global catalog | Sealed worker headers, individual task specifications, checked query chunks and prepared build partitions | `test_global_queries_not_read_by_independent_search_tasks`, `test_build_workers_read_only_their_partition`; 1M/four-shard/1,000-query benchmark |
| Chemistry preparation only used one CPU | Persistent bounded process batches through Loky; serial small-input path; cached fingerprint generator | `test_preparation_uses_bounded_process_workers`, `test_process_preparation_works_from_stdin_without_main_guard`; measured one/four-worker builds |
| SDF IDs and metadata changed type | Read original property strings without numeric conversion | `test_original_sdf_ids_remain_exact_strings` |
| Toolkit scientific changes did not invalidate caches | Separate preparation, operation and format versions in descriptors and compatibility checks | `test_pipeline_versions_invalidate_legacy_and_operation_caches`; legacy 0.1 data requires rebuilding in a new directory |
| Merge and verify accepted different dependency sets | Explicit required-unit manifests and shared validation, including zero-task query chunks and unexpected receipts | `test_merge_and_verify_require_every_prepared_query_dependency`, `test_unknown_fold_is_invalid_and_extra_receipts_are_rejected` |
| Benchmark flags could accept missing queries or swapped identities | Compare expected query coverage, ordered molecule IDs, integer counts and complete rejection sets | `test_benchmark_flags_reject_missing_queries_and_wrong_identities` |

The 780 random and constructed fingerprint candidate-completeness checks are
now permanent regression coverage in `test_random_and_constructed_fingerprint_candidates`.
Small-dataset scientific tests independently use exhaustive RDKit similarity.

## Additional release details

- Explicit ID mappings are validated, including header-only inputs. Automatic
  row IDs use a separate documented `id_column=None` mode.
- Atom maps are removed by default. Retaining maps is an explicit chemistry
  option; compatible settings remain mandatory on both sides.
- Completed staging reuses checked chemistry chunks without reparsing input.
  A damaged completion marker replays the reader while reusing valid chunks.
- Progress includes processed rows and completed staging chunks.
- `export_index` / `mol-sim compact` publishes a portable runtime index while
  preserving the original build. A regression deletes the original build and
  verifies search and filtering through the compact copy.
- `Result.open` and `mol-sim export` reopen/export disk-backed results.
- Contributor setup installs the documentation extra. The full API is available
  through a local MkDocs build and the private CI artifact.
- The Git installation URL uses HTTPS. A generated package description uses
  absolute links and PNG assets; a rendering regression checks it outside GitHub.
- The larger-query benchmark covers every match and rejection in a real
  four-shard workflow. Reports record code/input hashes and process-tree RSS.

## Second review

The follow-up review rechecked preparation order, invalid protected inputs,
identity and aliases, local/global ranking, cutoff arithmetic, worker input
validation, dependency completeness, resume publication and distribution files.
It added a regression for damaged staging markers and for the serial native
kernel mode used by the adapter and benchmark.

A preliminary direct-backend benchmark using FPSim2's internal four-worker
mode segfaulted once and passed on diagnostic rerun. Its cause remains unknown.
The package uses serial native kernels with bounded outer query concurrency;
the benchmark now explicitly matches that mode. The failed preliminary run is
not used as release performance evidence. See [performance notes](performance.md).

The completed checks and measured limits are in the [release review](release.md).
They cover the tested local workflows and CI platforms. They do not establish
whole-PubChem performance or production PBS/Slurm network-filesystem behavior.
