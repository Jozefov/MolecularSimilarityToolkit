"""Failures reproduced during the private pre-release review."""

import gzip
import importlib.util
import shutil
from pathlib import Path

import pytest
from rdkit import Chem

import molecular_similarity_toolkit as mst
from molecular_similarity_toolkit import index as index_module
from molecular_similarity_toolkit import operations, storage
from molecular_similarity_toolkit.resources import effective_workers
from molecular_similarity_toolkit.storage import atomic_json, batches, read_json, seal


def source():
    return mst.Source.from_records(
        [{"id": "a", "smiles": "CCO"}, {"id": "b", "smiles": "CCN"}]
    )


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("bad_first", [False, True])
def test_invalid_selected_sdf_is_never_silently_dropped(
    tmp_path, compressed, bad_first
):
    valid = (
        Chem.MolToMolBlock(Chem.MolFromSmiles("CCO")) + "\n>  <fold>\ntest\n\n$$$$\n"
    )
    bad = "broken\n  RDKit\n\nBAD COUNTS\n>  <fold>\ntest\n\n$$$$\n"
    path = tmp_path / ("ref.sdf.gz" if compressed else "ref.sdf")
    text = bad + valid if bad_first else valid + bad
    if compressed:
        with gzip.open(path, "wt") as f:
            f.write(text)
    else:
        path.write_text(text)
    protected = mst.Source(path, split_column="fold", splits=["test"])
    index = mst.build_index(protected, tmp_path / "ref")
    invalid = list(batches(index.invalid_path))
    assert len(invalid) == 1 and invalid[0]["error"] == "invalid SDF record"
    assert invalid[0]["source_row"] == (1 if bad_first else 2)
    pool = mst.build_index(source(), tmp_path / "pool")
    with pytest.raises(ValueError, match="invalid molecules"):
        mst.filter_pool(pool, protected, out=tmp_path / "filter")
    allowed = mst.filter_pool(
        pool, protected, out=tmp_path / "explicit", invalid_protected="quarantine"
    )
    assert allowed.summary["protected_invalid_count"] == 1


@pytest.mark.parametrize("parent", [False, True])
@pytest.mark.parametrize("remove_stereo", [False, True])
def test_ignored_isotopes_do_not_influence_tautomer_choice(
    tmp_path, parent, remove_stereo
):
    config = mst.ChemistryConfig(
        remove_isotopes=True,
        canonical_tautomer=True,
        normalized_parent=parent,
        remove_stereochemistry=remove_stereo,
    )
    pool = mst.build_index(
        mst.Source.from_records([{"smiles": "C[15NH]C(=N)N"}, {"smiles": "CNC(=N)N"}]),
        tmp_path / "pool",
        chemistry=config,
    )
    assert pool.size == 1
    protected = mst.build_index(
        mst.Source.from_records([{"smiles": "CNC(=N)N"}]),
        tmp_path / "protected",
        chemistry=config,
    )
    result = mst.filter_pool(pool, protected, out=tmp_path / "result")
    assert result.summary["rejected_count"] == 1
    assert next(result.batches("rejected")).to_pylist()[0]["similarity"] == 1


def test_original_sdf_ids_remain_exact_strings(tmp_path):
    path = tmp_path / "ids.sdf"
    m = Chem.MolFromSmiles("CCO")
    for key, value in [("id", "00123"), ("catalog", "00078"), ("exponent", "1e5")]:
        m.SetProp(key, value)
    with Chem.SDWriter(str(path)) as writer:
        writer.write(m)
    index = mst.build_index(
        mst.Source(path, id_column="id", metadata_columns=["catalog", "exponent"]),
        tmp_path / "index",
    )
    row = next(batches(index.aliases_path))
    assert row["original_id"] == "00123"
    assert __import__("json").loads(row["metadata"]) == {
        "catalog": "00078",
        "exponent": "1e5",
    }


def test_id_mapping_and_explicit_atom_map_policy(tmp_path):
    with pytest.raises(ValueError, match="missing mapped columns"):
        mst.build_index(
            mst.Source.from_records(
                [{"cid": "001", "smiles": "CCO"}], id_column="typo"
            ),
            tmp_path / "bad",
        )
    rows = [{"smiles": "CCO"}, {"smiles": "[CH3:1][CH2:2][OH:3]"}]
    default = mst.build_index(mst.Source.from_records(rows), tmp_path / "default")
    assert default.size == 1
    assert [r["original_id"] for r in batches(default.aliases_path)] == [
        "row-1",
        "row-2",
    ]
    hits = mst.search(
        default,
        mst.Source.from_records([{"smiles": "CCO"}]),
        exclude_identity=True,
        out=tmp_path / "excluded",
    )
    assert hits.summary["match_count"] == 0
    retained = mst.build_index(
        mst.Source.from_records(rows),
        tmp_path / "retained",
        chemistry=mst.ChemistryConfig(remove_atom_maps=False),
    )
    assert retained.size == 2


def test_completed_staging_does_not_reparse_input(tmp_path, monkeypatch):
    index = mst.build_index(source(), tmp_path / "index")

    def no_parse(*a, **kw):
        raise AssertionError("completed input replayed")

    monkeypatch.setattr(index_module, "selected_rows", no_parse)
    monkeypatch.setattr(
        index_module, "__version__", "99.0.0"
    )  # Audit-only version changes do not invalidate science.
    assert mst.build_index(source(), index.path).size == 2


def test_pipeline_versions_invalidate_legacy_and_operation_caches(
    tmp_path, monkeypatch
):
    index = mst.build_index(source(), tmp_path / "index")
    plan = mst.plan_search(index, source(), 1, tmp_path / "run")
    monkeypatch.setattr(operations, "OPERATION_VERSION", 99)
    with pytest.raises(ValueError, match="unsupported operation"):
        mst.run_task(plan, 0)
    manifest = read_json(index.path / "manifest.json")
    manifest["preparation_version"] = 1
    atomic_json(index.path / "manifest.json", seal(manifest))
    with pytest.raises(ValueError, match="rebuild legacy"):
        mst.MoleculeIndex.open(index.path)


def test_global_queries_not_read_by_independent_search_tasks(tmp_path, monkeypatch):
    pool = mst.build_index(source(), tmp_path / "pool", shard_size=1)
    queries = mst.Source.from_records(
        [{"id": str(i), "smiles": "CCO"} for i in range(12)]
    )
    plan = mst.plan_search(pool, queries, 1, tmp_path / "run", query_chunk_size=2)
    reads = []
    original = storage.digest

    def digest(path):
        reads.append(str(path))
        return original(path)

    monkeypatch.setattr(storage, "digest", digest)
    for task in read_json(plan)["tasks"]:
        mst.run_task(plan, task["id"])
    assert str(plan.parent / "units/queries/queries.parquet") not in reads
    assert str(pool.molecules_path) not in reads
    mst.merge(plan)
    assert reads.count(str(plan.parent / "units/queries/queries.parquet")) == 1

    def no_plan(*a, **kw):
        raise AssertionError("worker read complete plan")

    monkeypatch.setattr(operations, "_load_plan", no_plan)
    mst.run_task(plan, 0)


def test_build_workers_read_only_their_partition(tmp_path, monkeypatch):
    plan = mst.plan_build(source(), tmp_path / "index", shard_size=1)
    original = storage.digest
    reads = []

    def digest(path):
        reads.append(str(path))
        return original(path)

    monkeypatch.setattr(storage, "digest", digest)
    mst.run_task(plan, 0)
    assert str(plan.parent / "units/catalog/molecules.parquet") not in reads
    assert any("build-input-00000000" in p for p in reads)
    assert not any("build-input-00000001" in p for p in reads)


@pytest.mark.parametrize("empty", [False, True])
def test_merge_and_verify_require_every_prepared_query_dependency(tmp_path, empty):
    pool = mst.build_index(
        mst.Source.from_records([]) if empty else source(), tmp_path / "pool"
    )
    plan = mst.plan_search(pool, source(), 1, tmp_path / "run", query_chunk_size=1)
    mst.execute(plan)
    dependency = plan.parent / "units/query-chunk-00000000/queries.parquet"
    dependency.write_bytes(b"damaged")
    with pytest.raises(ValueError, match="corrupt"):
        mst.merge(plan)
    assert not mst.verify(plan.parent)["ok"]


def test_compact_index_is_self_contained(tmp_path):
    pool = mst.build_index(source(), tmp_path / "pool", shard_size=1)
    compact = mst.export_index(pool, tmp_path / "compact")
    assert mst.verify(compact.path)["ok"] and compact.manifest["compact"]
    assert not (compact.path / "build-plan.json").exists()
    shutil.rmtree(pool.path)
    result = mst.search(compact, source(), 1, tmp_path / "result")
    assert result.summary["match_count"] == 2
    reopened = mst.Result.open(result.path)
    assert reopened.summary == result.summary
    with pytest.raises(ValueError, match="immutable"):
        mst.build_index(source(), compact.path, shard_size=1)


def test_staging_status_exposes_completed_chunks(tmp_path, monkeypatch):
    monkeypatch.setenv("MST_FAILPOINT", "stage-0000-00000001:before_publish")
    with pytest.raises(RuntimeError):
        mst.build_index(source(), tmp_path / "index", chunk_size=1)
    value = mst.status(tmp_path / "index")
    assert (
        value["state"] == "staging"
        and value["staged_chunks"] == 1
        and value["selected_rows"] == 1
    )


def test_preparation_uses_bounded_process_workers(tmp_path):
    if effective_workers(2, 4) < 2:
        pytest.skip("one CPU allocation")
    records = mst.Source.from_records([{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}] * 16384)
    index = mst.build_index(records, tmp_path / "index", workers=2)
    progress = read_json(index.path / "progress.json")
    assert len(progress["worker_pids"]) == 2
    assert progress["preparation_cpu_seconds"] > 0
    assert index.size == 1 and index.manifest["alias_count"] == 16384


def test_benchmark_flags_reject_missing_queries_and_wrong_identities(tmp_path):
    path = Path(__file__).parents[1] / "benchmarks/benchmark.py"
    spec = importlib.util.spec_from_file_location("mst_benchmark", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    (tmp_path / "queries.tsv").write_text("id\tsmiles\n0\tCCO\n1\tCCN\n")
    module._write_validation(tmp_path, "package", {"0": [("a", 1, 1)]}, {"a"})
    module._write_validation(
        tmp_path, "direct", {"0": [("a", 1, 1)], "1": [("b", 1, 1)]}, {"b"}
    )
    checks = module.compare_results(tmp_path, ["package", "direct"], 1)
    assert not checks["package_query_coverage"]
    assert not checks["direct_matches_exact_identities"]
    assert not checks["direct_rejected_exact_identities"]


def test_changed_coordinator_plan_cannot_run_stale_worker_spec(tmp_path):
    index = mst.build_index(source(), tmp_path / "index")
    plan = mst.plan_search(index, source(), 1, tmp_path / "run")
    value = read_json(plan)
    value["k"] = 999
    atomic_json(plan, value)
    with pytest.raises(ValueError, match="plan changed"):
        mst.run_task(plan, 0)


def test_explicit_mapping_validated_for_header_only_sources(tmp_path):
    file = tmp_path / "empty.csv"
    file.write_text("smiles,id\n")
    with pytest.raises(ValueError, match="missing mapped columns"):
        mst.build_index(mst.Source(file, id_column="typo"), tmp_path / "bad")
    assert mst.build_index(mst.Source(file), tmp_path / "empty").size == 0


def test_process_preparation_works_from_stdin_without_main_guard(tmp_path):
    import subprocess
    import sys

    code = f"import molecular_similarity_toolkit as m; p=m.build_index(m.Source.from_records([{{'smiles':'CCO'}}]*2048),{str(tmp_path / 'stdin')!r},workers=2); assert p.manifest['alias_count']==2048"
    subprocess.run([sys.executable, "-c", code], check=True, timeout=60)


def test_unknown_fold_is_invalid_and_extra_receipts_are_rejected(tmp_path):
    pool = mst.build_index(source(), tmp_path / "pool")
    protected = mst.Source.from_records(
        [{"smiles": "CCO", "fold": "test"}, {"smiles": "CCN", "fold": ""}],
        split_column="fold",
        splits=["test"],
    )
    with pytest.raises(ValueError, match="invalid molecules"):
        mst.filter_pool(pool, protected, out=tmp_path / "bad")
    result = mst.search(pool, source(), 1, tmp_path / "run")
    shutil.copyfile(
        result.path / "receipts/task-00000000.json",
        result.path / "receipts/task-99999999.json",
    )
    assert not mst.verify(result.path)["ok"]
    with pytest.raises(ValueError, match="unexpected"):
        mst.merge(result.path / "plan.json")


def test_numeric_fold_selection_and_invalid_sequence_config(tmp_path):
    index = mst.build_index(
        mst.Source.from_records(
            [{"smiles": "CCO", "fold": 1}], split_column="fold", splits=[1]
        ),
        tmp_path / "index",
    )
    assert index.size == 1
    with pytest.raises(TypeError, match="sequences"):
        mst.Source.from_records([], split_column="fold", splits="test")
