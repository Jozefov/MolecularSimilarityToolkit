import pyarrow.parquet as pq
import pytest
from typer.testing import CliRunner

import molecular_similarity_toolkit as mst
from molecular_similarity_toolkit.cli import app
from molecular_similarity_toolkit.storage import atomic_json, read_json


def source():
    return mst.Source.from_records(
        [
            {"id": 0, "smiles": "CCO"},
            {"id": "x", "smiles": "CCN"},
            {"id": "x", "smiles": "CCC"},
        ]
    )


def rows(result, table=None):
    return [row for batch in result.batches(table) for row in batch.to_pylist()]


def test_build_plan_missing_and_retry(tmp_path):
    plan = mst.plan_build(source(), tmp_path / "index", shard_size=1, chunk_size=1)
    with pytest.raises(ValueError, match="incomplete"):
        mst.merge(plan)
    for task_id in [2, 0, 1]:
        mst.run_task(plan, task_id)
    index = mst.merge(plan)
    assert mst.verify(index.path)["ok"]
    assert mst.status(index.path)["state"] == "complete"
    aliases = pq.read_table(index.aliases_path).to_pylist()
    assert aliases[0]["original_id"] == "0"
    staged = index.path / "receipts/stage-0000-00000000.json"
    staged.unlink()
    assert not mst.verify(index.path)["ok"]
    with pytest.raises(ValueError, match="incomplete"):
        mst.MoleculeIndex.open(index.path)


@pytest.mark.parametrize(
    "point",
    [
        "task-00000000:before_publish",
        "task-00000000:after_publish",
        "merge:after_write",
        "merge:after_publish",
    ],
)
def test_operation_interruption_resume(tmp_path, monkeypatch, point):
    index = mst.build_index(source(), tmp_path / "index", shard_size=1)
    plan = mst.plan_search(index, source(), 2, tmp_path / "search", query_chunk_size=1)
    monkeypatch.setenv("MST_FAILPOINT", point)
    with pytest.raises(RuntimeError, match="interruption"):
        mst.execute(plan)
    monkeypatch.delenv("MST_FAILPOINT")
    result = mst.execute(plan)
    expected = mst.search(
        index, source(), 2, tmp_path / "expected", query_chunk_size=20, workers=1
    )
    assert rows(result) == rows(expected)
    assert mst.verify(result.path)["ok"]


def test_task_corruption_and_plan_change(tmp_path):
    index = mst.build_index(source(), tmp_path / "index", shard_size=2)
    plan = mst.plan_search(index, source(), 2, tmp_path / "search")
    result = mst.execute(plan)
    output = result.path / "units/task-00000000/matches.parquet"
    output.write_bytes(b"damaged")
    with pytest.raises(ValueError, match="corrupt"):
        mst.merge(plan)
    assert not mst.verify(result.path)["ok"]
    mst.run_task(plan, 0)
    assert mst.verify(mst.merge(plan).path)["ok"]
    value = read_json(plan)
    value["tasks"][0]["shard"] = "shard-00000001"
    atomic_json(plan, value)
    with pytest.raises(ValueError, match="manifest changed"):
        mst.merge(plan)


def test_cli_pipeline(tmp_path):
    runner = CliRunner()
    data = tmp_path / "molecules.csv"
    data.write_text("id,smiles\nfirst,CCO\nsecond,CCN\n")
    index = tmp_path / "index"
    for args in [
        ["build", str(data), str(index), "--shard-size", "1"],
        ["verify", str(index)],
        ["plan", "search", str(index), str(data), str(tmp_path / "search"), "--k", "1"],
    ]:
        result = runner.invoke(app, args)
        assert result.exit_code == 0, result.output
    plan = tmp_path / "search/plan.json"
    for i in range(len(read_json(plan)["tasks"])):
        result = runner.invoke(app, ["task", str(plan), str(i)])
        assert result.exit_code == 0, result.output
    result = runner.invoke(app, ["merge", str(plan)])
    assert result.exit_code == 0, result.output


@pytest.mark.parametrize(
    "point",
    [
        "stage-0000-00000000:before_publish",
        "shard-00000000:before_publish",
        "shard-00000000:after_publish",
    ],
)
def test_build_interruption_units(tmp_path, monkeypatch, point):
    monkeypatch.setenv("MST_FAILPOINT", point)
    with pytest.raises(RuntimeError, match="interruption"):
        mst.build_index(source(), tmp_path / "index", shard_size=1, chunk_size=1)
    monkeypatch.delenv("MST_FAILPOINT")
    result = mst.build_index(source(), tmp_path / "index", shard_size=1, chunk_size=1)
    assert result.size == 3
    assert mst.verify(result.path)["ok"]


def test_incompatible_indexes_and_tampered_manifest(tmp_path):
    pool = mst.build_index(source(), tmp_path / "pool")
    changed = mst.build_index(
        source(), tmp_path / "changed", fingerprint=mst.FingerprintConfig(radius=3)
    )
    with pytest.raises(ValueError, match="incompatible"):
        mst.search(pool, changed, out=tmp_path / "search")
    manifest = read_json(pool.path / "manifest.json")
    manifest["chemistry"]["remove_isotopes"] = True
    atomic_json(pool.path / "manifest.json", manifest)
    with pytest.raises(ValueError, match="manifest configuration"):
        mst.MoleculeIndex.open(pool.path)


def test_index_catalog_corruption_blocks_merge(tmp_path):
    pool = mst.build_index(source(), tmp_path / "pool")
    plan = mst.plan_filter(pool, source(), tmp_path / "filter")
    for task_id in range(len(read_json(plan)["tasks"])):
        mst.run_task(plan, task_id)
    pool.molecules_path.write_bytes(b"damaged")
    with pytest.raises(ValueError, match="catalog.*corrupt"):
        mst.merge(plan)


def test_execute_build_and_result_exports(tmp_path):
    plan = mst.plan_build(source(), tmp_path / "index", shard_size=2)
    index = mst.execute(plan)
    assert isinstance(index, mst.MoleculeIndex)
    result = mst.search(index, source(), 2, tmp_path / "search")
    exported = result.export(tmp_path / "matches.tsv")
    lines = exported.read_text().splitlines()
    assert "\t" in lines[0]
    assert len(lines) == result.summary["match_count"] + 1
    with pytest.raises(ValueError, match="unknown"):
        result.table_path("../../anything")
    with pytest.raises(ValueError, match="end in"):
        result.export(tmp_path / "matches.xlsx")
    with pytest.raises(ValueError, match="delimiter"):
        result.export(tmp_path / "matches.csv", delimiter="two")
