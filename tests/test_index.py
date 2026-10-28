import pytest

from molecular_similarity_toolkit.config import Source
from molecular_similarity_toolkit.index import build_index
from molecular_similarity_toolkit.storage import batches


def test_build_resume_aliases(tmp_path):
    source = Source.from_records(
        [
            {"id": "same", "smiles": "CCO", "x": 3},
            {"id": "same", "smiles": "OCC", "x": 4},
            {"smiles": "bad"},
        ],
        metadata_columns=("x",),
    )
    index = build_index(source, tmp_path / "index", shard_size=1, chunk_size=1)
    assert index.size == 1
    assert index.manifest["alias_count"] == 2
    assert index.manifest["invalid_count"] == 1
    assert len(list(batches(index.aliases_path))) == 2
    shard = index.path / "units/shard-00000000/fingerprints.h5"
    stamp = shard.stat().st_mtime_ns
    build_index(source, index.path, shard_size=1, chunk_size=1)
    assert shard.stat().st_mtime_ns == stamp
    with pytest.raises(ValueError, match="resume refused"):
        build_index(
            Source.from_records([{"smiles": "CC"}]),
            index.path,
            shard_size=1,
            chunk_size=1,
        )


def test_empty(tmp_path):
    index = build_index(Source.from_records([]), tmp_path / "empty")
    assert index.size == 0
    assert not index.manifest["shards"]


def test_interruption(tmp_path, monkeypatch):
    source = Source.from_records([{"smiles": "CCO"}, {"smiles": "CCC"}])
    monkeypatch.setenv("MST_FAILPOINT", "stage-0000-00000000:after_publish")
    with pytest.raises(RuntimeError, match="interruption"):
        build_index(source, tmp_path / "index", chunk_size=1, shard_size=1)
    monkeypatch.delenv("MST_FAILPOINT")
    index = build_index(source, tmp_path / "index", chunk_size=1, shard_size=1)
    assert index.size == 2
