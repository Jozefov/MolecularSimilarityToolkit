import gzip

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from rdkit import Chem, DataStructs

import molecular_similarity_toolkit as mst
from molecular_similarity_toolkit import chemistry
from molecular_similarity_toolkit.adapter import write_shard
from molecular_similarity_toolkit.readers import selected_rows
from molecular_similarity_toolkit.storage import batches


@pytest.mark.parametrize(
    "suffix",
    ["csv", "tsv", "csv.gz", "tsv.gz", "smi", "smi.gz", "parquet", "sdf", "sdf.gz"],
)
def test_formats(tmp_path, suffix):
    path = tmp_path / f"data.{suffix}"
    if suffix == "parquet":
        pq.write_table(pa.Table.from_pylist([{"id": "zero", "smiles": "CCO"}]), path)
    elif suffix.startswith("sdf"):
        mol = Chem.MolFromSmiles("CCO")
        mol.SetProp("_Name", "zero")
        text = Chem.MolToMolBlock(mol) + "\n$$$$\n"
        if suffix.endswith("gz"):
            with gzip.open(path, "wt") as handle:
                handle.write(text)
        else:
            path.write_text(text)
    else:
        text = (
            "CCO zero\n"
            if suffix.startswith("smi")
            else (
                "id\tsmiles\nzero\tCCO\n"
                if suffix.startswith("tsv")
                else "id,smiles\nzero,CCO\n"
            )
        )
        if suffix.endswith("gz"):
            with gzip.open(path, "wt") as handle:
                handle.write(text)
        else:
            path.write_text(text)
    rows = list(selected_rows(mst.Source(path)))
    assert rows[0][1]["smiles"] == "CCO"
    assert rows[0][1]["id"] == "zero"


def test_split_and_mapping(tmp_path):
    source = mst.Source.from_records(
        [{"smiles": "CCO", "split": "train"}, {"smiles": "CCN", "split": "test"}],
        split_column="split",
        splits=("train",),
    )
    index = mst.build_index(source, tmp_path / "index")
    assert index.size == 1
    assert next(batches(index.aliases_path))["source_row"] == 1
    with pytest.raises(ValueError, match="missing mapped columns"):
        mst.build_index(
            mst.Source.from_records(
                [{"smiles": "CCO", "split": "train"}],
                split_column="typo",
                splits=("train",),
            ),
            tmp_path / "bad",
        )


def test_zero_fp_and_overflow(tmp_path, monkeypatch):
    class Generator:
        def GetFingerprint(self, mol):
            return DataStructs.ExplicitBitVect(2048)

    monkeypatch.setattr(
        chemistry,
        "_generator",
        lambda config: Generator(),
    )
    with pytest.raises(ValueError, match="zero-bit"):
        chemistry.prepare("CCO", mst.ChemistryConfig(), mst.FingerprintConfig())
    with pytest.raises(ValueError, match="shard-local"):
        write_shard(
            tmp_path / "bad.h5",
            [{"local_id": 2**32, "fingerprint": b"\x01" * 256}],
            mst.FingerprintConfig(),
        )


def test_config_and_input_resume_refusals(tmp_path):
    path = tmp_path / "source.csv"
    path.write_text("smiles\nCCO\n")
    source = mst.Source(path)
    mst.build_index(source, tmp_path / "index")
    path.write_text("smiles\nCCN\n")
    with pytest.raises(ValueError, match="resume refused"):
        mst.build_index(source, tmp_path / "index")
    path.write_text("smiles\nCCO\n")
    with pytest.raises(ValueError, match="resume refused"):
        mst.build_index(
            source,
            tmp_path / "index",
            chemistry=mst.ChemistryConfig(remove_isotopes=True),
        )


def test_fingerprint_reuse_and_corrupt_chunk_repair(tmp_path, monkeypatch):
    from molecular_similarity_toolkit import index as index_module

    source = mst.Source.from_records([{"smiles": "CCO"}, {"smiles": "CCC"}])
    built = mst.build_index(source, tmp_path / "index", chunk_size=1)
    original = index_module.prepare
    calls = []

    def wrapped(*args):
        calls.append(1)
        return original(*args)

    monkeypatch.setattr(index_module, "prepare", wrapped)
    mst.build_index(source, built.path, chunk_size=1)
    assert not calls
    (built.path / "units/stage-0000-00000000/rows.parquet").write_bytes(b"bad")
    assert not mst.verify(built.path)["ok"]
    mst.build_index(source, built.path, chunk_size=1)
    assert len(calls) == 1
    assert mst.verify(built.path)["ok"]
