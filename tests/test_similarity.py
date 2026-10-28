from fractions import Fraction

import pytest
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

import molecular_similarity_toolkit as mst
from molecular_similarity_toolkit.storage import batches

SMILES = [
    "CCO",
    "CCCO",
    "CCCCO",
    "CCN",
    "CCC",
    "c1ccccc1",
    "Cc1ccccc1",
    "CCOC",
    "CC(=O)O",
    "CCCCCC",
    "C[C@H](O)F",
    "C[C@@H](O)F",
    "OCC",
    "[13CH3]CO",
]


def rdkit_counts(left, right):
    # Independent exhaustive oracle: generate fingerprints directly with RDKit,
    # never call toolkit preparation/counting or its FPSim2 adapter.
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, includeChirality=False
    )
    a, b = [generator.GetFingerprint(Chem.MolFromSmiles(s)) for s in (left, right)]
    numerator = (a & b).GetNumOnBits()
    denominator = (a | b).GetNumOnBits()
    assert numerator / denominator == DataStructs.BulkTanimotoSimilarity(a, [b])[0]
    return numerator, denominator


def rows(result, table=None):
    return [row for batch in result.batches(table) for row in batch.to_pylist()]


@pytest.mark.parametrize("shard_size", [1, 4, 100])
@pytest.mark.parametrize("workers", [1, 3])
def test_exact_topk(tmp_path, shard_size, workers):
    index = mst.build_index(
        mst.Source.from_records(
            [{"id": str(i), "smiles": s} for i, s in enumerate(SMILES)]
        ),
        tmp_path / "index",
        shard_size=shard_size,
    )
    molecules = list(batches(index.molecules_path))
    queries = mst.Source.from_records(
        [{"id": "repeated", "smiles": "CCO"}, {"id": "repeated", "smiles": "CCCO"}]
    )
    result = mst.search(
        index, queries, 5, tmp_path / "result", workers=workers, exclude_identity=True
    )
    query_index = mst.MoleculeIndex.open(result.path / "queries")
    qm = {r["molecule_id"]: r for r in batches(query_index.molecules_path)}
    actual = rows(result)
    assert len(actual) == 10
    for query_id in {row["query_id"] for row in actual}:
        hits = [r for r in actual if r["query_id"] == query_id]
        query = qm[hits[0]["query_molecule_id"]]
        expected = []
        for molecule in molecules:
            if molecule["molecule_id"] == query["molecule_id"]:
                continue
            a, b = rdkit_counts(query["canonical_smiles"], molecule["canonical_smiles"])
            expected.append((Fraction(a, b), molecule["molecule_id"]))
        expected.sort(key=lambda x: (-x[0], x[1]))
        assert [
            (Fraction(r["intersection"], r["union"]), r["molecule_id"]) for r in hits
        ] == expected[:5]
    assert result.summary["match_count"] == 10


@pytest.mark.parametrize("comparison", [">", ">="])
@pytest.mark.parametrize("threshold", ["0", "0.8", "1"])
def test_filter_exact(tmp_path, comparison, threshold):
    pool = mst.build_index(
        mst.Source.from_records([{"smiles": s} for s in SMILES] + [{"smiles": "bad"}]),
        tmp_path / "pool",
        shard_size=4,
    )
    protected = mst.build_index(
        mst.Source.from_records([{"smiles": "CCO"}, {"smiles": "CCCO"}]),
        tmp_path / "protected",
    )
    result = mst.filter_pool(
        pool, protected, threshold, comparison, tmp_path / "result"
    )
    expected = {}
    for mol in batches(pool.molecules_path):
        for query in batches(protected.molecules_path):
            a, b = rdkit_counts(mol["canonical_smiles"], query["canonical_smiles"])
            hit = (
                Fraction(a, b) > Fraction(threshold)
                if comparison == ">"
                else Fraction(a, b) >= Fraction(threshold)
            )
            if hit:
                expected.setdefault(mol["molecule_id"], query["molecule_id"])
    assert {
        r["molecule_id"]: r["protected_id"] for r in rows(result, "rejected")
    } == expected
    assert (
        result.summary["eligible_count"] + result.summary["rejected_count"] == pool.size
    )
    assert len(rows(result, "invalid")) == 1


def test_empty_and_invalid_protected(tmp_path):
    empty = mst.build_index(mst.Source.from_records([]), tmp_path / "empty")
    pool = mst.build_index(
        mst.Source.from_records([{"smiles": "CCO"}]), tmp_path / "pool"
    )
    assert (
        rows(
            mst.search(
                empty,
                mst.Source.from_records([{"smiles": "CCO"}]),
                out=tmp_path / "search",
            )
        )
        == []
    )
    with pytest.raises(ValueError, match="empty"):
        mst.filter_pool(pool, empty, out=tmp_path / "filter")
    with pytest.raises(ValueError, match="invalid"):
        mst.filter_pool(
            pool,
            mst.Source.from_records([{"smiles": "CCO"}, {"smiles": "bad"}]),
            out=tmp_path / "invalid",
        )


def test_k_larger_than_pool(tmp_path):
    index = mst.build_index(
        mst.Source.from_records([{"smiles": "CCO"}, {"smiles": "CCN"}]),
        tmp_path / "index",
    )
    assert (
        len(
            rows(
                mst.search(
                    index,
                    mst.Source.from_records([{"smiles": "CCO"}]),
                    99,
                    tmp_path / "search",
                )
            )
        )
        == 2
    )


def test_exact_eighty_percent_boundary(tmp_path):
    short = "O=C(O)C(=O)CC(O)C(O)COP(=O)(O)O"
    long = "O=C(O)C(=O)CC(O)C(O)C(O)COP(=O)(O)O"
    pool = mst.build_index(
        mst.Source.from_records([{"smiles": short}]), tmp_path / "pool"
    )
    protected = mst.Source.from_records([{"smiles": long}])
    strict = mst.filter_pool(pool, protected, "0.8", ">", tmp_path / "strict")
    inclusive = mst.filter_pool(pool, protected, "0.8", ">=", tmp_path / "inclusive")
    assert strict.summary["rejected_count"] == 0
    hit = rows(inclusive, "rejected")[0]
    assert (hit["intersection"], hit["union"]) == (24, 30)


def test_ties_and_fingerprint_identity_distinction(tmp_path):
    source = mst.Source.from_records([{"smiles": "C" * n} for n in range(8, 15)])
    pool = mst.build_index(source, tmp_path / "pool", shard_size=1)
    query = mst.Source.from_records([{"smiles": "C" * 10}])
    result = mst.search(pool, query, 2, tmp_path / "result", exclude_identity=True)
    hits = rows(result)
    assert all(hit["similarity"] == 1 for hit in hits)
    query_id = hits[0]["query_molecule_id"]
    expected = sorted(
        row["molecule_id"]
        for row in batches(pool.molecules_path)
        if row["molecule_id"] != query_id
    )[:2]
    assert [row["molecule_id"] for row in hits] == expected


def test_filter_parallel_equivalence_and_early_stop(tmp_path, monkeypatch):
    from molecular_similarity_toolkit.adapter import ShardEngine

    pool = mst.build_index(
        mst.Source.from_records([{"smiles": "C" * 10}]), tmp_path / "pool"
    )
    protected = mst.Source.from_records([{"smiles": "C" * n} for n in range(8, 50)])
    original = ShardEngine.threshold_candidates
    calls = []

    def count(self, fp, threshold):
        calls.append(1)
        return original(self, fp, threshold)

    monkeypatch.setattr(ShardEngine, "threshold_candidates", count)
    parallel = mst.filter_pool(pool, protected, out=tmp_path / "parallel", workers=4)
    assert len(calls) <= 4
    serial = mst.filter_pool(pool, protected, out=tmp_path / "serial", workers=1)
    assert rows(parallel, "rejected") == rows(serial, "rejected")
