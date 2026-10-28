import pytest
from molecular_similarity_toolkit.config import ChemistryConfig, FingerprintConfig, Source
from molecular_similarity_toolkit.chemistry import prepare, intersection_union


def canonical(s, **options):
    return prepare(s, ChemistryConfig(**options), FingerprintConfig())[0]


def test_default_identity_policy():
    assert canonical("C[C@H](O)F") == canonical("C[C@@H](O)F")
    assert canonical("CC(=O)O") != canonical("CC(=O)[O-]")
    assert canonical("CCO.[Na+]") != canonical("CCO")
    assert canonical("[13CH3]CO") != canonical("CCO")
    assert canonical("CC(=O)C") != canonical("CC(O)=C")
    assert canonical("[13CH3]CO", remove_isotopes=True) == canonical("CCO")
    assert canonical("CC(=O)[O-].[Na+]", normalized_parent=True) == canonical("CC(=O)O")
    assert canonical("CC(O)=C", canonical_tautomer=True) == canonical("CC(=O)C", canonical_tautomer=True)
    assert canonical("C[C@H](O)F", remove_stereochemistry=False) != canonical("C[C@@H](O)F", remove_stereochemistry=False)


def test_validation():
    with pytest.raises(ValueError):
        prepare("CC", ChemistryConfig(), FingerprintConfig(use_chirality=True))
    with pytest.raises(ValueError):
        prepare("", ChemistryConfig(), FingerprintConfig())
    fp = prepare("CCO", ChemistryConfig(), FingerprintConfig())[2]
    a, b = intersection_union(fp, fp)
    assert a == b > 0
    assert Source.from_records([{"smiles": "CC"}]).to_dict()["_records"]
