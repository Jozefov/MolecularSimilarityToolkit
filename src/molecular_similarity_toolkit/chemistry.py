"""Single structure/fingerprint implementation used by all pipeline stages."""

import hashlib
from functools import lru_cache

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize

from .config import ChemistryConfig, FingerprintConfig, validate_configs

RDKIT_VERSION = rdBase.rdkitVersion


class MoleculeError(ValueError):
    """A molecule cannot be parsed or has an empty fingerprint."""


def prepare(
    smiles: str, chemistry: ChemistryConfig, fingerprint: FingerprintConfig
) -> tuple[str, str, bytes, int]:
    """Return canonical identity, stable molecule ID, packed fingerprint and popcount."""
    validate_configs(chemistry, fingerprint)
    if not isinstance(smiles, str) or not smiles.strip():
        raise MoleculeError("missing SMILES")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not mol.GetNumAtoms():
        raise MoleculeError("invalid or empty molecule")
    # Remove ignored annotations before any canonical choice can depend on them.
    if chemistry.remove_isotopes:
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)
    if chemistry.remove_atom_maps:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    if chemistry.remove_stereochemistry:
        Chem.RemoveStereochemistry(mol)
    if chemistry.normalized_parent:
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
    if chemistry.canonical_tautomer:
        mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
    Chem.SanitizeMol(mol)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    generator = _generator(fingerprint)
    fp = generator.GetFingerprint(mol)
    count = fp.GetNumOnBits()
    if not count:
        raise MoleculeError("zero-bit fingerprint")
    return (
        canonical,
        hashlib.sha256(canonical.encode()).hexdigest(),
        DataStructs.BitVectToBinaryText(fp),
        count,
    )


def intersection_union(left: bytes, right: bytes) -> tuple[int, int]:
    """Exact binary Tanimoto numerator and denominator."""
    a, b = int.from_bytes(left, "little"), int.from_bytes(right, "little")
    return (a & b).bit_count(), (a | b).bit_count()


@lru_cache(maxsize=16)
def _generator(config: FingerprintConfig):
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=config.radius,
        fpSize=config.n_bits,
        includeChirality=config.use_chirality,
    )
