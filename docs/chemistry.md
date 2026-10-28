# Chemistry rules

A similarity result depends on both molecular preparation and fingerprint
settings. Store these settings with your experiment and reuse the same index
for repeat runs.

## Default: stereo-free

```python
chemistry = mst.ChemistryConfig()
fingerprint = mst.FingerprintConfig(radius=2, n_bits=2048, use_chirality=False)
```

| Property | Default behavior |
|---|---|
| Stereochemistry | Removed before canonicalization and fingerprinting |
| Formal charge | Kept; a charged and neutral form can have different IDs |
| Salts and disconnected fragments | Kept |
| Tautomer form | Kept |
| Isotope labels | Kept in canonical identity |
| Fingerprint | Binary Morgan, radius 2, 2,048 bits, chirality off |

The canonical SMILES is hashed with SHA-256 to make the stable molecule ID.
All original source memberships survive in aliases. Deduplication does not use
fingerprints or a truncated InChIKey.

Keeping a chemical feature in canonical identity does not mean every fingerprint
can distinguish it. Morgan fingerprints can collide, and isotope labels do not
guarantee different bit patterns. Different molecule IDs may therefore have
similarity 1.

## Normalized parent

```python
chemistry = mst.ChemistryConfig(normalized_parent=True)
```

This adds RDKit cleanup, parent-fragment selection, and uncharging, in that order.
It is useful when your definition of a molecule should combine salts and some
charge forms. It changes identity and scores, so build a new index. Tautomer
canonicalization is still off unless `canonical_tautomer=True` is selected.

## Stereo and isotopes are separate

```python
chemistry = mst.ChemistryConfig(remove_stereochemistry=False, remove_isotopes=True)
fingerprint = mst.FingerprintConfig(use_chirality=True)
```

Chiral fingerprints require retained stereochemistry. Isotope removal is an
independent option. The default differs from older pipelines that called
`MolToSmiles(..., isomericSmiles=False)`, because that older setting also drops
isotope labels.

## Boundary behavior

Tanimoto is the number of shared on-bits divided by the number of on-bits in
either molecule. Intersection and union are integers.

- `comparison='>'`, `threshold='0.8'`: remove scores strictly greater than 4/5.
- `comparison='>='`, `threshold='0.8'`: also remove scores equal to 4/5.

You may pass a decimal string to state the intended decimal boundary explicitly.
The backend searches a slightly lower float32 boundary to collect candidates;
integer cross-multiplication makes the final decision. Top-k boundary candidates
are rescored with integers, and equal scores are ordered by stable molecule ID.

## Compatibility and invalid molecules

Indexes must use identical chemistry settings, fingerprint settings, and RDKit
versions for comparison. The package refuses a mismatch. Dependency versions
and resolved settings are recorded in manifests.

Invalid source molecules and zero-bit fingerprints go into a quarantine table
with source-row provenance and a reason. Invalid protected molecules fail the
job by default. `invalid_protected='quarantine'` is an explicit opt-out; if you
use it, protection only covers the remaining valid reference molecules. An
empty protected set is always an error.
