<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/assets/logo-dark.svg">
  <img alt="MolecularSimilarityToolkit: find close matches and build protected training pools" src="docs/assets/logo.svg" width="960">
</picture>

You have a small set of molecules for testing your model and a much bigger
collection for training. This tool helps you remove molecules that look too
similar to the test set. It can also find the closest molecules to any query.

- **Build a training pool:** remove close matches to your validation and test molecules.
- **Find candidates:** search your whole collection and put the closest matches first.

Bring your own files: PubChem, ZINC, a benchmark, or any other molecule collection.
Build an index once and reuse it. Jobs can resume after an interruption and run
as separate tasks on a cluster.

<picture>
  <source media="(max-width: 600px) and (prefers-color-scheme: dark)" srcset="docs/assets/workflow-mobile-dark.svg">
  <source media="(max-width: 600px)" srcset="docs/assets/workflow-mobile.svg">
  <source media="(prefers-color-scheme: dark)" srcset="docs/assets/workflow-dark.svg">
  <img alt="Real molecule examples: aspirin returns ranked matches scoring 1.00, 0.67 and 0.61. Protecting aspirin at a cutoff above 0.8 removes the 1.00 match and keeps the other two." src="docs/assets/workflow.svg" width="1600">
</picture>

## Install

Python **3.11–3.13** on Linux x86-64 or macOS arm64. This is a private release
candidate. The repository remains private until release review is approved.

```bash
pip install 'git+https://github.com/Jozefov/MolecularSimilarityToolkit.git'
```

The HTTPS URL is the public installation URL prepared for release. While the
repository is private, authenticate Git or install from an authorized local
checkout with `pip install .`. Nothing is published on PyPI yet.

## Build a protected training pool

```python
import molecular_similarity_toolkit as mst

pool = mst.build_index(
    mst.Source('pool.tsv', name='pool', id_column='id', smiles_column='smiles'),
    out='indexes/pool',
)
protected = mst.build_index(
    mst.Source('benchmark.tsv', name='benchmark',
               split_column='fold', splits=['val', 'test']),
    out='indexes/protected',
)
clean = mst.filter_pool(
    pool, protected=protected, threshold=0.8, comparison='>',
    out='runs/clean',
)
clean.export('eligible.tsv')
```

The example protects validation and test. Training molecules are not protected.
You choose the files and folds to protect; a reference file without a split
selection protects every molecule in that file.

## Find the closest candidates

```python
hits = mst.search(
    pool, queries=mst.Source('queries.tsv', name='queries'),
    k=100, out='runs/candidates',
)
hits.export('candidates.tsv')
```

Matches are ranked from closest to furthest. An identical molecule can be a
match; use `exclude_identity=True` to leave it out. The result files stay on disk.
Use `hits.batches()` to read a little at a time. Later, reopen a run with
`mst.Result.open('runs/candidates')`. The CLI also exports tables:
`mol-sim export runs/candidates candidates.tsv`.

Try the complete tiny example from a checkout:

```bash
uv sync --locked --group dev
uv run python examples/quickstart.py
uv run mol-sim --help
```

## What does “similar” mean?

A **fingerprint** is a row of on/off bits describing parts of a molecule. A
**Tanimoto score** compares two fingerprints. A score of 1 means the bits are
identical; 0 means no bits overlap. Different molecules can sometimes have the
same fingerprint.

The default fingerprint is Morgan, radius 2, with 2,048 bits. The default
chemistry removes stereochemistry and keeps charges, salts, isotope labels,
and tautomer forms. Atom-map labels are removed by default. You can change
these rules, but both sides of a comparison
must use the same settings.

With `comparison='>'` and `threshold=0.8`:

| Closest protected score | Action |
|---|---|
| 0.67 | Keep |
| Exactly 0.80 | Keep |
| 1.00 | Remove |

Every removed molecule gets a **witness**: a protected molecule that caused its
removal, with the score. A witness is not necessarily the closest match.

Candidate output includes the query ID, molecule ID, rank, and score. Separate
alias tables link molecule IDs back to all original source rows, even when the
same molecule appears in several files.

## What this guarantees

Search is exact for the selected fingerprint and chemistry settings. Filtering
keeps molecules whose scores against the protected set meet your rule. It does
not check every other kind of data leakage, such as repeated measurements,
shared experimental sources, or labels copied between datasets.

The package uses packed fingerprints and CPU search from
[FPSim2](https://github.com/chembl/FPSim2), and molecular preparation from
[RDKit](https://www.rdkit.org/). Performance depends on your molecules, cutoff,
queries, disk, and hardware. See the measured benchmark notes before estimating
a large run.

## Learn more

- [Python API](docs/api.md) and [input/output tables](docs/data.md)
- [Chemistry and fingerprint rules](docs/chemistry.md)
- [Restart, progress, and portable compact indexes](docs/restart.md)
- [Cluster jobs: Slurm and PBS](docs/hpc.md)
- [Performance and benchmark commands](docs/performance.md)
- [Release checks](docs/release.md), [contributing](CONTRIBUTING.md), and [MIT license](LICENSE)

GPU execution, automatic database downloads, new train/val/test splits, and
approximate search are outside this first release.
