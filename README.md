# CellSimulator-Structures

A standalone, version-controlled catalog of experimental macromolecular structures for [CellSimulator](https://github.com/jgoers64/CellSimulator).

## Five starter proteins

| Protein | PDB ID | Organism |
| --- | --- | --- |
| Insulin | [2HIU](https://www.rcsb.org/structure/2HIU) | Human |
| Ubiquitin | [1UBQ](https://www.rcsb.org/structure/1UBQ) | Human |
| Lysozyme | [1LZ1](https://www.rcsb.org/structure/1LZ1) | Human |
| Myoglobin (mutant) | [3RGK](https://www.rcsb.org/structure/3RGK) | Human |
| Hemoglobin | [4HHB](https://www.rcsb.org/structure/4HHB) | Human |

Displayed names are concise (Insulin, Ubiquitin, Lysozyme, Myoglobin (mutant), Deoxyhemoglobin). All five entries are sourced from *Homo sapiens*. The available human myoglobin X-ray structure (3RGK) contains experimentally documented mutations; it is not labeled as wild type. Each entry has its own metadata file with a pinned PDB identifier, source URL, experimental method, organism, and rendering notes. The catalog **does not invent atomic coordinates**.

## Use with CellSimulator

This is a separate Git repository. A sibling checkout is recommended for development; a nested checkout is also supported. Run commands from this repository and supply the actual simulator target explicitly.

For sibling checkouts:

```powershell
python catalog.py list
python catalog.py validate
python catalog.py install all --target ..\CellSimulator
```

For a nested checkout, use `--target ..`. CellSimulator's optional `CELLSIM_STRUCTURES_ROOT` setting can identify either layout without relying on the terminal's working directory.

To install only insulin from the nested layout:

```powershell
python catalog.py install insulin --target ..
```

The installer also removes obsolete animal-protein JSON only when it recognizes its own previously generated files. Re-run `install all` after updating the repository to replace older chicken and whale models. The installer downloads the original coordinate files from [RCSB PDB](https://www.rcsb.org/) on first use and caches them in the Git-ignored `cache/` folder. It converts them into **view-only** JSON files under `CellSimulator/molecules/proteins/`, which the existing CellSimulator Molecules browser can load. CellSimulator owns the runtime catalog and resolves available structures through its shared library interface.

The local files are intended for **visual inspection**, not validated molecular dynamics. No missing hydrogen atoms, force-field parameters, or solvent molecules are fabricated. Display bonds are inferred for visual purposes only. NMR structures use their first deposited model.

For a download without installation:

```powershell
python catalog.py fetch all
```

## Repository layout

```text
catalog.json
catalog.py
structures/proteins/<protein>/metadata.json
tests/test_catalog.py
cache/                         # generated, ignored by Git
```

## Tests

```powershell
python -m unittest discover -s tests -v
```

Source PDB files remain unmodified in the cache so each visual model can be traced back to its experimental input.

## Dataset compatibility

`dataset_manifest.json` identifies this catalog release, its schema, and the SHA-256 of `catalog.json`. CellSimulator records the compatible release in `data_contracts.json` and reports the configured checkout in `python main.py --diagnose`.

A dataset version is a compatibility and provenance record. It does not mean its coordinates are ready for molecular dynamics. Source identity, assembly/model selection, processing and simulation preparation remain separate records. Catalog changes must update the manifest hash and dataset version together. Nuclear-library work in a draft branch is not part of this starter release.

Authored explanations belong in CellSimulator-Knowledge. BioSimGlossaryAgent is retired and is not a dependency.
