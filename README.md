# CellSimulator-Structures

A standalone, version-controlled catalog of experimental macromolecular structures for [CellSimulator](https://github.com/jgoers64/CellSimulator).

## Five starter proteins

| Protein | PDB ID | Organism |
| --- | --- | --- |
| Insulin | [2HIU](https://www.rcsb.org/structure/2HIU) | Human |
| Ubiquitin | [1UBQ](https://www.rcsb.org/structure/1UBQ) | Human |
| Lysozyme | [1LYZ](https://www.rcsb.org/structure/1LYZ) | Chicken |
| Myoglobin | [1MBN](https://www.rcsb.org/structure/1MBN) | Sperm whale |
| Hemoglobin | [4HHB](https://www.rcsb.org/structure/4HHB) | Human |

Each entry has its own metadata file with a pinned PDB identifier, source URL, experimental method, organism, and rendering notes. The catalog **does not invent atomic coordinates**.

## Use with CellSimulator

Clone this repo **next to** the `CellSimulator` folder, and run these commands from `CellSimulator-Structures`:

```powershell
python catalog.py list
python catalog.py validate
python catalog.py install all --target ..\CellSimulator
```

To install only insulin:

```powershell
python catalog.py install insulin --target ..\CellSimulator
```

The installer downloads the original coordinate files from [RCSB PDB](https://www.rcsb.org/) on first use and caches them in the Git-ignored `cache/` folder. It converts them into **view-only** JSON files under `CellSimulator/molecules/proteins/`, which the existing CellSimulator Molecules browser can load. The Macromolecules tab is not yet connected directly to this repository.

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
