# CellSimulator-Structures

A standalone, version-controlled catalog of experimental macromolecular structures for [CellSimulator](https://github.com/jgoers64/CellSimulator).

## Human nuclear collection

The repository also contains a source-aware nuclear library:

- **7,575 nuclear-associated protein gene records**, with evidence levels and verified UniProt source cross-references.
- **86 other molecule or family records**, preserving unresolved identities and unknown G1 quantities.
- **70 bundled chemical coordinate assets**, including 20 free amino-acid representatives. These use explicitly labeled **computed CCD ideal coordinates**, with exact dictionary chemistry and charge states.
- **13 pinned macromolecular sources**, including human nucleosomes, PCNA, NPM1, Ku–DNA, polymerases, RNA-processing assemblies, a lamin fragment, and an integrative nuclear-pore model.
- A separate **1,794-group HeLa abundance reference**, never assigned automatically to the default nucleus.

An inventory record is not a complete 3D structure. Cross-references still need appropriate sequence, assembly, state, and local-quality review. Unknown RNA sequences and lipid species receive no invented geometry. Read [the nuclear collection guide](docs/NUCLEAR_COLLECTION.md) and [the structure reliability policy](docs/STRUCTURE_RELIABILITY.md).

```powershell
python nucleus.py summary
python nucleus.py list --query PCNA
python nucleus.py assets
python -m pip install -r requirements-structures.txt
python nucleus.py fetch pdb-1vym
```

Install the PCNA example into a neighboring application checkout:

```powershell
python nucleus.py install pdb-1vym --target ../CellSimulator
```

Use `--target ..` when this repository is nested inside `CellSimulator`. Installation writes view-only files under `molecules/nucleus/`; it does not populate or change the nucleus scene. Macro sources download only when explicitly fetched or installed; searching the inventory never downloads thousands of atomic models. The original five-protein workflow below remains unchanged.

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

Clone this repo **inside** your `CellSimulator` folder or next to it. Run these commands from `CellSimulator-Structures`.

If cloned inside `CellSimulator` (your current layout):

```powershell
python catalog.py list
python catalog.py validate
python catalog.py install all --target ..
```

If the repositories are next to one another, use `--target ..\CellSimulator` instead.

To install only insulin from the nested layout:

```powershell
python catalog.py install insulin --target ..
```

The installer also removes obsolete animal-protein JSON only when it recognizes its own previously generated files. Re-run `install all` after updating the repository to replace older chicken and whale models. The installer downloads the original coordinate files from [RCSB PDB](https://www.rcsb.org/) on first use and caches them in the Git-ignored `cache/` folder. It converts them into **view-only** JSON files under `CellSimulator/molecules/proteins/`, which the existing CellSimulator Molecules browser can load. The Macromolecules tab is not yet connected directly to this repository.

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

Run tests only when explicitly requested. Starting the application or browsing the
collection must not run them. Nuclear-library commits use `[skip ci]` while the
current no-automatic-test instruction applies; existing test definitions remain intact.

```powershell
python -m unittest discover -s tests -v
```

Source PDB files remain unmodified in the cache so each visual model can be traced back to its experimental input.
