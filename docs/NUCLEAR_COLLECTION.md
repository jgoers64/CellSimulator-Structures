# Human nuclear collection

## Scope

This collection connects a broad human nuclear inventory with coordinate assets suitable for progressive molecular visualization in BioSim. It extends the structures repository without replacing the existing five starter proteins or their `catalog.py` workflow.

The inventory contains **7,575 nuclear-associated gene entries** and **86 other molecule, candidate, RNA-family, or lipid-family entries**. These are not 7,661 experimentally resolved structures, and they are not a claim that every entry is present in one nucleus. Some annotations are uncertain or state-specific. A protein gene, protein group, chemical species, molecular complex, and rendered instance are different counting units.

Default quantities for the generic diploid G1 environment remain unknown. The separate HeLa-derived protein-group estimates are research references tied to their original experimental context; they are not assigned to the default nucleus. A blank or `null` quantity means unknown, not absent.

Inventory sources:

- [Human Protein Atlas downloadable data](https://www.proteinatlas.org/about/download)
- [Human Protein Atlas subcellular data definitions](https://www.proteinatlas.org/humanproteome/subcellular/data#locations)
- [Itzhak et al., Global, quantitative and dynamic mapping of protein subcellular localization](https://doi.org/10.7554/eLife.16950)

Human Protein Atlas data retain their [CC BY 4.0 attribution requirements](https://www.proteinatlas.org/about/licence). Relevant reference: [Thul et al., A subcellular map of the human proteome](https://doi.org/10.1126/science.aal3321).

## Catalog records and coordinate assets

A **catalog record** identifies something of biological interest and preserves its evidence. A **coordinate asset** supplies an actual representation that can be fetched or installed. Only entries connected to a supported asset can become a molecular structure through the importer.

Use the collection manifest and the command-line asset listing for current availability and counts. Do not infer availability from the number of inventory entries. RNA and lipid families still require individual sequence or chemical-species selection; an entire family cannot be represented accurately by one arbitrary set of atoms.

The collection supports different evidence classes:

- Experimental macromolecular coordinate models, with a specified source and assembly.
- Integrative models, with their mixed evidence and resolution made explicit.
- Defined small molecules represented by bundled JSON assets, with their coordinate-generation method and chemical state recorded.
- Catalog-only entries for which no suitable coordinate asset is supplied.

An asset's presence is not a claim of completed scientific review. Read its provenance, coverage, review status, and limitations. The standards and interpretation are described in [Structure reliability](STRUCTURE_RELIABILITY.md).

## Inspect the collection

Run these commands from the `CellSimulator-Structures` directory:

```powershell
python nucleus.py summary
python nucleus.py list --query PCNA
python nucleus.py assets
```

`summary` reports collection coverage, `list` searches the inventory, and `assets` identifies supported coordinate assets. A search result without an asset remains a catalog record.

These commands inspect the catalog. They do not populate a nucleus with every listed molecule or determine a cell-specific abundance distribution.

The inventory is the identity/evidence layer; current coordinate availability is
joined from the asset manifests. A linked macromolecular assembly can represent
only a fragment or a component of a protein entry. It is not automatically a
full-length, isolated model of that protein.

## Predicted protein models

Many protein entries have a verified UniProt cross-reference to AlphaFold. A
cross-reference is not a downloaded or scientifically reviewed model. Predictions
are an explicit choice and retain their own evidence class:

```powershell
python nucleus.py resolve ERH --predicted
```

The resolver uses the provider's actual metadata and download URLs. It checks
the requested accession, human taxonomy when supplied, sequence bounds, and the
agreement between API sequence and coordinate sequence. It pins the model
version, source checksum, local confidence, and available PAE evidence. Multiple
accessions or fragment models require `--accession` or `--model-id` selection.
Relative-domain confidence remains unreviewed if PAE cannot be cached within the
size budget. No confidence number is treated as an accuracy percentage.

After resolution, `python nucleus.py assets` includes the exact cached model ID;
pass that returned ID to `export` or `install`. Prediction caches are ignored by
Git. Neither a search nor `fetch all` discovers predictions for every gene.
`fetch all` acts on the curated assets and any predictions already explicitly
selected in that local cache.

## Dependencies and source data

Curated macromolecular imports use source **PDBx/mmCIF** data and **Gemmi 0.7.5**. Install the importer dependency from the pinned requirements file:

```powershell
python -m pip install -r requirements-structures.txt
```

Raw downloaded sources are cached in Git-ignored storage. Repository metadata and source checksums make an import traceable without committing every raw structure download. Chemical assets are bundled as JSON; their metadata must distinguish idealized geometry from experimentally interpreted coordinates.

The original five-protein `catalog.py` commands remain available. The nuclear importer is an additional workflow, and installing its optional dependency does not change the starter catalog's source choices.

## Fetch or install a structure

For example, fetch the supported PCNA asset:

```powershell
python nucleus.py fetch pdb-1vym
```

Fetching retrieves the source or prepares the asset cache. It does not by itself install the structure into an application checkout.

If the repositories are next to each other:

```powershell
python nucleus.py install pdb-1vym --target ../CellSimulator
```

If this repository is nested inside the `CellSimulator` directory:

```powershell
python nucleus.py install pdb-1vym --target ..
```

The target must be an existing CellSimulator checkout. Installation makes the selected view-only asset available to the application's structure-loading workflow. It does not assign molecule counts, decide nuclear placement, or redesign the nucleus scene. In particular, a PCNA structure does not imply that a G1 scene contains active S-phase replication forks.

Use `python nucleus.py --help` and subcommand help for the exact supported batch options. Bulk operations apply to available coordinate assets, not to every inventory row. Large assemblies may require substantial downloads, conversion time, memory, or rendering resources. Import and view a small asset first when checking a new installation, then select larger assemblies deliberately.

The default conversion budget is 250,000 atoms, checked before expanding a
biological assembly. In this snapshot the full 7R5J pore assembly would contain
up to 4,937,064 source/model atoms, so full export or installation is refused at
the default budget. Its pinned source remains in the collection for future
chunked rendering. Raising `--max-atoms` is an explicit memory decision and does
not improve the scientific evidence or make the current renderer scalable.

The importer preserves source angstrom coordinates and the selected assembly
operations. It selects one residue-wide alternate conformation and omits source
waters, recording both choices. Other source ligands and modified residues are
retained. Within-residue bond templates come from the archive; consecutive
standard polymer links are separately disclosed display inferences. Crosslinks
and metal coordination are not reconstructed, so the bond graph is incomplete
and remains unsuitable for molecular dynamics.

## What installation does not establish

The imported structure is an experimental or modeled reference in a particular state. It is not an atomically complete living molecule unless the evidence actually supports that scope. Missing residues, hydrogens, cofactors, and flexible domains must not be silently filled with invented experimental coordinates.

A structure record does not specify a nuclear population recipe. Counts, compartment placement, interactions, and conformations depend on cell type and state. The inventory deliberately keeps those questions separate from source geometry. In particular:

- Protein-group abundance must not be copied to every member gene or accession.
- A complex and its bound constituents must not be counted as additional independent populations.
- Elemental ion measurements include bound pools and do not automatically specify free-ion concentrations.
- An RNA or lipid family is not a single chemical species with a universal geometry.
- A nucleosome template does not determine an entire chromosome's folding.
- Imported coordinates are not a molecular-dynamics parameterization.

Source-backed surfaces, backbones, and atomic detail can all be useful at different viewing scales. The renderer should load them selectively rather than instantiate every molecule or atom in the entire inventory simultaneously.

## Tests and reversibility

The collection changes are version controlled on a dedicated branch. Import commits use `[skip ci]` to respect the current instruction not to run the test suite automatically. Test execution remains an explicit later action; preparing a test does not mean it has passed.

When test execution is requested, run the repository suite from this directory:

```powershell
python -m unittest discover -s tests -v
```

This documentation is not a statement that the suite or a full scientific structure-validation workflow has already been run. Review the actual change report for completed checks and remaining limitations.
