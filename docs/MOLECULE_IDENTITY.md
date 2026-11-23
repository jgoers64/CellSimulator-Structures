# Reuse molecular identities across collections

The Human Proteome catalog remains the simulator's protein identity catalog.
The nuclear collection adds gene-localization evidence, source models, and
assembly relationships. It must not create another protein definition merely
because an existing protein also occurs in the nucleus.

## Three different records

| Record | Example | Reuse rule |
| --- | --- | --- |
| Protein identity | `human-proteome/P68431` | One existing UniProt protein record can have several source gene annotations. |
| Structural representation | PDB 2CV5, selected assembly and atom policy | Keep experimental/predicted models, source revisions, fragments and selections distinguishable. |
| Scene instance | A placed nucleosome in a chromatin patch | Many placements can share geometry; each placement is not a new protein species. |

The existing nucleus renderer already shares a nucleosome template. The nuclear
importer now recognizes that template instead of installing a second default
2CV5 entry. Its polymer-only coordinates and the fuller archive conversion are
**different representations**; reuse is not described as byte-for-byte equality.

The source inventory is still keyed by Ensembl gene because its observations
are gene-level evidence. Those rows are not standalone molecule definitions.
`structure_library.identity` derives unique verified UniProt references without
copying protein sequences or coordinates into a second protein catalog.

For the currently pinned inventory, 7,575 gene rows produce **7,394 unique verified
protein references**. There are also two unresolved accession labels and 140
gene rows without an accession. These are identifier-mapping counts, not numbers
of complete models or molecules in one nucleus. A row with several associations
can contribute to more than one reporting category.

For example, the ten H3C gene rows mapped to `P68431` retain ten sets of source
evidence, but share one `human-proteome/P68431` reference. Different UniProt
records are not merged because their names or sequences look alike.

## Inspect the user's actual downloaded catalog

From `CellSimulator-Structures`, when it is nested inside `CellSimulator`:

```powershell
python nucleus.py reconcile --target ..
python nucleus.py reconcile --target .. --json
python nucleus.py plan all --target ..
```

For neighboring repositories, use `--target ../CellSimulator` instead.

These commands are offline and read-only. Reconciliation reads the existing
`macromolecules/catalog/human_proteome.sqlite` using SQLite's read-only mode.
The full JSON report retains the gene evidence behind each canonical reference.

| Status | Meaning |
| --- | --- |
| `known` | Exact verified accession is already present in the inspected target database. |
| `reference_only` | The source mapping is verified, but this target database is missing or unreadable. |
| `missing_from_index` | Exact accession is absent from a readable target database. |
| Unresolved mapping | Source accession is missing, unverified, conflicting or unsupported. |

Missing local downloads are not evidence that the user never imported the Human
Proteome. No database is created by these commands. An indexed protein also does
not establish that its coordinates have been downloaded or that a model covers
the full protein.

## Installation decisions

`nucleus.py install` uses the same plan before fetching a source. The checks
cover existing molecule JSON, built-in structure metadata, current Human
Proteome cache sidecars and the bundled ChEBI pack index. A partial checkout,
unreadable entry, or malformed pack does not count as a clean catalog.

| Decision | Behavior |
| --- | --- |
| `reuse_protein_identity` | Return the existing Human Proteome reference; write no standalone nuclear protein JSON. |
| `reuse_exact_asset` | Reuse an already installed matching source/selection/geometry, including after a file rename. |
| `reuse_registered_asset` | Keep an intact registered asset with the same source checksum; do not claim exact equivalence to a new conversion. |
| `reuse_builtin_representation` | Keep the verified existing 2CV5 default; distinguish its atom selection from the fuller conversion. |
| `reuse_chemical_identity` | Keep the existing representation of a verified CCD state or the specifically checked legacy NIST water record. |
| `needs_review` | Stop this asset when state, identity, assembly selection, existing data integrity or overlap is unresolved. |
| `needs_proteome_index` | Use the application's Human Proteome workflow; do not work around it by installing another protein JSON. |
| `needs_chunking` | Leave an assembly above the configured atom budget uninstalled. |
| `install_new_asset` | No existing match was found in the inspected catalogs; a selected chemical/assembly representation may be added. |

Names, formulas and standard InChI descriptors can help find candidates. A name
or standard InChI match alone does not authorize merging stereoisomers, tautomers,
charge states or bound complexes. Chemical exports retain CCD identifiers,
formula, explicit formal charge and original descriptors. Existing legacy
records with insufficient identity metadata may need review instead of an
automatic match. The packed ChEBI index is not assumed to contain every synonym
or chemical descriptor: this is a conservative import guard, not proof that all
historical cross-database duplicates have already been eliminated.

Protein predictions with an existing canonical Human Proteome entry reuse that
entry. The nucleus importer does not replace the application's selected cached
model. An explicitly selected alternative remains available for `export` in the
Structures workflow. An isoform suffix is preserved and never stripped to force
a match with the current canonical-only simulator loader.

PDB source identity alone is not an exact structure key. The same PDB can contain
multiple assemblies, chains, models, unresolved regions and chemical states.
An existing source with a different or unverified atom selection is reported for
review. Older proteome coordinate caches without metadata sidecars also prevent
assuming that a PDB source is absent.

Reusing a protein identity does not replace an experimentally determined
nucleosome with independently positioned histone monomers. Assemblies keep
their source geometry. Component links derived from PDB entity metadata retain
their accession/isoform and explicit human biological source. These links are
not proof of full-length coverage, current accession validity, selected-assembly
membership or a particular modification state.

Repeated imports do not overwrite user files. Exact installed representations
are reused, conflicting contents under the same asset ID are rejected, and new
files are created without replacing an existing destination. Nuclear importer
writes are serialized with a target lock and the target is re-inspected before
creation. If a process is interrupted and leaves a lock, inspect that process
before manually removing the stale lock. Other importers and user edits are
outside this lock's scope.

## Counting and current scope

No scene populations are added by reconciliation or installation. Histones
contained in a nucleosome must not also be counted as a second population of
free histones. Free histones, when supported by evidence, would be separate
instances referencing the same identities and their specified states.

This change does not reorganize the simulator's navigator, rewrite its legacy
molecule catalog, deduplicate existing files destructively, or implement a
shared content-addressed cache for every protein download. In particular, the
existing per-accession proteome downloader can cache an entire PDB complex for
several accessions. Consolidating those caches requires chain/assembly-aware
work in the simulator and is deliberately separate from the nuclear importer.

## Verification

Focused regression tests are prepared in `tests/test_molecule_identity.py` and
`tests/test_nuclear_installation.py`. They must be run only when explicitly
requested, following the owner's existing no-automatic-tests instruction:

```powershell
python -m unittest discover -s tests -p "test_molecule_identity.py" -v
python -m unittest discover -s tests -p "test_nuclear_installation.py" -v
```

These tests use synthetic fixtures to check identity and file-handling behavior.
Passing them would not validate an experimental structure or nuclear abundance.

## Identifier references

- [UniProt accession numbers](https://www.uniprot.org/help/accession_numbers)
- [UniProt alternative products and isoform identifiers](https://www.uniprot.org/help/alternative_products)
- [wwPDB Chemical Component Dictionary and protonation variants](https://www.wwpdb.org/data/ccd)
- The pinned inventory and source snapshot in `collections/human-nucleus/`.
