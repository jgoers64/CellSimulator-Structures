# Structure reliability

## What this repository can verify

A structural model can be traceable, useful, and well supported without being a complete or perfectly accurate picture of a molecule in a living nucleus. Experimental and predicted structures both have limitations. Passing an importer check does not certify biological correctness, and a structure's reported resolution is not an accuracy percentage.

This repository keeps four questions separate:

1. **Identity:** Which molecule, species, sequence, construct, and biochemical state does the record describe?
2. **Structural evidence:** What supports its coordinates, and which regions are missing or uncertain?
3. **Import fidelity:** Did conversion preserve the selected source model, assembly, chemical identity, and physical scale?
4. **Nuclear context:** What supports its presence, quantity, and location in the selected cell type and cell-cycle state?

A verified source download answers neither the abundance question nor the placement question. A visually convincing scene also does not answer them.

Source: [RCSB, Assessing the Quality of 3D Structures](https://www.rcsb.org/docs/general-help/assessing-the-quality-of-3d-structures).

## Evidence classes

| Class | Meaning | Limits that must remain visible |
| --- | --- | --- |
| Experimental coordinate model | Coordinates interpreted from experimental data, such as diffraction, cryo-EM, or NMR measurements | Finite resolution, uncertain local regions, missing atoms, construct differences, and particular experimental conditions |
| Predicted model | Coordinates produced by a computational prediction method | Prediction confidence, sequence/model version, uncertain domains or interfaces, and absent biological partners |
| Integrative model | A structure assembled using multiple experimental and computational inputs | Different components can have different evidence and spatial resolution; the entire model must not be labeled as experimentally resolved atoms |
| Idealized chemical conformer | A computationally generated geometry for a defined chemical species | Correct connectivity or stereochemistry does not establish the pose of a free molecule in a nucleus |
| Schematic representation | An explicitly simplified educational shape or arrangement | It is not an experimentally determined molecular structure |
| Catalog-only record | An identity, family, or localization entry without an importable coordinate asset | **A catalog record is not a 3D structure.** A missing coordinate asset must remain missing |

Integrative structures can contain atomic components, coarse particles, multiple states, and ensembles. A PDB identifier alone does not establish that every part has experimental atomic support. See [RCSB, Integrative Structures](https://www.rcsb.org/docs/general-help/integrative-structures-on-rcsborg).

## Identity and assembly review

For a macromolecular asset, review the source accession and version, organism for each entity, polymer sequence, mutations, construct truncation, bound partners, and the selected biological assembly. Expression in a bacterial host does not by itself make a human protein sequence bacterial; record source organism and expression system separately when available.

The crystallographic asymmetric unit may contain only part of an oligomer or more than one copy of it. It must not be substituted automatically for the biologically relevant assembly. Record the chosen assembly and any symmetry operations used to generate it. If the biological assignment remains uncertain, preserve that uncertainty.

Assembly membership also matters for population accounting. A complex and its constituent chains are not independent extra populations. Domain constructs representing parts of one protein must not silently become several full proteins.

Source: [PDB-101, Biological Assemblies](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies).

## Coverage and missing coordinates

Record the sequence or construct that was deposited and the part actually modeled. Distinguish:

- A full biological sequence from a shortened experimental construct.
- Residues lacking coordinates from residues modeled with incomplete atoms.
- Deposited heavy atoms from an all-atom model that also contains hydrogens.
- Alternate conformers and their occupancies from independent copies of a molecule.
- A full atomic representation from a backbone or alpha-carbon-only trace.

Do not fill missing tails, loops, hydrogens, or domains with apparently measured coordinates. Any later reconstruction must be a separately identified modeling step. A missing measurement remains `null` or explicitly unavailable, never zero, perfect confidence, or an invented value. A low coverage fraction may still permit a useful view of the resolved domain, but it cannot support a claim that the displayed domain is the whole protein.

Source: [PDB-101, Missing Coordinates](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/missing-coordinates).

## Method-specific quality evidence

Archive validation reports assess model geometry, agreement with experimental data, or both. Preserve the original metrics and report link. A report being available or downloaded is different from someone reviewing its findings. Metadata inspection alone is not a full scientific review of maps, restraints, and local model quality.

| Method | Useful evidence to retain and review |
| --- | --- |
| X-ray crystallography | Resolution; R-work and R-free; clashes and geometric outliers; local density support such as RSRZ or RSCC; ligand fit; occupancy and alternate conformers |
| Cryo-EM | EMDB map accession; reported resolution and FSC evidence; map-model fit; local support, atom inclusion, and Q-scores where available; missing or weakly resolved regions |
| NMR | Ensemble identity and selected display member; chemical-shift data; distance/dihedral restraints and violations; well-defined versus uncertain regions |
| Predicted protein | Prediction method/version; target sequence; local confidence such as pLDDT; relative-domain or chain uncertainty such as PAE; relevant interface confidence when applicable |
| Integrative assembly | Input component models; restraints; representation scale; alternative states or ensemble; agreement with the particular input data |

No single resolution cutoff or invented composite score can establish suitability for every use. Review the region and question of interest. An overall surface can remain useful when a local side-chain interaction is poorly supported.

For AlphaFold models, pLDDT measures local prediction confidence. Even high-confidence domains can have uncertain relative positions, which requires examining PAE. **A pLDDT value of 90 does not mean “90% scientifically accurate.”** Predicted coordinates remain predicted even when confidence is high.

Sources:

- [RCSB, Assessing the Quality of 3D Structures](https://www.rcsb.org/docs/general-help/assessing-the-quality-of-3d-structures)
- [wwPDB, Validation Reports](https://www.wwpdb.org/validation/validation-reports)
- [Restraint Validation of Biomolecular Structures Determined by NMR in the Protein Data Bank](https://pmc.ncbi.nlm.nih.gov/articles/PMC11162339/)
- [EMBL-EBI, Understanding pLDDT](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/)
- [EMBL-EBI, Understanding PAE](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/pae-a-measure-of-global-confidence-in-alphafold-predictions/)

## Small molecules, ions, and chemical state

Check chemical identifiers, atom identity, bond orders, stereochemistry, formal charge, protonation, and tautomer where applicable. Keep formal chemical charge distinct from force-field partial charges and from a loader's placeholder charge.

The wwPDB Chemical Component Dictionary supplies chemical definitions and idealized coordinates. Those idealized coordinates are generated computationally. They must be labeled accordingly, even when the chemical identity is authoritative. The existence of a CCD component does not establish its abundance, predominant protonation state, or free conformation inside a nucleus.

An experimentally bound ligand pose is specific to its source complex. Review local experimental support and geometry when using that pose. Crystallization additives and buffer components must not automatically become nuclear constituents. Conversely, do not discard all heteroatom records: relevant cofactors, ions, and modified residues may use them.

Sources:

- [wwPDB, Chemical Component Dictionary](https://www.wwpdb.org/data/ccd)
- [RCSB, Ligand Structure Quality](https://www.rcsb.org/docs/general-help/ligand-structure-quality-in-pdb-structures)

## Import fidelity and reproducibility

Preserve the unmodified source where downloaded and record its checksum. Pin the source accession, relevant version or revision, assembly, and conversion method. A checksum establishes the identity of the downloaded bytes; it does not establish scientific accuracy.

Conversion checks should establish that coordinates are finite, atom identities and source selections are consistent, units are explicit, and transformations preserve physical distances. Source coordinates are conventionally in angstroms; conversion to nanometres uses 1 angstrom = 0.1 nanometres. Record any centering transform, selected ensemble member, omitted entities, reconstructed atoms, or inferred display bonds.

Coordinate, residue, and chain counts must be interpreted for the selected assembly and representation. An archive entry-level atom count may differ from the biological assembly count or the installed display count. Such a difference needs an explanation rather than an automatic scientific verdict.

Imported structures are for visualization. Coordinates alone do not supply a validated force field, chemical state, solvent environment, or realistic dynamics. Rendering detail and scientific evidence are separate properties.

## Review labels

Useful review outcomes include **not reviewed**, **reviewed with stated limitations**, and **quarantined with a reason**. These are workflow states, not universal accuracy grades. Keep provenance class, coverage, scientific review, import verification, and biological-context evidence separate.

Only report a check as completed when it was actually performed. Source metadata review, downloading an archive validation report, numerical import checks, and a detailed local map/model review are different activities.

Automated tests are opt-in for this project and run only when explicitly requested. Preparing tests or adding a record does not mean those tests have been executed. Application startup must not automatically run the test suite.
