![License](https://img.shields.io/github/license/pdbeurope/rellig) ![Documentation](https://github.com/PDBeurope/rellig/workflows/documentation/badge.svg)

# PDBe RelLig

## Relevant Ligands in PDB

With over 200,000 entries in the PDB, about 75% of these structures contain at least one ligand bound to a protein or nucleic acid. However, not all ligands are biologically relevant, some are present due to experimental necessities, such as aiding crystallisation or enabling cryoprotection, while others play biologically significant roles acting as cofactors, reactants or drugs. Unfortunately, the biological role of ligands present in PDB entries is not annotated in the PDB/mmCIF files.
PDBe RelLig is designed to bridge this gap by automatically annotating the ligand's functional role as the following:
 - cofactor
 - reactant
 - drug

### Installation

Create and activate a virtual environment and then install PDBe RelLig

```bash
$ pip install pdberellig
```

### Running the pipeline

There are three modes of pipelines:

### cofactors

```bash
$ pdberellig cofactors --cif <path_to_ligand_cif_file> --ligand_type <type_of_ligand> --out-dir <path_to_output>
```
**pipeline inputs**

```bash
Options:
  --cif TEXT                   path to input cif file  [required]
  --ligand-type [CCD|PRD|CLC]  type of ligand in the PDB  [required]
  --out-dir TEXT               path to output directory  [required]
  --help                       Show this message and exit.
```
**pipeline outputs**

<ligand_id>_cofactor_annotation.json - a json file containing ligand interacting proteins and similarity to template and representative molecules of cofactor classes .

Example

```json
{
    "HEM": {
        "template": {
            "id": "HEA",
            "similarity": 0.717
        },
        "representative": {
            "id": "HEA",
            "similarity": 0.717
        },
        "pdb_chains": [
            {
                "pdb_id": "3ks0",
                "auth_asym_id": "A",
                "struct_asym_id": "C",
                "uniprot_id": "P00175",
                "ec_number": "1.1.2.3"
            },
            {
                "pdb_id": "1ltd",
                "auth_asym_id": "A",
                "struct_asym_id": "A",
                "uniprot_id": "P00175",
                "ec_number": "1.1.2.3"
            },
        ]
    }
}
```
---

### reactants

```bash
$ pdberellig reactants --cif <path_to_ligand_cif_file> --ligand_type <type_of_ligand> --chebi-structure-file <csv_file_with_chebi_mol> --out-dir <path_to_output>
```
**pipeline inputs**
```bash
Options:
  --cif TEXT                     path to input cif file  [required]
  --ligand-type [CCD|PRD|CLC]    type of ligand in the PDB  [required]
  --chebi-structure-file TEXT    Path to the ChEBI SDF file  [required]
  --out-dir TEXT                 path to output directory  [required]
  --update-chebi                 Path to the ChEBI archive files
  --minimal-ligand-size INTEGER  Minimum ligand size.  [default: 5]
  --help                         Show this message and exit.
```

**pipeline outputs**

<ligand_id>_reactant_annotation.tsv - A tsv file containing ligand interacting proteins and similarity to the reaction participants present in the reactions catalysed the protein

Example

| pdb_id | auth_asym_id | struct_asym_id | uniprot_id | rhea_id | chebi_id | similarity |
| --- | --- | --- | --- | --- | --- | --- |
|1r3q | A | A | P06132 | 19865 | 57308 | 0.714
|1r3q | A | A | P06132 | 31239 | 62626 | 0.8
|1r3q | A | A | P06132 | 31239 | 62631 | 1.0
|1r3s | A | A | P06132 | 19865 | 57308 | 0.714
|1r3s | A | A | P06132 | 31239 | 62626 | 0.8
|1r3s | A | A | P06132 | 31239 | 62631 | 1.0
|1r3v | A | A | P06132 | 19865 | 57308 | 0.714
|1r3v | A | A | P06132 | 31239 | 62626 | 0.8
|1r3v | A | A | P06132 | 31239 | 62631 | 1.0
---


### drugs

```bash
$ pdberellig drugs --cif <path_to_ligand_cif_file> --ligand_type <type_of_ligand> --out-dir <path_to_output>
```
**pipeline inputs**

```bash
Options:
  --cif TEXT                   path to input cif file  [required]
  --out-dir TEXT               path to output directory  [required]
  --ligand-type [CCD|PRD|CLC]  type of ligand in the PDB  [required]
  --help                       Show this message and exit.
```
**pipeline outputs**

<ligand_id>_drug_annotation.tsv - a tsv file containing ligand interacting pharmacologically active drug-targets

Example

| pdb_id | auth_asym_id | struct_asym_id | uniprot_id | name | organism |
| --- | --- | --- | --- | --- | --- |
| 7n9g | A | A | P00519 | Tyrosine-protein kinase ABL1 | Humans |
| 3pyy | B | B | P00519 | Tyrosine-protein kinase ABL1 | Humans |
| 6npu | B | B | P00519 | Tyrosine-protein kinase ABL1 | Humans |
| 6npe | A | A | P00519 | Tyrosine-protein kinase ABL1 | Humans |
| 2hyy | A | A | P00519 | Tyrosine-protein kinase ABL1 | Humans |
| 6npu | A | A | P00519 | Tyrosine-protein kinase ABL1 | Humans |
| 6npv | B | B | P00519 | Tyrosine-protein kinase ABL1 | Humans |
| 6npe | B | B | P00519 | Tyrosine-protein kinase ABL1 | Humans |
