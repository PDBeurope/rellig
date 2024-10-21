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

```bash
Options:
  --cif TEXT                   path to input cif file  [required]
  --ligand-type [CCD|PRD|CLC]  type of ligand in the PDB  [required]
  --out-dir TEXT               path to output directory  [required]
  --help                       Show this message and exit.
```

<ligand_id>_cofactor_annotation.json - A json file containing ligand interacting PDB chains and similarity to template and representative molecules of cofactor classes .

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

<ligand_id>_reactant_annotation.tsv - A tsv file containing

---

### drugs

```bash
$ pdberellig drugs --cif <path_to_ligand_cif_file> --ligand_type <type_of_ligand> --out-dir <path_to_output>
```
```bash
Options:
  --cif TEXT                   path to input cif file  [required]
  --out-dir TEXT               path to output directory  [required]
  --ligand-type [CCD|PRD|CLC]  type of ligand in the PDB  [required]
  --help                       Show this message and exit.
```
