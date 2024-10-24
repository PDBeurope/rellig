
# Drugs

The PDBe RelLig pipeline identifies drug-like ligands in the PDB by utilizing drug-target annotations found in the ligand CIF files available from the [PDBe FTP](https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/). These CIF files are updated weekly and include drug-target data sourced from the [DrugBank database](https://go.drugbank.com/).


## Annotation process

The PDBe RelLig pipeline parses the input PDB ligand CIF file to retrieve any available drug-target annotations. It then uses the [PDBe API](https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/) to gather all proteins interacting with the ligand. If any of these proteins are pharmacologically active targets, the ligand is annotated as drug-like.

## Example

PDBe RelLig pipeline identifies STI (Imatinib) as drug-like molecule with drug-target Tyrosine-protein kinase ABL1 and write the output with ligand interacting proteins from the PDB

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
