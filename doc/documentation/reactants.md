```{eval-rst}
.. _reactants
```
# Reactants

The PDBe RelLig pipeline identifies reactant-like ligands in the PDB by comparing their 2D structural similarity to reactants in the [Rhea database](https://www.rhea-db.org/), a curated resource that uses the [ChEBI ontology](https://www.ebi.ac.uk/chebi/) to describe reaction participants and their structures.


## Annotation process

The pipeline analyzes input PDB ligands, then retrieves the proteins interacting with these ligands using the [PDBe API](https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/). For each protein, the corresponding catalysis reactions and reaction participants are obtained through the [Uniprot SPARQL](https://sparql.uniprot.org/) and [Rhea SPARQL](https://sparql.rhea-db.org/) endpoints, respectively. The PARITY method is used to compare the structural similarity of PDB ligands to reaction participants represented by the ChEBI ontology. Ligands with a minimum similarity score of 0.7 are classified as reactant-like, ensuring accurate identification of biologically relevant reactants (either the substrate or product of the reaction).

## Example

PDBe RelLig pipeline identifies 1CP (COPROPORPHYRINOGEN I) as reactant-like molecule and writes the output with interacting proteins from the PDB and similarities to reaction participants present in the reactions catalysed the the protein.

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
