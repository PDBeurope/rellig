```{eval-rst}
.. _reactants
```
# Reactants

The aim of the reactants pipeline is to map Rhea-annotated protein reactants to all chemical components (ligands, small molecules and monomers) in PDB entries sourced from PDBeChem.

[Rhea](https://www.rhea-db.org/) is an expert curated resource of biochemical reactions, it uses the [ChEBI](https://www.ebi.ac.uk/chebi/) (Chemical Entities of Biological Interest) ontology of small molecules to precisely describe reactions participants and their chemical structures. ChEBI compounds from Rhea constitute the first dataset to compare.

The second dataset to compare is the Chemical Component Dictionaries from [PDBeChem](https://www.ebi.ac.uk/pdbe-srv/pdbechem/).

This pipeline incorporates both annotated and similar reactants information into the PDBe infrastructure and provide weekly updates to the functional annotation of PDB entries.

## Database content

Reactants data is only loaded weekly to the Graph database, and not the Oracle.

Note that for ChEBI.csv, only ID's that have not been priorly loaded by external resource mappings are loaded to avoid clashes.

CSVs are produced for the current run (referred hereafter as current-run CVSs) in folder working_directory/date_run. The current-run CSVs are then sort uniq'd and **consolidated with Master CSV's** in folder working_directory/csv. On init mode, Master CSVs are simply sort uniq'd version of the current-run CSVs. On update mode, sort uniq'd version of current-run CSVs are combined with Master CSV's from the previous run. Master CSV's are then copied to the Neo4J loading folder.

#### (UniProt) <- [:REACTANT\_OF] - (Chemical\_Component)
**Chemical\_Component\_UniProt\_Rels.csv**  *:START_ID(ChemicalComponent),:END_ID(UniProt)*

#### (UniProt) - [:HAS\_REACTION] -> (Reaction).
**Reaction.csv** *ID:ID(Reaction),TEXT*
**UniProt\_Reaction\_Rels.csv** *:START_ID(UniProt),:END_ID(Reaction)*

#### (Reaction) - [:CONTAINS] -> (ChEBI)
**ChEBI.csv** *ID:ID(ChEBI)*
**Reaction\_ChEBI\_Rels.csv**  *:START_ID(Reaction),:END_ID(ChEBI)*

#### (Chemical\_Component) - [:SIMILAR_TO] -> (ChEBI)
**Chemical\_Component\_ChEBI\_Rels.csv**  *:START_ID(ChemicalComponent),:END_ID(ChEBI),SCORE*

## Pipeline workflow

1. Parse [Rhea db](https://www.rhea-db.org) content from UNIPROT db.
2. Extract [ChEBI](https://www.ebi.ac.uk/chebi/) compounds which are found in the Rhea db.
3. Analyze PDBeChem and Chebi compounds and fill queues for compounds with more than `N` heavy atoms.
4. Run comparisons either in init, update or file mode.

## Modes
### Update mode

1. First comparison is done between new PDBeChem compounds (ligands\_this\_week\_FILE) and all ChEBI compounds.
2. Second comparison is done between all new ChEBI compounds and all PDBeChem compounds.

### Init mode

1. Comparisons are made between all PDBeChem and ChEBI compounds

#### File mode

1. First comparison is done between compounds in the file and existing ChEBI compounds (in working_directory/chebi folder).
2. Second comparison is done between all new ChEBI compounds and compounds in the file.

## Example

PDB ID [3tz5](https://www.ebi.ac.uk/pdbe/entry/pdb/3tz5) (UniProt ID Q00972) is mapped to [Rhea ID 17301](https://www.rhea-db.org/reaction?id=17301).

The reaction for this protein is given by:
[3-methyl-2-oxobutanoate dehydrogenase]-L-serine + ATP = [3-methyl-2-oxobutanoate dehydrogenase]-O-phospho-L-serine + ADP + H(+)

Five ChEBI compounds are extracted for this protein:

- [15378](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:15378) H(+)
- [29999](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:29999) L-Serine
- [30616](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:30616) ATP
- [456216](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:456216) ADP
- [83421](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:83421) O-phospho-L-serine

H(+) and L-Serine do not make the cut-off for number of heavy atoms (5 atoms). The pipeline then searches for all ChEBI compounds similar to ATP, ADP and O-phospho-L-serine, and map them back to 3tz5. In Chemical_Component_UniProt_Rels.csv, these are shown to be comprised of a cumulative of 336 reactant-like PDBeChem compounds (note that this list also includes ADP and ATP themselves). To distinguish which compounds map to which ChEBI, refer to Chemical_Component_ChEBI_Rels.csv (for example, TPO maps to O-phospho-L-serine (83421)).

## Known problematic compounds

Some compounds time out or are not processed due to errors. This is not an exhaustive list, but during init run between 30422 CCD compounds with 3756 CHEBI compounds on 10 April 2020, these compounds had partial or full failure (reasons for failure unknown, but most likely because of their size):

Not processed entries in update: ['0TS', '1F2', '1RL', '1RM', '3Q8', '52G', '5WP', 'ARD', 'BL3', 'BY6', 'CE6', 'CGM', 'CWO', 'EL5', 'EM1', 'HFW', 'JAS', 'JRA', 'JS5', 'LHA', 'M1V', 'MBV', 'MJC', 'MT9', 'NRB', 'PVG', 'PVN', 'PXI', 'RAF', 'RBT', 'RIF', 'ROX', 'RPT', 'T17', 'TCG', 'TEL', 'TM6', 'VRB', 'ZIT', 'ZM3'].

Timeout error: ['B8B', 'GXZ', 'GXW', 'XCO', 'GX2', 'HFW', 'WNI', 'JGH'].
