```{eval-rst}
.. _cofactors
```
# Cofactors

The aim of the project is to incorporate cofactor information from [Cofactor database](https://www.ebi.ac.uk/thornton-srv/databases/CoFactor/) into PDBe infrastructure and provide weekly updates to the functional annotation of PDB entries.

The cofactor information is divided into 27 different classes. Each class contains a list of HET codes (ligands) that are members of the given cofactor class and a list of EC numbers which represent enzymes where these ligands act as cofactors.

There is also the [paper](https://doi.org/10.1093/bioinformatics/btz115) which describes the process and a few related things.

## Database content

Oracle database contains ``COFACTOR`` user that has access to the cofactor information in the following tables:

* ``COFACTOR`` - Contains just general description taken from the Cofactor database.
* ``COFACTOR_CLASS`` - Contains db mapping between cofactor class name and id
* ``COFACTOR_GROUP`` - No clue.
* ``COFACTOR_GROUP_REF`` - No clue.
* ``EC_COFACTOR`` - Mapping between cofactor class and EC numbers.
* ``HET_GROUPS`` - Mapping between cofactor class and HET codes. The parity score given is a similarity between a ligand and cofactor class representative not the template!!

## Annotation process

First, all the ligands are compared with templates for a given cofactor class using PARITY method implemented as a part of [pdbecccdutils](https://pdbe.gitdocs.ebi.ac.uk/ccdutils/) package. Ligands that can be part of polymers (that is `_chem_comp.type` !=  `non-polymer`) and those that can be subcomponents of other ligands are discarded.  If the similarity score is equal to, or above threshold defined for a given class, the ligand is further compared to the cofactor class representative. If the similarity score does not drop below defined threshold and the ligand is found in a PDB entry with EC number approved for the ligand class, the ligand is immediately marked as cofactor and CSV generating table `HET_GROUPS` is updated.

If the ligand's similarity score to cofactor class representative drops or the ligand is found in the PDB entry with EC number not found in the cofactor database a manual intervention is necessary. The process for manual intervention is described below. Also, cofactor classes, templates and their thresholds made by Abhik can be found in the table below:

| Class id     | Template     | Representative | Threshold |
|----------    |----------    |:--------------:|----------:|
|1          |TDPTHW     |[TDP](https://pdbe.org/chem/TDP)             |0.6        |
|2          |FAD        |[FAD](https://pdbe.org/chem/FAD)             |0.87       |
|3          |FMN        |[FMN](https://pdbe.org/chem/FMN)             |0.86       |
|4          |NAD2       |[NAD](https://pdbe.org/chem/NAD)             |0.68       |
|5          |PNS        |[PNS](https://pdbe.org/chem/PNS)             |1.0        |
|6          |COAM       |[COA](https://pdbe.org/chem/COA)             |0.67       |
|7          |PLP        |[PLP](https://pdbe.org/chem/PLP)             |0.78       |
|8          |0HG        |[GSH](https://pdbe.org/chem/GSH)             |0.42       |
|9          |BTN        |[BTN](https://pdbe.org/chem/BTN)             |0.58       |
|10         |FFO        |[FFO](https://pdbe.org/chem/FFO)             |0.94       |
|11         |B12        |[B12](https://pdbe.org/chem/B12)             |0.84       |
|12         |ASC        |[ASC](https://pdbe.org/chem/ASC)             |1.0        |
|13         |MQ7        |[MQ7](https://pdbe.org/chem/MQ7)             |1.0        |
|14         |UQ2        |[UQ1](https://pdbe.org/chem/UQ1)             |0.53       |
|15         |MSS        |[MGD](https://pdbe.org/chem/MGD)             |0.43       |
|16         |H4B        |[H4B](https://pdbe.org/chem/H4B)             |0.78       |
|17         |MDO        |[MDO](https://pdbe.org/chem/MDO)             |1.0        |
|18         |SAM3       |[SAM](https://pdbe.org/chem/SAM)             |0.6        |
|19         |F43        |[F43](https://pdbe.org/chem/F43)             |0.94       |
|20         |COM        |[COM](https://pdbe.org/chem/COM)             |1.0        |
|21         |TP7        |[TP7](https://pdbe.org/chem/TP7)             |1.0        |
|22         |HEA        |[HEA](https://pdbe.org/chem/HEA)             |0.52       |
|24         |DPM        |[DPM](https://pdbe.org/chem/DPM)             |0.97       |
|25         |PQQ        |[PQQ](https://pdbe.org/chem/PQQ)             |1.0        |
|26         |1TY        |[TPQ](https://pdbe.org/chem/TPQ)             |0.55       |
|27         |TRQ        |[TRQ](https://pdbe.org/chem/TRQ)             |0.88       |
|28         |LPA        |[LPA](https://pdbe.org/chem/LPA)             |1.0        |


Templates used for the cofactor similarity perception can be found as a part of ``relic`` package in the folder ``relic/data/cofactors/templates``.

## Manual annotation process

Whenever the pipeline cannot make an automatic guess on whether or not the ligand under questions is a cofactor a human intervention is necessary. This can happen in a certain number of cases:

During the process you may need to modify either of these files:

* ``/nfs/msd/release/data/cofactor/cofactors.backlog.json`` - This file contains all the ligands which needs to be either integrated in the database or scraped. After you process a ligand, please delete the entry from JSON file.

* ``/nfs/msd/release/data/cofactor/het_groups.csv`` - Modify this file to add het codes to the cofactor database. The Parity score represents a similarity of het code under question and cofactor class representative.

* ``/nfs/msd/release/data/cofactor/ec_cofactor.csv`` - Modify this file to add ec numbers to the cofactor database.

### EC number mismatch

The new ligand has a similarity to both template and representative, however, it is not found in the pdb entry with mapping to known EC classes for this cofactor class. Please verify that the ligand acts as a cofactor and add enzyme EC number to the database.

### Similarity score decreased

In certain cases the similarity score for ligand-representative is significantly lower than for ligand-template. These cases are generally interesting as it may mean that either representative or threshold was selected incorectly and needs special attention. Please contact ``Lukas Pravda (lpravda@ebi.ac.uk)`` who maintans the pipeline for further information.

## Examples

Cofactor class 1 - [TDP](https://pdbe.org/chem/TDP) cofactor class representative (framed) some other cofactor class representatives as identified by the pipeline: [8EL](https://pdbe.org/chem/8EL) (similarity: 1.0); [8EF](https://pdbe.org/chem/8EF) (similarity: 1.0); or [8FL](https://pdbe.org/chem/8FL) (similarity: 0.963);

<div align='center'>
    <img style="border: 2px solid black" alt="TDP" src='https://www.ebi.ac.uk/pdbe-srv/pdbechem/image/showNew?code=TDP&size=300' />
    <img alt="8EL" src='https://www.ebi.ac.uk/pdbe-srv/pdbechem/image/showNew?code=8EL&size=300' />
    <img alt="8EF" src='https://www.ebi.ac.uk/pdbe-srv/pdbechem/image/showNew?code=8EF&size=300' />
    <img alt="8FL" src='https://www.ebi.ac.uk/pdbe-srv/pdbechem/image/showNew?code=8FL&size=300' />

</div>
