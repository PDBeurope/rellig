```{eval-rst}
.. _cofactors
```
# Cofactors

The PDBe RelLig pipeline automatically identifies cofactor-like ligands in the PDB by comparing their 2D structural similarity to cofactor classes found in the [CoFactor database](https://www.ebi.ac.uk/thornton-srv/databases/CoFactor/). The CoFactor database contains 27 manually curated classes of organic enzyme cofactors and information about the associated enzymes, including their EC numbers. For each cofactor class, a representative small molecule was selected from the PDB based on its close structural match to the template molecule using [PARITY similarity](https://doi.org/10.1016/j.str.2018.02.009). A minimum similarity threshold was also defined, and the list of enzyme EC numbers was expanded using data from BRENDA.

## Annotation process

The pipeline analyzes input PDB ligands by calculating their similarity to template molecules in each cofactor class. If the similarity meets the minimum threshold for any cofactor class, the ligand is further compared to the representative molecule of that class. If the similarity remains above the threshold and the ligand interacts with a protein that has an EC number matching the cofactor class, the ligand is classified as [cofactor-like](https://doi.org/10.1093/bioinformatics/btz115).




## Data

Details of templates, representatives, thresholds and ec numbers used for the cofactor classes can be found in the folder ``pdberellig/data/cofactors``.

| Class id     | Template     | Representative | Threshold |
|----------    |----------    |:--------------:|----------:|
|1          |TDPTHW     |[TPP](https://pdbe.org/chem/TDP)             |0.6        |
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




## Examples

PDBe RelLig pipeline identifies the PDB Ligand 8FL as a cofactor-like molecule similar to TPP (Thiamine Diphosphate) and writes the outputs with ligand interacting proteins from the PDB with corresponding EC numbers.

```json
{
    "8FL": {
        "template": {
            "id": "TDPTHW",
            "similarity": 0.722
        },
        "representative": {
            "id": "TPP",
            "similarity": 0.963
        },
        "pdb_chains": [
            {
                "pdb_id": "5xvt",
                "auth_asym_id": "A",
                "struct_asym_id": "A",
                "uniprot_id": "P34736",
                "ec_number": "2.2.1.1"
            },
            {
                "pdb_id": "5xuf",
                "auth_asym_id": "A",
                "struct_asym_id": "A",
                "uniprot_id": "P34736",
                "ec_number": "2.2.1.1"
            },
            {
                "pdb_id": "8vzb",
                "auth_asym_id": "D",
                "struct_asym_id": "D",
                "uniprot_id": "A0A3C0TX30",
                "ec_number": "4.1.1.8"
            },
            {
                "pdb_id": "8vzb",
                "auth_asym_id": "A",
                "struct_asym_id": "A",
                "uniprot_id": "A0A3C0TX30",
                "ec_number": "4.1.1.8"
            },
            {
                "pdb_id": "8vza",
                "auth_asym_id": "C",
                "struct_asym_id": "C",
                "uniprot_id": "A0A3C0TX30",
                "ec_number": "4.1.1.8"
            },
            {
                "pdb_id": "8vza",
                "auth_asym_id": "B",
                "struct_asym_id": "B",
                "uniprot_id": "A0A3C0TX30",
                "ec_number": "4.1.1.8"
            },
            {
                "pdb_id": "8vzb",
                "auth_asym_id": "B",
                "struct_asym_id": "B",
                "uniprot_id": "A0A3C0TX30",
                "ec_number": "4.1.1.8"
            }
        ]
    }
}

```
