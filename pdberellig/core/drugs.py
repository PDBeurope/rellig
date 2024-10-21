#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2024 EMBL - European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on
# an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.

"""
Drugs pipeline data model
"""

import os
from collections import defaultdict

import pandas as pd
from gemmi import cif
from pdbeccdutils.core.component import Component

from pdberellig.helpers.utils import get_ligand_intx_chains, parse_ligand


class Drugs:
    """
    Drugs pipeline data model
    """

    def __init__(self, log, args):
        self.log = log
        self.args = args

    def process_entry(self):
        ligand = parse_ligand(self.args.cif, self.args.ligand_type)
        drugbank_targets = self.get_drugbank_targets(ligand)
        if drugbank_targets.empty:
            self.log.info("No target was found for {ligand.id} from DrugBank")
            return

        drug_targets = drugbank_targets[
            drugbank_targets["pharmacologically_active"] == "yes"
        ]
        if drug_targets.empty:
            self.log.info(
                "None of the targets from DrugBank found for {ligand.id} are pharmacologically active"
            )
            return

        # get ligand interacting PDB chains, uniprot ids and ec numbers
        ligand_intx_chains = get_ligand_intx_chains(ligand.id)
        if ligand_intx_chains.empty:
            self.log.warn(f"No interacting PDB chain was found for {ligand.id}")
            return

        ligand_drug_targets = pd.merge(
            ligand_intx_chains[
                ["pdb_id", "auth_asym_id", "struct_asym_id", "uniprot_id"]
            ],
            drug_targets[["uniprot_id", "name", "organism"]],
            on="uniprot_id",
            how="inner",
        )
        if ligand_drug_targets.empty:
            self.log.info(
                f"None of the pharmacologically active targets of {ligand.id} was found in the PDB"
            )

        ligand_drug_targets_file = os.path.join(
            self.args.out_dir, f"{ligand.id}_drug_annotation.tsv"
        )

        self.log.info(f"Writing drug annotations to {ligand_drug_targets_file}")
        ligand_drug_targets.to_csv(ligand_drug_targets_file, sep="\t", index=False)

    def get_drugbank_targets(self, component: Component) -> pd.DataFrame:
        """
        Returns target information from DrugBank

        Args:
            component: Component object of ligand
        """
        cif_block = component.ccd_cif_block
        if (
            "_pdbe_chem_comp_drugbank_targets."
            not in cif_block.get_mmcif_category_names()
        ):
            return
        items = ["name", "organism", "uniprot_id", "pharmacologically_active"]
        targets = cif_block.find("_pdbe_chem_comp_drugbank_targets.", items)
        targets_info = defaultdict(list)
        for row in targets:
            for item in items:
                targets_info[item].append(
                    cif.as_string(row[f"_pdbe_chem_comp_drugbank_targets.{item}"])
                )
        targets_df = pd.DataFrame.from_dict(targets_info)
        return targets_df
