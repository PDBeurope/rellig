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
Cofactors pipeline data model
"""

import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from multiprocessing import cpu_count
from typing import List, Union

import pandas as pd
from pdbeccdutils.core import ccd_reader

from pdberellig.conf import get_config, get_data_dir
from pdberellig.core.models import CofactorSim, CompareObj, Similarity
from pdberellig.helpers.utils import (
    get_ligand_intx_chains,
    init_rdkit_templates,
    parse_ligand,
)


class Cofactors:
    """Cofactors pipeline data model."""

    def __init__(self, log, args) -> None:
        self.log = log
        self.args = args
        self.templates = init_rdkit_templates(
            os.path.join(get_data_dir(), get_config("cofactor", "template_path"))
        )
        self.cofactor_details = self._get_cofactor_details()
        self.cofactor_ec = self._get_cofactor_ec()

    def process_entry(self) -> None:
        """Runs cofactor pipeline to check if a ligand
        acts a cofactor in the PDB

        * Calculates similarity to template molecules of cofactor classes
        * Calculates similarity to representative molecules of cofactor class
        * Checks if EC numbers associated with ligand interacting chains in the PDB
          has any overlap with EC numbers of cofactor classes

        """
        # parse ligand cif file
        component = parse_ligand(self.args.cif, self.args.ligand_type)
        ligand = CompareObj(component.id, component.mol_no_h)

        # get similarity of the ligand to cofactor templates
        cofactor_sim = self.get_similarity(ligand)
        if cofactor_sim:
            representative_score = round(
                cofactor_sim.representative_sim.result.similarity_score, 3
            )
            self.log.info(
                f"Possible new cofactor identified: {cofactor_sim.representative_sim.query_id} similar to"
                f" {cofactor_sim.representative_sim.target_id} representative with score {representative_score}."
            )

            # get ligand interacting PDB chains, uniprot ids and ec numbers
            ligand_intx_chains = get_ligand_intx_chains(cofactor_sim.query_id)
            if ligand_intx_chains.empty:
                self.log.warn(f"No interacting PDB chain was found for {ligand.id}")
                return
            cofactor_template_id = cofactor_sim.template_sim.target_id
            cofactor_id = self.cofactor_details[cofactor_template_id]["id"]
            ligand_cofactor_ec = pd.merge(
                self.cofactor_ec.loc[(self.cofactor_ec["COFACTOR_ID"] == cofactor_id)],
                ligand_intx_chains,
                left_on="EC_NO",
                right_on="ec_number",
            ).drop(columns=["EC_NO"])

            if not ligand_cofactor_ec.empty:
                self.log.info(f"""{cofactor_sim.query_id} is a cofactor-like
                            molecule similar to the cofactor template {cofactor_template_id}""")
                cofactor_results = {
                    cofactor_sim.query_id: {
                        "template": {
                            "id": cofactor_template_id,
                            "similarity": round(
                                cofactor_sim.template_sim.result.similarity_score, 3
                            ),
                        },
                        "representative": {
                            "id": cofactor_sim.representative_sim.target_id,
                            "similarity": round(
                                cofactor_sim.representative_sim.result.similarity_score,
                                3,
                            ),
                        },
                        "pdb_chains": ligand_cofactor_ec[
                            [
                                "pdb_id",
                                "auth_asym_id",
                                "struct_asym_id",
                                "uniprot_id",
                                "ec_number",
                            ]
                        ].to_dict("records"),
                    }
                }
                self._write_cofactor_results(cofactor_results)

    def _write_cofactor_results(self, cofactor_results: dict):
        """Writes the results of cofactor pipeline to a
        json file.

        Args:
            cofactor_results: results of cofactor pipeline as dictionary
        """
        ligand_id = list(cofactor_results.keys())[0]
        cofactor_results_path = os.path.join(
            self.args.out_dir, f"{ligand_id}_cofactor_annotation.json"
        )
        self.log.info(f"Writing cofactor annotations to {cofactor_results_path}")
        with open(cofactor_results_path, "w") as fh:
            json.dump(cofactor_results, fh, indent=4)

    def get_similarity(self, ligand: CompareObj) -> Union[CofactorSim, List]:
        """Returns the similarity of query molecule to template and
        representative molecules of cofactor class if it is above defined threshold

        Args:
            ligand: CompareObj of ligand
        """

        cofactor_sim = None
        with ThreadPoolExecutor(max_workers=cpu_count() - 1) as exec:
            future_to_result = {
                exec.submit(
                    template.similarity_to,
                    ligand,
                    self.cofactor_details[template.id]["threshold"] - 0.01,
                ): template.id
                for template in self.templates
            }

            for future in as_completed(future_to_result):
                template = future_to_result[future]
                try:
                    template_sim = future.result()
                    if not template_sim.result:
                        raise Exception(
                            f"Error occured in comparing {template_sim.target_id} "
                            f"to {template_sim.query_id}"
                        )

                    template_details = self.cofactor_details[template_sim.target_id]

                    if (
                        template_sim.result.similarity_score
                        < template_details["threshold"]
                    ):
                        continue
                    template_score = round(template_sim.result.similarity_score, 3)

                    self.log.info(
                        f"Possible new cofactor identified: {template_sim.query_id} similar to"
                        f" {template_sim.target_id} template with score {template_score}."
                    )

                    representative = self.get_representative(template_details)
                    representative_sim = self.get_representative_similarity(
                        ligand, representative
                    )
                    if (
                        representative_sim.result.similarity_score
                        < template_details["threshold"]
                    ):
                        continue

                    if not cofactor_sim:
                        cofactor_sim = CofactorSim(
                            ligand.id, template_sim, representative_sim
                        )
                    else:
                        prev_sim = (
                            cofactor_sim.template_sim.result.similarity_score
                            + cofactor_sim.representative_sim.result.similarity_score
                        ) / 2
                        curr_sim = (
                            template_sim.result.similarity_score
                            + representative_sim.result.similarity_score
                        ) / 2
                        if curr_sim > prev_sim:
                            cofactor_sim = CofactorSim(
                                ligand.id, template_sim, representative_sim
                            )

                except Exception as exc:
                    self.log.warn("%r generated an exception: %s" % (template, exc))

        return cofactor_sim

    def get_representative(self, template_details: dict) -> CompareObj:
        """Returns the representative molecule of a cofactor class

        Args:
            template_details: threshold and represeative details of cofactor class
        """

        representative_path = os.path.join(
            get_data_dir(),
            get_config("cofactor", "representative_path"),
            f"{template_details['representative']}.cif",
        )
        representative_component = ccd_reader.read_pdb_cif_file(
            representative_path
        ).component
        representative = CompareObj(
            representative_component.id, representative_component.mol_no_h
        )
        return representative

    def get_representative_similarity(
        self, query: CompareObj, representative: CompareObj
    ) -> Similarity:
        """Returns the similarity of query molecule to representative
        molecule of a cofactor class

        Args:
            query: query molecule
            representative: representative molecule
        """

        self.log.info(f"Running similarity to representative" f" {representative.id}")

        representative_sim = representative.similarity_to(query)
        if not representative_sim.result:
            raise Exception(
                f"Error occured in comparing {representative_sim.target_id} "
                f"to {representative_sim.query_id}"
            )

        return representative_sim

    def _get_cofactor_details(self) -> dict:
        """Returns the threshold and representative details of cofactor
        classes as dictionary"""
        path = os.path.join(
            os.path.join(get_data_dir(), get_config("cofactor", "details")),
        )
        with open(path) as f:
            obj = json.load(f)
            return {x["template"]: x for x in obj}

    @lru_cache(maxsize=None)
    def _get_cofactor_ec(self) -> pd.DataFrame:
        """Returns the EC numbers allowed for cofactor classes"""
        path = os.path.join(get_data_dir(), get_config("cofactor", "ec"))
        cofactor_ec = pd.read_csv(path, dtype={"EC_NO": str, "COFACTOR_ID": int})

        return cofactor_ec
