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

import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count
from typing import Union

import pandas as pd
from rdkit import Chem

from pdberellig.conf import get_config
from pdberellig.core.models import CompareObj
from pdberellig.helpers.utils import (
    download_chebi,
    get_ligand_intx_chains,
    parse_ligand,
    sparql_to_df,
)


class Reactants:
    """Reactants pipeline data model."""

    def __init__(self, log, args):
        self.log = log
        self.args = args

    def process_entry(self) -> None:
        """
        Runs the pipeline for reactant-like annotation of the ligand
        and writes the results to tsv files

        * parses ligand cif file
        * fetches all ligand interacting PDB chains and corresponding uniprot ids
        * calculates similarity to reactant participants
        """

        component = parse_ligand(self.args.cif, self.args.ligand_type)
        ligand_id = component.id
        if len(component.mol_no_h.GetAtoms()) < self.args.minimal_ligand_size:
            # ligand is too small.
            self.log.debug(f"""Number of atoms in {ligand_id} is less than
                           {self.args.minimal_ligand_size},
                           hence skipping""")
            return

        ligand = CompareObj(component.id, component.mol_no_h)

        # get ligand interacting PDB chains
        intx_chains = get_ligand_intx_chains(ligand.id)
        if intx_chains.empty:
            self.log.warn(f"No interacting PDB chain was found for {ligand.id}")
            return

        uniprot_ids = intx_chains.get("uniprot_id").to_list()

        # get similarity to reaction participants
        reactants_sim = self.get_reactant_annotation(ligand, uniprot_ids)
        if not reactants_sim.empty:
            intx_chains_reac_sim = pd.merge(
                intx_chains[["pdb_id", "auth_asym_id", "struct_asym_id", "uniprot_id"]],
                reactants_sim,
                on="uniprot_id",
                how="inner",
            )
            reactants_sim_file = os.path.join(
                self.args.out_dir, f"{ligand.id}_reactant_annotation.tsv"
            )
            self.log.info(f"Writing reactant annotations to {reactants_sim_file}")
            intx_chains_reac_sim.to_csv(reactants_sim_file, sep="\t", index=False)

    def get_reactant_annotation(
        self, ligand: CompareObj, uniprot_ids: list[str]
    ) -> pd.DataFrame:
        """
        Returns PARITY similarity of the ligand to all reaction participants corresponding to the
        input list of uniprot ids
        * Fetches rhea_ids corresponding to all the reactions participated
        by the list of proteins (uniprot_ids)
        * Fetches ChEBI ids of reaction participants present in the list of
        reactions
        * Calculates PARITY similarity of input ligand to the list of ChEBI molecules

        Args:
            ligand: CompareObj of ligand
            uniprot_ids: list of uniprot_ids corresponding to the ligand interacting PDB chains
        """

        # get reactions corresponding to the list of uniprot_ids
        reactions = self.get_reactions(uniprot_ids)
        if reactions.empty:
            self.log.warn(f"No reaction was fetched from Uniprot for {uniprot_ids}")
            return pd.DataFrame()

        rhea_ids = reactions.get("rhea_id").to_list()

        # get ChEBI ids of all the reaction participants
        reaction_participants_df = self.get_reaction_participants(rhea_ids)
        if reaction_participants_df.empty:
            self.log.warn(
                f"No reaction participants was fetched from Rhea for {rhea_ids}"
            )
            return pd.DataFrame()

        chebi_ids = reaction_participants_df.get("chebi_id").to_list()

        # add the structures of reaction particpants to templates
        templates = self.parse_chebi(chebi_ids)

        # calculate similarity of input ligand to all the reaction particiapnts
        chebi_similarities = self.get_similarities(ligand, templates)
        chebi_similarities_df = pd.DataFrame.from_dict(chebi_similarities)
        if chebi_similarities_df.empty:
            self.log.info(f"No similar reaction participant was found for {ligand.id}")
            return pd.DataFrame()

        reaction_chebi = pd.merge(
            reactions, reaction_participants_df, on="rhea_id", how="inner"
        )

        reaction_chebi_sim = pd.merge(
            reaction_chebi, chebi_similarities_df, on="chebi_id", how="inner"
        )

        return reaction_chebi_sim

    def parse_chebi(self, chebi_ids: list[str]) -> list[CompareObj]:
        """
        Parse ChEBI mol files from chebi_structure_file
        using RDKit and return as a list of CompareObj. If the
        update option is enabled new structure file is downloaded from ChEBI FTP.

        Args:
         chebi_ids: list of chebi ids

        Returns:
            a list of CompareObj

        """

        templates = []
        if not self.args.update_chebi:
            chebi_structure_file = self.args.chebi_structure_file
        else:
            chebi_structure_file = download_chebi(self.args.out_dir)

        self.chebi = pd.read_csv(chebi_structure_file, dtype=str)
        self.chebi = self.chebi.loc[
            (self.chebi["COMPOUND_ID"].isin(chebi_ids))
            & (self.chebi["TYPE"] == "mol")
            & (self.chebi["DEFAULT_STRUCTURE"] == "Y"),
            ["COMPOUND_ID", "STRUCTURE"],
        ]
        for _, row in self.chebi.iterrows():
            try:
                chebi_mol = Chem.MolFromMolBlock(row["STRUCTURE"])
                chebi_mol_no_h = Chem.RemoveHs(chebi_mol)
                if len(chebi_mol_no_h.GetAtoms()) < self.args.minimal_ligand_size:
                    # chebi is too small.
                    self.log.debug(f"""Number of atoms in {row["COMPOUND_ID"]} is less
                                    than {self.args.minimal_ligand_size},
                                    hence skipping""")
                else:
                    templates.append(CompareObj(row["COMPOUND_ID"], chebi_mol_no_h))

            except Exception:
                self.log.warn(f"Couldn't parse {row['COMPOUND_ID']} using RDKit")

        return templates

    def get_similarities(
        self, query: CompareObj, templates: list[CompareObj]
    ) -> dict[str, list[Union[str, float]]]:
        """
        Returns PARITY similarity of query molecules to ChEBI
        structures present in the list of templates

        Args:
            query: query molecule as a CompareObj
            templates: list of COmpareObj of template molecules to compare
        """
        chebi_similarities = defaultdict(list)
        threshold = float(get_config("main", "reactants_threshold"))
        with ThreadPoolExecutor(max_workers=cpu_count() - 1) as exec:
            future_to_result = {
                exec.submit(template.similarity_to, query, threshold): template.id
                for template in templates
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
                    if template_sim.result.similarity_score >= threshold:
                        chebi_similarities["chebi_id"].append(template_sim.target_id)
                        chebi_similarities["similarity"].append(
                            round(template_sim.result.similarity_score, 3)
                        )

                except Exception as exc:
                    self.log.warn("%r generated an exception: %s" % (template, exc))

        return chebi_similarities

    def get_reactions(self, uniprot_ids: list[str]) -> pd.DataFrame:
        """
        Fetches rhea_ids of all the reactions corresponding to the input
        list of uniprot_ids using Uniprot sparql endpoint

        Args:
            uniprot_ids: a list of uniprot_ids

        Returns:
            a datarame of uniprot_ids and rhea_ids
        """

        # Refer https://sparql.uniprot.org/.well-known/sparql-examples/ for examples of
        # saprql queries in uniprot
        sparql_uniprot_url = "https://sparql.uniprot.org/sparql"
        proteins = " ".join(
            [f"(uniprotkb:{uniprot_id})" for uniprot_id in uniprot_ids if uniprot_id]
        )

        reaction_query = (
            """
            #endpoint: https://sparql.uniprot.org/sparql
            #query: retrieve all reactions corresponding to a list of proteins.

            PREFIX up: <http://purl.uniprot.org/core/>
            PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>

            SELECT DISTINCT
            ?uniprot_id
            ?rhea_id

            WHERE {
            VALUES (?protein) {
                """
            + proteins
            + """
            }

            ?protein up:annotation ?annotation .
            ?annotation up:catalyticActivity ?ca .
            ?ca up:catalyzedReaction ?reaction .
            BIND(SUBSTR(STR(?protein),33) AS ?uniprot_id) .
            BIND(SUBSTR(STR(?reaction),24) AS ?rhea_id) .

            }
            """
        )

        reactions = sparql_to_df(reaction_query, sparql_uniprot_url)

        return reactions

    def get_reaction_participants(self, rhea_ids: list[str]) -> pd.DataFrame:
        """
        Fetches ChEBI ids of all the reaction participants
        corresponding to the input list of rhea_ids using
        Rhea sparql endpoint

        Args:
            rhea_ids: a list of rhea_ids

        Returns:
            a dataframe of rhea_ids and chebi_ids
        """
        sparql_rhea_url = "https://sparql.rhea-db.org/sparql"
        reaction_uris = " ".join([f"(rh:{rhea_id})" for rhea_id in rhea_ids])
        chebi_query = (
            """
            #endpoint: https://sparql.rhea-db.org/sparql
            #query: retrieve all chebi small molecules participating in the reactions

            PREFIX rh:<http://rdf.rhea-db.org/>

            SELECT DISTINCT
            ?rhea_id
            ?chebi_id

            WHERE {
            VALUES (?reaction) {
                """
            + reaction_uris
            + """
            }

            ?reaction rdfs:subClassOf rh:Reaction .
            ?reaction rh:status rh:Approved .
            ?reaction rh:accession ?rhea_accession .
            ?reaction rh:side/rh:contains/rh:compound ?compound .
            ?compound rdfs:subClassOf rh:SmallMolecule .
            ?compound rh:accession ?chebi_accession .
            BIND(SUBSTR(STR(?rhea_accession),6) AS ?rhea_id) .
            BIND(SUBSTR(STR(?chebi_accession),7) AS ?chebi_id) .
            }
            """
        )
        reaction_participants = sparql_to_df(chebi_query, sparql_rhea_url)
        return reaction_participants
