#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2019 EMBL - European Bioinformatics Institute
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
Various utilities used in the pipeline
"""

import importlib.metadata
import logging
import os
import sys
from collections import defaultdict

import pandas as pd
import pdbeccdutils
import rdkit
import requests
from pdbeccdutils.core import ccd_reader, clc_reader, prd_reader
from pdbeccdutils.core.component import Component
from rdkit import Chem
from requests.adapters import HTTPAdapter
from SPARQLWrapper import JSON, SPARQLWrapper
from urllib3 import Retry

from pdberellig.conf import get_config
from pdberellig.core.models import CompareObj


def setup_log(stage, mode):
    """Set up application log.

    Args:
        stage (str): Stage of the logger in the pipeline hierarchy.
        mode (str): Mode of the application.
    """
    frm = "[%(asctime)-15s] %(name)s:%(levelname)s:%(message)s"

    log = logging.getLogger(mode)
    log.setLevel(logging.DEBUG)

    stream = logging.StreamHandler(sys.stdout)
    stream.setLevel(logging.INFO)
    stream.setFormatter(logging.Formatter(frm))
    log.addHandler(stream)

    log.info(
        f"PDBe RelLig {stage} started | {mode}    v.{importlib.metadata.version('pdberellig')}"
    )
    log.info("Used packages:")
    log.info(f"  rdkit {rdkit.__version__}")
    log.info(f"  pdbeccdutils {importlib.metadata.version('pdbeccdutils')}")

    return log


def init_rdkit_templates(path) -> list[CompareObj]:
    """Returns list of templates for cofactor classes.

    Args:
        path (str): Path to the directory with templates
    """

    templates = [
        CompareObj(x.split(".")[0], Chem.MolFromMolFile(os.path.join(path, x)))
        for x in os.listdir(path)
    ]

    return templates


def get_ids_to_process_from_file(file_name: str):
    """Read in all the items to process from a file.

    Args:
        file_name (str): Path to the file.

    Returns:
        list of str: List of items to be processed.
    """
    with open(file_name, "r") as fp:
        return [x for x in fp.read().splitlines() if x]


def requests_retry_session(
    retries=4,
    backoff_factor=1,
    status_forcelist=(429, 500, 502, 503, 504),
    session=None,
):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )

    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


def get_ligand_intx_chains(ligand_id: str) -> pd.DataFrame:
    """Returns the details of ligand interacting
    PDB chains including chain numbers,uniprot ids
    and ec number by calling PDBe API

    Args:
        ligand_id: ligand identifier from PDB
    """
    api_base_url = get_config("main", "api_base_url")
    intx_chain_ec_url = f"{api_base_url}/pdb/compound/interacting_chains/{ligand_id}"
    try:
        response = requests_retry_session().get(
            intx_chain_ec_url, headers={"Content-Type": "application/json"}
        )
        response.raise_for_status()
        ligand_intx_chains = response.json()
        ligand_intx_chains_df = pd.DataFrame.from_records(
            ligand_intx_chains.get(ligand_id, {})
        )
        if "ec_number" in ligand_intx_chains_df.columns:
            ligand_intx_chains_df["ec_number"] = ligand_intx_chains_df[
                "ec_number"
            ].str.replace(".-", "")
    except requests.exceptions.HTTPError as err:
        if err.response.status_code == 404:
            ligand_intx_chains_df = pd.DataFrame()

    return ligand_intx_chains_df


def sparql_to_df(query, sparql_url: str) -> pd.DataFrame:
    """
    Returns the result of input sparql query as
    a data frame

    Args:
        sparql_url: input sparql query
    """
    sparql = SPARQLWrapper(sparql_url)
    sparql.method = "POST"
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = defaultdict(list)
    res = sparql.query().convert()
    for bindings in res["results"]["bindings"]:
        for key, value in bindings.items():
            results[key].append(value["value"])

    results_df = pd.DataFrame.from_dict(results)

    return results_df


def download_chebi(out_folder: str) -> str:
    """
    Downloads the structures.csv.gz file from
    ChEBI FTP and saves to the destination folder

    Args:
     out_folder: path to the output folder
    """
    chebi_structures_url = get_config("reactant", "chebi_structure_file")
    response = requests_retry_session().get(chebi_structures_url)
    response.raise_for_status()
    chebi_structure_file = os.path.join(
        out_folder, os.path.basename(chebi_structures_url)
    )

    with open(chebi_structure_file, "wb") as fh:
        fh.write(response.content)

    return chebi_structure_file


def parse_ligand(cif_path: str, ligand_type: str) -> Component:
    """
    Parse cif file of ligand and returns Component object

    Args:
        cif_path: path to ligand cif file
        ligand_type: type of ligand (CCD, PRD, CLC)
    """

    if ligand_type == "CCD":
        component = ccd_reader.read_pdb_cif_file(cif_path).component
    elif ligand_type == "CLC":
        component = clc_reader.read_clc_cif_file(cif_path).component
    elif ligand_type == "PRD":
        component = prd_reader.read_pdb_cif_file(cif_path).component
    else:
        raise ValueError("ligand_type should be either CCD, PRD or CLC")

    return component
