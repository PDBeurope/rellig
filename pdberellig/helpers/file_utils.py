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


import gzip
import os
import re
from datetime import datetime

RE_EC_NO = re.compile(r"EC=[\d\.-]+|$")  # /$ returns '' if match cannot be found
RE_CHEBI_ID = re.compile(r"CHEBI\:\d+")
RE_RHEA_ID = re.compile(r"RHEA\:\d+|$")
RE_DATE_MODIFIED = re.compile(r"^> <Last Modified>")
RE_ENTRY_END = re.compile(r"^\$\$\$\$")
RE_CHEBI_HOLDER = re.compile(r"^> <ChEBI ID>")


def get_chebi_iterator(chebi_file: str):
    """Returns an iterator over contents
    of the ChEBI_complete.sdf.gz file

    Args:
        chebi_file (str): Path to the chebi file

    Returns:

    """

    chebi_id = None
    chebi_data = []
    last_modified = None

    with gzip.open(chebi_file, "rt") as fp:
        for line in fp:
            chebi_data.append(line)
            if RE_CHEBI_HOLDER.match(line):
                line = next(fp)
                chebi_data.append(line)
                chebi_id = line.strip()

            if RE_DATE_MODIFIED.match(line):
                line = next(fp)
                chebi_data.append(line)
                last_modified = datetime.strptime(line.strip(), "%d %b %Y")

            if RE_ENTRY_END.match(line):
                chebi_file = "".join(chebi_data)
                chebi_data = []
                yield (chebi_id, last_modified, chebi_file)


def get_cif_path(base_dir: str, ccd: str):
    """Get the path to the cif file for the given ccd.

    Args:
        base_dir (str): Path to the base directory of the ccd files."""

    return os.path.join(base_dir, ccd[0], ccd, f"{ccd}.cif")
