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
Similarities pipeline data model
"""

from pdbeccdutils.core import ccd_reader

from pdbe_relic.conf import get_config
# from pdbe_relic.core.compare import compare
from pdbe_relic.helpers.file_utils import get_cif_path


class Similarities:
    """
    Similarities pipeline data model
    """

    def __init__(self, log, args):
        self.log = log
        self.args = args

    def run(self, source_ccd, target_ccd):
        self.log.info(f"Calculating similarities between {source_ccd} and {target_ccd}")
        src_ccd_path = get_cif_path(self.args.chem_comp_base_dir, source_ccd)
        target_ccd_path = get_cif_path(self.args.chem_comp_base_dir, target_ccd)
        src = ccd_reader.read_pdb_cif_file(src_ccd_path).component.mol_no_h
        target = ccd_reader.read_pdb_cif_file(target_ccd_path).component.mol_no_h
        threshold = float(get_config("main", "similarities_threshold"))

        result = compare((source_ccd, src), (target_ccd, target), threshold)

        if not result.result:
            return None

        if result.result.similarity_score >= threshold:
            score = round(result.result.similarity_score, 3)
            data_to_write = (
                source_ccd,
                target_ccd,
                ";".join(
                    [
                        f'{src.GetAtomWithIdx(k).GetProp("name")}:'
                        f'{target.GetAtomWithIdx(v).GetProp("name")}'
                        for k, v in result.result.mapping.items()
                    ]
                ),
                str(score),
            )
            self.log.info(
                f"{source_ccd} is similar to {target_ccd} with "
                f"similarity score {score}."
            )
            return ",".join(data_to_write)

        return None
