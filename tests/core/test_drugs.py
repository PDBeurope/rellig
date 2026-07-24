from types import SimpleNamespace

import pandas as pd
import pytest

from pdberellig.core import drugs as drugs_module
from pdberellig.core.drugs import Drugs


class FakeCifBlock:
    """Stands in for gemmi's cif block used by Component.ccd_cif_block."""

    def __init__(self, categories, rows=None):
        self._categories = categories
        self._rows = rows or []

    def get_mmcif_category_names(self):
        return self._categories

    def find(self, category, items):
        return self._rows


def make_ligand(ligand_id, cif_block):
    return SimpleNamespace(id=ligand_id, ccd_cif_block=cif_block)


@pytest.fixture(autouse=True)
def identity_as_string(monkeypatch):
    """Real gemmi cif values are CIF-quoted strings; rows in these tests are
    already plain strings, so as_string just needs to pass them through."""
    monkeypatch.setattr(drugs_module.cif, "as_string", lambda value: value)


def make_drugs(out_dir, logger, ligand_cif="ignored.cif", ligand_type="CCD"):
    return Drugs(ligand_cif, ligand_type, str(out_dir), logger)


class TestGetDrugbankTargets:
    def test_returns_empty_dataframe_when_category_absent(self, tmp_path, logger):
        drugs = make_drugs(tmp_path, logger)
        component = make_ligand("LIG", FakeCifBlock(categories=[]))

        result = drugs.get_drugbank_targets(component)

        assert isinstance(result, pd.DataFrame)
        assert result.empty

    def test_returns_targets_when_category_present(self, tmp_path, logger):
        rows = [
            {
                "_pdbe_chem_comp_drugbank_targets.name": "Target A",
                "_pdbe_chem_comp_drugbank_targets.organism": "Human",
                "_pdbe_chem_comp_drugbank_targets.uniprot_id": "P12345",
                "_pdbe_chem_comp_drugbank_targets.pharmacologically_active": "yes",
            }
        ]
        cif_block = FakeCifBlock(
            categories=["_pdbe_chem_comp_drugbank_targets."], rows=rows
        )
        drugs = make_drugs(tmp_path, logger)
        component = make_ligand("LIG", cif_block)

        result = drugs.get_drugbank_targets(component)

        assert list(result["uniprot_id"]) == ["P12345"]
        assert list(result["pharmacologically_active"]) == ["yes"]


class TestProcessEntry:
    def test_no_drugbank_targets_writes_nothing(self, tmp_path, logger, monkeypatch):
        component = make_ligand("LIG", FakeCifBlock(categories=[]))
        monkeypatch.setattr(drugs_module, "parse_ligand", lambda *a, **k: component)

        drugs = make_drugs(tmp_path, logger)
        drugs.process_entry()

        assert list(tmp_path.iterdir()) == []
        logger.info.assert_called_once()

    def test_no_pharmacologically_active_targets_writes_nothing(
        self, tmp_path, logger, monkeypatch
    ):
        rows = [
            {
                "_pdbe_chem_comp_drugbank_targets.name": "Target A",
                "_pdbe_chem_comp_drugbank_targets.organism": "Human",
                "_pdbe_chem_comp_drugbank_targets.uniprot_id": "P12345",
                "_pdbe_chem_comp_drugbank_targets.pharmacologically_active": "no",
            }
        ]
        cif_block = FakeCifBlock(
            categories=["_pdbe_chem_comp_drugbank_targets."], rows=rows
        )
        component = make_ligand("LIG", cif_block)
        monkeypatch.setattr(drugs_module, "parse_ligand", lambda *a, **k: component)

        drugs = make_drugs(tmp_path, logger)
        drugs.process_entry()

        assert list(tmp_path.iterdir()) == []

    def test_no_interacting_chains_writes_nothing(self, tmp_path, logger, monkeypatch):
        rows = [
            {
                "_pdbe_chem_comp_drugbank_targets.name": "Target A",
                "_pdbe_chem_comp_drugbank_targets.organism": "Human",
                "_pdbe_chem_comp_drugbank_targets.uniprot_id": "P12345",
                "_pdbe_chem_comp_drugbank_targets.pharmacologically_active": "yes",
            }
        ]
        cif_block = FakeCifBlock(
            categories=["_pdbe_chem_comp_drugbank_targets."], rows=rows
        )
        component = make_ligand("LIG", cif_block)
        monkeypatch.setattr(drugs_module, "parse_ligand", lambda *a, **k: component)
        monkeypatch.setattr(
            drugs_module, "get_ligand_intx_chains", lambda ligand_id: pd.DataFrame()
        )

        drugs = make_drugs(tmp_path, logger)
        drugs.process_entry()

        assert list(tmp_path.iterdir()) == []
        logger.warn.assert_called_once()

    def test_full_pipeline_writes_matched_targets(self, tmp_path, logger, monkeypatch):
        rows = [
            {
                "_pdbe_chem_comp_drugbank_targets.name": "Target A",
                "_pdbe_chem_comp_drugbank_targets.organism": "Human",
                "_pdbe_chem_comp_drugbank_targets.uniprot_id": "P12345",
                "_pdbe_chem_comp_drugbank_targets.pharmacologically_active": "yes",
            },
            {
                "_pdbe_chem_comp_drugbank_targets.name": "Target B",
                "_pdbe_chem_comp_drugbank_targets.organism": "Mouse",
                "_pdbe_chem_comp_drugbank_targets.uniprot_id": "Q99999",
                "_pdbe_chem_comp_drugbank_targets.pharmacologically_active": "no",
            },
        ]
        cif_block = FakeCifBlock(
            categories=["_pdbe_chem_comp_drugbank_targets."], rows=rows
        )
        component = make_ligand("LIG", cif_block)
        monkeypatch.setattr(drugs_module, "parse_ligand", lambda *a, **k: component)

        intx_chains = pd.DataFrame(
            {
                "pdb_id": ["1abc"],
                "auth_asym_id": ["A"],
                "struct_asym_id": ["A"],
                "uniprot_id": ["P12345"],
            }
        )
        monkeypatch.setattr(
            drugs_module, "get_ligand_intx_chains", lambda ligand_id: intx_chains
        )

        drugs = make_drugs(tmp_path, logger)
        drugs.process_entry()

        out_file = tmp_path / "LIG_drug_annotation.tsv"
        assert out_file.exists()
        result = pd.read_csv(out_file, sep="\t")
        assert list(result["uniprot_id"]) == ["P12345"]
        assert list(result["name"]) == ["Target A"]
        # the non-pharmacologically-active target must not leak into the output
        assert "Q99999" not in result["uniprot_id"].to_list()
