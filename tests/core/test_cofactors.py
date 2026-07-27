from types import SimpleNamespace

import pandas as pd
import pytest

from pdberellig.core import cofactors as cofactors_module
from pdberellig.core.cofactors import Cofactors
from pdberellig.core.models import CofactorSim, CompareObj, Similarity


def make_cofactors(out_dir, logger, ligand_cif="ignored.cif", ligand_type="CCD"):
    return Cofactors(ligand_cif, ligand_type, str(out_dir), logger)


def make_template(template_id, similarity_score):
    """A CompareObj whose similarity_to is stubbed instead of running RDKit/PARITY."""
    template = CompareObj(template_id, None)
    template.similarity_to = lambda other, threshold=0.01: Similarity(
        template_id, other.id, SimpleNamespace(similarity_score=similarity_score)
    )
    return template


class TestGetSimilarity:
    def test_returns_none_when_no_template_meets_its_threshold(self, tmp_path, logger):
        cofactors = make_cofactors(tmp_path, logger)
        ligand = CompareObj("LIG", None)
        templates = [make_template("FAD", 0.5)]
        cofactor_details = {
            "FAD": {"threshold": 0.87, "representative": "FAD", "id": 2}
        }

        result = cofactors.get_similarity(
            ligand, templates=templates, cofactor_details=cofactor_details
        )

        assert result is None

    def test_returns_match_when_template_and_representative_clear_threshold(
        self, tmp_path, logger, monkeypatch
    ):
        cofactors = make_cofactors(tmp_path, logger)
        ligand = CompareObj("LIG", None)
        templates = [make_template("FAD", 0.9)]
        cofactor_details = {
            "FAD": {"threshold": 0.87, "representative": "FAD", "id": 2}
        }

        representative = make_template("FAD", 0.9)
        monkeypatch.setattr(
            cofactors, "get_representative", lambda details: representative
        )

        result = cofactors.get_similarity(
            ligand, templates=templates, cofactor_details=cofactor_details
        )

        assert isinstance(result, CofactorSim)
        assert result.query_id == "LIG"
        assert result.template_sim.target_id == "FAD"
        assert result.representative_sim.result.similarity_score == 0.9

    def test_skips_when_representative_similarity_below_threshold(
        self, tmp_path, logger, monkeypatch
    ):
        cofactors = make_cofactors(tmp_path, logger)
        ligand = CompareObj("LIG", None)
        templates = [make_template("FAD", 0.9)]
        cofactor_details = {
            "FAD": {"threshold": 0.87, "representative": "FAD", "id": 2}
        }

        representative = make_template("FAD", 0.2)
        monkeypatch.setattr(
            cofactors, "get_representative", lambda details: representative
        )

        result = cofactors.get_similarity(
            ligand, templates=templates, cofactor_details=cofactor_details
        )

        assert result is None

    def test_picks_the_higher_combined_score_across_multiple_templates(
        self, tmp_path, logger, monkeypatch
    ):
        cofactors = make_cofactors(tmp_path, logger)
        ligand = CompareObj("LIG", None)
        templates = [make_template("FAD", 0.9), make_template("NAD2", 0.95)]
        # representative field is set equal to the template id to keep the
        # lookup in this test trivial; in production data they can differ.
        cofactor_details = {
            "FAD": {"threshold": 0.87, "representative": "FAD", "id": 2},
            "NAD2": {"threshold": 0.68, "representative": "NAD2", "id": 4},
        }
        representative_scores = {"FAD": 0.9, "NAD2": 0.95}
        monkeypatch.setattr(
            cofactors,
            "get_representative",
            lambda details: make_template(
                details["representative"],
                representative_scores[details["representative"]],
            ),
        )

        result = cofactors.get_similarity(
            ligand, templates=templates, cofactor_details=cofactor_details
        )

        assert result.template_sim.target_id == "NAD2"

    def test_continues_past_template_that_raises(self, tmp_path, logger):
        cofactors = make_cofactors(tmp_path, logger)
        ligand = CompareObj("LIG", None)

        broken = CompareObj("BROKEN", None)

        def raise_error(other, threshold):
            raise RuntimeError("boom")

        broken.similarity_to = raise_error
        # "good" scores below its own template threshold, so it is skipped via
        # the normal `continue` path without ever reaching get_representative -
        # keeping this test deterministic and free of real file I/O.
        good = make_template("FAD", 0.5)
        cofactor_details = {
            "BROKEN": {"threshold": 0.5, "representative": "FAD", "id": 2},
            "FAD": {"threshold": 0.87, "representative": "FAD", "id": 2},
        }

        result = cofactors.get_similarity(
            ligand, templates=[broken, good], cofactor_details=cofactor_details
        )

        # the exception on "broken" is swallowed and logged; get_similarity
        # itself must not raise.
        assert result is None
        assert logger.warn.called


class TestWriteCofactorResults:
    def test_writes_json_file(self, tmp_path, logger):
        cofactors = make_cofactors(tmp_path, logger)
        results = {
            "LIG": {
                "template": {"id": "FAD", "similarity": 0.9},
                "representative": {"id": "FAD", "similarity": 0.9},
                "pdb_chains": [],
            }
        }

        cofactors._write_cofactor_results(results)

        out_file = tmp_path / "LIG_cofactor_annotation.json"
        assert out_file.exists()
        import json

        with open(out_file) as fh:
            assert json.load(fh) == results


class TestProcessEntry:
    @pytest.fixture(autouse=True)
    def stub_lookup_tables(self, monkeypatch):
        cofactor_details = {
            "FAD": {"threshold": 0.87, "representative": "FAD", "id": 2}
        }
        cofactor_ec = pd.DataFrame(
            {"EC_NO": ["1.1.1.1"], "COFACTOR_ID": [2], "SOURCE": ["test"]}
        )
        monkeypatch.setattr(cofactors_module, "init_rdkit_templates", lambda: [])
        monkeypatch.setattr(
            cofactors_module, "get_cofactor_details", lambda: cofactor_details
        )
        monkeypatch.setattr(cofactors_module, "get_cofactor_ec", lambda: cofactor_ec)

    def test_writes_nothing_when_no_cofactor_similarity_found(
        self, tmp_path, logger, monkeypatch
    ):
        component = SimpleNamespace(id="LIG", mol_no_h=None)
        monkeypatch.setattr(cofactors_module, "parse_ligand", lambda *a, **k: component)

        cofactors = make_cofactors(tmp_path, logger)
        monkeypatch.setattr(cofactors, "get_similarity", lambda ligand: None)

        cofactors.process_entry()

        assert list(tmp_path.iterdir()) == []

    def test_writes_nothing_when_no_interacting_chains(
        self, tmp_path, logger, monkeypatch
    ):
        component = SimpleNamespace(id="LIG", mol_no_h=None)
        monkeypatch.setattr(cofactors_module, "parse_ligand", lambda *a, **k: component)

        sim = CofactorSim(
            "LIG",
            Similarity("FAD", "LIG", SimpleNamespace(similarity_score=0.9)),
            Similarity("FAD", "LIG", SimpleNamespace(similarity_score=0.9)),
        )
        cofactors = make_cofactors(tmp_path, logger)
        monkeypatch.setattr(cofactors, "get_similarity", lambda ligand: sim)
        monkeypatch.setattr(
            cofactors_module, "get_ligand_intx_chains", lambda ligand_id: pd.DataFrame()
        )

        cofactors.process_entry()

        assert list(tmp_path.iterdir()) == []
        logger.warn.assert_called_once()

    def test_writes_nothing_when_ec_numbers_do_not_overlap(
        self, tmp_path, logger, monkeypatch
    ):
        component = SimpleNamespace(id="LIG", mol_no_h=None)
        monkeypatch.setattr(cofactors_module, "parse_ligand", lambda *a, **k: component)

        sim = CofactorSim(
            "LIG",
            Similarity("FAD", "LIG", SimpleNamespace(similarity_score=0.9)),
            Similarity("FAD", "LIG", SimpleNamespace(similarity_score=0.9)),
        )
        cofactors = make_cofactors(tmp_path, logger)
        monkeypatch.setattr(cofactors, "get_similarity", lambda ligand: sim)

        intx_chains = pd.DataFrame(
            {
                "pdb_id": ["1abc"],
                "auth_asym_id": ["A"],
                "struct_asym_id": ["A"],
                "uniprot_id": ["P12345"],
                "ec_number": ["9.9.9.9"],  # does not overlap with cofactor_ec fixture
            }
        )
        monkeypatch.setattr(
            cofactors_module, "get_ligand_intx_chains", lambda ligand_id: intx_chains
        )

        cofactors.process_entry()

        assert list(tmp_path.iterdir()) == []

    def test_writes_results_when_ec_numbers_overlap(
        self, tmp_path, logger, monkeypatch
    ):
        component = SimpleNamespace(id="LIG", mol_no_h=None)
        monkeypatch.setattr(cofactors_module, "parse_ligand", lambda *a, **k: component)

        sim = CofactorSim(
            "LIG",
            Similarity("FAD", "LIG", SimpleNamespace(similarity_score=0.9)),
            Similarity("FAD", "LIG", SimpleNamespace(similarity_score=0.9)),
        )
        cofactors = make_cofactors(tmp_path, logger)
        monkeypatch.setattr(cofactors, "get_similarity", lambda ligand: sim)

        intx_chains = pd.DataFrame(
            {
                "pdb_id": ["1abc"],
                "auth_asym_id": ["A"],
                "struct_asym_id": ["A"],
                "uniprot_id": ["P12345"],
                "ec_number": ["1.1.1.1"],  # overlaps with cofactor_ec fixture
            }
        )
        monkeypatch.setattr(
            cofactors_module, "get_ligand_intx_chains", lambda ligand_id: intx_chains
        )

        cofactors.process_entry()

        out_file = tmp_path / "LIG_cofactor_annotation.json"
        assert out_file.exists()
