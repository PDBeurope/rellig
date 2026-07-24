from types import SimpleNamespace

import pandas as pd
import pytest
from rdkit import Chem

from pdberellig.core import reactants as reactants_module
from pdberellig.core.models import CompareObj, Similarity
from pdberellig.core.reactants import Reactants


class FakeMol:
    """Stands in for an RDKit Mol; only GetAtoms() is exercised by the pipeline."""

    def __init__(self, num_atoms):
        self._num_atoms = num_atoms

    def GetAtoms(self):
        return list(range(self._num_atoms))


def make_reactants(out_dir, logger, **kwargs):
    defaults = dict(
        ligand_cif="ignored.cif",
        ligand_type="CCD",
        chebi_structure_file=None,
        out_dir=str(out_dir),
        logger=logger,
        minimal_ligand_size=5,
        update_chebi=False,
    )
    defaults.update(kwargs)
    return Reactants(**defaults)


def make_template(template_id, similarity_score):
    template = CompareObj(template_id, None)
    template.similarity_to = lambda other, threshold=0.01: Similarity(
        template_id, other.id, SimpleNamespace(similarity_score=similarity_score)
    )
    return template


class TestProcessEntry:
    def test_skips_ligand_smaller_than_minimal_size(
        self, tmp_path, logger, monkeypatch
    ):
        component = SimpleNamespace(id="LIG", mol_no_h=FakeMol(3))
        monkeypatch.setattr(reactants_module, "parse_ligand", lambda *a, **k: component)

        reactants = make_reactants(tmp_path, logger, minimal_ligand_size=5)
        reactants.process_entry()

        assert list(tmp_path.iterdir()) == []
        logger.debug.assert_called_once()

    def test_writes_nothing_when_no_interacting_chains(
        self, tmp_path, logger, monkeypatch
    ):
        component = SimpleNamespace(id="LIG", mol_no_h=FakeMol(10))
        monkeypatch.setattr(reactants_module, "parse_ligand", lambda *a, **k: component)
        monkeypatch.setattr(
            reactants_module,
            "get_ligand_intx_chains",
            lambda ligand_id: pd.DataFrame(),
        )

        reactants = make_reactants(tmp_path, logger)
        reactants.process_entry()

        assert list(tmp_path.iterdir()) == []
        logger.warn.assert_called_once()

    def test_writes_nothing_when_reactant_annotation_is_empty(
        self, tmp_path, logger, monkeypatch
    ):
        component = SimpleNamespace(id="LIG", mol_no_h=FakeMol(10))
        monkeypatch.setattr(reactants_module, "parse_ligand", lambda *a, **k: component)
        intx_chains = pd.DataFrame(
            {
                "pdb_id": ["1abc"],
                "auth_asym_id": ["A"],
                "struct_asym_id": ["A"],
                "uniprot_id": ["P12345"],
            }
        )
        monkeypatch.setattr(
            reactants_module,
            "get_ligand_intx_chains",
            lambda ligand_id: intx_chains,
        )

        reactants = make_reactants(tmp_path, logger)
        monkeypatch.setattr(
            reactants,
            "get_reactant_annotation",
            lambda ligand, uniprot_ids: pd.DataFrame(),
        )
        reactants.process_entry()

        assert list(tmp_path.iterdir()) == []

    def test_writes_merged_annotation_on_success(self, tmp_path, logger, monkeypatch):
        component = SimpleNamespace(id="LIG", mol_no_h=FakeMol(10))
        monkeypatch.setattr(reactants_module, "parse_ligand", lambda *a, **k: component)
        intx_chains = pd.DataFrame(
            {
                "pdb_id": ["1abc"],
                "auth_asym_id": ["A"],
                "struct_asym_id": ["A"],
                "uniprot_id": ["P12345"],
            }
        )
        monkeypatch.setattr(
            reactants_module,
            "get_ligand_intx_chains",
            lambda ligand_id: intx_chains,
        )

        reactants_sim = pd.DataFrame(
            {
                "uniprot_id": ["P12345"],
                "rhea_id": ["12345"],
                "chebi_id": ["CHEBI:1"],
                "similarity": [0.8],
            }
        )
        reactants = make_reactants(tmp_path, logger)
        monkeypatch.setattr(
            reactants,
            "get_reactant_annotation",
            lambda ligand, uniprot_ids: reactants_sim,
        )
        reactants.process_entry()

        out_file = tmp_path / "LIG_reactant_annotation.tsv"
        assert out_file.exists()
        result = pd.read_csv(out_file, sep="\t")
        assert list(result["pdb_id"]) == ["1abc"]
        assert list(result["chebi_id"]) == ["CHEBI:1"]


class TestGetReactantAnnotation:
    def test_returns_empty_when_no_reactions_found(self, tmp_path, logger, monkeypatch):
        reactants = make_reactants(tmp_path, logger)
        monkeypatch.setattr(
            reactants, "get_reactions", lambda uniprot_ids: pd.DataFrame()
        )

        result = reactants.get_reactant_annotation(CompareObj("LIG", None), ["P12345"])

        assert result.empty
        logger.warn.assert_called_once()

    def test_returns_empty_when_no_reaction_participants_found(
        self, tmp_path, logger, monkeypatch
    ):
        reactants = make_reactants(tmp_path, logger)
        reactions = pd.DataFrame({"uniprot_id": ["P12345"], "rhea_id": ["12345"]})
        monkeypatch.setattr(reactants, "get_reactions", lambda uniprot_ids: reactions)
        monkeypatch.setattr(
            reactants, "get_reaction_participants", lambda rhea_ids: pd.DataFrame()
        )

        result = reactants.get_reactant_annotation(CompareObj("LIG", None), ["P12345"])

        assert result.empty
        logger.warn.assert_called_once()

    def test_returns_empty_when_no_similar_participant_found(
        self, tmp_path, logger, monkeypatch
    ):
        reactants = make_reactants(tmp_path, logger)
        reactions = pd.DataFrame({"uniprot_id": ["P12345"], "rhea_id": ["12345"]})
        participants = pd.DataFrame({"rhea_id": ["12345"], "chebi_id": ["CHEBI:1"]})
        monkeypatch.setattr(reactants, "get_reactions", lambda uniprot_ids: reactions)
        monkeypatch.setattr(
            reactants, "get_reaction_participants", lambda rhea_ids: participants
        )
        monkeypatch.setattr(reactants, "parse_chebi", lambda chebi_ids: [])
        monkeypatch.setattr(reactants, "get_similarities", lambda ligand, templates: {})

        result = reactants.get_reactant_annotation(CompareObj("LIG", None), ["P12345"])

        assert result.empty
        logger.info.assert_called_once()

    def test_merges_reactions_participants_and_similarities(
        self, tmp_path, logger, monkeypatch
    ):
        reactants = make_reactants(tmp_path, logger)
        reactions = pd.DataFrame({"uniprot_id": ["P12345"], "rhea_id": ["12345"]})
        participants = pd.DataFrame({"rhea_id": ["12345"], "chebi_id": ["CHEBI:1"]})
        monkeypatch.setattr(reactants, "get_reactions", lambda uniprot_ids: reactions)
        monkeypatch.setattr(
            reactants, "get_reaction_participants", lambda rhea_ids: participants
        )
        monkeypatch.setattr(reactants, "parse_chebi", lambda chebi_ids: [])
        monkeypatch.setattr(
            reactants,
            "get_similarities",
            lambda ligand, templates: {"chebi_id": ["CHEBI:1"], "similarity": [0.9]},
        )

        result = reactants.get_reactant_annotation(CompareObj("LIG", None), ["P12345"])

        assert list(result["uniprot_id"]) == ["P12345"]
        assert list(result["chebi_id"]) == ["CHEBI:1"]
        assert list(result["similarity"]) == [0.9]


class TestParseChebi:
    def _write_chebi_tsv(self, path, rows):
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    def test_filters_by_status_default_structure_and_size(self, tmp_path, logger):
        big_molfile = Chem.MolToMolBlock(
            Chem.MolFromSmiles("c1ccccc1")
        )  # benzene, 6 heavy atoms
        small_molfile = Chem.MolToMolBlock(
            Chem.MolFromSmiles("O")
        )  # water, 1 heavy atom

        chebi_file = tmp_path / "structures.tsv"
        self._write_chebi_tsv(
            chebi_file,
            {
                "compound_id": ["CHEBI:1", "CHEBI:2", "CHEBI:3", "CHEBI:4"],
                "status_id": ["1", "0", "1", "1"],  # CHEBI:2 excluded: not curated
                "default_structure": [
                    "true",
                    "true",
                    "false",
                    "true",
                ],  # CHEBI:3 excluded
                "molfile": [big_molfile, big_molfile, big_molfile, small_molfile],
            },
        )

        reactants = make_reactants(
            tmp_path,
            logger,
            chebi_structure_file=str(chebi_file),
            minimal_ligand_size=5,
        )

        templates = reactants.parse_chebi(["CHEBI:1", "CHEBI:2", "CHEBI:3", "CHEBI:4"])

        # only CHEBI:1 passes all filters: curated, default structure, big enough
        assert [t.id for t in templates] == ["CHEBI:1"]
        # CHEBI:4 was filtered out for being too small, and should be logged
        logger.debug.assert_called_once()

    def test_skips_unparseable_molfile(self, tmp_path, logger):
        chebi_file = tmp_path / "structures.tsv"
        self._write_chebi_tsv(
            chebi_file,
            {
                "compound_id": ["CHEBI:1"],
                "status_id": ["1"],
                "default_structure": ["true"],
                "molfile": ["not a valid molfile"],
            },
        )

        reactants = make_reactants(
            tmp_path, logger, chebi_structure_file=str(chebi_file)
        )

        templates = reactants.parse_chebi(["CHEBI:1"])

        assert templates == []
        logger.warn.assert_called_once()


class TestGetSimilarities:
    def test_keeps_only_templates_above_threshold(self, tmp_path, logger):
        reactants = make_reactants(tmp_path, logger)
        # main.reactants_threshold in conf.ini is 0.7
        templates = [make_template("CHEBI:1", 0.9), make_template("CHEBI:2", 0.3)]

        result = reactants.get_similarities(CompareObj("LIG", None), templates)

        assert result["chebi_id"] == ["CHEBI:1"]
        assert result["similarity"] == [0.9]
