import argparse
import os

import click

from pdberellig.conf import get_config
from pdberellig.core.cofactors import Cofactors
from pdberellig.core.drugs import Drugs
from pdberellig.core.reactants import Reactants
from pdberellig.helpers.utils import setup_log


@click.group("CLI", help="CLI for pdberellig pipelines.")
def main():
    """Application entry point."""
    return 0


@main.command(
    "cofactors",
    help=(
        "PDBe RelLig pipeline - cofactors mode\n"
        "=================\n\n"
        "This pipeline processes ligand CIF files from the PDB and "
        "annotates them as cofactor-like molecules"
    ),
)
@click.option("--cif", type=str, required=True, help="path to input cif file")
@click.option(
    "--ligand-type",
    type=click.Choice(["CCD", "PRD", "CLC"], case_sensitive=False),
    required=True,
    help="type of ligand in the PDB",
)
@click.option("--out-dir", type=str, required=True, help="path to output directory")
def cofactors(cif: str, ligand_type: str, out_dir: str):
    """Cofactors entry point."""

    log = setup_log("functional annotation pipeline", "cofactors")

    args = argparse.Namespace(cif=cif, ligand_type=ligand_type, out_dir=out_dir)
    cofactors = Cofactors(log, args)
    cofactors.process_entry()


@main.command(
    "reactants",
    help=(
        "PDBe RelLig pipeline - reactants mode\n"
        "=================\n\n"
        "This pipeline processes ligand CIF files from the PDB and "
        "annotates them as reactant-like molecule"
    ),
)
@click.option("--cif", type=str, required=True, help="path to input cif file")
@click.option(
    "--ligand-type",
    type=click.Choice(["CCD", "PRD", "CLC"], case_sensitive=False),
    required=True,
    help="type of ligand in the PDB",
)
@click.option(
    "--chebi-structure-file",
    type=str,
    required=True,
    help="Path to the ChEBI SDF file",
)
@click.option("--out-dir", type=str, required=True, help="path to output directory")
@click.option(
    "--update-chebi",
    is_flag=True,
    help="Flag to update ChEBI structures file",
)
@click.option(
    "--minimal-ligand-size",
    type=int,
    required=False,
    default=get_config("main", "minimal_ligand_size"),
    help="Minimum ligand size.",
    show_default=True,
)
def reactants(
    cif: str,
    ligand_type: str,
    chebi_structure_file: str,
    out_dir: str,
    update_chebi: str,
    minimal_ligand_size: int,
):
    """Reactants entry point."""

    log = setup_log("functional annotation pipeline", "reactants")

    if not os.path.isfile(chebi_structure_file):
        raise FileNotFoundError(
            f"Path to the ChEBI stuture file ({chebi_structure_file}) does not exist."
        )

    if not os.path.isdir(out_dir):
        raise FileNotFoundError(
            f"Path to the output directory ({out_dir}) does not exist."
        )

    args = argparse.Namespace(
        cif=cif,
        ligand_type=ligand_type,
        chebi_structure_file=chebi_structure_file,
        update_chebi=update_chebi,
        out_dir=out_dir,
        minimal_ligand_size=minimal_ligand_size,
    )
    reactants = Reactants(log, args)
    reactants.process_entry()


@main.command(
    "drugs",
    help="PDBe RelLig pipeline - drugs mode\n"
    "=================\n\n"
    "This pipeline parses ligand CIF files from the PDB and annotates "
    "them as drug-like molecule",
)
@click.option("--cif", type=str, required=True, help="path to input cif file")
@click.option("--out-dir", type=str, required=True, help="path to output directory")
@click.option(
    "--ligand-type",
    type=click.Choice(["CCD", "PRD", "CLC"], case_sensitive=False),
    required=True,
    help="type of ligand in the PDB",
)
def drugs(cif: str, ligand_type: str, out_dir: str):
    log = setup_log("functional annotation pipeline", "drugs")

    args = argparse.Namespace(
        cif=cif,
        ligand_type=ligand_type,
        out_dir=out_dir,
    )

    drugs = Drugs(log, args)
    drugs.process_entry()
