import argparse
import os

import click

from pdbe_relic.conf import get_config
from pdbe_relic.core.cofactors import Cofactors
from pdbe_relic.core.reactants import Reactants
from pdbe_relic.core.similarities import Similarities
from pdbe_relic.helpers.utils import setup_log


@click.group("CLI", help="CLI for relic pipelines.")
def main():
    """Application entry point."""
    return 0


@main.command(
    "cofactors",
    help=(
        "RELIc pipeline - cofactors mode\n"
        "=================\n\n"
        "This pipeline processes CCD CIF files and looks for a similarity "
        "to already established cofactors. Any matches are reported."
    ),
)
@click.option("--cif", type=str, required=True, help="path to input cif file")
@click.option("--ligand-type", 
              type=click.Choice(['CCD', 'PRD', 'CLC'], case_sensitive=False), 
              required=True,
              help="type of ligand in the PDB") 
@click.option("--out-dir", type=str, required=True, help="path to output directory")

def cofactors(cif: str, ligand_type: str, out_dir:str):
    """Cofactors entry point."""

    log = setup_log("functional annotation pipeline", "cofactors")

    args = argparse.Namespace(cif=cif, ligand_type=ligand_type, out_dir=out_dir)
    cofactors = Cofactors(log, args)
    cofactors.process_entry()


@main.command(
    "reactants",
    help=(
        "RELIc pipeline - reactants mode\n"
        "=================\n\n"
        "This pipeline processes Rhea database and establishes "
        "similarities to reactants defined by ChEBI and CCD CIF."
    ),
)
@click.option("--cif", 
              type=str, 
              required=True,
              help="path to input cif file")

@click.option("--ligand-type", 
              type=click.Choice(['CCD', 'PRD', 'CLC'], case_sensitive=False), 
              required=True,
              help="type of ligand in the PDB") 

@click.option(
    "--chebi-structure-file",
    type=str,
    required=True,
    help="Path to the ChEBI SDF file",
)

@click.option("--out-dir", 
              type=str, 
              required=True, 
              help="path to output directory")

@click.option(
    "--update-chebi",
    is_flag=True,
    help="Path to the ChEBI archive files",

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
    "similarities",
    help="RELIc pipeline - similarities mode\n"
    "=================\n\n"
    "This pipeline establishes atom level substructure matches "
    "between a pair of chem comps.",
)
@click.option("--source-ccd", type=str, required=True, help="Source CCD")
@click.option("--target-ccd", type=str, required=True, help="Target CCD")
@click.option(
    "-p",
    "--chem-comp-base-dir",
    type=str,
    required=True,
    help="Path to the CCD CIF base directory.",
)
def similarities(source_ccd: str, target_ccd: str, chem_comp_base_dir: str):
    log = setup_log("functional annotation pipeline", "similarities")

    if not os.path.isdir(chem_comp_base_dir):
        raise ValueError(
            f"Path to the CCD files ({chem_comp_base_dir}) does not exist."
        )

    args = argparse.Namespace(
        source_ccd=source_ccd,
        target_ccd=target_ccd,
        chem_comp_base_dir=chem_comp_base_dir,
    )

    similarities = Similarities(log, args)
    similarities.run(source_ccd, target_ccd)
