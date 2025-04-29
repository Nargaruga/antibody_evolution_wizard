import os
import subprocess
import tempfile

from pymol import CmdException

from antibody_evolution.mutation import Mutation
from antibody_evolution.residue import one_to_three


class AffinityComputationError(Exception):
    """Custom exception for affinity computation errors."""

    pass


class DDGComputationError(Exception):
    """Custom exception for DDG computation errors."""

    pass


def compute_affinity(
    molecule_file: str, antibody_chains: list[str], antigen_chains: list[str]
) -> float:
    """Compute the binding affinity of the given antibody-antigen complex.

    :param molecule_file: path to the PDB file containing the antibody-antigen complex
    :param antibody_chains: list of antibody chain identifiers (e.g., ["H", "L"])
    :param antigen_chains: list of antigen chain identifiers (e.g., ["A", "B"])

    :return the binding affinity, rounded to 2 decimal places
    :raises AffinityComputationError: if the Prodigy call fails or output is invalid
    """
    command = [
        "prodigy",
        molecule_file,
        "--selection",
        ",".join(antibody_chains),
        ",".join(antigen_chains),
        "--quiet",
    ]
    try:
        res = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
        )

        if not res.stdout or len(res.stdout.split()) != 2:
            raise AffinityComputationError(
                f"command {' '.join(command)} gave unexpected output: {res.stdout.strip()}"
            )

        return round(float(res.stdout.split()[1]), 2)

    except subprocess.CalledProcessError as e:
        raise AffinityComputationError(
            f"command {' '.join(command)} failed: {e.stderr}"
        )


def compute_ddg(
    molecule_file_path,
    chain,
    original_affinity,
    mutation: Mutation,
    antibody_chains,
    partner_chains,
    pymol_instance,
) -> float:
    """Compute the mutation's DDG and return it."""
    molecule="prova"
    pymol_instance.cmd.load(molecule_file_path, molecule)

    selection_string = (
        f"/{molecule}//{chain}/{mutation.start_residue.id}"
    )
    pymol_instance.cmd.select("tmp", selection_string)

    pymol_instance.cmd.wizard("mutagenesis")
    pymol_instance.cmd.do("refresh_wizard")
    try:
        pymol_instance.cmd.get_wizard().do_select("tmp")
    except CmdException as e:
        # Not sure why this happens. A molecule that causes this is 3L5X
        msg = f"Skipping mutation {mutation} due to error selecting residue: {e}"
        raise DDGComputationError(msg)

    pymol_instance.cmd.get_wizard().set_mode(one_to_three(mutation.target_resn))
    pymol_instance.cmd.frame(1)
    pymol_instance.cmd.get_wizard().apply()
    pymol_instance.cmd.set_wizard()

    mutated_molecule_file_handle, mutated_molecule_file_path = tempfile.mkstemp(
        suffix=".pdb"
    )
    pymol_instance.cmd.save(mutated_molecule_file_path, molecule)
    try:
        mutated_affinity = compute_affinity(
            mutated_molecule_file_path,
            antibody_chains,
            partner_chains,
        )
    except AffinityComputationError as e:
        msg = f"Skipping mutation {mutation} due to error computing new affinity: {e}"
        raise DDGComputationError(msg)
    finally:
        os.close(mutated_molecule_file_handle)
        os.remove(mutated_molecule_file_path)

    pymol_instance.cmd.delete("tmp")
    pymol_instance.cmd.delete("tmp_molecule")

    return round(mutated_affinity - original_affinity, 2)
