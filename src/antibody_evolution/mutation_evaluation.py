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

    print(f"Computing affinity betwen chains {antibody_chains} and {antigen_chains}")

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
            f"command {' '.join(command)} failed: {e}"
        ) from e


def apply_mutation(
    obj: str,
    chain: str,
    mutation: Mutation,
    pymol_cmd,
):
    """Apply the mutation to the given object in PyMOL."""

    print(f"Applying mutation {mutation} to {obj}")
    selection_string = f"/{obj}//{chain}/{mutation.start_residue.id}"
    pymol_cmd.select("tmp", selection_string)

    pymol_cmd.wizard("mutagenesis")
    pymol_cmd.do("refresh_wizard")
    pymol_cmd.get_wizard().do_select("tmp")

    pymol_cmd.get_wizard().set_mode(one_to_three(mutation.target_resn))
    pymol_cmd.frame(1)
    pymol_cmd.get_wizard().apply()
    pymol_cmd.set_wizard()

    pymol_cmd.delete("tmp")


def compute_ddg(
    molecule_file_path,
    chain,
    original_affinity,
    mutation: Mutation,
    antibody_chains,
    partner_chains,
    pymol_cmd,
) -> float:
    """Compute the mutation's DDG and return it."""

    molecule = "to_mutate"
    pymol_cmd.load(molecule_file_path, molecule)
    try:
        apply_mutation(
            molecule,
            chain,
            mutation,
            pymol_cmd,
        )
    except CmdException as e:
        # Not sure why this happens. A molecule that causes this is 3L5X
        raise DDGComputationError(
            f"Skipping mutation {mutation} due to error applying mutation."
        ) from e

    mutated_molecule_file_handle, mutated_molecule_file_path = tempfile.mkstemp(
        suffix=".pdb"
    )
    pymol_cmd.save(mutated_molecule_file_path, molecule)
    pymol_cmd.delete(molecule)

    try:
        mutated_affinity = compute_affinity(
            mutated_molecule_file_path,
            antibody_chains,
            partner_chains,
        )
        print(f"New affinity for {mutation} is {mutated_affinity}")
    except AffinityComputationError as e:
        raise DDGComputationError(
            f"Skipping mutation {mutation} due to error computing new affinity: {e}"
        ) from e
    finally:
        os.close(mutated_molecule_file_handle)
        os.remove(mutated_molecule_file_path)

    return round(mutated_affinity - original_affinity, 2)
