import os
import subprocess

from pymol import cmd

from .mutation import Mutation
from .residue import one_to_three


def compute_affinity(
    molecule_name: str, antibody_chain: str, antigen_chains: list[str]
) -> float:
    molecule_file = f"{molecule_name}_tmp.pdb"
    cmd.save(molecule_file, molecule_name)

    command = [
        "prodigy",
        molecule_file,
        "--selection",
        antibody_chain,
        ",".join(antigen_chains),
        "--quiet",
    ]
    print(f"Running command {' '.join(command)}")

    res = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    if res.stderr:
        print(f"Error: {res.stderr}")
        exit(-1)

    if res.stdout is None or len(res.stdout.split()) != 2:
        print(f"Error: could not parse Prodigy output. Got: {res.stdout}")
        exit(-1)

    os.remove(molecule_file)

    return round(float(res.stdout.split()[1]), 2)


def compute_ddg(
    molecule_name: str, mutation: Mutation, chain: str, partner_chains: list[str]
) -> float:
    original_affinity = compute_affinity(molecule_name, chain, partner_chains)
    print(f"Original affinity: {original_affinity}")

    cmd.wizard("mutagenesis")
    cmd.do("refresh_wizard")
    cmd.get_wizard().do_select(f"chain {chain} and resi {mutation.start_residue.id}")
    cmd.get_wizard().set_mode(one_to_three(mutation.target_resn))
    cmd.frame(1)
    cmd.get_wizard().apply()
    cmd.set_wizard()

    mutated_file_path = os.path.join(
        "pdbs",
        f"{molecule_name}_{chain}_{mutation.start_residue.name}{mutation.start_residue.id}{mutation.target_resn}.pdb",
    )
    cmd.save(mutated_file_path, molecule_name)
    cmd.delete("all")

    mutated_affinity = compute_affinity(mutated_file_path, chain, partner_chains)
    print(f"Mutated affinity: {mutated_affinity}")

    return mutated_affinity - original_affinity
