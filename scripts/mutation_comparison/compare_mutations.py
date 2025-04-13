from __future__ import annotations

import sys
import csv
import re
from dataclasses import dataclass
import tempfile


from pymol import cmd

from antibody_evolution.mutation_evaluation import compute_affinity
from antibody_evolution.mutation_suggestions import (
    is_residue_valid,
    get_mutation_suggestions,
)
from antibody_evolution.mutation import Mutation
from antibody_evolution.residue import Residue, one_to_three


@dataclass
class Experiment:
    """An experiment with a mutation."""

    mutation: Mutation
    partner_chains: list[str]
    ddg: float

    def serialize(self):
        """Serialize the experiment to a CSV file."""

        with open("experiments.csv", "a") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    self.mutation.start_residue.molecule,
                    self.mutation.start_residue.chain,
                    self.mutation.start_residue.id,
                    self.mutation.target_resn,
                    ",".join(self.partner_chains),
                    self.ddg,
                ]
            )


def parse_experiments(file_path):
    regex = re.compile(r"([A-Z])([0-9]+)([A-Z])")

    experiments: dict[str, dict[str, list[Experiment]]] = {}
    with open(file_path, "r") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            molecule = row[1]

            if molecule not in experiments:
                experiments[molecule] = {}

            mut_str = row[2]
            match = regex.match(mut_str)
            if not match:
                print(f"Invalid mutation string: {mut_str}")
                continue

            source_resn = match.group(1)
            source_resi = int(match.group(2))
            target_resn = match.group(3)
            ddg = float(row[5])
            chain = row[3]

            experiment = Experiment(
                Mutation(
                    start_residue=Residue(molecule, source_resn, source_resi, chain),
                    target_resn=target_resn,
                ),
                partner_chains=row[4].split(","),
                ddg=ddg,
            )

            if chain not in experiments[molecule]:
                experiments[molecule][chain] = []

            experiments[molecule][chain].append(experiment)

    return experiments


def suggest_experiments(verified_exps_by_mol: dict[str, dict[str, list[Experiment]]]):
    models = [
        "esm1b",
        # "esm1v1", "esm1v2", "esm1v3", "esm1v4", "esm1v5",
        "esm-msa",
    ]

    suggested_experiments: dict[str, dict[str, list[Experiment]]] = {}
    for molecule, verified_exps_by_chain in verified_exps_by_mol.items():
        if molecule not in suggested_experiments:
            suggested_experiments[molecule] = {}

        cmd.fetch(molecule)

        for chain, verified_experiments in verified_exps_by_chain.items():
            if chain not in suggested_experiments[molecule]:
                suggested_experiments[molecule][chain] = []

                annotated_residues = []
                cmd.iterate(
                    f"{molecule} and chain {chain} and name CA",
                    "annotated_residues.append((oneletter, resi))",
                    space=locals(),
                )

                sequence = []
                ids = []
                for oneletter, resi in annotated_residues:
                    sequence.append(oneletter)
                    ids.append(resi)

                suggestions = get_mutation_suggestions(
                    molecule,
                    "".join(sequence),
                    ids,
                    models,
                    chain,
                )

                for suggestion in suggestions:
                    if is_residue_valid(
                        molecule,
                        chain,
                        suggestion.mutation.start_residue.id,
                    ):
                        suggested_experiments[molecule][chain].append(
                            Experiment(
                                suggestion.mutation,
                                verified_experiments[0].partner_chains,
                                0.0,
                            )
                        )
                    else:
                        print(
                            f"Filtered out mutation for invalid residue: {suggestion}"
                        )

    return suggested_experiments


def compare_mutations(file_path):
    print(f"Parsing experiments from {file_path}")
    verified_exps_by_mol = parse_experiments(file_path)
    print(verified_exps_by_mol)

    print("Getting suggestions from Efficient Evolution...")
    suggested_exps_by_mol = suggest_experiments(verified_exps_by_mol)
    print(suggested_exps_by_mol)

    for molecule, suggested_exps_by_chain in suggested_exps_by_mol.items():
        cmd.delete("all")
        cmd.fetch(molecule)

        with (
            tempfile.NamedTemporaryFile(suffix=".pdb", delete=True) as molecule_file,
        ):
            cmd.save(molecule_file.name, molecule)

            for chain, experiments in suggested_exps_by_chain.items():
                if len(experiments) == 0:
                    continue

                # We assume partner chains are the same for all experiments in the same chain
                original_affinity = compute_affinity(
                    molecule_file.name, [chain], experiments[0].partner_chains
                )

                for experiment in experiments:
                    selection_string = f"{molecule} and chain {chain} and resi {experiment.mutation.start_residue.id}"
                    cmd.select("tmp", selection_string)

                    cmd.wizard("mutagenesis")
                    cmd.do("refresh_wizard")
                    cmd.get_wizard().do_select("tmp")
                    cmd.get_wizard().set_mode(
                        one_to_three(experiment.mutation.target_resn)
                    )
                    cmd.frame(1)
                    cmd.get_wizard().apply()
                    cmd.set_wizard()

                    with tempfile.NamedTemporaryFile(
                        suffix=".pdb", delete=True
                    ) as mutated_molecule_file:
                        cmd.save(mutated_molecule_file.name, molecule)
                        mutated_affinity = compute_affinity(
                            mutated_molecule_file.name,
                            [chain],
                            experiment.partner_chains,
                        )

                    ddg = round(mutated_affinity - original_affinity, 2)
                    experiment.ddg = ddg

    print("Verified experiments:")
    for molecule, verified_exps_by_chain in verified_exps_by_mol.items():
        print(f"  {molecule}:")
        for chain, experiments in verified_exps_by_chain.items():
            print(f"    {chain}:")
            for experiment in experiments:
                experiment.serialize()
                print(
                    f"      {experiment.mutation.to_string()} with DDG {experiment.ddg}"
                )

    print("Suggested experiments:")
    for molecule, suggested_exps_by_chain in suggested_exps_by_mol.items():
        print(f"  {molecule}:")
        for chain, experiments in suggested_exps_by_chain.items():
            print(f"    {chain}:")
            for experiment in experiments:
                print(
                    f"      {experiment.mutation.to_string()} with DDG {experiment.ddg}"
                )


def main():
    if len(sys.argv) < 2:
        print("Usage: python compare_mutations.py <path_to_file>")
        return

    file_path = sys.argv[1]
    compare_mutations(file_path)


if __name__ == "__main__":
    main()
