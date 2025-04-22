from __future__ import annotations

import os
import sys
import csv
import re
from dataclasses import dataclass
import tempfile

from pymol import cmd, CmdException

from antibody_evolution.mutation_evaluation import (
    compute_affinity,
    AffinityComputationError,
)
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
    ddg: float | None = None
    note: str | None = None


def serialize_experiments(
    experiments_by_molecule: dict[str, dict[str, list[Experiment]]], file_name: str
):
    """Serialize experiments to a CSV file."""

    with open(file_name, "w", newline="") as f:
        id = 0
        for molecule, experiments_by_chain in experiments_by_molecule.items():
            for chain, experiments in experiments_by_chain.items():
                for experiment in experiments:
                    writer = csv.writer(f)
                    start_resn = experiment.mutation.start_residue.name
                    start_resi = experiment.mutation.start_residue.id
                    target_resn = experiment.mutation.target_resn
                    writer.writerow(
                        [
                            id,
                            molecule,
                            f"{start_resn}{start_resi}{target_resn}",
                            chain,
                            ",".join(experiment.partner_chains),
                            experiment.ddg,
                            experiment.note,
                        ]
                    )
                    id += 1


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
            try:
                ddg = float(row[5])
            except ValueError:
                ddg = None
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


def compute_ddg(molecule, chain, original_affinity, experiment: Experiment) -> float:
    """Compute the mutation's DDG and return it."""

    selection_string = (
        f"{molecule} and chain {chain} and resi {experiment.mutation.start_residue.id}"
    )
    cmd.select("tmp", selection_string)

    cmd.wizard("mutagenesis")
    cmd.do("refresh_wizard")
    try:
        cmd.get_wizard().do_select("tmp")
    except CmdException as e:
        # Not sure why this happens. A molecule that causes this is 3L5X
        msg = f"Skipping mutation {experiment.mutation} due to error selecting residue: {e}"
        print(msg)
        experiment.note = msg
        return experiment
    cmd.get_wizard().set_mode(one_to_three(experiment.mutation.target_resn))
    cmd.frame(1)
    cmd.get_wizard().apply()
    cmd.set_wizard()

    mutated_molecule_file_handle, mutated_molecule_file_path = tempfile.mkstemp(
        suffix=".pdb"
    )
    cmd.save(mutated_molecule_file_path, molecule)
    try:
        mutated_affinity = compute_affinity(
            mutated_molecule_file_path,
            [chain],
            experiment.partner_chains,
        )
    except AffinityComputationError as e:
        msg = f"Skipping mutation {experiment.mutation} due to error computing new affinity: {e}"
        print(msg)
        experiment.note = msg
        return experiment
    finally:
        os.close(mutated_molecule_file_handle)
        os.remove(mutated_molecule_file_path)

    return round(mutated_affinity - original_affinity, 2)


def annotate_experiments(
    experiments_by_molecule: dict[str, dict[str, list[Experiment]]],
) -> dict[str, dict[str, list[Experiment]]]:
    """Annotate experiments with DDG values."""

    for molecule, experiments_by_chain in experiments_by_molecule.items():
        cmd.delete("all")
        cmd.fetch(molecule)

        molecule_file_handle, molecule_file_path = tempfile.mkstemp(suffix=".pdb")
        cmd.save(molecule_file_path, molecule)

        for chain, experiments in experiments_by_chain.items():
            if len(experiments) == 0:
                continue

            # We assume partner chains are the same for all experiments in the same chain
            try:
                original_affinity = compute_affinity(
                    molecule_file_path, [chain], experiments[0].partner_chains
                )
            except AffinityComputationError as e:
                msg = f"Skipping chain {chain} due to error computing original affinity: {e}"
                print(msg)
                for experiment in experiments:
                    experiment.note = msg
                continue

            for experiment in experiments:
                experiment.ddg = compute_ddg(
                    molecule,
                    chain,
                    original_affinity,
                    experiment,
                )

        os.close(molecule_file_handle)
        os.remove(molecule_file_path)

    return experiments_by_molecule


def suggest_experiments(verified_exps_by_mol: dict[str, dict[str, list[Experiment]]]):
    models = [
        "esm1b",
        "esm1v1",
        "esm1v2",
        "esm1v3",
        "esm1v4",
        "esm1v5",
        # "esm-msa",
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

                if len(sequence) == 0:
                    print(f"No residues found for {molecule} chain {chain}")
                    continue

                print(f"Generating suggestions for {molecule} chain {chain}...")
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
                        suggestion.mutation.start_residue,
                    ):
                        suggested_experiments[molecule][chain].append(
                            Experiment(
                                suggestion.mutation,
                                verified_experiments[0].partner_chains,
                            )
                        )
                    else:
                        print(
                            f"Filtered out mutation for invalid residue: {suggestion}"
                        )

    return suggested_experiments


def compare_mutations(
    verified_exps_by_molecule: dict[str, dict[str, list[Experiment]]],
    suggested_exps_by_molecule: dict[str, dict[str, list[Experiment]]],
):
    for molecule, verified_exps_by_chain in verified_exps_by_molecule.items():
        if molecule not in suggested_exps_by_molecule:
            print(f"Skipping {molecule} as it has no suggested experiments")
            continue

        for chain, verified_experiments in verified_exps_by_chain.items():
            if chain not in suggested_exps_by_molecule[molecule]:
                print(f"Skipping {chain} as it has no suggested experiments")
                continue

            print(
                f"---Comparing verified and suggested experiments for {molecule} chain {chain}---"
            )
            suggested_experiments = suggested_exps_by_molecule[molecule][chain]

            best_verified_mutation = None
            best_verified_ddg = float("inf")
            for experiment in verified_experiments:
                if experiment.ddg is not None and experiment.ddg < best_verified_ddg:
                    best_verified_ddg = experiment.ddg
                    best_verified_mutation = experiment.mutation
            print(f"{len(verified_experiments)} verified mutations.")
            print(
                f"Best verified mutation: {best_verified_mutation.to_string()} with DDG {best_verified_ddg}"
            )

            best_suggested_mutation = None
            best_suggested_ddg = float("inf")
            for experiment in suggested_experiments:
                if experiment.ddg is not None and experiment.ddg < best_suggested_ddg:
                    best_suggested_ddg = experiment.ddg
                    best_suggested_mutation = experiment.mutation
            print(f"{len(suggested_experiments)} suggested mutations.")
            print(
                f"Best suggested mutation: {best_suggested_mutation.to_string()} with DDG {best_suggested_ddg}"
            )

            print("-------")


def main():
    if len(sys.argv) <= 2:
        print(
            "Usage: python compare_mutations.py <compare|annotate> <path_to_file> [path_to_file2]"
        )
        return

    command = sys.argv[1]
    if command not in ["compare", "annotate"]:
        print("Invalid command. Use 'compare' or 'annotate'.")
        return

    file_path = sys.argv[2]
    print(f"Parsing verified experiments from {file_path}")
    verified_exps_by_molecule = parse_experiments(file_path)

    if len(sys.argv) > 3:
        print(f"Parsing generated suggestions from {sys.argv[3]}...")
        suggested_exps_by_molecule = parse_experiments(sys.argv[3])
    else:
        print("Getting suggestions from Efficient Evolution...")
        suggested_exps_by_molecule = suggest_experiments(verified_exps_by_molecule)

    if command == "annotate":
        print("Annotating verified experiments with DDG values...")
        verified_exps_by_molecule = annotate_experiments(verified_exps_by_molecule)
        serialize_experiments(
            verified_exps_by_molecule, "annotated_verified_experiments.csv"
        )

        print("Annotating generated experiments with DDG values...")
        suggested_exps_by_molecule = annotate_experiments(suggested_exps_by_molecule)
        serialize_experiments(
            suggested_exps_by_molecule, "annotated_generated_experiments.csv"
        )
    elif command == "compare":
        print("Comparing verified and suggested experiments...")
        compare_mutations(verified_exps_by_molecule, suggested_exps_by_molecule)

    # print("Verified experiments:")
    # for molecule, verified_exps_by_chain in verified_exps_by_molecule.items():
    #     print(f"  {molecule}:")
    #     for chain, experiments in verified_exps_by_chain.items():
    #         print(f"    {chain}:")
    #         for experiment in experiments:
    #             try:
    #                 print(
    #                     f"      {experiment.mutation.to_string()} with DDG {experiment.ddg}"
    #                 )
    #             except ValueError as e:
    #                 print(f"Error printing experiment {experiment}: {e}")

    # print("Suggested experiments:")
    # for molecule, suggested_exps_by_chain in suggested_exps_by_molecule.items():
    #     print(f"  {molecule}:")
    #     for chain, experiments in suggested_exps_by_chain.items():
    #         print(f"    {chain}:")
    #         for experiment in experiments:
    #             try:
    #                 print(
    #                     f"      {experiment.mutation.to_string()} with DDG {experiment.ddg}"
    #                 )
    #             except ValueError as e:
    #                 print(f"Error printing experiment {experiment}: {e}")


if __name__ == "__main__":
    main()
