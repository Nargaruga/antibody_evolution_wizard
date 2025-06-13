from __future__ import annotations

import os
import sys
import csv
import re
from dataclasses import dataclass
import tempfile

from pymol import cmd

from antibody_evolution.mutation_evaluation import (
    compute_affinity,
    compute_ddg,
    AffinityComputationError,
    DDGComputationError,
)
from antibody_evolution.mutation_suggestions import (
    get_mutation_suggestions,
)
from antibody_evolution.mutation import Mutation
from antibody_evolution.residue import Residue, oneletter_exists


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
            partner_chains = experiments[0].partner_chains

            try:
                original_affinity = compute_affinity(
                    molecule_file_path, [chain], partner_chains
                )
            except AffinityComputationError as e:
                msg = f"Skipping chain {chain} due to error computing original affinity: {e}"
                print(msg)
                for experiment in experiments:
                    experiment.ddg = None
                    experiment.note = msg
                continue

            for experiment in experiments:
                if experiment.mutation.start_residue.is_valid(cmd) and oneletter_exists(
                    experiment.mutation.target_resn
                ):
                    print(
                        f"Annotating {experiment.mutation.to_string()} in {molecule} chain {chain}..."
                    )
                    experiment = annotate_experiment(
                        molecule_file_path,
                        chain,
                        original_affinity,
                        experiment,
                        cmd,
                    )
                else:
                    msg = f"Filtered out mutation for invalid residue: {experiment}"
                    print(msg)
                    experiment.ddg = None
                    experiment.note = msg

        os.close(molecule_file_handle)
        os.remove(molecule_file_path)

    return experiments_by_molecule


def annotate_experiment(
    molecule_file_path, chain, original_affinity, experiment: Experiment, pymol_cmd
) -> Experiment:
    try:
        print("Computing DDG...")
        ddg = compute_ddg(
            molecule_file_path,
            chain,
            original_affinity,
            experiment.mutation,
            [chain],
            experiment.partner_chains,
            pymol_cmd,
        )

        print(f"DDG for {experiment.mutation.to_string()} is {ddg}")

        experiment.ddg = ddg
    except DDGComputationError as e:
        print(f"DDG computation for {experiment.mutation.to_string()} failed.")
        experiment.note = str(e)
    finally:
        return experiment


def suggest_experiments(
    verified_exps_by_mol: dict[str, dict[str, list[Experiment]]], file_name: str
):
    models = [
        "esm1b",
        "esm1v1",
        "esm1v2",
        "esm1v3",
        "esm1v4",
        "esm1v5",
        # "esm-msa",
    ]

    with open(file_name, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["molecule", "mutation", "partner_chains", "ddg", "note"])

        suggested_experiments: dict[str, dict[str, list[Experiment]]] = {}
        for molecule, verified_exps_by_chain in verified_exps_by_mol.items():
            print(f"Processing {molecule}...")
            cmd.delete("all")

            if molecule not in suggested_experiments:
                suggested_experiments[molecule] = {}

            cmd.fetch(molecule)
            chain_seqs = {}
            chain_resids = {}
            for chain in cmd.get_chains(molecule):
                chain_seqs[chain] = []
                chain_resids[chain] = []

                cmd.iterate(
                    f"{molecule} and chain {chain} and name CA",
                    "chain_seqs[chain].append(oneletter)",
                    space=locals(),
                )
                cmd.iterate(
                    f"{molecule} and chain {chain} and name CA",
                    "chain_resids[chain].append(resi)",
                    space=locals(),
                )

            molecule_file_handle, molecule_file_path = tempfile.mkstemp(suffix=".pdb")
            cmd.save(molecule_file_path, molecule)

            for chain, verified_experiments in verified_exps_by_chain.items():
                if len(verified_experiments) == 0:
                    print(
                        f"No verified experiments for {molecule} chain {chain}, skipping..."
                    )
                    continue

                if chain not in chain_seqs or chain not in chain_resids:
                    print(f"Chain {chain} not found in {molecule}, skipping...")
                    continue

                partner_chains = verified_experiments[0].partner_chains
                if len(partner_chains) == 0:
                    print(
                        f"No partner chains for {molecule} chain {chain}, skipping..."
                    )
                    continue

                print(f"Processing chain {chain} with partners {partner_chains}...")
                if chain not in suggested_experiments[molecule]:
                    suggested_experiments[molecule][chain] = []

                sequence = chain_seqs[chain]
                ids = chain_resids[chain]
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

                try:
                    original_affinity = compute_affinity(
                        molecule_file_path,
                        [chain],
                        partner_chains,
                    )
                    print(
                        f"Original affinity for {molecule} chain {chain} with partners {partner_chains} is {original_affinity}"
                    )
                except AffinityComputationError as e:
                    msg = f"Skipping chain {chain} due to error computing original affinity: {e}"
                    print(msg)
                    continue

                for suggestion in suggestions:
                    if suggestion.mutation.start_residue.is_valid(cmd):
                        experiment = annotate_experiment(
                            molecule_file_path,
                            chain,
                            original_affinity,
                            Experiment(
                                suggestion.mutation,
                                partner_chains,
                            ),
                            cmd,
                        )

                        print(
                            f"Annotated {experiment.mutation.to_string()} in {molecule} chain {chain} with DDG {experiment.ddg}"
                        )

                        writer.writerow(
                            [
                                molecule,
                                experiment.mutation.to_string(),
                                ",".join(experiment.partner_chains),
                                experiment.ddg if experiment.ddg is not None else "N/A",
                                experiment.note,
                            ]
                        )

                        suggested_experiments[molecule][chain].append(experiment)
                    else:
                        print(
                            f"Filtered out mutation for invalid residue: {suggestion}"
                        )

            os.close(molecule_file_handle)
            os.remove(molecule_file_path)

    return suggested_experiments


def print_summary(experiments: dict[str, list[Experiment]]):
    for chain, experiments in experiments.items():
        print(f"Processing chain {chain}:")
        n_improving_mutations = 0
        best_mutation = None
        lowest_ddg = float("inf")
        for experiment in experiments:
            if experiment.ddg < 0:
                n_improving_mutations += 1
            if experiment.ddg is not None and experiment.ddg < lowest_ddg:
                lowest_ddg = experiment.ddg
                best_mutation = experiment.mutation
        print(f"{len(experiments)} mutations.")
        print(f"Of those, {n_improving_mutations} have a DDG < 0.")
        print(f"Best mutation: {best_mutation.to_string()} with DDG {lowest_ddg}")


def main():
    if len(sys.argv) <= 1:
        print("Usage: python annotate_mutations.py <path_to_experiments>")
        return

    file_path = sys.argv[1]
    print(f"Parsing verified experiments from {file_path}")
    verified_exps_by_molecule = parse_experiments(file_path)

    print("Annotating experiments...")
    annotated_verified_exps_by_molecule = annotate_experiments(
        verified_exps_by_molecule
    )
    serialize_experiments(
        annotated_verified_exps_by_molecule, "annotated_verified_experiments.csv"
    )

    # print("Getting suggestions from Efficient Evolution...")
    # suggest_experiments(verified_exps_by_molecule, "annotated_verified_experiments.csv")


if __name__ == "__main__":
    main()
