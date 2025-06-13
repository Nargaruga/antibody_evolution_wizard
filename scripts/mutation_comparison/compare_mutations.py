from __future__ import annotations

import sys
import csv
import re
from dataclasses import dataclass
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import numpy as np

from antibody_evolution.mutation import Mutation
from antibody_evolution.residue import Residue, three_to_one

sns.set(style="whitegrid")


@dataclass
class Experiment:
    """An experiment with a mutation."""

    mutation: Mutation
    partner_chains: list[str]
    ddg: float | None = None
    note: str | None = None


# TODO: two distinct parsing functions are currently needed, but should be merged


def parse_verified_experiments(file_path):
    regex = re.compile(r"([A-Z])([0-9]+)([A-Z])")
    n_skipped = 0

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
                n_skipped += 1
                continue

            source_resn = match.group(1)
            source_resi = match.group(2)
            target_resn = match.group(3)
            try:
                ddg = float(row[5])
            except ValueError:
                print(f"Skipping {molecule}/{mut_str} due to invalid DDG value.")
                n_skipped += 1
                continue

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

    # remove molecules with no valid experiments
    for molecule in list(experiments.keys()):
        if not any(experiments[molecule].values()):
            del experiments[molecule]
            continue

        for chain in list(experiments[molecule].keys()):
            if not experiments[molecule][chain]:
                del experiments[molecule][chain]

    return experiments, n_skipped


def parse_suggested_experiments(file_path):
    regex = re.compile(r"([A-Z])\/([A-Z]+)([0-9]+[A-Z]?)->([A-Z]+)")
    n_skipped = 0

    experiments: dict[str, dict[str, list[Experiment]]] = {}
    with open(file_path, "r") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            molecule = row[0]

            if molecule not in experiments:
                experiments[molecule] = {}

            mut_str = row[1]
            match = regex.match(mut_str)
            if not match:
                print(f"Invalid mutation string: {mut_str}")
                n_skipped += 1
                continue

            chain = match.group(1)
            source_resn = three_to_one(match.group(2))
            source_resi = match.group(3)
            target_resn = three_to_one(match.group(4))
            try:
                ddg = float(row[3])
            except ValueError:
                print(f"Skipping {molecule}/{mut_str} due to invalid DDG value.")
                n_skipped += 1
                continue

            experiment = Experiment(
                Mutation(
                    start_residue=Residue(molecule, source_resn, source_resi, chain),
                    target_resn=target_resn,
                ),
                partner_chains=row[2].split(","),
                ddg=ddg,
            )

            if chain not in experiments[molecule]:
                experiments[molecule][chain] = []

            experiments[molecule][chain].append(experiment)

    return experiments, n_skipped


def get_summary(experiments: dict[str, dict[str, list[Experiment]]]):
    n_mutations = 0
    n_improving_mutations = 0
    n_molecules = 0
    n_improved_molecules = 0
    best_ddg_by_molecule = {}
    avg_ddg = 0.0

    for molecule, exps_by_chain in experiments.items():
        n_molecules += 1
        best_ddg_by_molecule[molecule] = float("inf")
        for chain, exps in exps_by_chain.items():
            n_mutations += len(exps)
            for exp in exps:
                if exp.ddg is not None:
                    if exp.ddg < 0:
                        avg_ddg += exp.ddg
                        n_improving_mutations += 1
                    if exp.ddg < best_ddg_by_molecule[molecule]:
                        best_ddg_by_molecule[molecule] = exp.ddg

        if best_ddg_by_molecule[molecule] < 0:
            n_improved_molecules += 1

    avg_ddg /= n_improving_mutations if n_improving_mutations > 0 else 1

    return {
        "n_mutations": n_mutations,
        "n_improving_mutations": n_improving_mutations,
        "n_molecules": n_molecules,
        "n_improved_molecules": n_improved_molecules,
        "best_ddg_by_molecule": best_ddg_by_molecule,
        "avg_ddg": avg_ddg,
    }


def get_matched_mutations(
    verified_exps: dict[str, dict[str, list[Experiment]]],
    suggested_exps: dict[str, dict[str, list[Experiment]]],
) -> list[Mutation]:
    """Get the mutations that are present in both verified and suggested experiments."""
    matched_mutations = []
    for molecule, verified_exps_by_chain in verified_exps.items():
        if molecule not in suggested_exps:
            continue

        for chain, verified_exps in verified_exps_by_chain.items():
            if chain not in suggested_exps[molecule]:
                continue

            for verified_exp in verified_exps:
                for suggested_exp in suggested_exps[molecule][chain]:
                    if (
                        verified_exp.mutation.start_residue
                        == suggested_exp.mutation.start_residue
                        and verified_exp.mutation.target_resn
                        == suggested_exp.mutation.target_resn
                    ):
                        matched_mutations.append(verified_exp.mutation)

    return matched_mutations


def plot_best_ddg(best_ddg_by_molecule):
    molecules = set()
    data_points = {"verified": [], "suggested": []}

    for molecule, ddg in best_ddg_by_molecule.items():
        molecules.add(molecule)
        print(f"{molecule}: {ddg['verified']:.2f} {ddg['suggested']:.2f}")
        data_points["verified"].append(ddg["verified"])
        data_points["suggested"].append(ddg["suggested"])

    min_ddg = min(
        min(data_points["verified"]),
        min(data_points["suggested"]),
    )
    max_ddg = max(
        max(data_points["verified"]),
        max(data_points["suggested"]),
    )
    print(f"{data_points['verified']}")
    print(f"{data_points['suggested']}")
    print(f"min_ddg: {min_ddg}, max_ddg: {max_ddg}")

    g = sns.scatterplot(
        x=data_points["verified"],
        y=data_points["suggested"],
        palette="deep",
    )

    lower_limit = min_ddg * 1.2
    upper_limit = max_ddg * 1.2

    # draw diagonal line
    x = np.linspace(lower_limit, upper_limit, 100)
    g.plot(x, x, color="black", linestyle="--")
    g.set_xlim(lower_limit, upper_limit)
    g.set_ylim(lower_limit, upper_limit)

    # g.set_xlabel("Best ΔΔG from Verified Mutations (kcal/mol)")
    # g.set_ylabel("Best ΔΔG from Suggested Mutations (kcal/mol)")
    g.set_xlabel("ΔΔG minimo mutazioni sperimentali (kcal/mol)")
    g.set_ylabel("ΔΔG minimo mutazioni suggerite (kcal/mol)")

    plt.savefig("best_ddg.png", dpi=300)

    # plt.show()


def compare_mutations(
    verified_exps_by_molecule: dict[str, dict[str, list[Experiment]]],
    suggested_exps_by_molecule: dict[str, dict[str, list[Experiment]]],
):
    # filter out molecules that are not in both verified and suggested experiments
    print(f"Length of verified_exps_by_molecule: {len(verified_exps_by_molecule)}")
    print(f"Length of suggested_exps_by_molecule: {len(suggested_exps_by_molecule)}")
    verified_exps_by_molecule = {
        molecule: verified_exps_by_chain
        for molecule, verified_exps_by_chain in verified_exps_by_molecule.items()
        if molecule in suggested_exps_by_molecule
    }
    suggested_exps_by_molecule = {
        molecule: suggested_exps_by_chain
        for molecule, suggested_exps_by_chain in suggested_exps_by_molecule.items()
        if molecule in verified_exps_by_molecule
    }
    print(
        f"Length of verified_exps_by_molecule after filtering: {len(verified_exps_by_molecule)}"
    )
    print(
        f"Length of suggested_exps_by_molecule after filtering: {len(suggested_exps_by_molecule)}"
    )

    verified_summary = get_summary(verified_exps_by_molecule)
    suggested_summary = get_summary(suggested_exps_by_molecule)

    matched_mutations = get_matched_mutations(
        verified_exps_by_molecule, suggested_exps_by_molecule
    )

    # fmt: off
    print("# molecules:")
    print(f"  verified: {verified_summary['n_molecules']}")
    print(f"  suggested: {suggested_summary['n_molecules']}")

    print("# mutations:")
    print(f"  verified: {verified_summary['n_mutations']}")
    print(f"  suggested: {suggested_summary['n_mutations']}")

    print(r"% matched mutations:")
    print(f"  {len(matched_mutations) / (suggested_summary['n_mutations'] + verified_summary['n_mutations']) * 100:.2f}%")

    print(r"% improving mutations:")
    print(f"  verified: {verified_summary['n_improving_mutations'] / verified_summary['n_mutations'] * 100:.2f}%")
    print(f"  suggested: {suggested_summary['n_improving_mutations'] / suggested_summary['n_mutations'] * 100:.2f}%")

    print("avg DDG of improving mutations:")
    print(f"  verified: {verified_summary['avg_ddg']:.2f}")
    print(f"  suggested: {suggested_summary['avg_ddg']:.2f}")

    print(r"% improved molecules:")
    print(f"  verified: {verified_summary['n_improved_molecules']}")
    print(f"  suggested: {suggested_summary['n_improved_molecules']}")

    # print("Best DDG by molecule (verified):")
    # for molecule, ddg in verified_summary["best_ddg_by_molecule"].items():
    #     print(f"  {molecule}: {ddg:.2f}")

    # print("Best DDG by molecule (suggested):")
    # for molecule, ddg in suggested_summary["best_ddg_by_molecule"].items():
    #     print(f"  {molecule}: {ddg:.2f}")
    # fmt: on

    return {
        "verified_summary": verified_summary,
        "suggested_summary": suggested_summary,
    }


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: python compare_mutations.py <path_to_verified_exps> <path_to_generated_exps>"
        )
        return

    verified_exps_path = sys.argv[1]
    print(f"Parsing verified experiments from {verified_exps_path}")
    verified_exps_by_molecule, n_verified_skipped = parse_verified_experiments(
        verified_exps_path
    )
    print(f"Skipped {n_verified_skipped} verified experiments due to missing DDG.")

    generated_exps_path = sys.argv[2]
    print(f"Parsing generated suggestions from {generated_exps_path}...")
    suggested_exps_by_molecule, n_suggested_skipped = parse_suggested_experiments(
        generated_exps_path
    )
    print(f"Skipped {n_suggested_skipped} suggested experiments due to missing DDG.")

    print("Comparing verified and suggested experiments...")
    res = compare_mutations(verified_exps_by_molecule, suggested_exps_by_molecule)

    verified_summary = res["verified_summary"]
    suggested_summary = res["suggested_summary"]

    best_ddg_by_molecule_verified = verified_summary["best_ddg_by_molecule"]
    best_ddg_by_molecule_suggested = suggested_summary["best_ddg_by_molecule"]

    best_ddg_by_molecule = {}
    for molecule, _ in best_ddg_by_molecule_verified.items():
        if molecule in best_ddg_by_molecule_suggested:
            best_ddg_by_molecule[molecule] = {
                "verified": best_ddg_by_molecule_verified[molecule],
                "suggested": best_ddg_by_molecule_suggested[molecule],
            }

    plot_best_ddg(best_ddg_by_molecule)

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
