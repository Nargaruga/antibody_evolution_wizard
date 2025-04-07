import sys
import os
import pathlib
import re
import csv

from pymol import cmd

from .mutation import Suggestion

class Experiment:
    partner_chains: list[str]
    mutation: Suggestion


def experiments_from_csv(file_path: str) -> dict[str, dict[str, list[Experiment]]]:
    regex = re.compile(r"([A-Z])([0-9]+)([A-Z])")

    print(f"Reading mutations from {file_path}")
    experimental_mutations: dict[str, dict[str, list[Experiment]]] = {}
    with open(file_path, "r") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            molecule = row[1]
            if molecule not in experimental_mutations:
                experimental_mutations[molecule] = {}

            mut_str = row[2]
            match = regex.match(mut_str)
            if not match:
                print(f"Invalid mutation string: {mut_str}")
                continue

            source_resn = match.group(1)
            source_resi = int(match.group(2))
            target_resn = match.group(3)
            ddg = float(row[5])
            mutation = Suggestion(
                
                target_resn=target_resn,
                ddg=ddg,
            )

            experiment = Experiment()
            experiment.partner_chains = row[4].split(",")
            experiment.mutation = mutation

            chain = row[3]
            if chain not in experimental_mutations[molecule]:
                experimental_mutations[molecule][chain] = []

            experimental_mutations[molecule][chain].append(experiment)

    return experimental_mutations

def main():
    input_file = sys.argv[1]
    experiments = experiments_from_csv(input_file)

    os.makedirs("pdbs", exist_ok=True)

    for molecule, experiments_by_chain in experiments.items():
        print(f"Molecule: {molecule}")

        cmd.fetch(molecule)
        cmd.save(os.path.join("pdbs", f"{molecule}.pdb"))

        for chain, experiments in experiments_by_chain.items():
            print(f"Chain: {chain}")
            for experiment in experiments:
                print(
                    f"{experiment.mutation} with partner chains {experiment.partner_chains}"
                )

                # Compute the ddG for the mutation
                ddg = compute_ddg(
                    os.path.join(
                        "pdbs",
                        f"{molecule}.pdb",
                    ),
                    experiment.mutation,
                    chain,
                    experiment.partner_chains,
                )
                print(f"Prodigy ddG: {ddg}")
                print(f"Experimental ddG: {experiment.mutation.ddg}")

    # Add a column to the CSV
    output_file = os.path.join(pathlib.Path(input_file).parent, "dataset_prodigy.csv")
    insert_prodigy_affinity_col(output_file)



if __name__ == "__main__":
    main()
