from __future__ import annotations

import subprocess
from dataclasses import dataclass

from .mutation import Mutation
from .residue import Residue


@dataclass
class Suggestion:
    mutation: Mutation
    occurrences: int

    def __str__(self):
        return f"{self.mutation.to_string()} - {self.occurrences} occurrences"

    @staticmethod
    def from_EE_output(line: str, molecule: str, chain: str) -> Suggestion:
        # Efficient Evolution outputs strings in the form
        # [start residue name][start residue id][target residue name] [occurrences]
        # example: E1M 2
        mut_str, count = line.split()
        start_resn = mut_str[0]
        start_resi = mut_str[1:-1]
        target = mut_str[-1]

        return Suggestion(
            mutation=Mutation(
                Residue(molecule, start_resn, int(start_resi), chain), target
            ),
            occurrences=int(count),
        )


def get_mutation_suggestions(
    molecule_name: str, sequence: str, models: list[str], antibody_chain: str
) -> list[Suggestion]:
    print(f"Running inference with models {models}")

    try:
        res = subprocess.run(
            [
                "conda",
                "run",
                "-n",
                "efficient-evolution",
                "recommend",
                sequence,
                "--model-names",
            ]
            + models,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Something went wrong while calling Efficient Evolution: {e}")
        exit(-1)

    suggestions = []
    for line in res.stdout.strip().splitlines():
        if line:
            suggestions.append(
                Suggestion.from_EE_output(line, molecule_name, antibody_chain)
            )

    return suggestions
