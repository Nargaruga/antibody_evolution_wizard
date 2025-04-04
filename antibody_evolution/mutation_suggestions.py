from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass

from .mutation import Mutation
from .residue import Residue


@dataclass
class Suggestion:
    """A suggestion for a mutation to apply to an antibody,
    annotated with the number of times such mutation was suggested."""

    mutation: Mutation
    occurrences: int

    def __str__(self):
        return f"{self.mutation.to_string()} - {self.occurrences} occurrences"

    @staticmethod
    def from_EE_output(line: str, molecule: str, chain: str) -> Suggestion:
        """Parse a line of Efficient Evolution output and return a Suggestion object."""

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
    molecule_name: str, sequence: str, models: list[str], chain: str
) -> list[Suggestion]:
    """Get mutation suggestions for the given sequence using the specified models."""

    if os.name == "nt":
        prefix = ["docker", "run", "--rm", "efficient-evolution"]
    else:
        prefix = ["conda", "run", "-n", "efficient-evolution", "recommend"]

    print(f"Running inference with models {models}")
    res = subprocess.run(
        prefix
        + [
            sequence,
            "--model-names",
        ]
        + models,
        check=True,
        capture_output=True,
        text=True,
    )

    suggestions = []
    for line in res.stdout.strip().splitlines():
        if line:
            suggestions.append(Suggestion.from_EE_output(line, molecule_name, chain))

    return suggestions
