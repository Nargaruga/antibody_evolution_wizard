from dataclasses import dataclass

from .residue import Residue, one_to_three


@dataclass
class Mutation:
    """An amino-acid mutation."""

    start_residue: Residue
    target_resn: str

    def to_string(self):
        """Get a string representation of the mutation."""
        return f"{self.start_residue.chain}/{one_to_three(self.start_residue.name)}{self.start_residue.id}->{one_to_three(self.target_resn)}"
