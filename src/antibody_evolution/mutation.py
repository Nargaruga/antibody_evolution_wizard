from dataclasses import dataclass

from .residue import Residue, one_to_three


@dataclass
class Mutation:
    """An amino-acid mutation."""

    start_residue: Residue
    target_resn: str

    def to_string(self):
        """Get a string representation of the mutation."""
        try:
            oneletter = one_to_three(self.start_residue.name)
        except ValueError:
            print(
                f"Invalid residue code: {self.start_residue.name}, using '?' placeholder"
            )
            oneletter = "?"

        return f"{self.start_residue.chain}/{oneletter}{self.start_residue.id}->{one_to_three(self.target_resn)}"
