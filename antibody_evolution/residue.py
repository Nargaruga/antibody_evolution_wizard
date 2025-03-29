from dataclasses import dataclass


# Map single letter amino acid codes to three letter codes
RESIDUE_NAME_MAPPING = {
    "C": "CYS",
    "D": "ASP",
    "S": "SER",
    "Q": "GLN",
    "K": "LYS",
    "I": "ILE",
    "P": "PRO",
    "T": "THR",
    "F": "PHE",
    "N": "ASN",
    "G": "GLY",
    "H": "HIS",
    "L": "LEU",
    "R": "ARG",
    "W": "TRP",
    "A": "ALA",
    "V": "VAL",
    "E": "GLU",
    "Y": "TYR",
    "M": "MET",
}


def one_to_three(oneletter):
    """Convert a single letter amino acid code to a three letter code."""
    return RESIDUE_NAME_MAPPING[oneletter]


@dataclass
class Residue:
    """A residue in a molecule."""

    molecule: str
    name: str
    id: int
    chain: str

    def __eq__(self, other):
        return (
            self.molecule == str(other.molecule)
            and self.name == str(other.name)
            and self.id == int(other.id)
            and self.chain == str(other.chain)
        )

    def get_selection_str(self):
        """Get the PyMOL selection string for the residue."""
        return f"{self.molecule} and resi {self.id} and chain {self.chain}"
