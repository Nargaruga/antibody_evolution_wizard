from dataclasses import dataclass

from pymol import cmd


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

    if oneletter not in RESIDUE_NAME_MAPPING:
        raise ValueError(f"Invalid single letter code: {oneletter}")

    return RESIDUE_NAME_MAPPING[oneletter]


@dataclass
class Residue:
    """A residue in a molecule."""

    molecule: str
    name: str
    id: str
    chain: str

    def __eq__(self, other):
        return (
            self.molecule == other.molecule
            and self.name == other.name
            and self.id == other.id
            and self.chain == other.chain
        )

    def get_selection_str(self):
        """Get the PyMOL selection string for the residue."""
        return f"{self.molecule} and resi {self.id} and chain {self.chain}"

    def is_valid(self, pymol_cmd=None) -> bool:
        try:
            one_to_three(self.name)
        except ValueError:
            print(f"Invalid residue code: {self.name}")
            return False

        if pymol_cmd is None:
            pymol_cmd = cmd

        pymol_cmd.select(
            "temp",
            f"byres ({self.molecule} and chain {self.chain} and resi {self.id}) and name CA",
        )

        if pymol_cmd.count_atoms("temp") <= 1:
            return False

        pymol_cmd.delete("temp")

        return True
