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

RESIDUE_NAME_REVERSE_MAPPING = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}


def one_to_three(oneletter):
    """Convert a single letter amino acid code to a three letter code."""

    if oneletter not in RESIDUE_NAME_MAPPING:
        raise ValueError(f"Invalid single letter code: {oneletter}")

    return RESIDUE_NAME_MAPPING[oneletter]


def three_to_one(threeletter):
    """Convert a three letter amino acid code to a single letter code."""

    if threeletter not in RESIDUE_NAME_REVERSE_MAPPING:
        raise ValueError(f"Invalid three letter code: {threeletter}")

    return RESIDUE_NAME_REVERSE_MAPPING[threeletter]


def resn_exists(resn):
    """Check if a residue name exists in the mapping."""
    return resn in RESIDUE_NAME_MAPPING.values()


def oneletter_exists(oneletter):
    """Check if a one-letter code exists in the mapping."""
    return oneletter in RESIDUE_NAME_MAPPING.keys()


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
        if not oneletter_exists(self.name):
            print(f"Invalid residue name: {self.name}")
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
