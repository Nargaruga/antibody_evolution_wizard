import os
from enum import IntEnum, auto
import threading
import subprocess
from pathlib import Path
from dataclasses import dataclass

from pymol.wizard import Wizard
from pymol import cmd

import yaml

from efficient_evolution.amis import reconstruct_multi_models


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
    segment: str

    def get_selection_str(self):
        """Get the PyMOL selection string for the residue."""
        return f"{self.molecule} and resi {self.id} and chain {self.chain} and segi {self.segment}"


@dataclass
class Mutation:
    """An amino-acid mutation."""

    molecule_name: str
    start_residue: Residue
    target: str
    occurrences: int

    def to_string(self):
        """Get a string representation of the mutation."""
        return f"{self.start_residue.segment}/{self.start_residue.chain}/{one_to_three(self.start_residue.name)}{self.start_residue.id}->{one_to_three(self.target)}"


class WizardState(IntEnum):
    """The possible states of the wizard."""

    INITIALIZING = auto()
    READY = auto()
    MOLECULE_SELECTED = auto()
    CHAIN_SELECTED = auto()
    RUNNING_INFERENCE = auto()
    SUGGESTIONS_GENERATED = auto()
    MUTATION_SELECTED = auto()
    MUTATION_APPLIED = auto()


class Antibody_evolution(Wizard):
    """A wizard for suggesting mutations to an antibody."""

    def __init__(self, _self=cmd):
        Wizard.__init__(self, _self)
        self.status = WizardState.INITIALIZING
        self.models = ["esm1b"]
        self.molecule = None
        self.chain = None
        self.mutations = {}
        self.selected_mutation = None
        self.binding_affinity = 0.0
        self.populate_molecule_choices()

        # Load the model in a separate thread
        initializer_thread = threading.Thread(target=self.load_model)
        initializer_thread.start()

    def load_model(self):
        # check which models were requested
        location = Path(__file__).parent
        with open(os.path.join(location, "config", "models.yml")) as f:
            contents = yaml.load(f, Loader=yaml.FullLoader)
            self.models = contents["models"]

        # TODO download the model
        self.status = WizardState.READY
        cmd.refresh_wizard()

    def get_prompt(self):  # type: ignore
        """Return the prompt for the current state of the wizard."""

        self.prompt = ["Binding affinity: %s" % self.binding_affinity]
        if self.status == WizardState.INITIALIZING:
            self.prompt.append("Initializing, please wait...")
        elif self.status == WizardState.READY:
            self.prompt.append("Select a molecule.")
        elif self.status == WizardState.MOLECULE_SELECTED:
            self.prompt.append("Select a chain.")
        elif self.status == WizardState.CHAIN_SELECTED:
            self.prompt.append(
                f"Run to generate mutation suggestions for {self.molecule}, chain {self.chain}."
            )
        elif self.status == WizardState.RUNNING_INFERENCE:
            self.prompt.append("Running inference, please wait...")
        elif self.status == WizardState.SUGGESTIONS_GENERATED:
            self.prompt.append("Select a mutation to apply.")
        elif (
            self.status == WizardState.MUTATION_SELECTED
            and self.selected_mutation is not None
        ):
            self.prompt.append(
                "Apply the mutation %s?" % self.selected_mutation.to_string()
            )
        elif self.status == WizardState.MUTATION_APPLIED:
            self.prompt.append("Mutation applied successfully.")

        return self.prompt

    def populate_molecule_choices(self):
        """Populate the menu with the available molecules in the session."""

        molecules = cmd.get_names("objects")
        self.menu["molecule"] = [[2, "Molecule", ""]]
        for m in molecules:
            self.menu["molecule"].append(
                [
                    1,
                    m,
                    'cmd.get_wizard().set_molecule("' + m + '")',
                ]
            )

    def populate_chain_choices(self):
        """Populate the menu with the available chains in the selected molecule."""

        if self.molecule is None:
            return

        chains = cmd.get_chains(self.molecule)
        self.menu["chain"] = [[2, "Chain", ""]]
        for c in chains:
            self.menu["chain"].append(
                [
                    1,
                    c,
                    'cmd.get_wizard().set_chain("' + c + '")',
                ]
            )

    def set_molecule(self, molecule):
        """Set the molecule to generate mutation suggestions for."""
        self.molecule = molecule
        self.status = WizardState.MOLECULE_SELECTED
        self.populate_chain_choices()
        cmd.refresh_wizard()

    def set_chain(self, chain):
        """Set the chain to generate mutation suggestions for."""
        self.chain = chain
        self.status = WizardState.CHAIN_SELECTED
        cmd.refresh_wizard()

    def run(self):
        """Run the wizard to generate suggestions for the selected molecule."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        residues = []

        def record_residue(molecule, oneletter, resi, chain, segi):
            """Record all the necessary information for each residue in the molecule."""
            residues.append(Residue(molecule, oneletter, resi, chain, segi))

        context = {
            "record_residue": record_residue,
            "molecule": self.molecule,
        }
        cmd.iterate(
            f"{self.molecule} and chain {self.chain} and name CA",
            "record_residue(molecule,oneletter,resi,chain,segi)",
            space=context,
        )

        self.status = WizardState.RUNNING_INFERENCE
        cmd.refresh_wizard()

        # Run inference on a separate thread
        worker_thread = threading.Thread(target=self.run_inference, args=[residues])
        worker_thread.start()

    def run_inference(self, residues):
        """Run the inference to generate mutation suggestions."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        sequence = ""
        for residue in residues:
            sequence += residue.name

        print(f"Running inference with models {self.models}")
        mutations_models = reconstruct_multi_models(
            sequence,
            self.models,
        )
        suggestions = []

        # Mutations k are in the form (before_idx, before_name, after_name)
        for k, v in sorted(mutations_models.items(), key=lambda item: -item[1]):
            suggestions.append(Mutation(self.molecule, residues[k[0]], k[2], v))

        self.populate_mutation_choices(suggestions)

        self.status = WizardState.SUGGESTIONS_GENERATED
        cmd.refresh_wizard()

    def populate_mutation_choices(self, suggestions):
        """Populate the menu with the generated mutation suggestions."""

        self.menu["mutations"] = [[2, "Mutations", ""]]
        suggestions.sort(key=lambda x: -x.occurrences)
        for mut in suggestions:
            self.mutations[mut.to_string()] = mut
            self.menu["mutations"].append(
                [
                    1,
                    f"{mut.to_string()} ({mut.occurrences})",
                    'cmd.get_wizard().set_mutation("' + mut.to_string() + '")',
                ]
            )

    def set_mutation(self, mutation_str):
        """Set the selected mutation to apply."""
        mutation = self.mutations[mutation_str]
        self.selected_mutation = mutation
        cmd.select("mutation", mutation.start_residue.get_selection_str())
        cmd.zoom("mutation", 2)

        self.status = WizardState.MUTATION_SELECTED
        cmd.refresh_wizard()

    def evaluate_binding_affinity(self):
        """Evaluate the binding affinity using Prodigy."""

        cmd.save(f"{self.molecule}.pdb", f"{self.molecule}")

        # TODO: specify chains
        res = subprocess.run(
            ["prodigy", f"{self.molecule}.pdb", "--pymol_selection", "--quiet"],
            capture_output=True,
        )
        self.parse_prodigy_output(res.stdout)

        os.remove(f"{self.molecule}.pdb")

        cmd.refresh_wizard()

    def parse_prodigy_output(self, output):
        """Parse the output of Prodigy to get the binding affinity."""
        print(output)

        if output is None or len(output.split()) != 2:
            print("Error: could not parse Prodigy output.")
            return

        self.binding_affinity = float(output.split()[1])

    def apply_mutation(self):
        """Apply the selected mutation to the molecule."""

        if self.selected_mutation is None:
            print("Please select a mutation.")
            return

        cmd.wizard("mutagenesis")
        cmd.do("refresh_wizard")
        cmd.get_wizard().set_mode("%s" % one_to_three(self.selected_mutation.target))
        selection = "/%s/%s/%s/%s" % (
            self.molecule,
            self.selected_mutation.start_residue.segment,
            self.selected_mutation.start_residue.chain,
            self.selected_mutation.start_residue.id,
        )
        cmd.select("sele", selection)
        cmd.get_wizard().do_select("sele")
        cmd.frame(str(1))
        cmd.get_wizard().apply()
        cmd.set_wizard()

        print(f"Applied mutation {self.selected_mutation.to_string()}.")

        # Remove the applied mutation from the list of suggestions
        del self.mutations[self.selected_mutation.to_string()]
        self.selected_mutation = None
        self.populate_mutation_choices(list(self.mutations.values()))

        self.status = WizardState.MUTATION_APPLIED
        cmd.refresh_wizard()

    def get_panel(self):
        """Return the menu panel for the wizard."""

        if self.molecule is None:
            molecule_label = "Choose molecule"
        else:
            molecule_label = self.molecule

        if self.mutations is None:
            mutations_label = "-"
        elif self.selected_mutation is None:
            mutations_label = "Select Mutation"
        else:
            mutations_label = self.selected_mutation.to_string()

        options = [
            [1, "Mutation Suggestions", ""],
        ]

        if self.status >= WizardState.READY:
            options.append(
                [3, molecule_label, "molecule"],
            )

        if self.status >= WizardState.MOLECULE_SELECTED:
            options.append(
                [
                    2,
                    "Evaluate Affinity",
                    "cmd.get_wizard().evaluate_binding_affinity()",
                ],
            )

            options.append(
                [3, self.chain, "chain"],
            )

        if self.status >= WizardState.CHAIN_SELECTED:
            options.append(
                [2, "Run Inference", "cmd.get_wizard().run()"],
            )

        if self.status >= WizardState.SUGGESTIONS_GENERATED:
            options.append(
                [3, mutations_label, "mutations"],
            )

        if self.status >= WizardState.MUTATION_SELECTED:
            options.append(
                [2, "Apply Mutation", "cmd.get_wizard().apply_mutation()"],
            )

        options.append([2, "Dismiss", "cmd.set_wizard()"])

        return options
