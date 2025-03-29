from __future__ import annotations

import os
from enum import IntEnum, auto
import threading
import subprocess
from pathlib import Path
from dataclasses import dataclass


from pymol.wizard import Wizard
from pymol.wizarding import WizardError
from pymol import cmd

import yaml

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


@dataclass
class Mutation:
    """An amino-acid mutation."""

    molecule_name: str
    start_residue: Residue
    target: str
    occurrences: int

    def to_string(self):
        """Get a string representation of the mutation."""
        return f"{self.start_residue.chain}/{one_to_three(self.start_residue.name)}{self.start_residue.id}->{one_to_three(self.target)}"

    @staticmethod
    def from_EE_output(line: str, molecule, chain: str) -> Mutation:
        # Efficient Evolution outputs strings in the form
        # [start residue name][start residue id][target residue name] [occurrences]
        # example: E1M 2
        mut_str, count = line.split()
        start_resn = mut_str[0]
        start_resi = mut_str[1:-1]
        target = mut_str[-1]

        return Mutation(
            molecule,
            Residue(molecule, start_resn, int(start_resi), chain),
            target,
            int(count),
        )


@dataclass
class HistoryEntry:
    """An entry in the history of the wizard."""

    binding_affinity: float
    mutations: list[Mutation]


class WizardState(IntEnum):
    """The possible states of the wizard."""

    INITIALIZING = auto()
    READY = auto()
    MOLECULE_SELECTED = auto()
    ANTIBODY_CHAIN_SELECTED = auto()
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
        self.antibody_chain = None
        self.antigen_chains = []
        self.mutations = {}
        self.selected_mutation = None
        self.history: list[HistoryEntry] = []
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

        self.prompt = []
        if self.status == WizardState.INITIALIZING:
            self.prompt.append("Initializing, please wait...")
        elif self.status == WizardState.READY:
            self.prompt.append("Select a molecule.")
        elif self.status == WizardState.MOLECULE_SELECTED:
            self.prompt.append("Select a chain.")
        elif self.status == WizardState.ANTIBODY_CHAIN_SELECTED:
            self.prompt.append(
                f"Run to generate mutation suggestions for {self.molecule}, chain {self.antibody_chain}."
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
            print("Please select a molecule.")
            return

        chains = cmd.get_chains(self.molecule)
        self.menu["antibody_chain"] = [[2, "Antibody Chain", ""]]
        self.menu["antigen_chains"] = [[2, "Antigen Chains", ""]]
        for c in chains:
            self.menu["antibody_chain"].append(
                [
                    1,
                    c,
                    'cmd.get_wizard().set_antibody_chain("' + c + '")',
                ]
            )

            self.menu["antigen_chains"].append(
                [
                    1,
                    c,
                    'cmd.get_wizard().set_antigen_chain("' + c + '")',
                ]
            )

    def set_molecule(self, molecule):
        """Set the molecule to generate mutation suggestions for."""
        self.molecule = molecule
        self.status = WizardState.MOLECULE_SELECTED
        self.populate_chain_choices()
        cmd.refresh_wizard()

    def set_antibody_chain(self, chain):
        """Set the antibody chain to generate mutation suggestions for."""
        self.antibody_chain = chain
        self.status = WizardState.ANTIBODY_CHAIN_SELECTED
        cmd.refresh_wizard()

    def set_antigen_chain(self, chain):
        """Set the antigen chain."""
        if chain in self.antigen_chains:
            self.antigen_chains.remove(chain)
        else:
            self.antigen_chains.append(chain)
        self.status = WizardState.ANTIBODY_CHAIN_SELECTED
        cmd.refresh_wizard()

    def find_mutation_for(self, sel: str):
        residues = []

        context = {
            "molecule": self.molecule,
            "residues": residues,
            "Residue": Residue,
        }
        cmd.iterate(
            f"{sel} and name CA",
            "residues.append(Residue(molecule,oneletter,resi,chain))",
            space=context,
        )

        if len(residues) == 0:
            raise WizardError("No residues selected.")

        residue = residues[0]
        for mutation_str, mutation in self.mutations.items():
            if mutation.start_residue == residue:
                return mutation_str

        return None

    def do_select(self, name: str):
        """
        Select a residue for mutation.

        :param selection: A PyMOL selection
        :type selection: string
        """
        try:
            mutation = self.find_mutation_for(name)
            if mutation is None:
                print("No mutation available the selected residue.")
            else:
                self.set_mutation(mutation)
        except WizardError as e:
            print(e)

        cmd.delete(name)

    def highlight_mutations(self):
        for mutation_str, mutation in self.mutations.items():
            r = mutation.start_residue
            cmd.color(
                "cyan",
                f"{self.molecule} and resi {r.id} and chain {r.chain} and state {cmd.count_states(self.molecule)}",
            )
            cmd.label(
                f"{r.get_selection_str()} and name CA and state {cmd.count_states(self.molecule)}",
                f'"{mutation_str}"',
            )

    def attach_affinity_label(self, affinity, state):
        if self.molecule is None:
            print("Please select a molecule.")
            return

        label_name = "big_label"

        if not cmd.get_object_list(label_name):
            cmd.create(label_name, "none")

        cmd.remove(f"{label_name} and state {state}")

        cmd.pseudoatom(
            label_name,
            pos=(cmd.get_coords(self.molecule)[0] - [0, 30, 0]).tolist(),
            label=f"Binding affinity: {affinity} kcal/mol",
            state=state,
        )
        cmd.set("label_size", 30, label_name)
        cmd.set("label_color", "yellow", label_name)
        cmd.set("float_labels", 1, label_name)

    def record_history(self, affinity):
        self.history.append(HistoryEntry(affinity, list(self.mutations.values())))
        print(f"Recorded history entry {len(self.history)}.")

    def undo(self):
        if self.molecule is None:
            print("Please select a molecule.")
            return

        if len(self.history) > 1:
            self.history.pop()
            self.binding_affinity = self.history[-1].binding_affinity
            self.mutations = {
                mut.to_string(): mut for mut in self.history[-1].mutations
            }
            self.selected_mutation = None
            cmd.delete_states(self.molecule, f"{cmd.count_states(self.molecule)}")
            cmd.delete_states("big_label", f"{cmd.count_states('big_label')}")
            self.populate_mutation_choices(list(self.mutations.values()))
            self.highlight_mutations()
            print("Undid last mutation")

    def run(self):
        """Run the wizard to generate suggestions for the selected molecule."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        self.status = WizardState.RUNNING_INFERENCE
        cmd.refresh_wizard()

        # Run inference on a separate thread
        worker_thread = threading.Thread(target=self.run_inference)
        worker_thread.start()

    def run_inference(self):
        """Run the inference to generate mutation suggestions."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        if self.antibody_chain is None:
            print("Please select a chain.")
            return

        fasta_str = cmd.get_fastastr(f"{self.molecule} and chain {self.antibody_chain}")
        sequence = fasta_str.split("\n")[1]

        print(f"Running inference with models {self.models}")
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
                + self.models,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"Something went wrong while calling Efficient Evolution: {e}")
            return

        mutations = []
        for line in res.stdout.strip().splitlines():
            if line:
                mutations.append(
                    Mutation.from_EE_output(line, self.molecule, self.antibody_chain)
                )
        self.populate_mutation_choices(mutations)
        self.highlight_mutations()

        if self.history == []:
            self.record_history(0.0)

        print("Select a mutation.")

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
        cmd.select("to_mutate", mutation.start_residue.get_selection_str())

        self.status = WizardState.MUTATION_SELECTED
        cmd.refresh_wizard()

    def update_binding_affinity(self):
        """Update the binding affinity for the current state of the molecule."""

        affinity = self.evaluate_binding_affinity()
        print(f"New affinity for state {cmd.get_state()}: {affinity} kcal/mol.")
        self.attach_affinity_label(affinity, cmd.get_state())

    def evaluate_binding_affinity(self):
        """Evaluate the binding affinity using Prodigy."""

        if self.antibody_chain is None or self.antigen_chains == []:
            print("Please select an antibody and antigen chain.")
            return

        print(f"Evaluating affinity for {self.molecule}, state {cmd.get_state()}.")
        cmd.save(f"{self.molecule}.pdb", f"{self.molecule}", state=cmd.get_state())

        res = subprocess.run(
            [
                "prodigy",
                f"{self.molecule}.pdb",
                "--selection",
                self.antibody_chain,
                ",".join(self.antigen_chains),
                "--quiet",
            ],
            capture_output=True,
            text=True,
        )

        os.remove(f"{self.molecule}.pdb")
        return self.parse_prodigy_output(res.stdout)

    def parse_prodigy_output(self, output):
        """Parse the output of Prodigy to get the binding affinity."""

        if output is None or len(output.split()) != 2:
            print(f"Error: could not parse Prodigy output. Got: {output}")
            return

        return float(output.split()[1])

    def apply_mutation(self):
        """Apply the selected mutation to the molecule."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        if self.selected_mutation is None:
            print("Please select a mutation.")
            return

        cmd.wizard("mutagenesis")
        cmd.do("refresh_wizard")
        cmd.get_wizard().set_mode("%s" % one_to_three(self.selected_mutation.target))

        # The mutation is applied to the last state
        cmd.create("last_state", self.molecule, source_state=-1, target_state=-1)
        cmd.select(
            "tmp",
            f"last_state and resi {self.selected_mutation.start_residue.id} and chain {self.selected_mutation.start_residue.chain}",
        )
        cmd.get_wizard().do_select("tmp")
        cmd.get_wizard().apply()
        cmd.join_states(self.molecule, "last_state", mode=-1)
        cmd.delete("last_state")
        cmd.set_wizard()

        cmd.frame(cmd.count_states(self.molecule))

        print(f"Applied mutation {self.selected_mutation.to_string()}.")

        # Remove the applied mutation from the list of suggestions
        del self.mutations[self.selected_mutation.to_string()]
        self.selected_mutation = None
        self.populate_mutation_choices(list(self.mutations.values()))

        cmd.label(f"{self.molecule}", "''")
        self.highlight_mutations()

        # Record the new binding affinity
        affinity = self.evaluate_binding_affinity()
        self.record_history(affinity)
        self.attach_affinity_label(affinity, cmd.count_states(self.molecule))

        cmd.refresh_wizard()

        self.status = WizardState.MUTATION_APPLIED
        cmd.refresh_wizard()

    def get_panel(self):
        """Return the menu panel for the wizard."""

        # Title
        options = [
            [1, "Antibody Evolution", ""],
        ]

        # Choose molecule
        if self.status >= WizardState.READY:
            if self.molecule is None:
                molecule_label = "Choose molecule"
            else:
                molecule_label = self.molecule

            options.append(
                [3, molecule_label, "molecule"],
            )

        # Add entries to select chains and evaluate affinity
        if self.status >= WizardState.MOLECULE_SELECTED:
            antibody_chain_label = "Antibody Chain: "
            if self.antibody_chain:
                antibody_chain_label += self.antibody_chain
            else:
                antibody_chain_label += "None"
            options.append(
                [3, antibody_chain_label, "antibody_chain"],
            )

            antigen_chains_label = "Antigen Chains: "
            if self.antigen_chains:
                antigen_chains_label += ", ".join(self.antigen_chains)
            else:
                antigen_chains_label += "None"
            options.append(
                [3, antigen_chains_label, "antigen_chains"],
            )

            options.append(
                [
                    2,
                    "Evaluate Affinity",
                    "cmd.get_wizard().update_binding_affinity()",
                ],
            )

        # Add button to generate suggestions
        if self.status >= WizardState.ANTIBODY_CHAIN_SELECTED:
            options.append(
                [2, "Run Inference", "cmd.get_wizard().run()"],
            )

        # Add entries to select and apply mutations
        if self.status >= WizardState.SUGGESTIONS_GENERATED or self.mutations:
            if self.mutations is None:
                mutations_label = "-"
            elif self.selected_mutation is None:
                mutations_label = "Select Mutation"
            else:
                mutations_label = self.selected_mutation.to_string()

            options.append(
                [3, mutations_label, "mutations"],
            )
        if self.status >= WizardState.MUTATION_SELECTED:
            options.append(
                [2, "Apply Mutation", "cmd.get_wizard().apply_mutation()"],
            )

        if len(self.history) > 1:
            options.append(
                [2, "Undo", "cmd.get_wizard().undo()"],
            )

        options.append([2, "Dismiss", "cmd.set_wizard()"])

        return options
