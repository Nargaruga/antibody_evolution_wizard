from __future__ import annotations

import os
from enum import IntEnum, auto
import threading
from pathlib import Path
from dataclasses import dataclass
import tempfile
import pkgutil


from pymol.wizard import Wizard
from pymol.wizarding import WizardError
from pymol import cmd
import pymol2

import yaml

from antibody_evolution.mutation_evaluation import (
    compute_affinity,
    AffinityComputationError,
)
from antibody_evolution.mutation_suggestions import (
    Suggestion,
    get_mutation_suggestions,
    is_residue_valid,
)
from antibody_evolution.residue import Residue, one_to_three


@dataclass
class HistoryEntry:
    """An entry in the history of the wizard."""

    binding_affinity: float
    suggestions: list[Suggestion]


class WizardInputState(IntEnum):
    INITIALIZING = auto()
    READY = auto()
    MOLECULE_SELECTED = auto()
    CHAIN_TO_MUTATE_SELECTED = auto()
    MUTATIONS_READY = auto()
    MUTATION_SELECTED = auto()


class WizardTaskState(IntEnum):
    IDLE = auto()
    GENERATING_SUGGESTIONS = auto()
    FINDING_BEST_MUTATION = auto()


class Evolution(Wizard):
    """A wizard for suggesting mutations to an antibody."""

    def __init__(self, _self=cmd):
        Wizard.__init__(self, _self)
        self.input_state = WizardInputState.INITIALIZING
        self.task_state = WizardTaskState.IDLE
        cmd.refresh_wizard()
        self.models = ["esm1b"]
        self.molecule = None
        self.chain_to_mutate = None
        self.antibody_chains = []
        self.antigen_chains = []
        self.suggestions: dict[str, Suggestion] = {}
        self.selected_suggestion = None
        self.history: list[HistoryEntry] = []
        self.populate_molecule_choices()

        # Load the model in a separate thread
        initializer_thread = threading.Thread(target=self.load_model)
        initializer_thread.start()

    def get_prompt(self):  # type: ignore
        """Return the prompt for the current state of the wizard."""

        prompt = []

        if self.input_state == WizardInputState.READY:
            prompt.append("Select a molecule.")
        elif self.input_state == WizardInputState.MOLECULE_SELECTED:
            prompt.append("Select a chain.")
        elif self.input_state == WizardInputState.CHAIN_TO_MUTATE_SELECTED:
            prompt.append(
                f"Run to generate mutation suggestions for {self.molecule}, chain {self.chain_to_mutate}."
            )
        elif self.input_state == WizardInputState.MUTATIONS_READY:
            prompt.append("Select a mutation to apply.")

        elif (
            self.input_state == WizardInputState.MUTATION_SELECTED
            and self.selected_suggestion is not None
        ):
            prompt.append(
                "Apply the mutation %s?" % self.selected_suggestion.mutation.to_string()
            )

        if self.task_state == WizardTaskState.GENERATING_SUGGESTIONS:
            prompt.append("Generating suggestions, please wait...")
        elif self.task_state == WizardTaskState.FINDING_BEST_MUTATION:
            prompt.append("Finding the best mutation, please wait...")

        return prompt

    def get_panel(self):  # type: ignore
        """Return the menu panel for the wizard."""

        # Title
        options = [
            [1, "Antibody Evolution", ""],
        ]

        # Choose molecule
        if self.input_state >= WizardInputState.READY:
            if self.molecule is None:
                molecule_label = "Choose molecule"
            else:
                molecule_label = self.molecule

            options.append(
                [3, molecule_label, "molecule"],
            )

        # Add entries to select chains and evaluate affinity
        if self.input_state >= WizardInputState.MOLECULE_SELECTED:
            chain_to_mutate_label = "Chain to Mutate: "
            if self.chain_to_mutate:
                chain_to_mutate_label += self.chain_to_mutate
            else:
                chain_to_mutate_label += "None"
            options.append(
                [3, chain_to_mutate_label, "chain_to_mutate"],
            )

            antibody_chains_label = "Antibody Chains: "
            if self.antibody_chains:
                antibody_chains_label += ", ".join(self.antibody_chains)
            else:
                antibody_chains_label += "None"
            options.append(
                [3, antibody_chains_label, "antibody_chains"],
            )

            antigen_chains_label = "Antigen Chains: "
            if self.antigen_chains:
                antigen_chains_label += ", ".join(self.antigen_chains)
            else:
                antigen_chains_label += "None"
            options.append(
                [3, antigen_chains_label, "antigen_chains"],
            )

            if self.antibody_chains and self.antigen_chains:
                options.append(
                    [
                        2,
                        "Evaluate Affinity",
                        "cmd.get_wizard().update_binding_affinity()",
                    ]
                )

        # Add button to generate suggestions
        if self.input_state >= WizardInputState.CHAIN_TO_MUTATE_SELECTED:
            options.append(
                [2, "Generate Suggestions", "cmd.get_wizard().run()"],
            )

        # Add entries to select and apply mutations
        if self.input_state >= WizardInputState.MUTATIONS_READY:
            if self.suggestions is None:
                mutations_label = "-"
            elif self.selected_suggestion is None:
                mutations_label = "Select Mutation"
            else:
                mutations_label = self.selected_suggestion.mutation.to_string()

            options.extend(
                [
                    [
                        2,
                        "Select Best Mutation",
                        "cmd.get_wizard().select_best_mutation()",
                    ],
                    [3, mutations_label, "mutations"],
                ],
            )
        if self.input_state >= WizardInputState.MUTATION_SELECTED:
            options.append(
                [2, "Apply Mutation", "cmd.get_wizard().apply_selected_mutation()"],
            )

        if len(self.history) > 1:
            options.append(
                [2, "Undo", "cmd.get_wizard().undo()"],
            )

        options.append([2, "Dismiss", "cmd.set_wizard()"])

        return options

    def update_input_state(self):
        """Update the state of the wizard based on the current inputs."""

        if self.molecule:
            self.input_state = WizardInputState.MOLECULE_SELECTED
        else:
            self.input_state = WizardInputState.READY

        if self.chain_to_mutate:
            self.input_state = WizardInputState.CHAIN_TO_MUTATE_SELECTED

        if self.suggestions:
            self.input_state = WizardInputState.MUTATIONS_READY

        if self.selected_suggestion:
            self.input_state = WizardInputState.MUTATION_SELECTED

        cmd.refresh_wizard()

    def load_model(self):
        raw_yaml = pkgutil.get_data(
            "antibody_evolution", os.path.join("config", "models.yaml")
        )
        if raw_yaml is None:
            raise WizardError("Could not load models file.")

        yaml_str = raw_yaml.decode("utf-8")
        models = yaml.safe_load(yaml_str)
        self.models = [model["name"] for model in models["models"] if model["use"]]

        self.input_state = WizardInputState.READY
        cmd.refresh_wizard()

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
        self.menu["chain_to_mutate"] = [[2, "Chain to Mutate", ""]]
        self.menu["antibody_chains"] = [[2, "Antibody Chains", ""]]
        self.menu["antigen_chains"] = [[2, "Antigen Chains", ""]]
        for c in chains:
            self.menu["chain_to_mutate"].append(
                [
                    1,
                    c,
                    'cmd.get_wizard().set_chain_to_mutate("' + c + '")',
                ]
            )

            self.menu["antibody_chains"].append(
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
        self.populate_chain_choices()
        self.update_input_state()

    def set_chain_to_mutate(self, chain):
        """Set the chain to mutate."""

        self.chain_to_mutate = chain
        self.update_input_state()

    def set_antibody_chain(self, chain):
        """Set the antibody chain."""

        if chain in self.antibody_chains:
            self.antibody_chains.remove(chain)
        else:
            self.antibody_chains.append(chain)

        self.update_input_state()

    def set_antigen_chain(self, chain):
        """Set the antigen chain."""

        if chain in self.antigen_chains:
            self.antigen_chains.remove(chain)
        else:
            self.antigen_chains.append(chain)

        self.update_input_state()

    def find_mutation_for(self, sel: str):
        """Find the mutation suggestion for the selected residue, if any."""

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
        for suggestion_str, suggestion in self.suggestions.items():
            if suggestion.mutation.start_residue == residue:
                return suggestion_str

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
                self.set_suggestion(mutation)
        except WizardError as e:
            print(e)

        cmd.delete(name)

    def highlight_mutations(self):
        """Highlight the residues that have been suggested for mutation."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        for mutation_str, suggestion in self.suggestions.items():
            start_residue = suggestion.mutation.start_residue
            cmd.color(
                "cyan",
                f"{self.molecule} and resi {start_residue.id} and chain {start_residue.chain} and state {cmd.count_states(self.molecule)}",
            )
            cmd.label(
                f"{start_residue.get_selection_str()} and name CA and state {cmd.count_states(self.molecule)}",
                f'"{mutation_str}"',
            )

    def attach_affinity_label(self, affinity, state):
        """Attach a label with the binding affinity to the molecule."""

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
        self.history.append(HistoryEntry(affinity, list(self.suggestions.values())))
        print(f"Recorded history entry {len(self.history)}.")

    def undo(self):
        """Undo the last mutation."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        if len(self.history) > 1:
            self.history.pop()
            print(f"Removed history entry {len(self.history) + 1}.")
            self.binding_affinity = self.history[-1].binding_affinity
            self.suggestions = {
                s.mutation.to_string(): s for s in self.history[-1].suggestions
            }
            self.selected_suggestion = None

            cmd.delete_states(self.molecule, f"{cmd.count_states(self.molecule)}")
            if "big_label" in cmd.get_names():
                cmd.delete_states("big_label", f"{cmd.count_states('big_label')}")
            self.populate_mutation_choices(list(self.suggestions.values()))
            self.highlight_mutations()
            print("Undid last mutation")

            self.input_state = WizardState.SUGGESTIONS_GENERATED
            cmd.refresh_wizard()
        else:
            print("Nothing to undo.")

    def run(self):
        """Run the wizard to generate suggestions for the selected molecule."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        def aux():
            self.input_state = WizardTaskState.GENERATING_SUGGESTIONS
            cmd.refresh_wizard()

            try:
                self.suggest_mutations()
                self.input_state = WizardInputState.MUTATIONS_READY
            except Exception as e:
                print(f"Error generating suggestions: {e}")
                return
            finally:
                self.task_state = WizardTaskState.IDLE
                cmd.refresh_wizard()

        # Generate suggestions on a separate thread
        worker_thread = threading.Thread(target=aux)
        worker_thread.start()

    def suggest_mutations(self):
        """Generate mutation suggestions."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        if self.chain_to_mutate is None:
            print("Please select a chain.")
            return

        annotated_residues = []
        cmd.iterate(
            f"{self.molecule} and chain {self.chain_to_mutate} and name CA",
            "annotated_residues.append((oneletter, resi))",
            space=locals(),
        )

        sequence = []
        ids = []
        for oneletter, resi in annotated_residues:
            sequence.append(oneletter)
            ids.append(resi)

        annotated_mutations = get_mutation_suggestions(
            self.molecule, "".join(sequence), ids, self.models, self.chain_to_mutate
        )

        filtered_mutations = []
        for suggestion in annotated_mutations:
            if is_residue_valid(
                self.molecule,
                self.chain_to_mutate,
                suggestion.mutation.start_residue.id,
            ):
                filtered_mutations.append(suggestion)
            else:
                print(f"Filtered out mutation for invalid residue: {suggestion}")

        self.populate_mutation_choices(filtered_mutations)
        self.highlight_mutations()

        if self.history == []:
            self.record_history(0.0)

        print("Select a mutation.")

    def populate_mutation_choices(self, suggestions: list[Suggestion]):
        """Populate the menu with the generated mutation suggestions."""

        self.menu["mutations"] = [[2, "Mutations", ""]]
        suggestions.sort(key=lambda x: -x.occurrences)
        for suggestion in suggestions:
            self.suggestions[suggestion.mutation.to_string()] = suggestion
            self.menu["mutations"].append(
                [
                    1,
                    f"{suggestion.mutation.to_string()} ({suggestion.occurrences})",
                    'cmd.get_wizard().set_suggestion("'
                    + suggestion.mutation.to_string()
                    + '")',
                ]
            )

    def set_suggestion(self, mutation_str):
        """Set the selected mutation to apply."""
        suggestion = self.suggestions[mutation_str]
        self.selected_suggestion = suggestion
        cmd.select("to_mutate", suggestion.mutation.start_residue.get_selection_str())

        self.update_input_state()

    def update_binding_affinity(self):
        """Update the binding affinity for the current state of the molecule."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        if self.antibody_chains is None or self.antigen_chains == []:
            print("Please select an antibody and antigen chain.")
            return

        prodigy_outfile_handle, prodigy_outfile_path = tempfile.mkstemp(suffix=".pdb")

        cmd.save(prodigy_outfile_path, self.molecule)
        try:
            affinity = compute_affinity(
                prodigy_outfile_path, self.antibody_chains, self.antigen_chains
            )
        except AffinityComputationError as e:
            print(f"Error computing affinity: {e}")
            return
        finally:
            os.close(prodigy_outfile_handle)
            os.remove(prodigy_outfile_path)

        print(f"New affinity for state {cmd.get_state()}: {affinity} kcal/mol.")
        self.attach_affinity_label(affinity, cmd.get_state())

    def apply_selected_mutation(self):
        """Apply the selected mutation to the molecule."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        if self.antibody_chains is None or self.antigen_chains == []:
            print("Please select an antibody and antigen chain.")
            return

        if self.selected_suggestion is None:
            print("Please select a mutation.")
            return

        # The mutation is applied to the last state
        cmd.create(
            "last_state",
            self.molecule,
            cmd.count_states(self.molecule),
        )

        cmd.wizard("mutagenesis")
        cmd.do("refresh_wizard")
        cmd.get_wizard().do_select(
            f"/last_state//{self.chain_to_mutate}/{self.selected_suggestion.mutation.start_residue.id}"
        )
        cmd.get_wizard().set_mode(
            one_to_three(self.selected_suggestion.mutation.target_resn)
        )
        cmd.frame(1)
        cmd.get_wizard().apply()
        cmd.set_wizard()

        cmd.join_states(self.molecule, "last_state", mode=0)
        cmd.delete("last_state")
        cmd.frame(cmd.count_states(self.molecule))

        print(f"Applied mutation {self.selected_suggestion.mutation.to_string()}.")

        # Remove the applied mutation from the list of suggestions
        del self.suggestions[self.selected_suggestion.mutation.to_string()]
        self.selected_suggestion = None
        self.populate_mutation_choices(list(self.suggestions.values()))

        cmd.label(f"{self.molecule}", "''")
        self.highlight_mutations()

        # Record the new binding affinity
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=True) as tmp_file:
            cmd.save(tmp_file.name, self.molecule)
            affinity = compute_affinity(
                tmp_file.name, self.antibody_chains, self.antigen_chains
            )
        self.record_history(affinity)
        self.attach_affinity_label(affinity, cmd.count_states(self.molecule))

        self.update_input_state()

    def select_best_mutation(self):
        """Select the mutation with the highest impact on the binding affinity."""
        if self.molecule is None:
            print("Please select a molecule.")
            return

        if self.antibody_chains is None or self.antigen_chains == []:
            print("Please select an antibody and antigen chain.")
            return

        self.task_state = WizardTaskState.FINDING_BEST_MUTATION

        ddgs = {}
        with (
            tempfile.NamedTemporaryFile(suffix=".pdb", delete=True) as molecule_file,
        ):
            # TODO use thread to avoid blocking gui
            cmd.save(molecule_file.name, self.molecule)
            original_affinity = compute_affinity(
                molecule_file.name, self.antibody_chains, self.antigen_chains
            )

            ddgs = {}
            with pymol2.PyMOL() as bg_pymol:
                # TODO I'd put this in separate functions, but passing the PyMOL instance object causes problems
                for mutation_str, suggestion in self.suggestions.items():
                    bg_pymol.cmd.load(molecule_file.name, self.molecule)

                    selection_string = f"{self.molecule} and chain {self.chain_to_mutate} and resi {suggestion.mutation.start_residue.id}"
                    bg_pymol.cmd.select("tmp", selection_string)

                    bg_pymol.cmd.wizard("mutagenesis")
                    bg_pymol.cmd.do("refresh_wizard")
                    bg_pymol.cmd.get_wizard().do_select("tmp")
                    bg_pymol.cmd.get_wizard().set_mode(
                        one_to_three(suggestion.mutation.target_resn)
                    )
                    bg_pymol.cmd.frame(1)
                    bg_pymol.cmd.get_wizard().apply()
                    bg_pymol.cmd.set_wizard()

                    with tempfile.NamedTemporaryFile(
                        suffix=".pdb", delete=True
                    ) as mutated_molecule_file:
                        bg_pymol.cmd.save(mutated_molecule_file.name, self.molecule)
                        mutated_affinity = compute_affinity(
                            mutated_molecule_file.name,
                            self.antibody_chains,
                            self.antigen_chains,
                        )

                    ddg = round(mutated_affinity - original_affinity, 2)

                    print(f"Computed ddG for {mutation_str}: {ddg}")
                    ddgs[mutation_str] = ddg

                    bg_pymol.cmd.delete(self.molecule)
                    bg_pymol.cmd.delete("tmp")

        best_suggestion = min(ddgs, key=ddgs.get)
        print(f"Best mutation: {best_suggestion} with ddG of {ddgs[best_suggestion]}.")
        self.set_suggestion(best_suggestion)

        self.task_state = WizardTaskState.IDLE
        self.update_input_state()
