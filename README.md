# Antibody Evolution Wizard
PyMOL wizard for performing ML-based antibody evolution through [brianhie/efficient-evolution](https://github.com/brianhie/efficient-evolution). Utilizes [haddocking/prodigy](https://github.com/haddocking/prodigy) for measuring the binding affinity of the loaded antibody-antigen complex and assessing the impact of mutations.

## Installation
The wizard can be installed with the [Wizard Installer](https://github.com/Nargaruga/pymol_wizard_installer) tool.

## Usage
The wizard can be accessed in one of two ways:
- by writing `wizard evolution` in the PyMOL console;
- by going to `Wizard->Antibody Evolution` in the external GUI (or the internal one, if in VR mode);

### Obtaining Mutations
After selecting the molecule from the drop-down menu, specify the chain you want to obtain mutations for and click `Generate Suggestions` to acquire them. You can then select the desired mutation from the list, where each mutation is accompanied by the consensus score, and apply it. The `Select Best Mutation` button automatically selects the mutation which brings the greatest reduction in binding free energy.

### Evaluating Binding Affinity
Evaluating the binding affinity requires you to select the participating chains in the appropriate drop-down menus. Doing so will reveal the `Evaluate Affinity` button: click it to obtain the binding free energy of the complex in kcal/mol.

### The History Feature
Each applied mutation adds a state to the selected object, creating a timeline of applied mutations. Each individual frame can be annotated with its binding affinity to track its changes through each round of evolution. Through the `Undo` button it is possible to revert the last mutation performed.

