import os
import subprocess


def compute_affinity(
    molecule_file: str, antibody_chain: str, antigen_chains: list[str]
) -> float:
    """Compute the binding affinity of the given antibody-antigen complex."""

    command = [
        "prodigy",
        molecule_file,
        "--selection",
        antibody_chain,
        ",".join(antigen_chains),
        "--quiet",
    ]
    print(f"Running command {' '.join(command)}")

    res = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    if res.stderr:
        print(f"Error: {res.stderr}")
        exit(-1)

    if res.stdout is None or len(res.stdout.split()) != 2:
        print(f"Error: could not parse Prodigy output. Got: {res.stdout}")
        exit(-1)

        try:
            os.remove(molecule_file)
        except FileNotFoundError:
            pass

    return round(float(res.stdout.split()[1]), 2)
