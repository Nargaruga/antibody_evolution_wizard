import subprocess


class AffinityComputationError(Exception):
    """Custom exception for affinity computation errors."""

    pass


def compute_affinity(
    molecule_file: str, antibody_chains: list[str], antigen_chains: list[str]
) -> float:
    """Compute the binding affinity of the given antibody-antigen complex.

    :param molecule_file: path to the PDB file containing the antibody-antigen complex
    :param antibody_chains: list of antibody chain identifiers (e.g., ["H", "L"])
    :param antigen_chains: list of antigen chain identifiers (e.g., ["A", "B"])

    :return the binding affinity, rounded to 2 decimal places
    :raises AffinityComputationError: if the Prodigy call fails or output is invalid
    """
    command = [
        "prodigy",
        molecule_file,
        "--selection",
        ",".join(antibody_chains),
        ",".join(antigen_chains),
        "--quiet",
    ]
    try:
        res = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
        )

        if not res.stdout or len(res.stdout.split()) != 2:
            raise AffinityComputationError(
                f"command {' '.join(command)} gave unexpected output: {res.stdout.strip()}"
            )

        return round(float(res.stdout.split()[1]), 2)

    except subprocess.CalledProcessError as e:
        raise AffinityComputationError(
            f"command {' '.join(command)} failed: {e.stderr}"
        )
