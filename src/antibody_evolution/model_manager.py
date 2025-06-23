import os
import subprocess

MODEL_NAME_MAPPING = {
    "esm1b": "esm1b_t33_650M_UR50S",
    "esm1v1": "esm1v_t33_650M_UR90S_1",
    "esm1v2": "esm1v_t33_650M_UR90S_2",
    "esm1v3": "esm1v_t33_650M_UR90S_3",
    "esm1v4": "esm1v_t33_650M_UR90S_4",
    "esm1v5": "esm1v_t33_650M_UR90S_5",
    "esm-msa": "esm_msa1_t12_100M_UR50S",
}


def map_model_name(model_name: str) -> str:
    """Map the model name to the correct format."""
    if model_name not in MODEL_NAME_MAPPING:
        raise ValueError(f"Invalid model name: {model_name}")
    return MODEL_NAME_MAPPING[model_name]


def get_powershell_prefix() -> str:
    if os.path.exists("C:\\Program Files\\PowerShell\\7\\pwsh.exe"):
        return ["pwsh", "-Command"]
    else:
        return ["powershell"]


def check_powershell_version(prefix: list[str], min_version: str) -> bool:
    """Check if the PowerShell version is greater than or equal to the specified version."""
    major = int(
        subprocess.run(
            prefix + ["$PSVersionTable.PSVersion.Major"],
            capture_output=True,
            text=True,
        ).stdout.strip()
    )

    minor = int(
        subprocess.run(
            prefix + ["$PSVersionTable.PSVersion.Minor"],
            capture_output=True,
            text=True,
        ).stdout.strip()
    )

    return float(f"{major}.{minor}") >= float(min_version)


def download_model(model: str):
    continue_flag = ""
    if os.name == "nt":
        powershell_min_version = "6.1"
        if check_powershell_version(get_powershell_prefix(), powershell_min_version):
            continue_flag = "-Resume"
        else:
            print(
                f"Note that resuming downloads is only supported in PowerShell version {powershell_min_version} or higher."
            )
            print("Starting download from scratch.")
    elif os.name == "posix":
        continue_flag = "-c"

    print(f"Downloading {model}...")

    checkpoints_dir = os.path.join(
        os.path.expanduser("~"),
        ".cache",
        "torch",
        "hub",
        "checkpoints",
    )

    os.makedirs(checkpoints_dir, exist_ok=True)

    output = subprocess.run(
        (get_powershell_prefix() if os.name == "nt" else [])
        + (["Invoke-WebRequest"] if os.name == "nt" else ["wget"])
        + [continue_flag]
        + [
            f"https://dl.fbaipublicfiles.com/fair-esm/models/{map_model_name(model)}.pt",
            ("-OutFile" if os.name == "nt" else "-P"),
            os.path.join(
                checkpoints_dir,
                f"{map_model_name(model)}.pt",
            ),
        ],
    )

    return output


def download_models(models):
    for model in models:
        download_model(model)


def is_downloaded(model):
    return os.path.exists(
        os.path.join(
            os.path.expanduser("~"),
            ".cache",
            "torch",
            "hub",
            "checkpoints",
            f"{map_model_name(model)}.pt",
        )
    )
