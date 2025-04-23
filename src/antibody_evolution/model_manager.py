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


def download_model(model: str):
    return subprocess.run(
        (["powershell.exe"] if os.name == "nt" else [])
        + [
            "wget",
            f"https://dl.fbaipublicfiles.com/fair-esm/models/{map_model_name(model)}.pt",
            ("-OutFile" if os.name == "nt" else "-P"),
            os.path.join(
                os.path.expanduser("~"),
                ".cache",
                "torch",
                "hub",
                "checkpoints",
                f"{map_model_name(model)}.pt",
            ),
        ]
    )


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
