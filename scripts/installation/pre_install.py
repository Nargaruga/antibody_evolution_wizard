import os
import sys
import subprocess


def download_model(model_name: str):
    subprocess.run(
        [
            "wget",
            f"https://dl.fbaipublicfiles.com/fair-esm/models/{model_name}.pt",
            "-P",
            os.path.join(
                os.path.expanduser("~"),
                ".cache",
                "torch",
                "hub",
                "checkpoints",
            ),
        ]
    )


def main():
    wizard_root = sys.argv[1]

    if os.name == "posix":
        try:
            subprocess.run(
                "conda list --name efficient-evolution",
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
                shell=True,
            ).check_returncode()
        except subprocess.CalledProcessError:
            subprocess.run(
                [
                    "conda",
                    "env",
                    "create",
                    "--name",
                    "efficient-evolution",
                    "-f",
                    "environment.yml",
                ],
                cwd=os.path.join(wizard_root, "ext", "efficient-evolution"),
                check=True,
            )

        subprocess.run(
            [
                "conda",
                "run",
                "-n",
                "efficient-evolution",
                "pip",
                "install",
                "-e",
                ".",
            ],
            cwd=os.path.join(wizard_root, "ext", "efficient-evolution"),
            check=True,
        )

        try:
            print(
                "You can choose to download the models now (may take a while) or have them downloaded automatically on first use. Download now? (Y/n)"
            )
            answer = input().strip().lower() or "y"
        except KeyboardInterrupt:
            print("Aborted by user.")
            exit(0)

        if answer == "y":
            models = [
                "esm1v_t33_650M_UR90S_1",
                "esm1v_t33_650M_UR90S_2",
                "esm1v_t33_650M_UR90S_3",
                "esm1v_t33_650M_UR90S_4",
                "esm1v_t33_650M_UR90S_5",
                "esm_msa1_t12_100M_UR50S",
                "esm1b_t33_650M_UR50S",
            ]

            for model in models:
                download_model(model)


if __name__ == "__main__":
    main()
