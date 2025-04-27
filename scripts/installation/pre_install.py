import os
import sys
import subprocess
from pathlib import Path


def main():
    wizard_root = sys.argv[1]

    plugin_name = "evolution_settings_plugin"

    if os.name == "nt":
        Path(os.path.join(wizard_root, "checkpoints")).mkdir(exist_ok=True)

        print("Building Efficient-Evolution image...")
        subprocess.run(
            [
                "docker",
                "build",
                "-f",
                os.path.join("docker", "efficient_evolution.Dockerfile"),
                "--tag",
                "efficient-evolution:latest",
                os.path.join("ext", "efficient-evolution"),
            ],
            cwd=wizard_root,
            check=True,
        )
    else:
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

    if os.name == "nt":
        subprocess.run(
            [
                "powershell.exe",
                "Compress-Archive",
                "-Path",
                plugin_name,
                "-DestinationPath",
                "plugin.zip",
                "-Force",
            ],
            cwd=wizard_root,
            check=True,
        )
    else:
        subprocess.run(
            [
                "zip",
                "-r",
                "plugin.zip",
                plugin_name,
            ],
            cwd=wizard_root,
            check=True,
        )


if __name__ == "__main__":
    main()
