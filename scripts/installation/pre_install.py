import os
import sys
import subprocess


def main():
    wizard_root = sys.argv[1]

    if os.name == "nt":
        prefix = ["powershell.exe"]
    else:
        prefix = []

    subprocess.run(
        prefix
        + [
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
        prefix
        + [
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


if __name__ == "__main__":
    main()
