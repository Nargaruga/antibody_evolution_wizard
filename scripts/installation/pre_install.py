import os
import sys
import subprocess


def main():
    wizard_root = sys.argv[1]

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
    )


if __name__ == "__main__":
    main()
