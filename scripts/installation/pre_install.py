import os
import sys
import subprocess


def main():
    wizard_root = sys.argv[1]
    env_name = sys.argv[2]

    subprocess.run(
        ["conda", "run", "-n", env_name, "pip", "install", "-e", "."],
        cwd=os.path.join(wizard_root, "ext", "efficient-evolution"),
    )


if __name__ == "__main__":
    main()
