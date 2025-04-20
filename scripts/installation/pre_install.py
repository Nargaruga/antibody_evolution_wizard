import os
import sys
import subprocess


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

        if os.name == "nt":
            subprocess.run(
                [
                    "Compress-Archive",
                    "-Path",
                    "wizard_settings_plugin",
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
                    "wizard_settings_plugin",
                ],
                cwd=wizard_root,
                check=True,
            )


if __name__ == "__main__":
    main()
