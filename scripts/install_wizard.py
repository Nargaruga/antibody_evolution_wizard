import os
from pathlib import Path
import subprocess
import shutil
import re
import stat
import json
from util.constants import PYMOL_PYTHON_VERSION, DEFAULT_ENV_NAME, WIZARD_NAME

wizard_root = Path(__file__).parent.parent


def get_env_file(prefix):
    envs_dir = os.path.join(prefix, "envs")
    return (
        os.path.join(envs_dir, "windows_env.yml")
        if os.name == "nt"
        else os.path.join(envs_dir, "linux_env.yml")
    )


def install_openvr(clone_dir, conda_base_path, env_name):
    """Clone, build and install OpenVR."""

    print("Installing OpenVR...")

    env_dir = os.path.join(conda_base_path, "envs", env_name)

    if not os.path.exists(os.path.join(clone_dir, "openvr")):
        subprocess.run(
            [
                "git",
                "clone",
                "-b",
                "v1.0.17",
                "git@github.com:ValveSoftware/openvr.git",
                os.path.join(clone_dir, "openvr"),
            ],
            check=True,
        )

    if os.name == "posix":
        subprocess.run(
            [
                "cmake",
                "-S",
                ".",
                "-B",
                "build",
                "-DCMAKE_BUILD_TYPE=Release",
            ],
            cwd=os.path.join(clone_dir, "openvr"),
            check=True,
        )

        subprocess.run(
            ["cmake", "--build", "build", "--config", "Release"],
            cwd=os.path.join(clone_dir, "openvr"),
            check=True,
        )

        subprocess.run(
            ["sudo", "make", "install"],
            cwd=os.path.join(clone_dir, "openvr", "build"),
            check=True,
        )
    elif os.name == "nt":
        subprocess.run(
            [
                "cmake",
                "-S",
                ".",
                "-B",
                "build",
                f"-DCMAKE_INSTALL_PREFIX={env_dir}",
                "-DBUILD_SHARED=1",
            ],
            cwd=os.path.join(clone_dir, "openvr"),
            check=True,
        )

        subprocess.run(
            [
                "cmake",
                "--build",
                "build",
                "--config",
                "Release",
                "--target",
                "install",
            ],
            cwd=os.path.join(clone_dir, "openvr"),
            check=True,
        )

        # Copy the openvr.h header
        shutil.copy(
            os.path.join(clone_dir, "openvr", "headers", "openvr.h"),
            os.path.join(env_dir, "include"),
        )

        # Rename the .lib file
        shutil.move(
            os.path.join(
                env_dir,
                "Lib",
                "openvr_api64.lib",
            ),
            os.path.join(
                env_dir,
                "Lib",
                "openvr_api.lib",
            ),
        )

        # Move the .dll to the right directory
        shutil.move(
            os.path.join(
                env_dir,
                "Lib",
                "openvr_api64.dll",
            ),
            os.path.join(
                env_dir,
                "Library",
                "bin",
                "openvr_api64.dll",
            ),
        )


def install_pymol(clone_dir, conda_base_path, env_name):
    install_openvr(clone_dir, conda_base_path, env_name)

    if not os.path.exists(os.path.join(clone_dir, "pymol-open-source")):
        print("Installing PyMOL...")
        subprocess.run(
            [
                "git",
                "clone",
                "-b",
                "vr_atom_selection",
                "git@github.com:Nargaruga/pymol-open-source.git",
                os.path.join(clone_dir, "pymol-open-source"),
            ],
            check=True,
        )

    subprocess.run(
        f'conda run -n {env_name} pip install --config-settings openvr=True {os.path.join(clone_dir, "pymol-open-source")}',
        shell=True,
        check=True,
    )


def add_line_after(file, to_insert, pattern_to_insert, target_pattern):
    """Adds a line to the file after the specified point. If the line is already present, the file is unchanged."""

    with open(file, "r") as f:
        contents = f.read()

    if re.search(pattern_to_insert, contents):
        print(f"Entry already exists in {file}, skipping...")
        return

    target = target_pattern.search(contents)
    if target is None:
        print(f"Could not find target in {file}")
        return

    target_end = target.end()

    contents = contents[:target_end] + f"{to_insert}" + contents[target_end:]

    with open(file, "w") as f:
        f.write(contents)


# Entry point
current_env = os.environ.get("CONDA_DEFAULT_ENV")
if current_env is None:
    print("Could not detect conda environment. Is conda installed?")
    exit(1)

print(
    f"You are currently about to install the plugin in the {current_env} environment. Do you wish to create a new conda environment instead? (Y/n)"
)
try:
    answer = input().strip().lower() or "y"
except KeyboardInterrupt:
    print("Aborted by user.")
    exit(0)

if answer == "y":
    print(
        f'Please enter the name of the new environment (leave empty for default, "{DEFAULT_ENV_NAME}"):'
    )
    try:
        new_env = input().strip() or DEFAULT_ENV_NAME
    except KeyboardInterrupt:
        print("Aborted by user.")
        exit(0)

    if new_env == "":
        new_env = DEFAULT_ENV_NAME

    try:
        subprocess.run(
            f"conda list --name {new_env}",
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
            shell=True,
        ).check_returncode()
        print(
            f"Environment {new_env} already exists. Do you wish to overwrite it, use it or abort? (o/u/A)"
        )
        answer = ""
        while answer not in ["o", "u", "a"]:
            try:
                answer = input().strip().lower() or "a"
            except KeyboardInterrupt:
                print("Aborted by user.")
                exit(0)
            if answer == "o":
                print(f"Overwriting existing environment {new_env}.")
                if new_env == current_env:
                    print(
                        "Cannot overwrite an active environment. Please deactivate it before retrying."
                    )
                    exit(1)
                try:
                    subprocess.run(
                        f"conda env remove -n {new_env} -y",
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.STDOUT,
                        shell=True,
                    )
                    subprocess.run(
                        f"conda env create -n {new_env} -f {get_env_file(wizard_root)}",
                        check=True,
                        shell=True,
                    )
                except subprocess.CalledProcessError as e:
                    print(
                        f"Something went wrong while overwriting the environment: {e}"
                    )
                    exit(1)
            elif answer == "u":
                print(f"Using existing environment {new_env}.")
                continue
            elif answer == "a":
                print("Aborted by user.")
                exit(0)
            else:
                print(
                    "Invalid input. Please enter 'o' (overwrite), 'u' (use) or 'a' (abort)."
                )
    except subprocess.CalledProcessError:
        print(f"Creating new environment {new_env}.")
        subprocess.run(
            f"conda env create -n {new_env} -f {get_env_file(wizard_root)}",
            check=True,
            shell=True,
        )

    current_env = new_env

try:
    conda_base_path = str(
        subprocess.check_output("conda info --base", shell=True), "utf-8"
    ).strip()
except subprocess.CalledProcessError:
    print("Failed to retrieve conda base path.")
    exit(1)

prefix = os.path.join(conda_base_path, "envs", current_env)
if prefix is None:
    print("Something went wrong while creating the new environment.")
    exit(1)

if os.name == "nt":
    pymol_dir = os.path.join(
        prefix,
        "Lib",
        "site-packages",
        "pymol",
    )
else:
    pymol_dir = os.path.join(
        prefix,
        "lib",
        f"python{PYMOL_PYTHON_VERSION}",
        "site-packages",
        "pymol",
    )

if not os.path.exists(pymol_dir):
    print(
        f"PyMOL is not installed in the {current_env} environment. Do you wish to install it? (Y/n)"
    )
    try:
        answer = input().strip().lower() or "y"
    except KeyboardInterrupt:
        print("Aborted by user.")
        exit(0)
    if answer == "y":
        clone_dir_path = os.path.join(".", "tmp")
        Path(clone_dir_path).mkdir(parents=True, exist_ok=True)
        try:
            install_pymol(clone_dir_path, conda_base_path, current_env)
        except subprocess.CalledProcessError as e:
            print(f"Failed to install PyMOL: {e}")
            exit(1)

print("Copying files...")
installed_wizard_dir = os.path.join(pymol_dir, "wizard")
try:
    shutil.copy(
        os.path.join(wizard_root, f"{WIZARD_NAME}.py"),
        os.path.join(installed_wizard_dir, f"{WIZARD_NAME}.py"),
    )
    shutil.copytree(
        os.path.join(wizard_root, "config"),
        os.path.join(installed_wizard_dir, "config"),
        dirs_exist_ok=True,
    )
except shutil.Error as e:
    print(f"Failed to copy files: {e}")
    exit(1)

print("Adding menu entries...")
# Edit the openvr wizard to add a menu item in the internal menu
openvr_wizard = os.path.join(installed_wizard_dir, "openvr.py")
openvr_entry = f'\n[1, "Antibody Evolution", "wizard {WIZARD_NAME}"],'
openvr_entry_pattern = re.compile(
    openvr_entry.replace("[", r"\[").replace("]", r"\]").replace('"', r'["\']')
)
openvr_target_pattern = re.compile(r'\[2, ["\']Wizard Menu["\'], ["\']["\']\],')
add_line_after(openvr_wizard, openvr_entry, openvr_entry_pattern, openvr_target_pattern)

# Edit the openvr wizard to add a menu item in the internal menu
gui_file = os.path.join(pymol_dir, "_gui.py")
external_entry = f'\n("command", "Antibody Evolution", "wizard {WIZARD_NAME}"),'
external_entry_pattern = re.compile(
    external_entry.replace("(", r"\(").replace(")", r"\)").replace('"', r'["\']')
)
external_target_pattern = re.compile(r'\(\s*["\']menu["\'],\s*["\']Wizard["\'],\s*\[')
add_line_after(
    gui_file, external_entry, external_entry_pattern, external_target_pattern
)

print("Done!")

# Record the environment used for the installation
installation_data_file = os.path.join(wizard_root, "installation_data.json")
os.makedirs(os.path.dirname(installation_data_file), exist_ok=True)
with open(installation_data_file, "w") as f:
    json.dump({"conda_env": current_env}, f)


def remove_readonly(func, path, _):
    """Clear the readonly bit and remove the file."""

    os.chmod(path, stat.S_IWRITE)
    func(path)


if os.path.exists(os.path.join(wizard_root, "tmp")):
    print("Do you wish to clear the installation files? (Y/n)")
    try:
        answer = input().strip().lower() or "y"
    except KeyboardInterrupt:
        print("Aborted by user.")
        exit(0)
    if answer == "y":
        try:
            shutil.rmtree(os.path.join(wizard_root, "tmp"), onerror=remove_readonly)
            print("Files removed.")
        except FileNotFoundError:
            print("No files to remove.")
            pass

print(f"Remember to activate the {current_env} conda environment before running PyMOL.")
