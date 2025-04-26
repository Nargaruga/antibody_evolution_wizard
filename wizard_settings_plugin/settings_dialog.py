import os
import pkgutil
from pathlib import Path
from importlib.resources import is_resource, path
from functools import partial
import yaml

from pymol.Qt import QtWidgets, QtCore, QtGui

import antibody_evolution
import antibody_evolution.model_manager
from .settings_gui import Ui_Form

# global reference to avoid garbage collection of our dialog
# dialog = None

MODELS_REL_PATH = os.path.join("config", "models.yaml")
ICONS_PATH = os.path.join(Path(__file__).resolve().parent, "icons")

live_threads = []


class PluginError(Exception):
    """Custom exception for configuration errors."""

    pass


class DownloadThread(QtCore.QThread):
    finished = QtCore.pyqtSignal(bool)

    def __init__(self, model_name):
        super().__init__()
        self.model_name = model_name

    def run(self):
        result = antibody_evolution.model_manager.download_model(self.model_name)
        self.finished.emit(result.returncode == 0)


def run_plugin_gui():
    dialog = QtWidgets.QDialog()
    form = Ui_Form()
    form.setupUi(dialog)
    load_model_list(form)

    form.apply_button.clicked.connect(lambda: save_model_list(form))
    form.cancel_button.clicked.connect(lambda: dialog.close())
    form.download_button.clicked.connect(lambda: start_download(form))
    form.stop_button.clicked.connect(lambda: stop_downloads(form))

    dialog.show()


def update_icon(item):
    icon = QtGui.QIcon(os.path.join(ICONS_PATH, "download.png"))
    if antibody_evolution.model_manager.is_downloaded(item.text()):
        icon = QtGui.QIcon(os.path.join(ICONS_PATH, "check.png"))

    item.setIcon(icon)


def download_finished(success, form, list_pos):
    """Handle the download finished signal."""

    if success:
        print("Download completed successfully.")
    else:
        print("Download failed.")
        if os.name == "nt":
            print(
                'WINDOWS NOTE: if you are seeing "Invoke-WebRequest: Object reference not set to an instance of an object.", it probably means you attempted to resume download of a fully-downloaded model. In that case, you can ignore this error.'
            )

    update_icon(form.model_list.item(list_pos))


def start_download(form):
    """Start the download of selected models."""

    selected_models = get_selected_models(form)

    if not selected_models:
        print("No models selected for download.")
        return

    for list_pos, model_name in selected_models:
        form.model_list.item(list_pos).setIcon(
            QtGui.QIcon(os.path.join(ICONS_PATH, "wait.png"))
        )

        thread = DownloadThread(model_name)
        thread.finished.connect(
            partial(download_finished, form=form, list_pos=list_pos)
        )
        thread.start()
        live_threads.append(thread)


def stop_downloads(form):
    """Stop all running download threads."""
    for thread in live_threads:
        thread.terminate()
        thread.wait()

    for i in range(form.model_list.count()):
        item = form.model_list.item(i)
        if item.checkState() == QtCore.Qt.Checked:
            update_icon(item)

    print("Stopped all downloads.")


def get_selected_models(form):
    """Get the list of selected models from the form."""
    selected_models = []
    for i in range(form.model_list.count()):
        item = form.model_list.item(i)
        if item.checkState() == QtCore.Qt.Checked:
            selected_models.append((i, item.text()))

    return selected_models


def load_model_list(form: Ui_Form):
    raw_yaml = pkgutil.get_data("antibody_evolution", MODELS_REL_PATH)
    if raw_yaml is None:
        raise PluginError("Could not load models.")

    yaml_str = raw_yaml.decode("utf-8")
    models = yaml.safe_load(yaml_str)["models"]

    for model in models:
        checkable_item = QtWidgets.QListWidgetItem(model["name"])
        checkable_item.setFlags(checkable_item.flags() | QtCore.Qt.ItemIsUserCheckable)
        checkable_item.setCheckState(
            QtCore.Qt.Checked if model["use"] else QtCore.Qt.Unchecked
        )

        update_icon(checkable_item)
        form.model_list.addItem(checkable_item)

    file_path = os.path.join(Path(antibody_evolution.__file__).parent, MODELS_REL_PATH)
    print(f"Loaded model list from {file_path}.")


def save_model_list(form):
    models = []
    for i in range(form.model_list.count()):
        item = form.model_list.item(i)
        model = {
            "name": item.text(),
            "use": item.checkState() == QtCore.Qt.Checked,
        }
        models.append(model)

    yaml_str = yaml.dump({"models": models}, default_flow_style=False)
    file_path = os.path.join(Path(antibody_evolution.__file__).parent, MODELS_REL_PATH)
    print(file_path)
    with open(file_path, "w") as f:
        f.write(yaml_str)

    print(f"Saved model list to {file_path}.")
