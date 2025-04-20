from pymol.plugins import addmenuitemqt

from .settings_dialog import run_plugin_gui

def __init_plugin__(app=None):
    addmenuitemqt('Antibody Evolution Wizard Settings', run_plugin_gui)
