from pymol.wizard import Wizard
from pymol import cmd
from general_functions import get_mouse_config


class Sphere(Wizard):
    def __init__(self):
        #print "New instance of flexsphere Wizard"
        Wizard.__init__(self)
        self.mouse_mode = get_mouse_config()

    def done(self):
        cmd.config_mouse(self.mouse_mode)
        cmd.set_wizard()

    def get_prompt(self):
        # TODO: change prompt based on os
        self.prompt = None
        self.prompt = ["Mouse: Press Shift + Mouse3 (Wheel Click) to move the sphere.",
                       "Trackpad: Press Command + Mouse click to move the sphere.",
                       "Use the scaler in the interface to edit its radius."]
        return self.prompt

    def get_panel(self):
        return [
            [1, 'Move sphere', ''],
            [2, 'Done', 'cmd.set_wizard().done()'],
        ]