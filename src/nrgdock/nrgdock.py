from PyQt5.QtCore import pyqtSignal, QThread
import os
from src.nrgdock.run_nrgdock import run_nrgdock


class WorkerThread(QThread):
    def __init__(self, form, ligand_set_folder_path, install_dir):
        super().__init__()
        self.form = form
        self.ligand_set_folder_path = ligand_set_folder_path
        self.install_dir = install_dir

    def run(self):
        temp_path = os.path.join(self.form.temp_line_edit.text(), 'NRGDock')
        run_nrgdock(self.form, temp_path, self.ligand_set_folder_path, self.install_dir)