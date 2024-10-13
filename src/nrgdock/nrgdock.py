from PyQt5.QtWidgets import QWidget, QPushButton, QVBoxLayout, QApplication
from PyQt5.QtCore import pyqtSignal, QThread
import os
from src.nrgdock.run_nrgdock import run_nrgdock


class WorkerThread(QThread):
    def __init__(self, form, temp_path, ligand_set_folder_path, install_dir):
        super().__init__()
        self.form = form
        self.temp_path = temp_path
        self.ligand_set_folder_path = ligand_set_folder_path
        self.install_dir = install_dir

    def run(self):
        run_nrgdock(self.form, self.temp_path, self.ligand_set_folder_path, self.install_dir)


class NRGDockRunner(QWidget):
    def __init__(self):
        super().__init__()

    def run_task(self, form, ligand_set_folder_path, install_dir):
        print('Starting the task...')
        temp_path = os.path.join(form.temp_line_edit.text(), 'NRGDock')
        self.worker = WorkerThread(form, temp_path, ligand_set_folder_path, install_dir)
        self.worker.start()