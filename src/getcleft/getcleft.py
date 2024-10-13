from PyQt5.QtWidgets import QWidget, QPushButton, QVBoxLayout, QApplication
from PyQt5.QtCore import pyqtSignal, QThread
from src.getcleft.run_getcleft import run_getcleft


class WorkerThread(QThread):
    def __init__(self, form, binary_folder_path, binary_suffix, temp_path, install_dir):
        super().__init__()
        self.form = form
        self.binary_folder_path = binary_folder_path
        self.binary_suffix = binary_suffix
        self.temp_path = temp_path
        self.install_dir = install_dir

    def run(self):
        run_getcleft(self.form, self.binary_folder_path, self.binary_suffix, self.temp_path, self.install_dir)


class GetCleftRunner(QWidget):
    def __init__(self):
        super().__init__()

    def run_task(self, form, binary_folder_path, binary_suffix, install_dir):
        print('Starting the task...')
        temp_path = form.temp_line_edit.text()
        self.worker = WorkerThread(
            form,
            binary_folder_path,
            binary_suffix,
            temp_path,
            install_dir
        )
        self.worker.start()
