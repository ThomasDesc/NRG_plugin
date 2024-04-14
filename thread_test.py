from pymol.Qt import QtCore
import time
from subprocess import run
import os


class FileUpdaterThread(QtCore.QThread):
    def __init__(self, folder_path):
        super(FileUpdaterThread, self).__init__()
        self.folder_path = folder_path
        self.running = True

    def run(self):
        while self.running:
            update_file_path = os.path.join(self.folder_path, ".update")
            if os.path.exists(update_file_path):
                os.remove(update_file_path)

    def stop(self):
        self.running = False


class WorkerThread(QtCore.QThread):
    def __init__(self, command, folder_path):
        super(WorkerThread, self).__init__()
        self.command = command
        self.folder_path = folder_path
        self.file_updater_thread = FileUpdaterThread(folder_path)

    finished = QtCore.pyqtSignal()

    def run(self):
        self.file_updater_thread.start()
        run(self.command, shell=True)
        self.file_updater_thread.stop()
        self.file_updater_thread.wait()
        self.finished.emit()
