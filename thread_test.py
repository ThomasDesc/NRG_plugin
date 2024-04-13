from pymol.Qt import QtCore
import time
from subprocess import run


class WorkerThread(QtCore.QThread):
    def __init__(self, command):
        super(WorkerThread, self).__init__()
        self.command = command
    finished = QtCore.pyqtSignal()

    def run(self):
        run(self.command, shell=True)
        self.finished.emit()
