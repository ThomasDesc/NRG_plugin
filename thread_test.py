from pymol.Qt import QtCore
import time
from subprocess import run


class WorkerThread(QtCore.QThread):
    def __init__(self, command):
        super(WorkerThread, self).__init__()
        print('yo')
        self.command = command
    finished = QtCore.pyqtSignal()

    def run(self):
        run(self.command, shell=True)
        self.finished.emit()


def yo(command):
    worker = WorkerThread(command)
    time.sleep(1)
    worker.start()
    worker.finished.connect(worker.quit)
    worker.finished.connect(lambda: print('finished done done'))
