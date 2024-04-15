from pymol.Qt import QtCore
import time
from pymol.Qt import QtGui
from pymol.Qt import QtWidgets
import general_functions
from subprocess import run
import os


class FileUpdaterThread(QtCore.QThread):
    generation_signal = QtCore.pyqtSignal(int)
    table_signal = QtCore.pyqtSignal(list)

    def __init__(self, folder_path, table_widget, hex_colour_list):
        super(FileUpdaterThread, self).__init__()
        self.folder_path = folder_path
        self.table_widget = table_widget
        self.hex_colour_list = hex_colour_list
        self.current_generation = 0
        self.update_file_path = os.path.join(self.folder_path, ".update")
        self.running = True

    def run(self):
        while self.running:
            if os.path.exists(self.update_file_path):
                self.read_update(self.table_widget, self.hex_colour_list)
                os.remove(self.update_file_path)
                self.current_generation += 1

    def colour_specific_cell(self, table_widget, data):
        self.table_signal.emit([table_widget, data])

    def read_update(self, table_widget, hex_colour_list, num_results=5):
        number_color_list = general_functions.create_number_list(num_results, len(hex_colour_list))
        with open(self.update_file_path, "r") as f:
            for line_counter, line in enumerate(f):
                if line_counter > 1 and self.current_generation % 10 == 0:
                    line = line.split()
                    top_number = int(line[0]) + 1
                    cf = line[-5]
                    fitness = line[-1]
                    rmsd = 'N/A'
                    data = (hex_colour_list[number_color_list[top_number - 1]], top_number, cf, fitness, rmsd)
                    self.generation_signal.emit(self.current_generation)
                    self.colour_specific_cell(table_widget, data)

    def stop(self):
        self.running = False


class WorkerThread(QtCore.QThread):
    generation_signal_received = QtCore.pyqtSignal(int)
    table_signal_received = QtCore.pyqtSignal(list)

    def __init__(self, command, folder_path, table_widget, hex_colour_list):
        super(WorkerThread, self).__init__()
        self.command = command
        self.folder_path = folder_path
        self.hex_colour_list = hex_colour_list
        self.file_updater_thread = FileUpdaterThread(folder_path, table_widget, hex_colour_list)
        self.file_updater_thread.generation_signal.connect(self.generation_signal_received)
        self.file_updater_thread.table_signal.connect(self.table_signal_received)

    finished = QtCore.pyqtSignal()

    def run(self):
        self.file_updater_thread.start()
        run(self.command, shell=True)
        self.file_updater_thread.stop()
        self.file_updater_thread.wait()
        self.finished.emit()
