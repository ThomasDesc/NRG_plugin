from pymol.Qt import QtCore
import time
from pymol.Qt import QtGui
from pymol.Qt import QtWidgets
import general_functions
from subprocess import run
import os


class FileUpdaterThread(QtCore.QThread):
    def __init__(self, folder_path, table_widget, progress_bar, hex_colour_list):
        super(FileUpdaterThread, self).__init__()
        self.folder_path = folder_path
        self.table_widget = table_widget
        self.progress_bar = progress_bar
        self.hex_colour_list = hex_colour_list
        self.update_file_path = os.path.join(self.folder_path, ".update")
        self.running = True

    def run(self):
        while self.running:
            if os.path.exists(self.update_file_path):
                self.read_update(self.table_widget, self.progress_bar, self.hex_colour_list)
                os.remove(self.update_file_path)

    def colour_specific_cell(self, table_widget, data):
        num_column = table_widget.columnCount()
        row = int(data[1]) - 1
        for column_counter, column in enumerate(range(num_column)):
            if column_counter == 0:
                item = QtWidgets.QTableWidgetItem()
                item.setBackground(QtGui.QColor(data[column_counter]))
            else:
                item = QtWidgets.QTableWidgetItem()
                item = QtWidgets.QTableWidgetItem(str(data[column_counter]))
            table_widget.setItem(row, column, item)

    def read_update(self, table_widget, progress_bar, hex_colour_list, num_results=5):
        number_color_list = general_functions.create_number_list(num_results, len(hex_colour_list))
        generation = 0
        with open(self.update_file_path, "r") as f:
            for line_counter, line in enumerate(f):
                if line.startswith("Generation"):
                    generation = int(line[11:].strip())
                    progress_bar.setValue(generation)
                if line_counter > 1 and generation % 10 == 0:
                    line = line.split()
                    top_number = int(line[0]) + 1
                    cf = line[-5]
                    fitness = line[-1]
                    rmsd = 'N/A'
                    data = (hex_colour_list[number_color_list[top_number - 1]], top_number, cf, fitness, rmsd)
                    self.colour_specific_cell(table_widget, data)

    def stop(self):
        self.running = False


class WorkerThread(QtCore.QThread):
    def __init__(self, command, folder_path, table_widget):
        super(WorkerThread, self).__init__()
        self.command = command
        self.folder_path = folder_path
        self.file_updater_thread = FileUpdaterThread(folder_path, table_widget)

    finished = QtCore.pyqtSignal()

    def run(self):
        self.file_updater_thread.start()
        run(self.command, shell=True)
        self.file_updater_thread.stop()
        self.file_updater_thread.wait()
        self.finished.emit()
