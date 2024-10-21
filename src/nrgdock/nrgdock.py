# from PyQt5.QtCore import QRunnable, QThreadPool, pyqtSignal, QObject

# class WorkerSignals(QObject):
#     finished = pyqtSignal()
#     progress = pyqtSignal(str)
#
# class WorkerRunnable(QObject):
#     def __init__(self, target_path, targets, overwrite, run_getcleft, deps_path):
#         super().__init__()
#         self.target_path = target_path
#         self.targets = targets
#         self.overwrite = overwrite
#         self.run_getcleft = run_getcleft
#         self.deps_path = deps_path
#         self.signals = WorkerSignals()
#
#     def run(self):
#         try:
#             print('here')
#             from src.nrgdock.process_target import main as process_target
#             process_target(self.target_path, self.targets, overwrite=self.overwrite,
#                            run_getcleft=self.run_getcleft, deps_path=self.deps_path)
#             self.signals.progress.emit("Processing completed successfully.")
#         except Exception as e:
#             print(e)
#             self.signals.progress.emit(f"Error occurred: {str(e)}")
#
#         self.signals.finished.emit()
#
# # Using QThreadPool
# def start_processing():
#     worker = WorkerRunnable('/Users/thomasdescoteaux/Documents/apo_holo/targets/holo',
#                             ['akt2'], overwrite=True, run_getcleft=False,
#                             deps_path='/Users/thomasdescoteaux/pymol-open-source-build/lib/python/pmg_tk/startup/NRGSuite_Qt/deps/nrgdock')
#     # Use a global or instance-level QThreadPool
#     thread_pool = QThreadPool.globalInstance()
#     thread_pool.start(worker)
import multiprocessing as mp
mp.set_start_method("fork")
from PyQt5.QtCore import QThread, pyqtSignal, QObject


class Runner(QObject):
    # Qt signal to communicate with the main thread
    progress = pyqtSignal(str)
    finished = pyqtSignal()

    def __init__(self, target_path, targets, overwrite, run_getcleft, deps_path):
        super().__init__()
        self.target_path = target_path
        self.targets = targets
        self.overwrite = overwrite
        self.run_getcleft = run_getcleft
        self.deps_path = deps_path
        self.run()

    def run(self):
        from src.nrgdock.process_target import main as process_target
        print('run')
        try:
            p = mp.Process(target=process_target, args=(self.target_path, self.targets, self.overwrite,
                               self.run_getcleft, None, self.deps_path))
            p.start()
        except Exception as e:
            print(f"An error occurred while starting the process: {e}")



def start_processing():
    # Create the worker thread
    from src.nrgdock.process_target import main as process_target
    import subprocess
    # process_target('/Users/thomasdescoteaux/Documents/apo_holo/targets/holo', ['akt2'], True, False, None, '/Users/thomasdescoteaux/pymol-open-source-build/lib/python/pmg_tk/startup/NRGSuite_Qt/deps/nrgdock')
    print('Running subprocess')
    subprocess.Popen(['python', '/Users/thomasdescoteaux/pymol-open-source-build/lib/python/pmg_tk/startup/NRGSuite_Qt/src/nrgdock/process_target.py', '-p', '/Users/thomasdescoteaux/Documents/apo_holo/targets/holo', '-t', 'akt2', '-o', '-d', '/Users/thomasdescoteaux/pymol-open-source-build/lib/python/pmg_tk/startup/NRGSuite_Qt/deps/nrgdock'])
    # runner_thread = QThread()
    # runner = Runner('/Users/thomasdescoteaux/Documents/apo_holo/targets/holo', ['akt2'], True, False, '/Users/thomasdescoteaux/pymol-open-source-build/lib/python/pmg_tk/startup/NRGSuite_Qt/deps/nrgdock')
    # # runner.msg_from_job.connect(handle_msg)
    # runner.moveToThread(runner_thread)


