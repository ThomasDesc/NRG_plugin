import threading
import time
import random

# Function that processes data in chunks
def process_data(data_chunk):
    sleep_time = random.randint(0, 10)
    print(f"Sleeping for {sleep_time} seconds.")
    time.sleep(sleep_time)  # Simulate a delay in processing
    for item in data_chunk:
        print(f"Processing: {item}")

# Function to run the processing in parallel, non-blocking version
def run_in_parallel_non_blocking(data_list, num_threads=4):
    # Splitting the list into chunks
    chunk_size = len(data_list) // num_threads
    chunks = [data_list[i:i + chunk_size] for i in range(0, len(data_list), chunk_size)]

    # Handling remainder if list length is not perfectly divisible
    if len(chunks) > num_threads:
        chunks[num_threads - 1].extend(chunks[num_threads])
        chunks = chunks[:num_threads]

    # Function that creates and runs threads in parallel
    def thread_worker():
        threads = []
        for i in range(num_threads):
            thread = threading.Thread(target=process_data, args=(chunks[i],))
            threads.append(thread)
            thread.start()

        # Optionally, you can wait for all threads to finish if needed (but without blocking PyMOL)
        for thread in threads:
            thread.join()

    # Running thread_worker in a background thread to avoid blocking the main thread (PyMOL)
    background_thread = threading.Thread(target=thread_worker)
    background_thread.daemon = True  # Daemon thread ensures it won't block PyMOL when closing
    background_thread.start()

# Example usage of the non-blocking parallel processing
def main():
    my_list = list(range(1, 100))  # Replace with your actual data list
    run_in_parallel_non_blocking(my_list, num_threads=4)
    print("Done")