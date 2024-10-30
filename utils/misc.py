import subprocess
import concurrent.futures
from tqdm import tqdm

from utils.Logger import *


def run_bash_command(bash_command):
    return subprocess.Popen(bash_command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")[:-1]


def get_files_list(path, redirector=None, regex=".*"):
    bash_command = f"ls {path} | grep -E \"{regex}\" | sort -n"
    if redirector is not None:
        bash_command = f"xrdfs {redirector} {bash_command}"
    files_list = run_bash_command(bash_command).split("\n")
    if redirector is not None:
        files_list = [f"{redirector}{file_name}" for file_name in files_list]
    else:
        files_list = [f"{path}/{file_name}" for file_name in files_list]

    return files_list


def process_in_parallel(files_list, process_function, n_workers, max_n_files=5000):
    files_list_batches = []
    for i in range(1 + len(files_list)//max_n_files):
        files_list_batches.append(files_list[i * max_n_files: (i+1) * max_n_files])

    results = []
    for i, files_list_batch in enumerate(files_list_batches):
        log.info(f"Processing batch {i+1}/{len(files_list_batches)}...")
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
            results += list(tqdm(executor.map(process_function, files_list_batch), total=len(files_list_batch)))
    
    return results