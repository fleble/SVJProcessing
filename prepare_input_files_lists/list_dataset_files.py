import argparse
import csv 
from importlib import import_module
from pathlib import Path

import uproot

from utils.misc import get_files_list, process_in_parallel
from utils.Logger import *


def __get_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-d", "--datasets",
        help="Comma-separated list of dataset names for which to produce files list",
        required=True,
        type=lambda x: x.split(","),
    )
    parser.add_argument(
        "-y", "--year",
        help="Data-taking year",
        required=True,
    )
    parser.add_argument(
        '-nano', '--nano_aod',
        help='Set True if input files are (PF)NanoAOD files', 
        default=False, 
        action='store_true',
    )
    parser.add_argument(
        "-c", "--config",
        help="Config file describing the location of the datasets",
        required=True,
    )
    parser.add_argument(
        "-o", "--output",
        help="Output directory where the dataset files list will be stored",
        required=True,
    )
    parser.add_argument(
        "-n", "--n_workers",
        help="Number of cores to use in parallel, default=%(default)s",
        type=int,
        default=6,
    )
    
    return parser.parse_args()


def __run_bash_command(bash_command):
    return subprocess.Popen(bash_command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")[:-1]


def __get_files_list_from_info_dict(info_dict,nano_aod):
    redirector = info_dict["redirector"]
    path = info_dict["path"]
    regex = info_dict["regex"]
    bash_command = f"xrdfs {redirector} ls {path} | grep -E \"{regex}\" | sort -n"
    files_list = __run_bash_command(bash_command).split("\n")
    files_list = [f"{redirector}{file_name}" for file_name in files_list]
    files_list = [file_name for file_name in files_list if "log" not in file_name]
    files_list = [file_name for file_name in files_list if "test.sh" not in file_name]
    if nano_aod:
        if "nano_data" in files_list[0].split("/")[-1]:
            files_list = sorted(files_list, key=lambda x: int(x.split("/")[-1].split("_")[2].replace(".root","")))
        elif "PFNanoSuper" in files_list[0].split("/")[-1]:
            files_list = sorted(files_list, key=lambda x: int(x.split("/")[-1].split("-")[2].replace(".root","")))
        elif "PFNanoAOD_SVJL_" in files_list[0].split("/")[-1]:
            files_list = sorted(files_list, key=lambda x: int(x.split("/")[-1].split("-")[6].replace(".root","")))
        elif "PFNanoAOD_SVJtaus_" in files_list[0].split("/")[-1]:
            files_list = sorted(files_list, key=lambda x: int(x.split("/")[-1].split("-")[6].replace(".root","")))
        elif "PFNANOAOD" in files_list[0].split("/")[-1]:
            files_list = sorted(files_list, key=lambda x: int(x.split("/")[-1].split("-")[2].replace(".root","")))
        elif "PFNanoAOD" in files_list[0].split("/")[-1]:
            files_list = sorted(files_list, key=lambda x: int(x.split("/")[-1].split("-")[2].replace(".root","")))
        else:
            files_list = sorted(files_list, key=lambda x: int(x.split("/")[-1].split("_")[0]))
    else:
        files_list = sorted(files_list, key=lambda x: int(x.split("/")[-1].split("_")[0]))
    return files_list


def __list_files(dataset_info,nano_aod):
    files_list = []
    for info_dict in dataset_info:
        files_list += __get_files_list_from_info_dict(info_dict, nano_aod)


    return files_list


def __get_number_of_events(file_name, nano_aod, tree_name=""):  
    file_ = uproot.open(file_name)
    if nano_aod:
        tree_name = "Events"
    else:
        tree_name = "TreeMaker2/PreSelection"
    events = file_[tree_name]
    f0 = events.keys()[0]
    return events[f0].num_entries


def __write_dataset_info(
        dataset,
        dataset_info,
        year,
        output_directory,
        n_workers,
        nano_aod=False,
    ):
    
    files_list = __list_files(dataset_info,nano_aod)

    output_directory_ = f"{output_directory}/files_list/{year}"
    Path(output_directory_).mkdir(parents=True, exist_ok=True)


    files_list_batches = []
    max_n_files = 5000
    for i in range(1 + len(files_list)//max_n_files):
        files_list_batches.append(files_list[i * max_n_files: (i+1) * max_n_files])

    number_of_events = []
    for i, files_list_batch in enumerate(files_list_batches):
        log.info(f"Processing batch {i+1}/{len(files_list_batches)}...")
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
            number_of_events += list(tqdm(executor.map(__get_number_of_events, files_list_batch, len(files_list_batch)*[nano_aod]), total=len(files_list_batch)))


    output_file_name = f"{output_directory_}/{dataset}.csv"
    header = ["file_name", "number_of_events"]
    with open(output_file_name, "w") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(header)
        for file_name, n_events in zip(files_list, number_of_events):
            writer.writerow([file_name, n_events])

    log.info(f"{output_file_name} was written.")


def main ():

    args = __get_arguments()
    datasets_info = import_module(args.config).datasets_info
    datasets = args.datasets

    for dataset in datasets:
        __write_dataset_info(
            dataset,
            datasets_info[args.year][dataset],
            args.year,
            args.output,
            args.n_workers,
            args.nano_aod,
        )


if __name__ == "__main__":
    main()

