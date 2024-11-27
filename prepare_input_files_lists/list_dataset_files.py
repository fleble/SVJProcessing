import argparse
import csv 
from importlib import import_module
from pathlib import Path
from functools import partial

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
        '-nano_scout','--nano_aod_scouting',
        help='Set True if input files are (PF)NanoAOD Scouting files', 
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


def __list_files(dataset_info, nano_aod):
    files_list = []
    for info_dict in dataset_info:
        files_list_ = get_files_list(
            path=info_dict["path"],
            redirector=info_dict["redirector"],
            regex=info_dict["regex"],
        )

        print("path: ", info_dict["path"])
        

        # TODO: Do something more clever should be done insteead of relying of the file name
        try:
            if nano_aod:
                if "nano_data" in files_list[0].split("/")[-1]:
                    files_list_ = sorted(files_list_, key=lambda x: int(x.split("/")[-1].split("_")[2].replace(".root","")))
                elif "PFNanoSuper" in files_list[0].split("/")[-1]:
                    files_list_ = sorted(files_list_, key=lambda x: int(x.split("/")[-1].split("-")[2].replace(".root","")))
                elif "PFNanoAOD_SVJL_" in files_list[0].split("/")[-1]:
                    files_list_ = sorted(files_list_, key=lambda x: int(x.split("/")[-1].split("-")[6].replace(".root","")))
                elif "PFNanoAOD_SVJtaus_" in files_list[0].split("/")[-1]:
                    files_list_ = sorted(files_list_, key=lambda x: int(x.split("/")[-1].split("-")[6].replace(".root","")))
                elif "PFNANOAOD" in files_list[0].split("/")[-1] and files_list[0].split("/")[-1].split("-")[2].replace(".root","").isdigit():
                    files_list_ = sorted(files_list_, key=lambda x: int(x.split("/")[-1].split("-")[2].replace(".root","")))
                elif "PFNANOAOD" in files_list[0].split("/")[-1] and files_list[0].split("/")[-1].split("-")[3].replace(".root","").isdigit():
                    files_list_ = sorted(files_list_, key=lambda x: int(x.split("/")[-1].split("-")[3].replace(".root","")))
                elif "PFNanoAOD" in files_list[0].split("/")[-1]:
                    files_list_ = sorted(files_list_, key=lambda x: int(x.split("/")[-1].split("-")[2].replace(".root","")))
                else:
                    files_list_ = sorted(files_list_, key=lambda x: int(x.split("/")[-1].split("_")[0]))
            else:
                files_list_ = sorted(files_list_, key=lambda x: int(x.split("/")[-1].split("_")[0]))
        except:
            pass

        files_list += files_list_

    return files_list


def __get_number_of_events(file_name, tree_name):
    
    file_ = uproot.open(file_name)
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
        nano_aod_scout=False,
    ):
    
    files_list = __list_files(dataset_info, nano_aod)
    output_directory_ = f"{output_directory}/files_list/{year}"
    Path(output_directory_).mkdir(parents=True, exist_ok=True)

    if nano_aod:
        tree_name = "Events"
    elif nano_aod_scout:
        tree_name = "mmtree/Events"
    else:
        tree_name = "TreeMaker2/PreSelection"

    process_function = partial(__get_number_of_events, tree_name=tree_name)
    number_of_events = process_in_parallel(
        files_list=files_list,
        process_function=process_function,
        n_workers=n_workers,
    )

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

    print("datasets: ", datasets)

    for dataset in datasets:
        __write_dataset_info(
            dataset,
            datasets_info[args.year][dataset],
            args.year,
            args.output,
            args.n_workers,
            args.nano_aod,
            args.nano_aod_scouting,
        )


if __name__ == "__main__":
    main()

