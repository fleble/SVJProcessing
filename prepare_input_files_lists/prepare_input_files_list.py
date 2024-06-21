import os
import argparse
import csv
from pathlib import Path

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
        "-s", "--selection",
        help="Name of the selection applied (as defined in previous step). "
             "If None, selection efficiency will be taken as 1.",
        required=False,
    )
    parser.add_argument(
        "-y", "--year",
        help="Data-taking year",
        required=True,
    )
    parser.add_argument(
        "-i", "--input",
        help="Input directory where the dataset files were stored",
        required=True,
    )
    parser.add_argument(
        "-o", "--output",
        help="Output directory where the dataset files list will be stored",
        required=True,
    )
    parser.add_argument(
        "-m", "--max_events",
        help="Maximum number of events per file after selection, default=%(default)s.",
        default=50000,
        type=int,
    )
    
    return parser.parse_args()


def __prepare_input_files_list(
        dataset,
        selection,
        year,
        max_events,
        input_directory,
        output_directory,
    ):

    input_files_list = f"{input_directory}/files_list/{year}/{dataset}.csv"
    if not os.path.exists(input_files_list):
        log.critical(f"Input files list {input_files_list} does not exist!")
        exit(1)

    selection_efficiency_file = f"{input_directory}/selections/{year}/{selection}/{dataset}.txt"
    if not os.path.exists(selection_efficiency_file):
        log.critical(f"Input selection file {selection_efficiency_file} does not exist!")
        exit(1)
    with open(selection_efficiency_file) as file:
        lines = file.readlines()
        selection_efficiency = float(lines[0][:-1])
 
    output_files_list = []
    with open(input_files_list) as f:
        csv_reader = csv.reader(f)
        sum = 0
        output_files_list.append([])
        next(csv_reader)
        for line in csv_reader:
            file_name = line[0]
            n_events = int(line[1])
            n_selected_events = n_events * selection_efficiency
            if sum + n_selected_events <= max_events:
                output_files_list[-1].append(file_name)
            else:
                output_files_list.append([file_name])
                sum = 0
            sum += n_selected_events
        
    output_directory_ = f"{output_directory}/skim_input_files_list/{year}/{selection}/{dataset}"
    Path(output_directory_).mkdir(parents=True, exist_ok=True)
    
    n_input_files = 0
    for x in output_files_list: n_input_files += len(x)
    n_output_files = len(output_files_list)
    log.info(f"Preparing skim files list, merging {n_input_files} input files "
             f"into {n_output_files} output skimmed files")

    for idx, files_list in enumerate(output_files_list):
        output_file_name = f"{output_directory_}/part-{idx}.txt"
        with open(output_file_name, "w") as output_file:
            for file_name in files_list:
                output_file.write(f"{file_name}\n")
        log.info(f"{output_file_name} was written.")


def main ():

    args = __get_arguments()
    datasets = args.datasets

    for dataset in datasets:
        args.primary_dataset = dataset
        __prepare_input_files_list(
            dataset,
            args.selection,
            args.year,
            args.max_events,
            args.input,
            args.output,
        )


if __name__ == "__main__":
    main()

