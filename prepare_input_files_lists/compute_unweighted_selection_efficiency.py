import os
import argparse
import csv
from pathlib import Path
import math

from skimmer.skim import add_coffea_args, skim
from utils.Logger import *


def __get_arguments():
    parser = argparse.ArgumentParser(description="Compute *unweighted* selection efficiency")

    parser.add_argument(
        "-d", "--datasets",
        help="Comma-separated list of dataset names for which to produce files list",
        required=True,
        type=lambda x: x.split(","),
    )
    parser.add_argument(
        "-s", "--selection",
        help="Name of the selection applied (short hand for the selection config file name)",
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
        "-precision", "--precision",
        help="Minimum relative uncertainty on the efficiency calculation, default=%(default)s.",
        default=5,
        type=int,
    )
    
    add_coffea_args(parser)

    return parser.parse_args()


def __read_files_list_csv_file(files_list):
    input_files_names = []
    number_of_events = []
    with open(files_list) as f:
        csv_reader = csv.reader(f)
        next(csv_reader)
        for line in csv_reader:
            file_name = line[0]
            n_events = int(line[1])
            input_files_names.append(file_name)
            number_of_events.append(n_events)

    return input_files_names, number_of_events
 

def __compute_uncertainty(n_initial, n_final):
    k = n_final
    n = n_initial
    eff = k / n
    return math.sqrt(eff * (1 - eff) / n)


def __compute_effiency(
        input_files_names,
        number_of_events,
        precision,
        coffea_args,
    ):

    relative_uncertainty = 100
    n_input_files = 1
    while relative_uncertainty > precision:
        log.info(f"Computing (unweighted) efficiency for {n_input_files} file(s)...")
        accumulator = skim(
            input_files_list=input_files_names[:n_input_files],
            coffea_args=coffea_args,
        )

        n_initial = sum(number_of_events[:n_input_files])
        n_final = len(accumulator["events"].value)
        efficiency = n_final / n_initial
        log.info(f"Efficiency: {efficiency}")
        if efficiency == 0:
            relative_uncertainty = 100
        else:
            uncertainty = __compute_uncertainty(n_initial, n_final)
            relative_uncertainty = 100 * (uncertainty / efficiency)

        # If relative uncertainty less than required precision,
        # increase the number of events to exceed the target
        log.info(f"Relative uncertainty: {relative_uncertainty:.2f}%")
        if relative_uncertainty > precision:
            if n_input_files < len(input_files_names):
                log.info(f"Will increase number of files to reach target "
                         f"precision of {precision}%")
            else:
                log.info(f"Cannot reach target precision of {precision}% "
                         f"as maximum number of files as been reached!")
                return efficiency

            if efficiency == 0:
                n_input_files *= 10
            else:
                uncertainty_target = precision * efficiency / 100
                n_events_target = efficiency * (1 - efficiency) / (uncertainty_target ** 2)
        
                n_events = 0
                current_n_input_files = n_input_files
                for i_input_files in range(len(input_files_names)):
                    n_events += number_of_events[i_input_files]
                    if n_events > n_events_target:
                        n_input_files = math.ceil((i_input_files + 1) * 1.1)  # +1 to exceed precision
                        break
                n_input_files = i_input_files + 1
 
    return efficiency


def __compute_unweighted_selection_efficiency(
        dataset,
        selection,
        precision,
        coffea_args,
        input_directory,
        output_directory,
    ):

    input_files_list = f"{input_directory}/files_list/{coffea_args.year}/{dataset}.csv"
    if not os.path.exists(input_files_list):
        log.critical(f"Input files list {input_files_list} does not exist!")
        exit(1)

    input_files_names, number_of_events = __read_files_list_csv_file(input_files_list)

    efficiency = __compute_effiency(
        input_files_names,
        number_of_events,
        precision,
        coffea_args,
    )
       
    output_directory_ = f"{output_directory}/selections/{coffea_args.year}/{selection}"
    Path(output_directory_).mkdir(parents=True, exist_ok=True)

    output_file_name = f"{output_directory_}/{dataset}.txt"
    with open(output_file_name, "w") as output_file:
        output_file.write(f"{efficiency}\n")
    
    log.info(f"{output_file_name} was written.")


def main ():

    args = __get_arguments()
    datasets = args.datasets

    for dataset in datasets:
        __compute_unweighted_selection_efficiency(
            dataset,
            args.selection,
            args.precision,
            args,
            args.input,
            args.output,
        )


if __name__ == "__main__":
    main()

