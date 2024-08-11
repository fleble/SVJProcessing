import argparse
from functools import partial

import pandas as pd
import numpy as np
import uproot

from utils.misc import get_files_list, process_in_parallel
from utils.coffea.n_tree_maker_schema import BaseSchema, NTreeMakerSchema
from coffea.nanoevents import NanoEventsFactory
from skimmer import skimmer_utils
from utils.Logger import *


def __get_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i", "--input_files_list",
        help="Path to CSV file with list of all files to skim",
        required=True,
    )
    parser.add_argument(
        "-o", "--output_directory",
        help="Output directory with all files",
        required=True,
    )
    parser.add_argument(
        '-nano', '--nano_aod',
        help='Set True if input files are (PF)NanoAOD files', 
        default=False, 
        action='store_true',
    )
    parser.add_argument(
        "-n", "--n_workers",
        help="Number of worker nodes (default=%(default)s)",
        default=4,
        type=int,
    )
    parser.add_argument(
        "-t", "--tree_name",
        help="To override the default event tree name in input files",
    )

    return parser.parse_args()


def __get_number_of_events_from_event_tree(file_name, tree_name, schema):
    events = NanoEventsFactory.from_root(
        file=file_name,
        treepath=tree_name,
        schemaclass=schema,
    ).events()
    
    return skimmer_utils.get_number_of_events(events)


def __get_number_of_events_from_cutflow(file_name):
    cut_flow = uproot.open(f"{file_name}:CutFlow").arrays()
    return cut_flow["Initial"][0]


def __compare_n_significant_digits(a, b, n=8):
    exp_a = np.floor(np.log10(a))
    exp_b = np.floor(np.log10(b))
    a_prime = int(a * 10**(n-exp_a-1))
    b_prime = int(b * 10**(n-exp_b-1))
    return a_prime == b_prime


def main():
    args = __get_arguments()
    
    output_redirector = "/".join(args.output_directory.split("/")[:3]) + "/"
    output_path = "/".join(args.output_directory.split("/")[3:]) 
    output_files = get_files_list(
        path=output_path,
        redirector=output_redirector,
    )
    output_files = sorted(output_files, key=lambda x: int(x.split("/")[-1].split("-")[-1].split(".")[0]))

    log.info(f"Input files list: {args.input_files_list}")
    log.info("Found the following output files:")
    for x in output_files:
        log.info(f"\t{x}")

    if args.nano_aod:
        schema = BaseSchema
    else:
        schema = NTreeMakerSchema

    if args.tree_name:
        input_tree_name = args.tree_name
    else:
        if args.nano_aod:
            input_tree_name = "Events"
        else:
            input_tree_name = "TreeMaker2/PreSelection"

    events_file_0 = NanoEventsFactory.from_root(
        file=output_files[0],
        treepath="Events",
        schemaclass=schema,
        entry_stop=10,
    ).events()
    is_data = skimmer_utils.is_data(events_file_0)

    df = pd.read_csv(args.input_files_list)

    if is_data:
        n_events_input = np.sum(df["number_of_events"])

    else:
        process_function = partial(__get_number_of_events_from_event_tree, tree_name=input_tree_name, schema=schema)
        files_list = df["file_name"].to_list()
        n_events_per_file = process_in_parallel(
            files_list=files_list,
            process_function=process_function,
            n_workers=args.n_workers,
        )
        n_events_input = np.sum(n_events_per_file)
    
    n_events_per_file = process_in_parallel(
        files_list=output_files,
        process_function=__get_number_of_events_from_cutflow,
        n_workers=args.n_workers,
    )
    n_events_output = np.sum(n_events_per_file)
 
    if is_data:
        fmt = "%d"
    else:
        fmt = "%.4f"
    if n_events_input == n_events_output:
        log.info("Same number of initial events in input and output!")
        log.info("N events: " + fmt %(n_events_input))
    elif 1 - n_events_output / n_events_input < 0.01:
        log.warning("Subpercent difference between initial number of events in input and output!")
        log.warning("N events input: " + fmt %(n_events_input))
        log.warning("N events output: " + fmt %(n_events_output))
        for n in range(8, 4, -1):
            if __compare_n_significant_digits(n_events_input, n_events_output, n=n):
                log.info(f"This seems to be due to numerical rounding as numbers match up to {n} significant digits!")
                break
    else:
        log.error("Large difference between initial number of events in input and output!")
        log.error("N events input: " + fmt %(n_events_input))
        log.error("N events output: " + fmt %(n_events_output))


if __name__ == "__main__":
    main()
