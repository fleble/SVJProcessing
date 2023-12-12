import argparse
from importlib import import_module

import awkward as ak
from coffea import processor

import utils.uprootUtilities as uprootUtl
from utils.coffea.akArrayAccumulator import AkArrayAccumulator
from utils.coffea.NTreeMakerSchema import NTreeMakerSchema
from skimmer import skimmerUtils


class Skimmer(processor.ProcessorABC):

    def __init__(self, process_function):
        self.accumulator = AkArrayAccumulator()
        self.converted_branches = []
        self.process_function = process_function

        
    def process(self, events):

        cut_flow = {}
        skimmerUtils.update_cut_flow(cut_flow, "Initial", events)

        events, cut_flow = self.process_function(events, cut_flow)
        skimmerUtils.update_cut_flow(cut_flow, "Final", events)

        accumulator = {
            "events": AkArrayAccumulator(ak.copy(events)),
            "cut_flow": cut_flow,
        }

        return accumulator


    def postprocess(self, accumulator):
        return super().postprocess(accumulator)


def __get_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i", "--input_file_names",
        help="Input file names",
        required=True,
    )
    parser.add_argument(
        "-o", "--output_file_name",
        help="Output file name",
        required=True,
    )
    parser.add_argument(
        "-m", "--process_module_name",
        help="Process module name, e.g. analysisConfigs.config",
        required=True,
    )
    parser.add_argument(
        "-y", "--year",
        help="Data-taking year",
        required=True,
    )

    return parser.parse_args()


def main():

    args = __get_arguments()
    input_file_names = args.input_file_names.split(",")
    process_module = import_module(args.process_module_name)
    process_function = lambda x, y: process_module.process(x, y, year=args.year)

    # TODO: replace by read arguments
    executor = processor.iterative_executor
    executor_args = {
        "schema": NTreeMakerSchema,
        "workers": 1,
    }
    chunk_size = 400
    max_chunks = None

    # Calculate new branches
    accumulator = processor.run_uproot_job(
        {"fileset": input_file_names},
        treename="TreeMaker2/PreSelection",
        processor_instance=Skimmer(process_function),
        executor=executor,
        executor_args=executor_args,
        chunksize=chunk_size,
        maxchunks=max_chunks,
        )


    # Making output ROOT file
    trees = {
        "CutFlow": accumulator["cut_flow"]
    }

    uprootUtl.write_tree_maker_root_file(
        output_file_name=args.output_file_name,
        events=accumulator["events"].value,
        trees=trees,
    )


if __name__ == "__main__":
    main()
