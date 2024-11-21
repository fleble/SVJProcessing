import argparse
from importlib import import_module

import awkward as ak
from coffea import processor
from coffea.nanoevents import BaseSchema

import utils.uproot_utilities as uproot_utl
from utils.coffea.ak_array_accumulator import AkArrayAccumulator
from utils.coffea.dict_accumulator import DictAccumulator
from utils.coffea.n_tree_maker_schema import NTreeMakerSchema
from utils.coffea.job_submission_helper import get_executor, get_executor_args
from skimmer import skimmer_utils
from utils.Logger import *


class Skimmer(processor.ProcessorABC):

    def __init__(
            self,
            process_function,
            year,
            is_mc=False,
            variation=None,
            weight_variations=[],
            nano_aod=False,
            pfnano_corr_file=None,
        ):
        self.process_function = process_function
        self.variation = variation
        self.weight_variations = weight_variations
        self.nano_aod = nano_aod
        self.pfnano_corr_file = pfnano_corr_file
        self.is_mc = "mc" if is_mc else "data"
        self.year = year

    def process(self, events):

        cut_flow = {}
        skimmer_utils.update_cut_flow(cut_flow, "Initial", events)

        if skimmer_utils.is_mc(events):
            # Calculate and store the weight variations
            if "scale" in self.weight_variations:
                events, sumw_scale_up, sumw_scale_down = skimmer_utils.apply_scale_variations(events,is_nano=self.nano_aod)
                skimmer_utils.update_cut_flow(cut_flow, "InitialScaleUp", sumw=sumw_scale_up)
                skimmer_utils.update_cut_flow(cut_flow, "InitialScaleDown", sumw=sumw_scale_down)
            
            
            if "pdf" in self.weight_variations:
                events, sumw_pdf_up, sumw_pdf_down = skimmer_utils.apply_pdf_variations(events,is_nano=self.nano_aod)
                skimmer_utils.update_cut_flow(cut_flow, "InitialPDFUp", sumw=sumw_pdf_up)
                skimmer_utils.update_cut_flow(cut_flow, "InitialPDFDown", sumw=sumw_pdf_down)

        if self.nano_aod:
            events = skimmer_utils.apply_variation_pfnano(events, variation=self.variation, year=self.year, run=self.is_mc, pfnano_sys_file=self.pfnano_corr_file)
        else:
            events = skimmer_utils.apply_variation(events, self.variation)

        events, cut_flow = self.process_function(events, cut_flow)
        skimmer_utils.update_cut_flow(cut_flow, "Final", events)

        accumulator = {
            "events": AkArrayAccumulator(ak.copy(events)),
            "cut_flow": DictAccumulator(cut_flow.copy()),
        }

        return accumulator


    def postprocess(self, accumulator):
        return super().postprocess(accumulator)


def add_coffea_args(parser):
    parser.add_argument(
        "-p", "--process_module_name",
        help="Process module name, e.g. analysis_configs.t_channel_pre_selection",
        required=True,
    )
    parser.add_argument(
        "-y", "--year",
        help="Data-taking year",
        required=True,
    )
    parser.add_argument(
        "-mc", "--is_mc",
        help='Set True if input files MC files', 
        default=False, 
        action='store_true',
    )
    parser.add_argument(
        "-pd", "--primary_dataset",
        help="If data, name of the primary dataset",
        default="",
    )
    parser.add_argument(
        "-c", "--chunk_size",
        help="Size of the data chunks (default=%(default)s)",
        default=10000,
        type=int,
    )
    parser.add_argument(
        "-m", "--max_chunks",
        help="Maximum number of chunks to process, no flag means no maximum",
        type=int,
    )
    parser.add_argument(
        "-n", "--n_workers",
        help="Number of worker nodes (default=%(default)s)",
        default=6,
        type=int,
    )
    parser.add_argument(
        "-nc", "--cores",
        help="Number of cores request for HTCondor or SLURM (default=%(default)s)",
        default=2,
        type=int,
    )
    parser.add_argument(
        "-mem", "--memory",
        help="Memory request for HTCondor or SLURM (default=%(default)s)",
        default="4GB",
    )
    parser.add_argument(
        "-t", "--walltime",
        help="Time requested for HTCondor or SLURM (default=%(default)s)",
        default="00:30:00",
    )
    parser.add_argument(
        "-q", "--queue",
        help="Queue for HTCondor or SLURM (default=%(default)s)",
        default="00:30:00",
    )
    parser.add_argument(
        "-disk", "--disk",
        help="Disk space request for HTCondor (default=%(default)s)",
        default="100MB",
    )
    parser.add_argument(
        "--skip_bad_files",
        help="Skip bad files",
        action="store_true",
    )
    parser.add_argument(
        "-e", "--executor_name",
        choices=[
            "iterative",
            "futures",
            "dask/slurm",
            "dask/lpccondor",
            "dask/etpcondor",
        ],
        default="futures",
        help="The type of executor to use (default=%(default)s)"
    )
    parser.add_argument(
        '-skim_source', '--skim_source',       
        help='Set True if input files are skim files', 
        default=False, 
        action='store_true',
    )
    parser.add_argument(
        '-nano', '--nano_aod',
        help='Set True if input files are (PF)NanoAOD files', 
        default=False, 
        action='store_true',
    )
    parser.add_argument(
        '-nano_scout', '--nano_aod_scouting',
        help='Set True if input files are (PF)NanoAOD files for scouting', 
        default=False, 
        action='store_true',
    )
    parser.add_argument(
        '-xsec', '--cross_section',
        help='If cross-section not in file, adding it', 
        type=float,
        action='store',
    )

    parser.add_argument(
        '-corrfile', '--pfnano_corr_file',
        help='Precompiled corrections file for PFNanoAOD', 
        type=str,
        action='store',
        default=None,
    )

    parser.add_argument(
        '-pn_tagger', '--pn_tagger',       
        help='Add particleNet jet tagger score', 
        default=False, 
        action='store_true',
    )
    parser.add_argument(
        "-port", "--port",
        help="Port for dask distributed computation (default=%(default)s)",
        type=int,
        default=8787,
    )
    parser.add_argument(
        "-var", "--variation",
        help="What systematic variation to compute (choice=%(choices)s)",
        type=str,
        choices=["jec_up", "jec_down", "jer_up", "jer_down", "unclEn_up", "unclEn_down"],
        default=None,
    )
    parser.add_argument(
        "-wvar", "--weight_variations",
        help="What systematic weight variations to compute (choice=%(choices)s). "
             "Can choose several weights variations.",
        type=str,
        nargs="*",
        choices=["scale", "pdf"],
        default=[],
    )


def __get_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i", "--input_files",
        help="Input files, either comma-separated list of files or "
             "txt file with one file per line. Empty lines and lines "
             "starting with # will be ignored.",
        required=True,
    )
    parser.add_argument(
        "-o", "--output_file_name",
        help="Output file name",
        required=True,
    )

    add_coffea_args(parser)

    return parser.parse_args()


def __get_input_files(input_files_arg):
    if input_files_arg.endswith(".txt"):
        with open(input_files_arg) as file:
            lines = file.readlines()
            input_file_names = []
            for line in lines:
                if line.startswith("#"): continue
                if len(line) < 3: continue
                input_file_names.append(line[:-1])
    else:
        input_file_names = input_files_arg.split(",")
    
    return input_file_names


def __prepare_cut_flow_tree(cut_flow_dict):
    cut_flow_tree = {
        key: [value] for key, value in cut_flow_dict.items()
    }
    return cut_flow_tree


def __prepare_uproot_job_kwargs_from_coffea_args(args):

    process_module = import_module(args.process_module_name)
    process_function = lambda x, y: process_module.process(
        x,
        y,
        year=args.year,
        primary_dataset=args.primary_dataset,
        pn_tagger=args.pn_tagger,
    )

    executor = get_executor(args.executor_name)
    if args.nano_aod or args.nano_aod_scouting:
        executor_args = {
            "schema": BaseSchema,
        }
    else:
        executor_args = {
            "schema": NTreeMakerSchema,
        }
    executor_args.update(get_executor_args(
        executor_name=args.executor_name,
        n_workers=args.n_workers,
        cores=args.cores,
        memory=args.memory,
        disk=args.disk,
        time=args.walltime,
        partition=args.queue,
        skip_bad_files=args.skip_bad_files,
        port=args.port,
    ))

    if args.skim_source or args.nano_aod:
        treename = "Events"
    elif args.nano_aod_scouting:
        treename = "mmtree/Events"
    else:
        treename = "TreeMaker2/PreSelection"

    uproot_job_kwargs = {
        "treename": treename,
        "processor_instance": Skimmer(
            process_function,
            args.year,
            args.is_mc,
            args.variation,
            args.weight_variations,
            args.nano_aod,
            args.pfnano_corr_file,
        ),
        "executor": executor,
        "executor_args": executor_args,
        "chunksize": args.chunk_size,
        "maxchunks": args.max_chunks,
    }

    return uproot_job_kwargs


def skim(
        input_files_list,
        coffea_args,
    ):

    uproot_job_kwargs = __prepare_uproot_job_kwargs_from_coffea_args(coffea_args)

    # Skim files
    accumulator = processor.run_uproot_job(
        {"fileset": input_files_list},
        **uproot_job_kwargs
    )

    return accumulator


def main():

    args = __get_arguments()
    input_file_names = __get_input_files(args.input_files)

    accumulator = skim(input_file_names, args)

    # Making output ROOT file
    cut_flow_tree = __prepare_cut_flow_tree(accumulator["cut_flow"].value)
    if args.skim_source:
        # use original values from the skim's cutFlow
        cut_flow_tree = skimmer_utils.get_cut_flow_from_skims(args.input_files, cut_flow_tree)

    trees = {
        "CutFlow": cut_flow_tree
    }
   
    events = accumulator["events"].value
    if len(events) == 0:
        log.warning("No events passed selection")
        return

    if args.cross_section:
        # Add cross-section in the custom PFNanoAOD way
        trees["Metadata"] = {
            "GenCrossSection": [args.cross_section],
        }
        # But also in the same way as TreeMaker, which is better
        events["CrossSection"] = ak.Array([args.cross_section for i in range(len(events))])
 
    if args.nano_aod or args.nano_aod_scouting:
        uproot_utl.write_nano_aod_root_file(
            output_file_name=args.output_file_name,
            events=events,
            trees=trees,
        )
    else:
        uproot_utl.write_tree_maker_root_file(
            output_file_name=args.output_file_name,
            events=events,
            trees=trees,
        )


if __name__ == "__main__":
    main()

