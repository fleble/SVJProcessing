# SVJProcessing

## Overview

This repo is meant to produce skim from NanoAOD-like NTuples by performing event-selection, branch addition and merging of several input files into a single output skim file.

Input files to process are configurable, see for instance [dataset_configs/t_channel_datasets_paths.py](https://github.com/fleble/SVJProcessing/blob/main/dataset_configs/t_channel_datasets_paths.py).   
The event-selection and branch addition are configurable, see for instance [analysis_configs/t_channel_pre_selection.py](https://github.com/fleble/SVJProcessing/blob/main/analysis_configs/t_channel_pre_selection.py).    
The number of input files to merge is automatically determined depending on the selection efficiency and the target maximum number of events in the output skim.

For the moment only TreeMaker NTuples are supported, but only minimal changes are required to process any nanoAOD-like NTuples format.


## Installation and environment setup

Install the repo:
```
git clone git@github.com:fleble/SVJProcessing.git
```
At each new login do:
```
source setup.sh
```

The environment to run the code is borrowed from the [t-channel analysis framework](https://github.com/cms-svj/t-channel_Analysis/tree/master). Follow instructions there to install the repo and initialize the singularity container.


## Preparing list of input files to skim

Fill in the dataset config, see instructions in [analysis_configs/t_channel_pre_selection.py](https://github.com/fleble/SVJProcessing/blob/main/analysis_configs/t_channel_pre_selection.py). Each process / "datatet" must only list files with same cross-section! E.g. different QCD bins are treated as different datasets.    
Prepare the skim preparation config, see e.g. [analysis_configs/t_channel_pre_selection.py](https://github.com/fleble/SVJProcessing/blob/main/analysis_configs/t_channel_pre_selection.py)

The preparation on the input files proceeds in 3 steps:     
1- Fetching the number of events in each input file    
2- Computing the expected unweighted selection efficiency    
3- Preparing the list of N files to skim into one output skim file such that each output skim file has at most N events.    

The goal is that:    
1- the number of skim files is reduced compared to the number of initial input files to reduce I/O operations     
2- all output skim files have the same size and can be processed with similar memory request

If the cut efficiency is very low, it can take time to estimate the selection efficiency.

An example bash script to perform these steps is provided in [prepare_input_files_lists/prepare_input_files_list.sh](https://github.com/fleble/SVJProcessing/blob/main/prepare_input_files_lists/prepare_input_files_list.sh).


## Making skims

An example bash script to produce skims is provided in [skimmer/make_skims.sh](https://github.com/fleble/SVJProcessing/blob/main/skimmer/make_skims.sh).

You can run locally on the interactive nodes where you are logged in, or distributed, using the LPCConderCluster (only at LPC!) or the SLURMCluster at any site with a SLURM batch submission system:    
* If you run locally, you need to adjust the number of workers to not use too many and bother other users on that login node!
* If you run distributed, you need to adjust the memory request: the more you ask, the more your jobs will queue, but if you do not request enough the jobs will run out of memory and crash. For 50k events per file and TreeMaker NTuples containing PF candidates, 4 GB of memory should be enough!

For processes with high efficiency, the bottleneck is the queuing time... It is much faster to run locally requesting one node if the queuing time is non-zero...
