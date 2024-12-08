# SVJProcessing

## Overview

This repo is meant to produce skim from NanoAOD-like NTuples by performing event-selection, branch addition and merging of several input files into a single output skim file.

Input files to process are configurable, see for instance [dataset_configs/t_channel_datasets_paths.py](https://github.com/fleble/SVJProcessing/blob/main/dataset_configs/t_channel_datasets_paths.py).   
The event-selection and branch addition are configurable, see for instance [analysis_configs/t_channel_pre_selection.py](https://github.com/fleble/SVJProcessing/blob/main/analysis_configs/t_channel_pre_selection.py).    
The number of input files to merge is automatically determined depending on the selection efficiency and the target maximum number of events in the output skim.

For the moment only TreeMaker and NanoAOD NTuples are supported, as well as job submission with HTCondor at LPC, and with SLURM.


## Installation and environment setup

Install the repo:
```
git clone git@github.com:fleble/SVJProcessing.git
```
At each new login do:
```
source setup.sh
```

The environment to run the code at the LPC is borrowed from the [FNAL t-channel analysis framework](https://github.com/cms-svj/t-channel_Analysis/tree/master). Follow instructions there to install the repo and initialize the singularity container.

The environment to run the code at the PSI T3 is borrowed from the [ETH t-channel analysis framework](https://github.com/eth-svj/SVJanalysis). Follow instructions [here](https://github.com/jniedzie/SVJanalysis_wiki/wiki/Creating-SVJ-virtual-environment) to install the virtual environment.

To check your environment, run the automated tests (see below).


## Tests

Automated tests were implemented to check the execution of the code (input files list preparation and skimming).
To run those tests, do:
```
source setup.sh
cd tests/
./run_tests.sh
```


## Preparing list of input files to skim

#### Creating dataset config

Fill in the dataset config, see example in [analysis_configs/t_channel_pre_selection.py](https://github.com/fleble/SVJProcessing/blob/main/analysis_configs/t_channel_pre_selection.py). Each process / "datatet" must only list files with same cross section! E.g. different QCD bins are treated as different datasets. 

Example to decribe the location of 2016 JetHT data, stored as TreeMaker NTuples at FNAL:
```python
datasets_info["2016"] = {   # First key is the year
    "JetHT": [  # Seco<nd key if the name of the dataset
        {
            # Specify the location on the GRID or use None if files are available locally
            "redirector": "root://cmseos.fnal.gov/",
            # Path to the files
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016F-UL2016-v2/JetHT/", 
            # Can use a regular expression here to select only some files within that directory, e.g. useful o all your signal samples are together in the same directory
            "regex": "",
        },
        # Can add different locations adding as many of this dictionary to the list
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016G-UL2016-v2/JetHT/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016H-UL2016-v2/JetHT/",
            "regex": "",
        },
    ]
}
```

#### Creating skim config

Prepare the skim config, defining the selections to apply to produce the skims, see e.g. [analysis_configs/t_channel_pre_selection.py](https://github.com/fleble/SVJProcessing/blob/main/analysis_configs/t_channel_pre_selection.py).
To first test your config, you can use the following example script: [skimmer/examples.sh](https://github.com/fleble/SVJProcessing/blob/main/skimmer/example.sh). See instructions in the [next section](#making-skims).

The preparation on the input files proceeds in 3 steps:     
1- Fetching the number of events in each input file    
2- Computing the expected unweighted selection efficiency    
3- Preparing the list of N files to skim into one output skim file such that each output skim file has at most N events.    

The goal is that:    
1- the number of skim files is reduced compared to the number of initial input files to reduce I/O operations     
2- all output skim files have the same size and can be processed with similar memory request

If the cut efficiency is very low, it can take time to estimate the selection efficiency. In that case, you can just create the efficiency file by hand, located at `<dataset_directory>/selections/<year>/<selection_name>/<dataset_names>.txt` (replace the terms enclosed in `<>` by the relevant values), and write `1e-9` inside.

If the maximum number of events per output file is too high (`-m/--max_events` flag), then you will need to request a large amount of memory for the diferent chunks to be aggregated at the end of the jobs. A typical good number is:
* 50k for MC (tested on TreeMaker samples with PF candidates)
* 5k for data (tested on TreeMaker samples with PF candidates)

An example bash script to perform these steps is provided in [prepare_input_files_lists/prepare_input_files_list_t_channel.sh](https://github.com/fleble/SVJProcessing/blob/main/prepare_input_files_lists/prepare_input_files_list_t_channel.sh).    
It shows how to run the following scripts:
* list_dataset_files.py
* compute_unweighted_selection_efficiency.py
* prepare_input_files_list.py
For the usage of these scripts, do `python <script.py> -h`.

#### Important note when running on data
When running on data, in order to handle the overlap between the different primary datasets, the primary dataset name must be provided to the selection config. This is done by naming the subdirectory holding the list of files just as the primary dataset.


## Making skims
### Introduction

An example bash script to produce skims is provided in [skimmer/examples.sh](https://github.com/fleble/SVJProcessing/blob/main/skimmer/example.sh).

To see the full usage of `skimmer/skim.py`, do:
```python
cd skimmer
python skim.py -h
```

Some explanation of the flags relative to the job submission:
* -e/--executor_name: Where to run the job, use futures to run locally, dask/slurm to send jobs with SLURM, and dask/lpccondor to send jobs with HTCondor at LPC.
* -n/--n_workers: Number of threads used locally, or number of jobs if running on a batch system.
* -c/--chunk_size: Number of events processed in each chunk. If you set the memory request too low or the chunksize too large, then you will run out of memory while processing.
* -nc/--cores: Number of cores requested per job.
* -mem/--memory: Memory requested per job.
* -t/--walltime: The maximum time of the job (for SLURM only).
* -q/--queue: The parition name for SLURM jobs.

Note if you are running on PFNano for offline analysis or scouting, set the relative flags ```-nano``` or ```-nano_scout``` as options to ```skim.py``` as shown in the examples ```make_skims_pfnano.sh``` and ```make_skims_s_channel_scouting.sh```. In those examples can be found also the flags to set the x-section of a given input dataset which will be used to compute the weights.

You can run locally on the interactive nodes where you are logged in, or distributed, using the LPCConderCluster (only at LPC!) or the SLURMCluster at any site with a SLURM batch submission system:    
* If you run locally, you need to adjust the number of workers to not use too many and bother other users on that login node!
* If you run distributed, you need to adjust the memory request: the more you ask, the more your jobs will queue, but if you do not request enough the jobs will run out of memory and crash. For 50k events per file and TreeMaker NTuples containing PF candidates, 4 GB of memory should be enough! When running ParticleNet inference with 50k events per file, 6 GB of memory is needed.

For processes with high efficiency, the bottleneck is the queuing time... It is much faster to run locally requesting one node if the queuing time is non-zero...

If the total number of workers is too high, the OS limit of maximum number of files opened by one process is reached. At FNAL, this limit seems to be around 200 workers. It is recommended to ask for a maximum of 150 workers at once. It seems that the code terminates well if asking for more workers, but workers will all crash throwing an OSError...

### Running systematics
When running on ```Treemaker``` ntuples the systematics variations need not to be recomputed, and the code fetches the necessary information form the input files to create the varied skimmed datasets. An example of how to run the systematics variations for ```Tremaker``` can be found ```make_skims_t_channel.sh```.

When running on ```(PF)NanoAOD``` ntuples, the systematics  need to be recomputed. This is handled in the framework using [coffea](https://coffeateam.github.io/coffea/) tools and [correctionlib](https://github.com/cms-nanoAOD/correctionlib). There are two main steps to run the skims with the systematics variations:
1. Produce a ```.coffea``` file with the whole set of compiled corrections (corrections factory). This can be done by running the code ```skimmer/make_compiled_corrections_file_pfnano.sh```. The bash script will run ```corrections_pfnano_coffea.py``` which will create a ```data``` directory in the repository with inside all the needed ```.txt``` files to be used to create the corrections factory, and it will output in the same directory also the ```.coffa``` file. It is possible to include in the corrections factory specific corrections, such as JME ones, with flags like ```-jme``` or include all possible ones implemented in ```corrections_pfnano_coffea.py```  with ```-all```. Suported arguments for ```corrections_pfnano_coffea.py```: ```-all``` (include all implemented corrections), ```-jme``` (include only JME corrections)
2. Run the skimmer code. An example can be found in ```make_skims_pfnano_sys.sh```. There you will need to specify the full path to the corrections factory file in ```pfnano_corrections_file```. In order to add ```pdfs``` and ```scales``` variations you can set also ```add_weights_variations=1```. When running variations it is important to set the ```-mc``` flag.

#### Tips for running at LPC

At LPC, to avoid the code to stop when closing the terminal, run the code in a `screen` session:
```
screen -S <name_of_session>
<set up the environment>
cd skimmer
./make_skims.sh
Ctrl + A + D   # to exit the screen session
```
Then you can close your terminal, disconnect from internet, the code will run in the screen session. Remember the machine on which you logged in, you need to log in on the same machine to access your screen session again.

To access the screen session again, log in on the same machine and do:
```
screen -ls  # to find the number of the screen session
screen -r <screen_session_number>
```
Then either exit the screen session again (`Ctrl + A` followed by `D`) to continue to run the code, kill the code running (`Ctrl + C`) to run something different, or kill the screen session (`Ctrl +D`).


#### Tips for running at the PSI T3

At the PSI T3, to avoid the code to stop when closing the terminal, run the code with `nohup`:
```
nohup ./make_skims.sh > log.log 2>&1 &
```

## Checking skims completion

To check that no event is missing in the output skims compared to the input NTuples list, adapt and execute the bash script [skimmer/check_number_of_events.sh](https://github.com/fleble/SVJProcessing/blob/main/skimmer/check_number_of_events.sh).     
The code is designed to check if you have processed the entire input files list, not partially.


## ParticleNet Jet Tagger Score

By default, the skimming code does not inference on the ParticleNet jet tagger. To add ParticleNet jet tagger score to the skims, use `analysis_configs.t_channel_pre_selection.py` as the module and add the `-pn_tagger` flag to the skimming command. The inferencing will be done on the triton server, and the network scores will be saved to the skims. This part of the skimming code can be done using either TreeMaker NTuples or skims as input files. To run the code over skims, add the `-skim_source` flag.
