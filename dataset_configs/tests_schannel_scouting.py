###################################  README  ###################################
#
# Is called "dataset" a set of files corresponding to the same physics process.
# The object `datasets_info` describes the location of the different datasets.
# Its structure is the following:
#    * Keys are year
#    * Values are the "datasets_per_year_info"
# The structure of the `datasets_per_year_info` is the following:
#    * Keys are dataset names
#    * Values are the "dataset_info" defining which files belong to the dataset
# 
# The "dataset_info" has the following structure. It is a list of dict, which
# has 2 keys:
#    * "redirector": The XRootD redirector to the remote storage element
#    * "path": The path to the directory at which the files are located
#    * "regex": The regex to apply to select some files from that directory. 
#               The regex must be "" if no regex is applied.
#
################################################################################


years = ["2018"]

datasets_info = {
    year: {} for year in years
}

signal_models = [
    "mMed-1500GeV_mDark-20GeV_rinv-0.3_alpha-peak_13TeV"
]

for year in years:
    datasets_info[year].update({
        signal_model: [
            {

                "redirector": "root://storage01.lcg.cscs.ch:1096//",
                "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/scouting/truth_study/SVJ_std2_UL2018_scouting_truth_study/SVJ_{signal_model}/",
                "regex": f"PFNano_s-channel_{signal_model.replace('GeV','')}-pythia8_n-1000_part-[1-9].root",

            },
        ]
        for signal_model in signal_models
    })
