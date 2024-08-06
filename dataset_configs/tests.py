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


years = ["2016", "2016APV", "2017", "2018"]

datasets_info = {
    year: {} for year in years
}

signal_models = [
    "t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

for year in years:
    datasets_info[year].update({
        signal_model: [
            {
                "redirector": "root://cmseos.fnal.gov/",
                "path": f"/store/user/lpcdarkqcd/tchannel_UL/{year}/Full/PrivateSamples/SVJ_UL{year}_{signal_model}_13TeV-madgraphMLM-pythia8_n-1000",
                "regex": "/[0-9]_RA2AnalysisTree.root",
            },
        ]
        for signal_model in signal_models
    })

