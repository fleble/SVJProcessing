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

# TODO: Add signals for all years
# TODO: Add all bkgs for all years

datasets_info = {
    "2016": {},
    "2017": {},
    "2018": {},
}

signal_models = [
    "t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1",
    "t-channel_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

qcd_bins = [
    "QCD_Pt_170to300",
    "QCD_Pt_300to470",
    "QCD_Pt_470to600",
    "QCD_Pt_600to800",
    "QCD_Pt_800to1000",
    "QCD_Pt_1000to1400",
    "QCD_Pt_1400to1800",
    "QCD_Pt_1800to2400",
    "QCD_Pt_2400to3200",
    "QCD_Pt_3200toInf",
]

ttjets_bins = [
    "TTJets",
    "TTJets_HT-600to800",
    "TTJets_HT-800to1200",
    "TTJets_HT-1200to2500",
    "TTJets_HT-2500toInf",
]


datasets_info["2018"].update({
    signal_model: [
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": f"/store/user/lpcdarkqcd/tchannel_UL/2018/Full_11142022/PrivateSamples/SVJ_UL2018_{signal_model}_13TeV-madgraphMLM-pythia8_n-1000",
            "regex": "",
        }
    ]
    for signal_model in signal_models
})

datasets_info["2018"].update({
    qcd_bin: [
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": f"/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Summer20UL18/{qcd_bin}_TuneCP5_13TeV_pythia8/",
            "regex": "",
        }
    ]
    for qcd_bin in qcd_bins
})

datasets_info["2018"].update({
    ttjets_bin: [
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": f"/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Summer20UL18/{ttjets_bin}_TuneCP5_13TeV-madgraphMLM-pythia8/",
            "regex": "",
        }
    ]
    for ttjets_bin in ttjets_bins
})

