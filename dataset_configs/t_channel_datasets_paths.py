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
    "TTJets_SingleLeptFromT",
    "TTJets_SingleLeptFromTbar",
    "TTJets_DiLept",
    "TTJets_SingleLeptFromT_genMET-150", 
    "TTJets_SingleLeptFromTbar_genMET-150", 
    "TTJets_DiLept_genMET-150",
]

wjets_bins = [
    "WJetsToLNu",
    "WJetsToLNu_HT-70To100",
    "WJetsToLNu_HT-70To100",
    "WJetsToLNu_HT-100To200",
    "WJetsToLNu_HT-200To400",
    "WJetsToLNu_HT-400To600",
    "WJetsToLNu_HT-600To800",
    "WJetsToLNu_HT-800To1200",
    "WJetsToLNu_HT-1200To2500",
    "WJetsToLNu_HT-2500ToInf",
    "WJetsToQQ_HT-200to400",
    "WJetsToQQ_HT-400to600",
    "WJetsToQQ_HT-600to800",
    "WJetsToQQ_HT-800toInf",
]

zjets_bins = [
    "ZJetsToNuNu_HT-100To200",
    "ZJetsToNuNu_HT-200To400",
    "ZJetsToNuNu_HT-400To600",
    "ZJetsToNuNu_HT-600To800",
    "ZJetsToNuNu_HT-800To1200",
    "ZJetsToNuNu_HT-1200To2500",
    "ZJetsToNuNu_HT-2500ToInf",
]

single_top_bins = [
    "ST_s-channel_4f_hadronicDecays",
    "ST_s-channel_4f_leptonDecays",
    "ST_t-channel_antitop_5f_InclusiveDecays",
    "ST_t-channel_top_5f_InclusiveDecays",
    "ST_tW_top_5f_inclusiveDecays",
    "ST_tW_antitop_5f_inclusiveDecays",
    "tZq_ll_4f_ckm_NLO",
]

diboson_bins = [
    "ZZTo2Q2Nu",
    "WZTo2Q2Nu",
    "WZTo1L1Nu2Q",
    "WWTo1L1Nu2Q",
]

background_bins = (
    qcd_bins
    + ttjets_bins
    + wjets_bins
    + zjets_bins
    + single_top_bins
    + diboson_bins
)


for year in years:
    datasets_info[year].update({
        signal_model: [
            {
                "redirector": "root://cmseos.fnal.gov/",
                "path": f"/store/user/lpcdarkqcd/tchannel_UL/{year}/Full/PrivateSamples/SVJ_UL{year}_{signal_model}_13TeV-madgraphMLM-pythia8_n-1000",
                "regex": "",
            },
        ]
        for signal_model in signal_models
    })


for year in years:
    year_tag = year[2:]

    for bin in background_bins:
        if "QCD_" in bin:
              suffix = "TuneCP5_13TeV_pythia8"
        elif "TTJets_" in bin or "WJets" in bin or "ZJets" in bin:
            suffix = "TuneCP5_13TeV-madgraphMLM-pythia8"
        elif "ST_s-channel" in bin or "tZq_" in bin:
            suffix = "TuneCP5_13TeV-amcatnlo-pythia8"
        elif "ST_t-channel" in bin or "ST_tW" in bin:
            suffix = "TuneCP5_13TeV-powheg-pythia8"
        elif "WWTo" in bin or "ZZTo" in bin or "WZTo":
            suffix = "TuneCP5_13TeV-amcatnloFXFX-pythia8"
        else:
            print(f"Unknown background {bin}")
            exit(1)
    
        datasets_info[year].update({
            bin: [
                {
                    "redirector": "root://cmseos.fnal.gov/",
                    "path": f"/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Summer20UL{year_tag}/{bin}_{suffix}/",
                    "regex": "",
                }
            ]
        })

# Data
# EGamma 2018
datasets_info["2018"].update({
    "EGamma": [
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": f"/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018A-UL2018-v1/EGamma/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": f"/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018B-UL2018-v1/EGamma/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": f"/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018C-UL2018-v1/EGamma/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": f"/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018D-UL2018-v2/EGamma/",
            "regex": "",
        },
    ]
})

# SingleMuon 2018
datasets_info["2018"].update({
    "SingleMuon": [
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": f"/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018A-UL2018-v3/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": f"/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018B-UL2018-v2/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": f"/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018C-UL2018-v2/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": f"/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018D-UL2018-v3/SingleMuon/",
            "regex": "",
        },
    ]
})



# To check the content of the dataset config dict
# import json
# print(json.dumps(datasets_info["2018"], indent=4))

