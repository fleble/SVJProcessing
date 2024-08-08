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

#mMed-5000GeV_mDark-32GeV_rinv-0.3
signal_models_svjl = [
    #SVJL signal models
    "mMed-1500GeV_mDark-8GeV_rinv-0.3",
    "mMed-1500GeV_mDark-16GeV_rinv-0.3",
    "mMed-1500GeV_mDark-32GeV_rinv-0.3",
    #"mMed-1500GeV_mDark-64GeV_rinv-0.3",
    "mMed-1500GeV_mDark-8GeV_rinv-0.5",
    "mMed-1500GeV_mDark-16GeV_rinv-0.5",
    "mMed-1500GeV_mDark-32GeV_rinv-0.5",
    #"mMed-1500GeVmDark-64GeV_rinv-0.5",
    "mMed-1500GeV_mDark-8GeV_rinv-0.7",
    "mMed-1500GeV_mDark-16GeV_rinv-0.7",
    "mMed-1500GeV_mDark-32GeV_rinv-0.7",
    #"mMed-1500GeV_mDark-64GeV_rinv-0.7",
    #repeat for 2000,3000,4000,5000
    "mMed-2000GeV_mDark-8GeV_rinv-0.3",
    "mMed-2000GeV_mDark-16GeV_rinv-0.3",
    "mMed-2000GeV_mDark-32GeV_rinv-0.3",
    #"mMed-2000GeV_mDark-64GeV_rinv-0.3",
    "mMed-2000GeV_mDark-8GeV_rinv-0.5",
    "mMed-2000GeV_mDark-16GeV_rinv-0.5",
    "mMed-2000GeV_mDark-32GeV_rinv-0.5",
    #"mMed-2000GeV_mDark-64GeV_rinv-0.5",
    "mMed-2000GeV_mDark-8GeV_rinv-0.7",
    "mMed-2000GeV_mDark-16GeV_rinv-0.7",
    "mMed-2000GeV_mDark-32GeV_rinv-0.7",
    #"mMed-2000GeV_mDark-64GeV_rinv-0.7",
    #3000
    "mMed-3000GeV_mDark-8GeV_rinv-0.3",
    "mMed-3000GeV_mDark-16GeV_rinv-0.3",
    "mMed-3000GeV_mDark-32GeV_rinv-0.3",
    #"mMed-3000GeV_mDark-64GeV_rinv-0.3",
    "mMed-3000GeV_mDark-8GeV_rinv-0.5",
    "mMed-3000GeV_mDark-16GeV_rinv-0.5",
    "mMed-3000GeV_mDark-32GeV_rinv-0.5",
    #"mMed-3000GeV_mDark-64GeV_rinv-0.5",
    "mMed-3000GeV_mDark-8GeV_rinv-0.7",
    "mMed-3000GeV_mDark-16GeV_rinv-0.7",
    "mMed-3000GeV_mDark-32GeV_rinv-0.7",
    #"mMed-3000GeV_mDark-64GeV_rinv-0.7",
    #4000
    "mMed-4000GeV_mDark-8GeV_rinv-0.3",
    "mMed-4000GeV_mDark-16GeV_rinv-0.3",
    "mMed-4000GeV_mDark-32GeV_rinv-0.3",
    #"mMed-4000GeV_mDark-64GeV_rinv-0.3",
    "mMed-4000GeV_mDark-8GeV_rinv-0.5",
    "mMed-4000GeV_mDark-16GeV_rinv-0.5",
    "mMed-4000GeV_mDark-32GeV_rinv-0.5",
    #"mMed-4000GeV_mDark-64GeV_rinv-0.5",
    "mMed-4000GeV_mDark-8GeV_rinv-0.7",
    "mMed-4000GeV_mDark-16GeV_rinv-0.7",
    "mMed-4000GeV_mDark-32GeV_rinv-0.7",
    #"mMed-4000GeV_mDark-64GeV_rinv-0.7",
    #5000
    "mMed-5000GeV_mDark-8GeV_rinv-0.3",
    "mMed-5000GeV_mDark-16GeV_rinv-0.3",
    "mMed-5000GeV_mDark-32GeV_rinv-0.3",
    #"mMed-5000GeV_mDark-64GeV_rinv-0.3",
    "mMed-5000GeV_mDark-8GeV_rinv-0.5",
    "mMed-5000GeV_mDark-16GeV_rinv-0.5",
    "mMed-5000GeV_mDark-32GeV_rinv-0.5",
    #"mMed-5000GeV_mDark-64GeV_rinv-0.5",
    "mMed-5000GeV_mDark-8GeV_rinv-0.7",
    "mMed-5000GeV_mDark-16GeV_rinv-0.7",
    "mMed-5000GeV_mDark-32GeV_rinv-0.7",
    #"mMed-5000GeV_mDark-64GeV_rinv-0.7",
]

signal_models_svjtau = [
  #SVJtau signal models
  "mMed-1500GeV_mDark-6.4GeV_rinv-0.3_alpha-0.3",
  "mMed-1500GeV_mDark-6.4GeV_rinv-0.3_alpha-0.5",
  "mMed-1500GeV_mDark-6.4GeV_rinv-0.3_alpha-0.7",
  "mMed-1500GeV_mDark-8GeV_rinv-0.3_alpha-0.3",
  "mMed-1500GeV_mDark-8GeV_rinv-0.3_alpha-0.5",
  "mMed-1500GeV_mDark-8GeV_rinv-0.3_alpha-0.7",

  "mMed-2000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.3",
  "mMed-2000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.5",
  "mMed-2000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.7",
  "mMed-2000GeV_mDark-8GeV_rinv-0.3_alpha-0.3",
  "mMed-2000GeV_mDark-8GeV_rinv-0.3_alpha-0.5",
  "mMed-2000GeV_mDark-8GeV_rinv-0.3_alpha-0.7",

  "mMed-3000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.3",
  "mMed-3000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.5",
  "mMed-3000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.7",
  "mMed-3000GeV_mDark-8GeV_rinv-0.3_alpha-0.3",
  "mMed-3000GeV_mDark-8GeV_rinv-0.3_alpha-0.5",
  "mMed-3000GeV_mDark-8GeV_rinv-0.3_alpha-0.7",
  
  "mMed-4000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.3",
  "mMed-4000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.5",
  "mMed-4000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.7",
  "mMed-4000GeV_mDark-8GeV_rinv-0.3_alpha-0.3",
  "mMed-4000GeV_mDark-8GeV_rinv-0.3_alpha-0.5",
  "mMed-4000GeV_mDark-8GeV_rinv-0.3_alpha-0.7",

  "mMed-5000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.3",
  "mMed-5000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.5",
  "mMed-5000GeV_mDark-6.4GeV_rinv-0.3_alpha-0.7",
  "mMed-5000GeV_mDark-8GeV_rinv-0.3_alpha-0.3",
  "mMed-5000GeV_mDark-8GeV_rinv-0.3_alpha-0.5",
  "mMed-5000GeV_mDark-8GeV_rinv-0.3_alpha-0.7",  


]


qcd_bins = [
    #"QCD_Pt_170to300",
    #"QCD_Pt_300to470",
    "qcd_pt-470to600GeV",
    "qcd_pt-600to800GeV",
    "qcd_pt-800to1000GeV",
    "qcd_pt-1000to1400GeV",
    "qcd_pt-1400to1800GeV",
    "qcd_pt-1800to2400GeV",
    "qcd_pt-2400to3200GeV",
    "qcd_pt-3200toInfGeV",
]

ttjets_bins = [
    "TTJets_ht-600to800GeV",
    "TTJets_ht-800to1200GeV",
    "TTJets_ht-1200to2500GeV",
    "TTJets_ht-2500toInfGeV",
]

wjets_bins = [
    "WJetsToLNu_ht-100to200GeV",
    "WJetsToLNu_ht-200to400GeV",
    "WJetsToLNu_ht-400to600GeV",
    "WJetsToLNu_ht-600to800GeV",
    "WJetsToLNu_ht-800to1200GeV",
    "WJetsToLNu_ht-1200to2500GeV",
    "WJetsToLNu_ht-2500toInfGeV",
]

zjets_bins = [
    "ZJetsToNuNu_ht-100to200GeV",
    "ZJetsToNuNu_ht-200to400GeV",
    "ZJetsToNuNu_ht-400to600GeV",
    "ZJetsToNuNu_ht-600to800GeV",
    "ZJetsToNuNu_ht-800to1200GeV",
    "ZJetsToNuNu_ht-1200to2500GeV",
    "ZJetsToNuNu_ht-2500toInfGeV",
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
    #+ single_top_bins
    #+ diboson_bins
)

#SVJL model
for year in years:
    datasets_info[year].update({
        signal_model: [
            {
                "redirector": "root://t3se01.psi.ch:1094//",
                "path": f"/store/t3groups/ethz-susy/darkshowers/samples/signal/PFNanoUL/UL2018/SVJleptons_mDark_scan/PFNanoAOD/",
                "regex": f"PFNanoAOD_SVJL_{signal_model}",
            },
        ]
        for signal_model in signal_models_svjl
    })

#SVJtau model
for year in years:
    datasets_info[year].update({
        signal_model: [
            {
                "redirector": "root://t3se01.psi.ch:1094//",
                "path": f"/store/t3groups/ethz-susy/darkshowers/samples/signal/PFNanoUL/UL2018/SVJtaus/PFNanoAOD/",
                "regex": f"PFNanoAOD_SVJtaus_{signal_model}",
            },
        ]
        for signal_model in signal_models_svjtau
    })


for year in years:
    year_tag = year[2:]

    for bin in background_bins:
        if "qcd_" in bin:
              suffix = ""
        elif "TTJets_" in bin or "WJetsToLNu_" in bin or "ZJets" in bin:
            suffix = ""
        elif "ST_s-channel" in bin or "tZq_" in bin:
            suffix = "TuneCP5_13TeV-amcatnlo-pythia8"
        elif "ST_t-channel" in bin or "ST_tW" in bin:
            suffix = "TuneCP5_13TeV-powheg-pythia8"
        elif "WWTo" in bin or "ZZTo" in bin or "WZTo":
            suffix = "TuneCP5_13TeV-amcatnloFXFX-pythia8"
        else:
            print(f"Unknown background {bin}")
            exit(1)
    
        if "qcd_" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://t3se01.psi.ch:1094//",
                        "path": f"/store/user/jniedzie/svjets/cmssw/qcd/PFNano_preselectionBranches/{bin}/",
                        "regex": "",
                    }
                ]
            })

        
        elif "TTJets" in bin and "600to800GeV" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://storage01.lcg.cscs.ch:1096/",
                        "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/fleble/SVJ/TTJets/samples/UL2018/main_chain/PFNANOAOD/",
                        "regex": "PFNANOAOD_TTJets_ht-600to800GeV",
                    }
                ]
            })

        elif "TTJets" in bin and "800to1200GeV" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://storage01.lcg.cscs.ch:1096/",
                        "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/fleble/SVJ/TTJets/samples/UL2018/main_chain/PFNANOAOD/",
                        "regex": "PFNANOAOD_TTJets_ht-800to1200GeV",
                    }
                ]
            })

        elif "TTJets" in bin and "1200to2500GeV" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://storage01.lcg.cscs.ch:1096/",
                        "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/fleble/SVJ/TTJets/samples/UL2018/main_chain/PFNANOAOD/",
                        "regex": "PFNANOAOD_TTJets_ht-1200to2500GeV",
                    }
                ]
            })

        elif "TTJets" in bin and "2500toInfGeV" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://storage01.lcg.cscs.ch:1096/",
                        "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/fleble/SVJ/TTJets/samples/UL2018/main_chain/PFNANOAOD/",
                        "regex": "PFNANOAOD_TTJets_ht-2500toInfGeV",
                    }
                ]
            })

        elif "WJetsToLNu_ht-100to200GeV" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://t3se01.psi.ch:1094//",
                        "path": f"/store/t3groups/ethz-susy/darkshowers/samples/backgrounds/PFNano_cms/UL2018/wjets_new/PFNanoAOD/{bin}/",
                        "regex": "",
                    }
                ]
            })

        elif "WJetsToLNu_ht-200to400GeV" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://t3se01.psi.ch:1094//",
                        "path": f"/store/t3groups/ethz-susy/darkshowers/samples/backgrounds/PFNano_cms/UL2018/wjets_new/PFNanoAOD/{bin}/",
                        "regex": "",
                    }
                ]
            })

        elif "WJetsToLNu_ht-400to600GeV" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://t3se01.psi.ch:1094//",
                        "path": f"/store/t3groups/ethz-susy/darkshowers/samples/backgrounds/PFNano_cms/UL2018/wjets_new/PFNanoAOD/{bin}/",
                        "regex": "",
                    }
                ]
            })

        elif "WJetsToLNu_ht-600to800GeV" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://t3se01.psi.ch:1094//",
                        "path": f"/store/t3groups/ethz-susy/darkshowers/samples/backgrounds/PFNano_cms/UL2018/wjets_new/PFNanoAOD/{bin}/",
                        "regex": "",
                    }
                ]
            })

        elif "WJetsToLNu_ht-800to1200GeV" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://t3se01.psi.ch:1094//",
                        "path": f"/store/t3groups/ethz-susy/darkshowers/samples/backgrounds/PFNano_cms/UL2018/wjets_new/PFNanoAOD/{bin}/",
                        "regex": "",
                    }
                ]
            })

        elif "WJetsToLNu_ht-1200to2500GeV" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://t3se01.psi.ch:1094//",
                        "path": f"/store/t3groups/ethz-susy/darkshowers/samples/backgrounds/PFNano_cms/UL2018/wjets_new/PFNanoAOD/{bin}/",
                        "regex": "",
                    }
                ]
            })

        elif "WJetsToLNu_ht-2500toInfGeV" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://t3se01.psi.ch:1094//",
                        "path": f"/store/t3groups/ethz-susy/darkshowers/samples/backgrounds/PFNano_cms/UL2018/wjets_new/PFNanoAOD/{bin}",
                        "regex": "",
                    }
                ]
            })

        elif "ZJetsToNuNu_" in bin:
            datasets_info[year].update({
                bin: [
                    {
                        "redirector": "root://t3se01.psi.ch:1094//",
                        "path": f"/store/t3groups/ethz-susy/darkshowers/samples/backgrounds/PFNano_cms/UL2018/zjets_new/PFNanoAOD/{bin}/",
                        "regex": "",
                    }
                ]
            })

       

# Data
# SingleElectron 2016
datasets_info["2016"].update({
    "SingleElectron": [
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016B-UL2016_HIPM-ver2-v2/SingleElectron/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016C-UL2016_HIPM-v2/SingleElectron/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016D-UL2016_HIPM-v2/SingleElectron/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016E-UL2016_HIPM-v5/SingleElectron/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016F-UL2016_HIPM-v2/SingleElectron/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016G-UL2016-v2/SingleElectron/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016H-UL2016-v2/SingleElectron/",
            "regex": "",
        },
    ]
})


# SingleElectron 2017
datasets_info["2017"].update({
    "SingleElectron": [
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2017B-UL2017-v1/SingleElectron/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2017C-UL2017-v1/SingleElectron/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2017D-UL2017-v1/SingleElectron/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2017E-UL2017-v1/SingleElectron/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2017F-UL2017-v1/SingleElectron/",
            "regex": "",
        },
    ]
})



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


# SingleMuon 2016
datasets_info["2016"].update({
    "SingleMuon": [
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016B-UL2016_HIPM-ver2-v2/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016C-UL2016_HIPM-v2/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016D-UL2016_HIPM-v2/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016E-UL2016_HIPM-v2/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016F-UL2016_HIPM-v2/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016G-UL2016-v2/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2016H-UL2016-v2/SingleMuon/",
            "regex": "",
        },
    ]
})


# SingleMuon 2017
datasets_info["2017"].update({
    "SingleMuon": [
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2017B-UL2017-v1/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2017C-UL2017-v1/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2017D-UL2017-v1/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2017E-UL2017-v1/SingleMuon/",
            "regex": "",
        },
        {
            "redirector": "root://cmseos.fnal.gov/",
            "path": "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2017F-UL2017-v1/SingleMuon/",
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


# JetHT 2018
datasets_info["2018"].update({
    "JetHT": [
        ##########Run2018A
        {
            "redirector": "root://storage01.lcg.cscs.ch:1096/",
            "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/kadatta/PFNano/106x_v02/JetHT/Run2018A-UL2018_MiniAODv2-v1_PFNanov2pt2/221108_162639/0000/",
            "regex": "",
        },
        {
            "redirector": "root://storage01.lcg.cscs.ch:1096/",
            "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/kadatta/PFNano/106x_v02/JetHT/Run2018A-UL2018_MiniAODv2-v1_PFNanov2pt2/221108_162639/0001/",
            "regex": "",
        },
        ##########Run2018B
        {
            "redirector": "root://storage01.lcg.cscs.ch:1096/",
            "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/kadatta/PFNano/106x_v02/JetHT/Run2018B-UL2018_MiniAODv2-v1_PFNanov2pt2/221108_162412/0000/",
            "regex": "",
        },
        ##########Run2018C
        {
            "redirector": "root://storage01.lcg.cscs.ch:1096/",
            "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/kadatta/PFNano/106x_v02/JetHT/Run2018C-UL2018_MiniAODv2-v1_PFNanov2pt2/221108_163047/0000/",
            "regex": "",
        },
        ##########Run2018D
        {
            "redirector": "root://storage01.lcg.cscs.ch:1096/",
            "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/kadatta/PFNano/106x_v02/JetHT/Run2018D-UL2018_MiniAODv2-v2_PFNanov2pt2/221108_163418/0000/",
            "regex": "",
        },
        {
            "redirector": "root://storage01.lcg.cscs.ch:1096/",
            "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/kadatta/PFNano/106x_v02/JetHT/Run2018D-UL2018_MiniAODv2-v2_PFNanov2pt2/221108_163418/0001/",
            "regex": "",
        },
        {
            "redirector": "root://storage01.lcg.cscs.ch:1096/",
            "path": f"/pnfs/lcg.cscs.ch/cms/trivcat/store/user/kadatta/PFNano/106x_v02/JetHT/Run2018D-UL2018_MiniAODv2-v2_PFNanov2pt2/221108_163418/0002/",
            "regex": "",
        },
    ]
})



# To check the content of the dataset config dict
# import json
# print(json.dumps(datasets_info["2018"], indent=4))

