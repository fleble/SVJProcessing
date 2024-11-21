from coffea.jetmet_tools.JetResolution import JetResolution
from coffea.lookup_tools import extractor
import os

#define path for data directory
path_current_file = os.path.dirname(os.path.abspath(__file__))
out_path_corrections = path_current_file.replace("utils","data")


jer_names_dict_pt = { 

    "2016preVFPmc": {
        "Jet": "Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs",
        "FatJet": "Summer20UL16APV_JRV3_MC_PtResolution_AK8PFPuppi",
    },
    "2016postVFPmc": {
        "Jet": "Summer20UL16_JRV3_MC_PtResolution_AK4PFchs",
        "FatJet": "Summer20UL16_JRV3_MC_PtResolution_AK8PFPuppi",
    },
    "2017mc": {
        "Jet": "Summer19UL17_JRV3_MC_PtResolution_AK4PFchs",
        "FatJet": "Summer19UL17_JRV3_MC_PtResolution_AK8PFPuppi",
    },
    "2018mc": {
        "Jet": "Summer19UL18_JRV2_MC_PtResolution_AK4PFchs",
        "FatJet": "Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi",
    },

}

jer_names_dict_eta = { 

    "2016preVFPmc": {
        "Jet": "Summer20UL16APV_JRV3_MC_EtaResolution_AK4PFchs",
        "FatJet": "Summer20UL16APV_JRV3_MC_EtaResolution_AK8PFPuppi",
    },
    "2016postVFPmc": {
        "Jet": "Summer20UL16_JRV3_MC_EtaResolution_AK4PFchs",
        "FatJet": "Summer20UL16_JRV3_MC_EtaResolution_AK8PFPuppi",
    },
    "2017mc": {
        "Jet": "Summer19UL17_JRV3_MC_EtaResolution_AK4PFchs",
        "FatJet": "Summer19UL17_JRV3_MC_EtaResolution_AK8PFPuppi",
    },
    "2018mc": {
        "Jet": "Summer19UL18_JRV2_MC_EtaResolution_AK4PFchs",
        "FatJet": "Summer19UL18_JRV2_MC_EtaResolution_AK8PFPuppi",
    },

}

jer_names_dict_phi = { 

    "2016preVFPmc": {
        "Jet": "Summer20UL16APV_JRV3_MC_PhiResolution_AK4PFchs",
        "FatJet": "Summer20UL16APV_JRV3_MC_PhiResolution_AK8PFPuppi",
    },
    "2016postVFPmc": {
        "Jet": "Summer20UL16_JRV3_MC_PhiResolution_AK4PFchs",
        "FatJet": "Summer20UL16_JRV3_MC_PhiResolution_AK8PFPuppi",
    },
    "2017mc": {
        "Jet": "Summer19UL17_JRV3_MC_PhiResolution_AK4PFchs",
        "FatJet": "Summer19UL17_JRV3_MC_PhiResolution_AK8PFPuppi",
    },
    "2018mc": {
        "Jet": "Summer19UL18_JRV2_MC_PhiResolution_AK4PFchs",
        "FatJet": "Summer19UL18_JRV2_MC_PhiResolution_AK8PFPuppi",
    },

}


#gen matching criteria based on: https://gitlab.cern.ch/cms-analysis/general/CMSJMECalculators/-/blob/main/src/JMESystematicsCalculators.cc?ref_type=heads#L220 

def jetmet_evaluator(files):
    extract = extractor()
    directory=f'{out_path_corrections}/jme/jerc'
    for filename in files:
        extract.add_weight_sets([f"* * {directory+'/'+filename}"])
    extract.finalize()

    return extract.make_evaluator()


def fetch_resolution(correction_key,jet_coll, var):
    if var == "pt":
        jer_names = [jer_names_dict_pt[correction_key][jet_coll]]
        files = [f"{jer_names[0]}.jr.txt"]
    elif var == "eta":
        jer_names = [jer_names_dict_eta[correction_key][jet_coll]]
        files = [f"{jer_names[0]}.txt"]
    elif var == "phi":
        jer_names = [jer_names_dict_phi[correction_key][jet_coll]]
        files = [f"{jer_names[0]}.txt"]
    
    
    evaluator = jetmet_evaluator(files)
    reso_obj = JetResolution(**{name: evaluator[name] for name in jer_names})
    return reso_obj


def check_pt_resolution(reco_pt, gen_pt, resolution):
    m_genMatch_dPtmax = 0.3
    return abs(gen_pt - reco_pt) < m_genMatch_dPtmax * resolution * reco_pt


def get_jets_resolution(reco_pt, reco_eta, reco_rho, correction_key,jet_coll,var):
    reso_obj = fetch_resolution(correction_key,jet_coll,var)
    resolution = reso_obj.getResolution(JetEta=reco_eta, Rho=reco_rho, JetPt=reco_pt)
    return resolution