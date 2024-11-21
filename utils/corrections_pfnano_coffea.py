#import relevant libraries
import correctionlib
import os
import awkward as ak
import argparse

from coffea.lookup_tools import extractor
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory
from coffea.util import save

from datetime import datetime
import requests
from bs4 import BeautifulSoup

now = datetime.now()


#define path for data directory
path_current_file = os.path.dirname(os.path.abspath(__file__))
out_path_corrections = path_current_file.replace("utils","data")



def __get_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-jme', '--include_jme_corrections',
        help='Set True if input files are (PF)NanoAOD files', 
        default=False, 
        action='store_true',
    )

    parser.add_argument(
        '-all', '--include_all_corrections',
        help='Set True if input files are (PF)NanoAOD files', 
        default=False, 
        action='store_true',
    )
    

    return parser.parse_args()


def download_raw_files(url,output_dir):
    
    print(f"Downloading raw files from {url}...")

    # Send a GET request to the URL
    response = requests.get(url)
    
    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.text, 'html.parser')
    
    # Find all links on the page
    links = soup.find_all('a', href=True)

    # Create a directory to store the downloaded files
    os.makedirs(output_dir, exist_ok=True)
    
    # Iterate through the links and download raw files
    for link in links:
        href = link['href']
        if href.endswith('.txt'):
            file_name = href.split('/')[-1]
            file_path = os.path.join(output_dir, file_name)
            
            #check if file_path already exists
            if os.path.exists(file_path):
                print(f"{file_name} already exists in {output_dir}")
                continue

            # Construct the raw file URL
            raw_url = f"https://raw.githubusercontent.com{href.replace('blob','refs/heads')}"

            print(f"Downloading {file_name}...")
            os.system(f"wget {raw_url} -P {output_dir}")

            print(f"{file_name} downloaded successfully.")



def fetch_and_save_files_corrections(pog,observable,year=None):

    path_saved_corrections = ""
    
    if pog == 'jme':
        
        #path where to save the corrections
        path_saved_corrections = f"{out_path_corrections}/jme/"
        
        if observable == 'met':
            
            #check if file met.json.gz is already in the directory
            if os.path.exists(f"{path_saved_corrections}/met/{year}_UL/met.json.gz"):
                print(f"File met.json.gz already exists in {path_saved_corrections}/met/{year}_UL")
                return path_saved_corrections
                
            #get the file met.json.gz from https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME, save it in out_path_corrections/jme/year_UL
            os.system(f"wget https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/raw/master/POG/JME/{year}_UL/met.json.gz -P {path_saved_corrections}/met/{year}_UL")

        if observable == 'jet':

            download_raw_files(url="https://github.com/mcremone/decaf/tree/UL/analysis/data/jerc", output_dir=path_saved_corrections + "/jerc")


    return path_saved_corrections


def XY_MET_Correction(year, npv, run, pt, phi, isData):
    
    npv = ak.where((npv>200),ak.full_like(npv,200),npv)
    pt  = ak.where((pt>1000.),ak.full_like(pt,1000.),pt)

    if '2016preVFP' in year:
        #run = ak.where((run<271036),ak.full_like(run,271036),run)
        run = ak.where((run>278771),ak.full_like(run,278771),run)
    if '2016postVFP' in year:
        #run = ak.where((run<271036),ak.full_like(run,271036),run)
        run = ak.where((run>284045),ak.full_like(run,284045),run)
    if '2017' in year:
        #run = ak.where((run<294927),ak.full_like(run,294927),run)
        run = ak.where((run>306463),ak.full_like(run,306463),run)
    if '2018' in year:
        #run = ak.where((run<314472),ak.full_like(run,314472),run)
        run = ak.where((run>325274),ak.full_like(run,325274),run)
        
    
    path_saved_corrections = fetch_and_save_files_corrections(pog='jme',observable='met',year=year)
    evaluator = correctionlib.CorrectionSet.from_file(f'{path_saved_corrections}/met/met.json.gz')

    if isData:
        corrected_pt = evaluator['pt_metphicorr_pfmet_data'].evaluate(pt,phi,npv,run)
        corrected_phi = evaluator['phi_metphicorr_pfmet_data'].evaluate(pt,phi,npv,run)

    if not isData:
        corrected_pt = evaluator['pt_metphicorr_pfmet_mc'].evaluate(pt,phi,npv,run)
        corrected_phi = evaluator['phi_metphicorr_pfmet_mc'].evaluate(pt,phi,npv,run)

    return corrected_pt, corrected_phi



def jet_factory_factory(files,jec_name_map):
    ext = extractor()
    directory=f'{out_path_corrections}/jme/jerc'
    for filename in files:
        ext.add_weight_sets([f"* * {directory+'/'+filename}"])
    ext.finalize()
    jec_stack = JECStack(ext.make_evaluator())
    return CorrectedJetsFactory(jec_name_map, jec_stack)


def build_jet_and_corrections():

    #fetch the corrections for the jets
    fetch_and_save_files_corrections(pog='jme',observable='jet')

    jec_name_map = {
    'JetPt': 'pt',
    'JetMass': 'mass',
    'JetEta': 'eta',
    'JetA': 'area',
    'ptGenJet': 'pt_gen',
    'ptRaw': 'pt_raw',
    'massRaw': 'mass_raw',
    'Rho': 'event_rho',
    'METpt': 'pt',
    'METphi': 'phi',
    'JetPhi': 'phi',
    'UnClusteredEnergyDeltaX': 'MetUnclustEnUpDeltaX',
    'UnClusteredEnergyDeltaY': 'MetUnclustEnUpDeltaY',
    }


    jet_factory = {
        "2016preVFPmc": jet_factory_factory(
            files=[
                "Summer19UL16APV_V7_MC_L1FastJet_AK4PFchs.jec.txt",
                "Summer19UL16APV_V7_MC_L2Relative_AK4PFchs.jec.txt",
                "Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.junc.txt",
                "Summer19UL16APV_V7_MC_Uncertainty_AK4PFchs.junc.txt",
                "Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.jr.txt",
                "Summer20UL16APV_JRV3_MC_SF_AK4PFchs.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2016preVFPmcNOJER": jet_factory_factory(
            files=[
                "Summer19UL16APV_V7_MC_L1FastJet_AK4PFchs.jec.txt",
                "Summer19UL16APV_V7_MC_L2Relative_AK4PFchs.jec.txt",
                "Summer19UL16APV_V7_MC_Uncertainty_AK4PFchs.junc.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2016preVFPmcNOJEC": jet_factory_factory(
            files=[
                "Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.jr.txt",
                "Summer20UL16APV_JRV3_MC_SF_AK4PFchs.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2016postVFPmc": jet_factory_factory(
            files=[
                "Summer19UL16_V7_MC_L1FastJet_AK4PFchs.jec.txt",
                "Summer19UL16_V7_MC_L2Relative_AK4PFchs.jec.txt",
                "Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.junc.txt",
                "Summer19UL16_V7_MC_Uncertainty_AK4PFchs.junc.txt",
                "Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.jr.txt",
                "Summer20UL16_JRV3_MC_SF_AK4PFchs.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2016postVFPmcNOJER": jet_factory_factory(
            files=[
                "Summer19UL16_V7_MC_L1FastJet_AK4PFchs.jec.txt",
                "Summer19UL16_V7_MC_L2Relative_AK4PFchs.jec.txt",
                "Summer19UL16_V7_MC_Uncertainty_AK4PFchs.junc.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2016postVFPmcNOJEC": jet_factory_factory(
            files=[
                "Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.jr.txt",
                "Summer20UL16_JRV3_MC_SF_AK4PFchs.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2017mc": jet_factory_factory(
            files=[
                "Summer19UL17_V5_MC_L1FastJet_AK4PFchs.jec.txt",
                "Summer19UL17_V5_MC_L2Relative_AK4PFchs.jec.txt",
                "Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.junc.txt",
                "Summer19UL17_V5_MC_Uncertainty_AK4PFchs.junc.txt",
                "Summer19UL17_JRV3_MC_PtResolution_AK4PFchs.jr.txt",
                "Summer19UL17_JRV3_MC_SF_AK4PFchs.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2017mcNOJER": jet_factory_factory(
            files=[
                "Summer19UL17_V5_MC_L1FastJet_AK4PFchs.jec.txt",
                "Summer19UL17_V5_MC_L2Relative_AK4PFchs.jec.txt",
                "Summer19UL17_V5_MC_Uncertainty_AK4PFchs.junc.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2017mcNOJEC": jet_factory_factory(
            files=[
                "Summer19UL17_JRV3_MC_PtResolution_AK4PFchs.jr.txt",
                "Summer19UL17_JRV3_MC_SF_AK4PFchs.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2018mc": jet_factory_factory(
            files=[
                "Summer19UL18_V5_MC_L1FastJet_AK4PFchs.jec.txt",
                "Summer19UL18_V5_MC_L2Relative_AK4PFchs.jec.txt",
                "Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.junc.txt",
                "Summer19UL18_V5_MC_Uncertainty_AK4PFchs.junc.txt",
                "Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.jr.txt",
                "Summer19UL18_JRV2_MC_SF_AK4PFchs.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2018mcNOJER": jet_factory_factory(
            files=[
                "Summer19UL18_V5_MC_L1FastJet_AK4PFchs.jec.txt",
                "Summer19UL18_V5_MC_L2Relative_AK4PFchs.jec.txt",
                "Summer19UL18_V5_MC_Uncertainty_AK4PFchs.junc.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2018mcNOJEC": jet_factory_factory(
            files=[
                "Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.jr.txt",
                "Summer19UL18_JRV2_MC_SF_AK4PFchs.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
    }

    
    fatjet_factory = {
        "2016preVFPmc": jet_factory_factory(
            files=[
                "Summer19UL16APV_V7_MC_L1FastJet_AK8PFPuppi.jec.txt",
                "Summer19UL16APV_V7_MC_L2Relative_AK8PFPuppi.jec.txt",
                "Summer19UL16APV_V7_MC_UncertaintySources_AK8PFPuppi.junc.txt",
                "Summer19UL16APV_V7_MC_Uncertainty_AK8PFPuppi.junc.txt",
                "Summer20UL16APV_JRV3_MC_PtResolution_AK8PFPuppi.jr.txt",
                "Summer20UL16APV_JRV3_MC_SF_AK8PFPuppi.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2016preVFPmcNOJER": jet_factory_factory(
            files=[
                "Summer19UL16APV_V7_MC_L1FastJet_AK8PFPuppi.jec.txt",
                "Summer19UL16APV_V7_MC_L2Relative_AK8PFPuppi.jec.txt",
                "Summer19UL16APV_V7_MC_Uncertainty_AK8PFPuppi.junc.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2016preVFPmcNOJEC": jet_factory_factory(
            files=[
                "Summer20UL16APV_JRV3_MC_PtResolution_AK8PFPuppi.jr.txt",
                "Summer20UL16APV_JRV3_MC_SF_AK8PFPuppi.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2016postVFPmc": jet_factory_factory(
            files=[
                "Summer19UL16_V7_MC_L1FastJet_AK8PFPuppi.jec.txt",
                "Summer19UL16_V7_MC_L2Relative_AK8PFPuppi.jec.txt",
                "Summer19UL16_V7_MC_UncertaintySources_AK8PFPuppi.junc.txt",
                "Summer19UL16_V7_MC_Uncertainty_AK8PFPuppi.junc.txt",
                "Summer20UL16_JRV3_MC_PtResolution_AK8PFPuppi.jr.txt",
                "Summer20UL16_JRV3_MC_SF_AK8PFPuppi.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2016postVFPmcNOJER": jet_factory_factory(
            files=[
                "Summer19UL16_V7_MC_L1FastJet_AK8PFPuppi.jec.txt",
                "Summer19UL16_V7_MC_L2Relative_AK8PFPuppi.jec.txt",
                "Summer19UL16_V7_MC_Uncertainty_AK8PFPuppi.junc.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2016postVFPmcNOJEC": jet_factory_factory(
            files=[
                "Summer20UL16_JRV3_MC_PtResolution_AK8PFPuppi.jr.txt",
                "Summer20UL16_JRV3_MC_SF_AK8PFPuppi.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2017mc": jet_factory_factory(
            files=[
                "Summer19UL17_V5_MC_L1FastJet_AK8PFPuppi.jec.txt",
                "Summer19UL17_V5_MC_L2Relative_AK8PFPuppi.jec.txt",
                "Summer19UL17_V5_MC_UncertaintySources_AK8PFPuppi.junc.txt",
                "Summer19UL17_V5_MC_Uncertainty_AK8PFPuppi.junc.txt",
                "Summer19UL17_JRV3_MC_PtResolution_AK8PFPuppi.jr.txt",
                "Summer19UL17_JRV3_MC_SF_AK8PFPuppi.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2017mcNOJER": jet_factory_factory(
            files=[
                "Summer19UL17_V5_MC_L1FastJet_AK8PFPuppi.jec.txt",
                "Summer19UL17_V5_MC_L2Relative_AK8PFPuppi.jec.txt",
                "Summer19UL17_V5_MC_Uncertainty_AK8PFPuppi.junc.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2017mcNOJEC": jet_factory_factory(
            files=[
                "Summer19UL17_JRV3_MC_PtResolution_AK8PFPuppi.jr.txt",
                "Summer19UL17_JRV3_MC_SF_AK8PFPuppi.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2018mc": jet_factory_factory(
            files=[
                "Summer19UL18_V5_MC_L1FastJet_AK8PFPuppi.jec.txt",
                "Summer19UL18_V5_MC_L2Relative_AK8PFPuppi.jec.txt",
                "Summer19UL18_V5_MC_UncertaintySources_AK8PFPuppi.junc.txt",
                "Summer19UL18_V5_MC_Uncertainty_AK8PFPuppi.junc.txt",
                "Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi.jr.txt",
                "Summer19UL18_JRV2_MC_SF_AK8PFPuppi.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2018mcNOJER": jet_factory_factory(
            files=[
                "Summer19UL18_V5_MC_L1FastJet_AK8PFPuppi.jec.txt",
                "Summer19UL18_V5_MC_L2Relative_AK8PFPuppi.jec.txt",
                "Summer19UL18_V5_MC_Uncertainty_AK8PFPuppi.junc.txt",
            ],
            jec_name_map=jec_name_map,
        ),
        "2018mcNOJEC": jet_factory_factory(
            files=[
                "Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi.jr.txt",
                "Summer19UL18_JRV2_MC_SF_AK8PFPuppi.jersf.txt",
            ],
            jec_name_map=jec_name_map,
        ),
    }
    
    met_factory = CorrectedMETFactory(jec_name_map)

    return jet_factory, fatjet_factory, met_factory


def build_corrections(corrections_to_include_dict):

    corrections = {}
    corrections = {
        'jet_factory':              corrections_to_include_dict["jet_factory"],
        'fatjet_factory':           corrections_to_include_dict["fatjet_factory"],
        'met_factory':              corrections_to_include_dict["met_factory"],
        #'get_met_xy_correction':    corrections_to_include_dict["met_xy_correction"],
    }

    return corrections

def main():

    args = __get_arguments()

    #get actual date and time
    date_time = now.strftime("%Y-%m-%d_%H-%M-%S")

    #build corrections based on the arguments, if include_jme_corrections is True
    corrections_to_include_dict = {}

    print("Building corrections...")

    if args.include_jme_corrections or args.include_all_corrections:
        print("Including JME corrections...")

        jet_factory, fatjet_factory, met_factory = build_jet_and_corrections()

        corrections_to_include_dict["jet_factory"] = jet_factory
        corrections_to_include_dict["fatjet_factory"] = fatjet_factory
        corrections_to_include_dict["met_factory"] = met_factory
        #corrections_to_include_dict["met_xy_correction"] = XY_MET_Correction    


    corrections = build_corrections(corrections_to_include_dict)

    print("Corrections built successfully !")

    #build a string with the arguments
    arguments = ""

    if args.include_all_corrections:
        arguments += f"_all_corr"
        
    if args.include_jme_corrections and not args.include_all_corrections:
            arguments += f"_jme_corr"

    print(f"Saving corrections to {out_path_corrections}/corrections_{date_time}{arguments}.coffea")
    save(corrections, f'{out_path_corrections}/corrections_{date_time}{arguments}.coffea')
    print(f"Corrections saved successfully.")

if __name__ == "__main__":

    main()