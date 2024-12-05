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
        '-pu', '--include_pu_corrections',
        help='Set True if input files are (PF)NanoAOD files', 
        default=False, 
        action='store_true',
    )

    #parser.add_argument(
    #    '-leptons', '--include_leptons_sf',
    #    help='Set True if input files are (PF)NanoAOD files',
    #    default=False,
    #    action='store_true',
    #)


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
        if href.endswith('.txt') or href.endswith('.gz'):
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
            download_raw_files(url=f"https://github.com/mcremone/decaf/tree/UL/analysis/data/JetMETCorr/{year}_UL", output_dir=path_saved_corrections + f"/met/{year}_UL")

        if observable == 'jet':

            #get the file jerc from https://github.com/mcremone/decaf/tree/UL/analysis/data/jerc
            download_raw_files(url="https://github.com/mcremone/decaf/tree/UL/analysis/data/jerc", output_dir=path_saved_corrections + "/jerc")


    if pog == 'lum':

        #path where to save the corrections
        path_saved_corrections = f"{out_path_corrections}/lum/"

        #if f"{path_saved_corrections}/{year}_UL/ is not already created, create it
        os.makedirs(f"{path_saved_corrections}/{year}_UL", exist_ok=True)

        if observable == 'pu':

            #check if file puWeights.json.gz is already in the directory
            if os.path.exists(f"{path_saved_corrections}/{year}_UL/puWeights.json.gz"):
                print(f"File puWeights.json.gz already exists in {path_saved_corrections}/{year}_UL")
                return path_saved_corrections


            #get the file puWeights.json.gz from https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/LUM, save it in out_path_corrections/lum/year_UL
            download_raw_files(url=f"https://github.com/mcremone/decaf/tree/UL/analysis/data/PUweight/{year}_UL", output_dir=path_saved_corrections + f"/{year}_UL/")



    #if pog == 'egamma':
    #
    #    #path where to save the corrections
    #    path_saved_corrections = f"{out_path_corrections}/EGammaSF/"
    #
    #    #if f"{path_saved_corrections}/{year}_UL/ is not already created, create it
    #    os.makedirs(f"{path_saved_corrections}/{year}_UL", exist_ok=True)
    #
    #    #check if file electron.json.gz is already in the directory
    #    if os.path.exists(f"{path_saved_corrections}/{year}_UL/electron.json.gz"):
    #        print(f"File electron.json.gz already exists in {path_saved_corrections}/{year}_UL")
    #        return path_saved_corrections
    #
    #    #get the file electron.json.gz from https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM, save it in out_path_corrections/EGammaSF/year_UL
    #    download_raw_files(url=f"


    return path_saved_corrections


####
# Electron ID scale factor
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaSFJSON
# jsonPOG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####
'''
def get_ele_loose_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file(f'{path}/EGammaSF/'+year+'_UL/electron.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["UL-Electron-ID-SF"].evaluate(year, "sf", "Loose", flateta, flatpt)

    return ak.unflatten(weight, counts=counts)

def get_ele_tight_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file(f'{path}/EGammaSF/'+year+'_UL/electron.json.gz')

    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator["UL-Electron-ID-SF"].evaluate(year, "sf", "Tight", flateta, flatpt)
    
    return ak.unflatten(weight, counts=counts)

def get_ele_reco_sf_below20(year, eta, pt):
    ele_reco_files_below20 = {
        '2016postVFP': f"{path}/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root:EGamma_SF2D",
        '2016preVFP': f"{path}/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2016preVFP.root:EGamma_SF2D",
        '2017': f"{path}/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2017.root:EGamma_SF2D",
        '2018': f"{path}/ElectronRecoSF/egammaEffi_ptBelow20.txt_EGM2D_UL2018.root:EGamma_SF2D"
    }

    corr = convert.from_uproot_THx(ele_reco_files_below20[year])
    evaluator = corr.to_evaluator()
    
    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    eta = ak.where((eta<2.399), ak.full_like(eta,-2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    pt = ak.where((pt<10.), ak.full_like(pt,10.), pt)
    pt = ak.where((pt>19.99), ak.full_like(pt,19.99), pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator.evaluate(flateta, flatpt)
    return ak.unflatten(weight, counts=counts)
    #get_ele_reco_err_below20[year]=lookup_tools.dense_lookup.dense_lookup(ele_reco_hist.variances() ** 0.5, ele_reco_hist.axes)


def get_ele_reco_sf_above20(year, eta, pt):
    ele_reco_files_above20 = {
        '2016postVFP': f"{path}/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root:EGamma_SF2D",
        '2016preVFP': f"{path}/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root:EGamma_SF2D",
        '2017': f"{path}/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2017.root:EGamma_SF2D",
        '2018': f"{path}/ElectronRecoSF/egammaEffi_ptAbove20.txt_EGM2D_UL2018.root:EGamma_SF2D"
    }
    
    corr = convert.from_uproot_THx(ele_reco_files_above20[year])
    evaluator = corr.to_evaluator()
    
    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    eta = ak.where((eta<2.399), ak.full_like(eta,-2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)
    
    pt = ak.where((pt<20.), ak.full_like(pt,20.), pt)
    pt = ak.where((pt>499.99), ak.full_like(pt,499.99), pt)
    flatpt = ak.flatten(pt)
    
    weight = evaluator.evaluate(flateta, flatpt)
    return ak.unflatten(weight, counts=counts)
    #get_ele_reco_err_above20[year]=lookup_tools.dense_lookup.dense_lookup(ele_reco_hist.variances() ** 0.05, ele_reco_hist.axes)
    

####
# Muon ID scale factor
# https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2018n?topic=MuonUL2018
# jsonPOG: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/MUO
# /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration
####

def get_mu_loose_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file(f'{path}/MuonSF/'+year+'_UL/muon_Z.json.gz')

    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    if year == '2018':
        weight = evaluator["NUM_LooseID_DEN_TrackerMuons"].evaluate(year+'_UL', flateta, flatpt, "sf")
    else:
        weight = evaluator["NUM_LooseID_DEN_genTracks"].evaluate(year+'_UL', flateta, flatpt, "sf")

    return ak.unflatten(weight, counts=counts)

def get_mu_tight_id_sf (year, eta, pt):
    evaluator = correctionlib.CorrectionSet.from_file(f'{path}/MuonSF/'+year+'_UL/muon_Z.json.gz')
    
    eta = ak.where((eta>2.399), ak.full_like(eta,2.399), eta)
    flateta, counts = ak.flatten(eta), ak.num(eta)

    pt  = ak.where((pt<15.),ak.full_like(pt,15.),pt)
    flatpt = ak.flatten(pt)
    
    if year == '2018':
        weight = evaluator["NUM_TightID_DEN_TrackerMuons"].evaluate(year+'_UL', flateta, flatpt, "sf")
    else:
        weight = evaluator["NUM_TightID_DEN_genTracks"].evaluate(year+'_UL', flateta, flatpt, "sf")
    
    return ak.unflatten(weight, counts=counts)


###
# Muon scale and resolution (i.e. Rochester)
# https://twiki.cern.ch/twiki/bin/view/CMS/RochcorMuon
###

#tag = 'roccor.Run2.v5'
#get_mu_rochester_sf = {}
#for year in ['2016postVFP', '2016preVFP', '2017','2018']:
#    if '2016postVFP' in year: 
#        fname = f'{path}/{tag}/RoccoR2016bUL.txt'
#    elif '2016preVFP' in year:  
#        fname = f'{path}/{tag}/RoccoR2016aUL.txt'
#    else:
#        fname = f'{path}/{tag}/RoccoR{year}UL.txt'
#    sfs = lookup_tools.txt_converters.convert_rochester_file(fname,loaduncs=True)
#    get_mu_rochester_sf[year] = lookup_tools.rochester_lookup.rochester_lookup(sfs)
'''



####
# PU weight
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/LUM
####
def get_pu_weight(year, trueint):
    correction = {'2018': 'Collisions18_UltraLegacy_goldenJSON',
                  '2017': 'Collisions17_UltraLegacy_goldenJSON',
                  '2016preVFP': 'Collisions16_UltraLegacy_goldenJSON',
                  '2016postVFP':'Collisions16_UltraLegacy_goldenJSON'}
    evaluator = correctionlib.CorrectionSet.from_file(f'{out_path_corrections}/lum/'+year+'_UL/puWeights.json.gz')
    weight_nom = evaluator[correction[year]].evaluate(trueint, 'nominal')
    weight_up = evaluator[correction[year]].evaluate(trueint, 'up')
    weight_down = evaluator[correction[year]].evaluate(trueint, 'down')

    return weight_nom, weight_up, weight_down



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
        
    
    evaluator = correctionlib.CorrectionSet.from_file(f'{out_path_corrections}/jme/met/{year}_UL/met.json.gz')

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


def fetch_lum_corrections():

    for year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
        _ = fetch_and_save_files_corrections(pog="lum",observable="pu",year=year)

    return 


def fetch_xy_met_corrections():
    
        #fetch the corrections for the met
        for year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
          _ = fetch_and_save_files_corrections(pog='jme',observable='met',year=year)
    
        return

def fetch_leptons_corrections():

    #fetch the corrections for the electrons
    for year in ["2016preVFP", "2016postVFP", "2017", "2018"]:
        fetch_and_save_files_corrections(pog='egamma',observable='ele',year=year)
        fetch_and_save_files_corrections(pog='muon',observable='mu',year=year)

    #fetch the corrections for the muons
    fetch_and_save_files_corrections(pog='lum',observable='mu',year='2016preVFP')
    fetch_and_save_files_corrections(pog='lum',observable='mu',year='2016postVFP')
    fetch_and_save_files_corrections(pog='lum',observable='mu',year='2017')
    fetch_and_save_files_corrections(pog='lum',observable='mu',year='2018')

    return



def build_corrections(corrections_to_include_dict):

    corrections = {}
    for key, value in corrections_to_include_dict.items():
        corrections[key] = value


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

        fetch_xy_met_corrections()
        corrections_to_include_dict["get_met_xy_correction"] = XY_MET_Correction    

    #if args.include_leptons_sf or args.include_all_corrections:
    #    print("Including lepton scale factors...")

    #    fetch_leptons_corrections()

        #TODO: add mva id sf for electrons
    #    corrections_to_include_dict["get_ele_loose_id_sf"] = get_ele_loose_id_sf
    #    corrections_to_include_dict["get_ele_tight_id_sf"] = get_ele_tight_id_sf


    #    corrections_to_include_dict["get_mu_loose_id_sf"] = get_mu_loose_id_sf
    #    corrections_to_include_dict["get_mu_tight_id_sf"] = get_mu_tight_id_sf

        #corrections_to_include_dict["get_mu_rochester_sf"] = get_mu_rochester_sf


    if args.include_pu_corrections or args.include_all_corrections:
        print("Including PU corrections...")

        fetch_lum_corrections()

        corrections_to_include_dict["get_pu_weight"] = get_pu_weight

    corrections = build_corrections(corrections_to_include_dict)

    print("Corrections built successfully !")

    #build a string with the arguments
    arguments = ""

    if args.include_all_corrections:
        arguments += "_all_corr"
        
    if args.include_jme_corrections and not args.include_all_corrections:
            arguments += "_jme"

    if args.include_pu_corrections and not args.include_all_corrections:
            arguments += "_pu"

    if not args.include_all_corrections:
            arguments += "_corr"

    print(f"Saving corrections to {out_path_corrections}/corrections_{date_time}{arguments}.coffea")
    save(corrections, f'{out_path_corrections}/corrections_{date_time}{arguments}.coffea')
    print(f"Corrections saved successfully.")

if __name__ == "__main__":

    main()