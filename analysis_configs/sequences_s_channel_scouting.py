import numpy as np
import awkward as ak

from skimmer import skimmer_utils
from utils.variables_computation import event_variables as event_vars
from utils.variables_computation import jet_variables as jet_vars
from utils.data.triggers import primary_dataset_triggers
from utils.tree_maker.triggers import trigger_table
from analysis_configs import objects_definition_s_channel_scouting as obj

from utils.Logger import *


def __unflatten_to_event_level(self, *args):
    arrays = [ak.unflatten(jet_variable, self.n_jets, axis=0)
                for jet_variable in args]
    if len(arrays) == 1:
        return arrays[0]
    else:
        return arrays

def apply_good_ak8_jet_filter(events):
    #analysis_jets = events.FatJet[obj.is_analysis_ak8_jet(events.FatJet)]
    #good_ak8_jets_filter = ak.all(obj.is_good_ak8_jet(analysis_jets), axis=1)
    events = add_branches_for_ak8_jet_id(events)
    events = add_jet_id_branch(events)
    ak8_jets = ak.zip(
        {
            "pt": events.FatJet_pt,
            "mass": events.FatJet_mass,
            "eta": events.FatJet_eta,
            "phi": events.FatJet_phi,
            "id": events.FatJet_jetId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    
    analysis_jets = ak8_jets[obj.is_analysis_ak8_jet(ak8_jets)]
    good_ak8_jets_filter = ak.all(obj.is_good_ak8_jet(analysis_jets), axis=1)
    events = events[good_ak8_jets_filter]
    return events


def add_branches_for_ak8_jet_id(events):
        
        all_variables = ["chHEF", "neHEF", "electron_energy_fraction", "muon_energy_fraction", "photon_energy_fraction"]

        if which == "all":
            which = all_variables

        for variable_name in which:
            function = getattr(jet_vars, f"calculate_{variable_name}")
            variable = function(jet_pf_cands_per_jet, flat_jets)
            branch_name = f"FatJet_{variable_name}"
            events[branch_name] = __unflatten_to_event_level(variable)

def add_jet_id_branch(events):
    #CZZ: here custom implementation of id function from offline analysis, to be substituted with variable in ntuplizer
    # passID
    CHM = events.FatJet_chargedhadron_multiplicity +events.FatJet_electron_multiplicity + events.FatJet_muon_multiplicity
    passID = ((abs(events.FatJet_eta)<=2.6) 
                & (events.FatJet_electron_energy_fraction<0.8) 
                   & (CHM>0) & (events.FatJet_chHEF>0) 
                   & (events.FatJet_multiplicity>1) & (events.FatJet_photon_energy_fraction<0.9) 
                   & (events.FatJet_muon_energy_fraction<0.8) & (events.FatJet_neHEF<0.9) )

    events["FatJet_jetId"] = passID
    return events

def add_good_ak8_jet_branch(events):
    ak8_jets = ak.zip(
        {
            "pt": events.FatJet_pt,
            "mass": events.FatJet_mass,
            "eta": events.FatJet_eta,
            "phi": events.FatJet_phi,
            "id": events.FatJet_jetId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    
    is_good_analysis_ak8_jet = (
        obj.is_analysis_ak8_jet(ak8_jets)
        & obj.is_good_ak8_jet(ak8_jets)
    )

    #add new branch to the events
    events["FatJet_isGood"] = is_good_analysis_ak8_jet

    return events



def add_analysis_branches(events):

    # Event variables
    good_jets_ak8_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=events.FatJet_pt[events.FatJet_isGood],
        eta=events.FatJet_eta[events.FatJet_isGood],
        phi=events.FatJet_phi[events.FatJet_isGood],
        mass=events.FatJet_mass[events.FatJet_isGood],
    )
    met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET_pt,
            phi=events.MET_phi,
    )

    #add rt variable
    mt = event_vars.calculate_transverse_mass(good_jets_ak8_lv, met)
    rt = met.pt / mt
    events["RTFatJet"] = rt

    #add mt variable
    events["MT01FatJetMET"] = mt

    #add deltaeta the two leading jets
    events["DeltaEtaJ0J1FatJet"] = event_vars.calculate_delta_eta(good_jets_ak8_lv)
    
    #add minimum delta phi between the MET and the two leading jets
    events["DeltaPhiMinFatJetMET"] = event_vars.calculate_delta_phi_min(good_jets_ak8_lv, met)
    
    return events


def remove_collections(events):
    events = events[[x for x in events.fields if x != "JetsAK15"]]
    events = events[[x for x in events.fields if x != "GenJetsAK15"]]
    return events
