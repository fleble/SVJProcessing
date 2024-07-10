import numpy as np
import awkward as ak

from skimmer import skimmer_utils
from utils.variables_computation import event_variables as event_vars
from utils.data.triggers import primary_dataset_triggers
from utils.tree_maker.triggers import trigger_table
from analysis_configs import objects_definition_s_channel_leptons as obj

from utils.Logger import *


def remove_primary_dataset_overlap(events, year, primary_dataset):

    def trigger_list_to_index_list(year, primary_dataset):
        return [trigger_table[trigger] for trigger in primary_dataset_triggers[primary_dataset][year]]

    def is_passing_primary_dataset_triggers(events, year, primary_dataset):
        trigger_indices = trigger_list_to_index_list(year, primary_dataset)
        primary_dataset_trigger_mask = np.zeros(len(events.TriggerPass[0]), dtype=int)
        primary_dataset_trigger_mask[trigger_indices] = 1
        primary_dataset_trigger_mask = np.tile(primary_dataset_trigger_mask, [len(events), 1])
        return ak.any(events.TriggerPass * primary_dataset_trigger_mask, axis=1)

    pass_jetht_triggers = is_passing_primary_dataset_triggers(events, year, "JetHT")
    pass_met_triggers = is_passing_primary_dataset_triggers(events, year, "MET")

    if primary_dataset == "JetHT":
        overlap_mask = np.ones(len(events), dtype=bool)
    elif primary_dataset == "MET":
        overlap_mask = np.bitwise_not(pass_jetht_triggers)
    elif primary_dataset == "HTMHT":
        overlap_mask = np.bitwise_not(pass_jetht_triggers | pass_met_triggers)
    else:
        log.critical(f"Invalid primary dataset name {primary_dataset}!")
        exit(1)
    
    return events[overlap_mask]


def remove_single_lepton_primary_dataset_overlap(events, year, primary_dataset):

    def trigger_list_to_index_list(year, primary_dataset):
        return [trigger_table[trigger] for trigger in primary_dataset_triggers[primary_dataset][year]]

    def is_passing_primary_dataset_triggers(events, year, primary_dataset):
        trigger_indices = trigger_list_to_index_list(year, primary_dataset)
        primary_dataset_trigger_mask = np.zeros(len(events.TriggerPass[0]), dtype=int)
        primary_dataset_trigger_mask[trigger_indices] = 1
        primary_dataset_trigger_mask = np.tile(primary_dataset_trigger_mask, [len(events), 1])
        return ak.any(events.TriggerPass * primary_dataset_trigger_mask, axis=1)

    pass_single_muon_triggers = is_passing_primary_dataset_triggers(events, year, "SingleMuon")
    if year in ["2016", "2016APV", "2017"]:
        single_electron_trigger = "SingleElectron"
    else:
        single_electron_trigger = "EGamma"

    if primary_dataset == "SingleMuon":
        overlap_mask = np.ones(len(events), dtype=bool)
    elif primary_dataset == single_electron_trigger:
        overlap_mask = np.bitwise_not(pass_single_muon_triggers)
    else:
        log.critical(f"Invalid primary dataset name {primary_dataset}!")
        exit(1)
    
    return events[overlap_mask]


def apply_good_ak8_jet_filter(events):
    #analysis_jets = events.FatJet[obj.is_analysis_ak8_jet(events.FatJet)]
    #good_ak8_jets_filter = ak.all(obj.is_good_ak8_jet(analysis_jets), axis=1)
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


def apply_lepton_veto(
        events,
        electron_extra_condition=None,
        muon_extra_condition=None,
    ):

    electron_condition = obj.is_veto_electron(events.Electrons)
    if electron_extra_condition is not None:
        electron_condition = electron_condition & electron_extra_condition

    muon_condition = obj.is_veto_muon(events.Muons)
    if muon_extra_condition is not None:
        muon_condition = muon_condition & muon_extra_condition

    veto_electrons = events.Electrons[electron_condition]
    veto_muons = events.Muons[muon_condition]
    n_veto_electrons = ak.count(veto_electrons.pt, axis=1)
    n_veto_muons = ak.count(veto_muons.pt, axis=1)
    n_veto_leptons = n_veto_electrons + n_veto_muons
    events = events[n_veto_leptons == 0]
    return events


def require_n_veto_leptons(events, n):
    is_veto_electron = obj.is_veto_electron(events.Electrons)
    is_veto_muon = obj.is_veto_muon(events.Muons)
    veto_electrons = events.Electrons[is_veto_electron]
    veto_muons = events.Muons[is_veto_muon]
    n_veto_electrons = ak.count(veto_electrons.pt, axis=1)
    n_veto_muons = ak.count(veto_muons.pt, axis=1)
    n_veto_leptons = n_veto_electrons + n_veto_muons
    events = events[n_veto_leptons == n]
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
