import numpy as np
import awkward as ak

from skimmer import skimmer_utils
from utils.variables_computation import event_variables, jet_variables
from utils.inference_particlenet import run_jet_tagger
from utils.data.triggers import primary_dataset_triggers
from utils.tree_maker.triggers import trigger_table
from analysis_configs import objects_definition as obj

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


def apply_good_jet_filter(events):
    analysis_jets = events.JetsAK8[obj.is_analysis_jet(events.JetsAK8)]
    good_jets_filter = ak.all(obj.is_good_jet(analysis_jets), axis=1)
    events = events[good_jets_filter]
    return events


def add_good_jet_branch(events):
    is_good_analysis_jet = (
        obj.is_analysis_jet(events.JetsAK8)
        & obj.is_good_jet(events.JetsAK8)
    )
    events["JetsAK8"] = ak.with_field(
        events["JetsAK8"],
        is_good_analysis_jet,
        "isGood",
    )
    return events


def apply_lepton_veto(events):
    veto_electrons = events.Electrons[obj.is_veto_electron(events.Electrons)]
    veto_muons = events.Muons[obj.is_veto_muon(events.Muons)]
    n_veto_electrons = ak.count(veto_electrons.pt, axis=1)
    n_veto_muons = ak.count(veto_muons.pt, axis=1)
    n_veto_leptons = n_veto_electrons + n_veto_muons
    events = events[n_veto_leptons == 0]
    return events


def add_analysis_branches(events):

    # Jets AK8 variables
    new_branches = {}
    jets_ak8 = events.JetsAK8
    jets_ak8_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=jets_ak8.pt,
        eta=jets_ak8.eta,
        phi=jets_ak8.phi,
        mass=jets_ak8.mass,
    )
    met_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=events.MET,
        phi=events.METPhi,
    )

    # Kinematics
    new_branches["mass"] = jets_ak8_lv.mass
    new_branches["deltaPhiMET"] = jet_variables.calculate_delta_phi_with_met(jets_ak8_lv, met_lv)
    new_branches["LundJetPlaneZ"] = jet_variables.calculate_lund_jet_plane_z_with_met(jets_ak8_lv, met_lv)
    new_branches["MTMET"] = jet_variables.calculate_invariant_mass_with_met(jets_ak8_lv, met_lv)

    # Dark jet branches
    if skimmer_utils.is_mc(events):
        gen_jets_ak8 = events.GenJetsAK8
        if "hvCategory" in gen_jets_ak8.fields:  # If signal samples
            gen_jet_index = ak.mask(jets_ak8.genIndex, jets_ak8.genIndex >= 0)
            hv_category = gen_jets_ak8.hvCategory[gen_jet_index]
            hv_category = ak.fill_none(hv_category, -9999)
            new_branches["hvCategory"] = hv_category

            new_branches["isDarkJetTightNoMix"] = (
                (hv_category == 1)
                | (hv_category == 3)
                | (hv_category == 5)
                | (hv_category == 9)
            )

            new_branches["isDarkJetTight"] = (
                (hv_category != -9999)
                & (hv_category != 0)
                & (hv_category < 16)
            )

            new_branches["isDarkJetMedium"] = (
                (hv_category != -9999)
                & (hv_category != 0)
                & (hv_category != 16)
                & (hv_category != 17)
            )

            new_branches["isDarkJetLoose"] = (
                (hv_category != -9999)
                & (hv_category != 0)
                & (hv_category != 16)
            )

    for branch_name, branch in new_branches.items():
        events["JetsAK8"] = ak.with_field(
            events["JetsAK8"],
            branch,
            branch_name,
        )


    # Event variables
    nan_value = 0.  # Natural choice for missing values for LJP variables and delta eta / phi!
    good_jets_ak8 = events.JetsAK8[events.JetsAK8.isGood]
    good_jets_ak8_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=good_jets_ak8.pt,
        eta=good_jets_ak8.eta,
        phi=good_jets_ak8.phi,
        mass=good_jets_ak8.mass,
    )

    n_jets_max = 4
    for index_0 in range(n_jets_max):
        for index_1 in range(index_0+1, n_jets_max):
            delta_eta = event_variables.calculate_delta_eta(
                physics_objects=good_jets_ak8_lv,
                indices=(index_0, index_1),
                absolute_value=False,
                nan_value=None,
            )
            delta_phi = event_variables.calculate_delta_phi(
                physics_objects=good_jets_ak8_lv,
                indices=(index_0, index_1),
                absolute_value=False,
                nan_value=None,
            )
            delta_r = event_variables.calculate_delta_r(
                physics_objects=good_jets_ak8_lv,
                indices=(index_0, index_1),
                nan_value=nan_value,
            )
            dijet_mass = event_variables.calculate_invariant_mass(
                physics_objects=good_jets_ak8_lv,
                indices=(index_0, index_1),
                nan_value=nan_value,
            )
            lund_jet_plane_z = event_variables.calculate_lund_jet_plane_z(
                physics_objects=good_jets_ak8_lv,
                indices=(index_0, index_1),
                nan_value=nan_value,
            )
            delta_eta_abs = abs(delta_eta)
            delta_phi_abs = abs(delta_phi)
            delta_eta = ak.fill_none(delta_eta, nan_value)
            delta_phi = ak.fill_none(delta_phi, nan_value)
            delta_eta_abs = ak.fill_none(delta_eta_abs, nan_value)
            delta_phi_abs = ak.fill_none(delta_phi_abs, nan_value)

            events[f"DeltaEta{index_0}{index_1}GoodJetsAK8"] = delta_eta
            events[f"DeltaPhi{index_0}{index_1}GoodJetsAK8"] = delta_phi
            events[f"DeltaR{index_0}{index_1}GoodJetsAK8"] = delta_r
            events[f"DeltaEtaAbs{index_0}{index_1}GoodJetsAK8"] = delta_eta_abs
            events[f"DeltaPhiAbs{index_0}{index_1}GoodJetsAK8"] = delta_phi_abs
            events[f"DijetMass{index_0}{index_1}GoodJetsAK8"] = dijet_mass
            events[f"LundJetPlaneZ{index_0}{index_1}GoodJetsAK8"] = lund_jet_plane_z
   
    events["DeltaPhiMinGoodJetsAK8"] = ak.min(abs(good_jets_ak8.deltaPhiMET), axis=1)
    events["ST"] = events.MET + events.HT
    events["ATLASDeltaPhiMinMax"] = event_variables.calculate_atlas_delta_phi_max_min(
        jets=good_jets_ak8_lv,
        met=met_lv,
        nan_value=nan_value,
    )
    events["ATLASPtBalance"] = event_variables.calculate_atlas_momentum_balance(
        jets=good_jets_ak8_lv,
        met=met_lv,
        nan_value=nan_value,
    )

    return events


def add_particle_net_tagger(events):
    jets_ak8 = events.JetsAK8
    events["JetsAK8"] = ak.with_field(
        events["JetsAK8"],
        run_jet_tagger(events, jets_ak8),
        "pNetJetTaggerScore",
    )
    return events


def remove_collections(events):
    events = events[[x for x in events.fields if x != "JetsAK15"]]
    events = events[[x for x in events.fields if x != "GenJetsAK15"]]
    return events
