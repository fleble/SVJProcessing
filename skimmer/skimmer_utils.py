import awkward as ak
import numpy as np
import numba as nb
from coffea.nanoevents.methods import vector
import uproot
import pickle  

from utils.awkward_array_utilities import as_type
from utils.tree_maker.triggers import trigger_table as trigger_table_treemaker
from utils.systematics import calc_jec_variation, calc_jer_variation
from utils.Logger import *
from analysis_configs import objects_definition as obj

# Needed so that ak.zip({"pt": [...], "eta": [...], "phi": [...], "mass": [...]},
#                         with_name="PtEtaPhiMLorentzVector")
# is understood as a PtEtaPhiMLorentzVector from coffea.nanoevents.methods.vector_
ak.behavior.update(vector.behavior)


def update_cut_flow(cut_flow, cut_name, events):
    """Update cut flow table in a coffea accumulator.

    Args:
        cut_flow (dict[str, float])
        cut_name (str): the name of the cut to appear in the cut flow tree
        events (EventsFromAkArray)
        use_raw_events (bool)
    """

    if cut_name in cut_flow.keys():
        cut_flow[cut_name] += get_number_of_events(events)
    else:
        cut_flow[cut_name] = get_number_of_events(events)


def apply_trigger_cut(events, trigger_list):
    """Filter events using an or of all triggers.

    Args:
        events (ak.Array)
        trigger_list (list[str])

    Returns:
        ak.Array
    """

    #check if TriggerPass is inside events or not, and use it if it is (for treemaker skims)
    if is_tree_maker(events):
        trigger_filter = ak.zeros_like(events.EvtNum, dtype=bool)
        for trigger_name in trigger_list:
            trigger_index = trigger_table_treemaker[trigger_name]
            trigger_branch = events.TriggerPass[:, trigger_index]
            trigger_filter = trigger_filter | (trigger_branch == 1)
    else:
        for idx,trigger_name in enumerate(trigger_list):
            trigger_branch = getattr(events, trigger_name)
            if idx == 0:
                trigger_filter = (trigger_branch == 1)
            else:
                trigger_filter = trigger_filter | (trigger_branch == 1)
 
    events = events[trigger_filter]

    return events


def apply_met_filters_cut(events, met_filter_names):
    """MET filters cuts.

    Args:
        events (ak.Array)
        met_filter_names (list[str])

    Returns:
        ak.Array
    """

    for met_filter_name in met_filter_names:
        if met_filter_name not in events.fields:
            branch_name = "Flag_" + met_filter_name
        else:
            branch_name = met_filter_name

        met_filter = getattr(events, branch_name)
        events = events[met_filter>0]
            

    return events


def make_pt_eta_phi_mass_lorentz_vector(pt, eta=None, phi=None, mass=None):
    """Take pt, eta, phi, mass awkward arrays and return the corresponding PtEtaPhiMLorentzVector.

    eta and mass can be None, e.g. for the MET.
    """

    if eta is None:
        eta = ak.zeros_like(pt)
    if mass is None:
        mass = ak.zeros_like(pt)

    vec = ak.zip(
        {
            "pt": pt,
            "eta": eta,
            "phi": phi,
            "mass": mass,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    return vec

def make_pt_eta_phi_energy_lorentz_vector(pt, eta=None, phi=None, energy=None):
    """Take pt, eta, phi, mass awkward arrays and return the corresponding PtEtaPhiMLorentzVector.

    eta and mass can be None, e.g. for the MET.
    """

    if eta is None:
        eta = ak.zeros_like(pt)

    vec = ak.zip(
        {
            "pt": pt,
            "eta": eta,
            "phi": phi,
            "energy": energy,
        },
        with_name="PtEtaPhiELorentzVector",
    )

    return vec


def is_tree_maker(events):
    return "TriggerPass" in events.fields


def is_mc(events):
    if is_tree_maker(events):
        return "Weight" in events.fields
    else:
        return "genWeight" in events.fields


def is_data(events):
    return not is_mc(events)


def jet_var_i(var,i,pad_value=np.Inf):
    padded_var = ak.fill_none(ak.pad_none(var,i+1),pad_value)
    return padded_var[:,i]


def get_phi_spike_filter(hot_spots_dict,var_name,j_eta_i,j_phi_i,rad):
    hot_etas_i, hot_phis_i = hot_spots_dict[var_name]
    hot_etas_i_reshaped = np.reshape(hot_etas_i,(len(hot_etas_i),1))
    hot_phis_i_reshaped = np.reshape(hot_phis_i,(len(hot_phis_i),1))
    j_eta_i_reshaped = np.broadcast_to(list(j_eta_i),(len(hot_etas_i),len(j_eta_i)))
    j_phi_i_reshaped = np.broadcast_to(list(j_phi_i),(len(hot_phis_i),len(j_phi_i)))
    return np.prod((j_eta_i_reshaped - hot_etas_i_reshaped)**2 + (j_phi_i_reshaped - hot_phis_i_reshaped)**2 > rad, axis=0, dtype=bool)


def apply_phi_spike_filter(events, year, hot_spots_pkl, channel="t", jet_eta_branch_name="Jet_eta", jet_phi_branch_name="Jet_phi"):
    with open(hot_spots_pkl,"rb") as infile:
        phi_spike_hot_spots = pickle.load(infile)
    rad = 0.028816*0.35 # the factor of 0.35 was optimized from the signal vs. background sensitivity study for s-channel
    hot_spots_dict = phi_spike_hot_spots[year]
    jets_eta = getattr(events, jet_eta_branch_name)
    jets_phi = getattr(events, jet_phi_branch_name)
    if channel == "t":
        pass1 = get_phi_spike_filter(hot_spots_dict,"j1Phivsj1Eta",jetVar_i(jets_eta,0),jetVar_i(jets_phi,0),rad)
        pass2 = get_phi_spike_filter(hot_spots_dict,"j2Phivsj2Eta",jetVar_i(jets_eta,1),jetVar_i(jets_phi,1),rad)
        pass3 = get_phi_spike_filter(hot_spots_dict,"j3Phivsj3Eta",jetVar_i(jets_eta,2),jetVar_i(jets_phi,2),rad)
        pass4 = get_phi_spike_filter(hot_spots_dict,"j4Phivsj4Eta",jetVar_i(jets_eta,3),jetVar_i(jets_phi,3),rad)
    elif channel == "s":
        pass1 = get_phi_spike_filter(hot_spots_dict,"j1Phivsj1Eta",jetVar_i(jets_eta,0),jetVar_i(jets_phi,0),rad)
        pass2 = get_phi_spike_filter(hot_spots_dict,"j2Phivsj2Eta",jetVar_i(jets_eta,1),jetVar_i(jets_phi,1),rad)
    events = events[pass1 & pass2 & pass3 & pass4]
    return events


def apply_hem_veto(events):

    if is_tree_maker(events):
        jets = events.Jets
        electrons = events.Electrons
        muons = events.Muons

    else:
        jets = ak.zip({
            "pt": events["Jet_pt"],
            "eta": events["Jet_eta"],
            "phi": events["Jet_phi"],
        })
        electrons = ak.zip({
            "pt":  events["Electron_pt"],
            "eta": events["Electron_eta"],
            "phi": events["Electron_phi"],
        })
        muons = ak.zip({
            "pt":  events["Muon_pt"],
            "eta": events["Muon_eta"],
            "phi": events["Muon_phi"],
        })

    jet_condition = obj.is_good_ak4_jet(jets)
    electron_condition = obj.is_veto_electron(electrons)
    muon_condition = obj.is_veto_muon(muons)
    good_ak4_jets = jets[jet_condition]
    veto_electrons = electrons[electron_condition]
    veto_muons = muons[muon_condition]

    jet_hem_condition = (
        (good_ak4_jets.eta > -3.05)
        & (good_ak4_jets.eta < -1.35)
        & (good_ak4_jets.phi > -1.62)
        & (good_ak4_jets.phi < -0.82)
    )
    electron_hem_condition = (
        (veto_electrons.eta > -3.05)
        & (veto_electrons.eta < -1.35)
        & (veto_electrons.phi > -1.62)
        & (veto_electrons.phi < -0.82)
    )
    muon_hem_condition = (
        (veto_muons.eta > -3.05)
        & (veto_muons.eta < -1.35)
        & (veto_muons.phi > -1.62)
        & (veto_muons.phi < -0.82)
    )
    veto = (
        ((ak.num(jets) > 0) & ak.any(jet_hem_condition, axis=1))
        | ((ak.num(muons) > 0) & ak.any(muon_hem_condition, axis=1))
        | ((ak.num(electrons) > 0) & ak.any(electron_hem_condition, axis=1))
    )
    hem_filter = ~veto

    if is_mc(events):
        events = events[hem_filter]
    else:
        run_number = events.RunNum if is_tree_maker(events) else events.run
        hem_filter_data = (
            ((run_number >= 319077) & hem_filter)
            | (run_number < 319077)
        )
        events = events[hem_filter_data]

    return events


def is_clean(
        cleaned_collection,
        cleaning_collection,
        radius,
    ):

    cleaned_objects = make_pt_eta_phi_mass_lorentz_vector(
        pt=cleaned_collection.pt,
        eta=cleaned_collection.eta,
        phi=cleaned_collection.phi,
        mass=cleaned_collection.mass,
    )
    cleaning_objects = make_pt_eta_phi_mass_lorentz_vector(
        pt=cleaning_collection.pt,
        eta=cleaning_collection.eta,
        phi=cleaning_collection.phi,
        mass=cleaning_collection.mass,
    )

    
    cleaned_objects_, cleaning_objects_ = ak.unzip(ak.cartesian([cleaned_objects, cleaning_objects], axis=-1, nested=True))
    min_delta_r = ak.min(cleaned_objects_.delta_r(cleaning_objects_), axis=-1)
    min_delta_r = ak.fill_none(min_delta_r, 2*radius)
    filter = (min_delta_r > radius)
    filter = as_type(filter, bool)

    return filter


def collections_matching(
        collection1,
        collection2,
    ):
    """Return the mapping of closest object in collection1 to each object of collection2.
    
    Use as:
    mapping = collections_matching(collection1, collection2)
    matched_collection2 = collection2[mapping]

    The matching is a simple closest delta R with repetition!
    There can be several times the same object from collection2 matched.

    Args:
        events (EventsFromAkArray)
        collection1 (ak.Array): ak array with fields "pt", "eta", "phi", "mass"
        collection2 (ak.Array): ak array with fields "pt", "eta", "phi", "mass" 
    """

    collection1_4vector = make_pt_eta_phi_mass_lorentz_vector(
        pt=collection1.pt,
        eta=collection1.eta,
        phi=collection1.phi,
        mass=collection1.mass,
    )
    collection2_4vector = make_pt_eta_phi_mass_lorentz_vector(
        pt=collection2.pt,
        eta=collection2.eta,
        phi=collection2.phi,
        mass=collection2.mass,
    )

    collection1_, collection2_ = ak.unzip(ak.cartesian([collection1_4vector, collection2_4vector], axis=-1, nested=True))
    mapping = ak.argmin(collection1_.delta_r(collection2_), axis=-1)
    mapping = as_type(mapping, int)

    return mapping


def get_b_tagging_wp(year):
    """Tight b-tagging working points."""

    if year == "2016":
        return 0.6502
    elif year == "2016APV":
        return 0.6377
    elif year == "2017":
        return 0.7476
    elif year == "2018":
        return 0.7100


def get_number_of_events(events):
    if is_mc(events):
        if is_tree_maker(events):
            return ak.sum(events.Weight)
        else:
            return ak.sum(events.genWeight)
    else:
        return len(events)


def get_cut_flow_from_skims(input_file, cut_flow_tree):
    f = uproot.open(input_file)
    cut_flow = f["CutFlow"].arrays(cut_flow_tree.keys(),  library="pd")
    return cut_flow.to_dict("list")


def apply_variation(events, variation):
    if variation is None: return events
    elif variation in ["jec_up", "jec_down", "jer_up", "jer_down"]:
        # Calculate the varied jet kinematics
        if "jec" in variation:
            pt_var, eta_var, phi_var, energy_var, permutation = calc_jec_variation(
                events.JetsAK8.pt,
                events.JetsAK8.eta,
                events.JetsAK8.phi,
                events.JetsAK8.energy,
                events.JetsAK8.jerFactor,
                events.JetsAK8.jecUnc,
                events.JetsAK8.origIndex,
                events.JetsAK8JECup.o if "_up" in variation else events.JetsAK8JECdown.o,
                events.JetsAK8JECup.j if "_up" in variation else events.JetsAK8JECdown.j,
            )
        else:
            pt_var, eta_var, phi_var, energy_var, permutation = calc_jer_variation(
                events.JetsAK8.pt,
                events.JetsAK8.eta,
                events.JetsAK8.phi,
                events.JetsAK8.energy,
                events.JetsAK8.jerFactor,
                events.JetsAK8.origIndex,
                events.JetsAK8JERup.o if "_up" in variation else events.JetsAK8JERdown.o,
                events.JetsAK8.jerFactorUp if "_up" in variation else events.JetsAK8.jerFactorDown,
            )

        corrected_jets = make_pt_eta_phi_energy_lorentz_vector(
            pt=pt_var,
            eta=eta_var,
            phi=phi_var,
            energy=energy_var,
        )

        events["JetsAK8"] = ak.with_field(
            events["JetsAK8"],
            corrected_jets.pt,
            "pt",
        )
        events["JetsAK8"] = ak.with_field(
            events["JetsAK8"],
            corrected_jets.eta,
            "eta",
        )
        events["JetsAK8"] = ak.with_field(
            events["JetsAK8"],
            corrected_jets.phi,
            "phi",
        )
        events["JetsAK8"] = ak.with_field(
            events["JetsAK8"],
            corrected_jets.energy,
            "energy",
        )
        events["JetsAK8"] = ak.with_name(
            events["JetsAK8"],
            "PtEtaPhiELorentzVector",
        )
        for k in events["JetsAK8"].fields:
            if k in ["pt", "eta", "phi", "energy"]: continue

            events["JetsAK8"] = ak.with_field(
                events["JetsAK8"],
                events["JetsAK8"][k][permutation],
                k,
            )

        return events