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

# Needed so that ak.zip({"pt": [...], "eta": [...], "phi": [...], "mass": [...]},
#                         with_name="PtEtaPhiMLorentzVector")
# is understood as a PtEtaPhiMLorentzVector from coffea.nanoevents.methods.vector_
ak.behavior.update(vector.behavior)


def update_cut_flow(cut_flow, cut_name, events=None, sumw=None):
    """Update cut flow table in a coffea accumulator.

    Args:
        cut_flow (dict[str, float])
        cut_name (str): the name of the cut to appear in the cut flow tree
        events (ak.Array): the events array from which to compute the number of events. 
            None only if `sumw` is not None.
        sumw (float): The sum of weights to add to the cut flow. 
            None if `events` is not None.
    """

    if sumw is not None:
        n_events = sumw
    else:
        n_events = get_number_of_events(events)

    if cut_name in cut_flow.keys():
        cut_flow[cut_name] += n_events
    else:
        cut_flow[cut_name] = n_events


def add_variations_to_cutflow(cut_flow, var_name, nominal, up, down):
    """Add normalization factors for nominal/up/down variations

    Args:
        cut_flow (dict[str, float])
        var_name (str): the name of the cut to appear in the cut flow tree
        nominal (float)
        up (float)
        down (float)
    """

    cut_flow[f"sumw_{var_name}_nominal"] = nominal
    cut_flow[f"sumw_{var_name}_up"] = up
    cut_flow[f"sumw_{var_name}_down"] = down


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


def __jet_var_i(var,i,pad_value=np.Inf):
    padded_var = ak.fill_none(ak.pad_none(var,i+1),pad_value)
    return padded_var[:,i]


def __get_phi_spike_filter(hot_spots_dict,var_name,j_eta_i,j_phi_i,rad):
    hot_etas_i, hot_phis_i = hot_spots_dict[var_name]
    hot_etas_i_reshaped = np.reshape(hot_etas_i,(len(hot_etas_i),1))
    hot_phis_i_reshaped = np.reshape(hot_phis_i,(len(hot_phis_i),1))
    j_eta_i_reshaped = np.broadcast_to(list(j_eta_i),(len(hot_etas_i),len(j_eta_i)))
    j_phi_i_reshaped = np.broadcast_to(list(j_phi_i),(len(hot_phis_i),len(j_phi_i)))
    return np.prod((j_eta_i_reshaped - hot_etas_i_reshaped)**2 + (j_phi_i_reshaped - hot_phis_i_reshaped)**2 > rad, axis=0, dtype=bool)


def apply_phi_spike_filter(
        events,
        year,
        hot_spots_pkl,
        n_jets,
        jets_eta,
        jets_phi,
    ):

    if year == "2016APV": year = "2016"
    with open(hot_spots_pkl,"rb") as infile:
        phi_spike_hot_spots = pickle.load(infile)
    rad = 0.028816*0.35 # the factor of 0.35 was optimized from the signal vs. background sensitivity study for s-channel
    hot_spots_dict = phi_spike_hot_spots[year]
    conditions = np.ones(len(events), dtype=bool)
    for i in range(n_jets):
        conditions &= __get_phi_spike_filter(
            hot_spots_dict,
            f"j{i+1}Phivsj{i+1}Eta",
            __jet_var_i(jets_eta, i),
            __jet_var_i(jets_phi, i),
            rad,
        )
    events = events[conditions]
    return events


def get_hem_veto_filter(*objects_list):

    eta_min = -3.05
    eta_max = -1.35
    phi_min = -1.62
    phi_max = -0.82

    veto = ak.zeros_like(range(len(objects_list[0])), dtype=bool)
    for objects in objects_list:
        hem_veto = (
            (objects.eta > eta_min)
            & (objects.eta < eta_max)
            & (objects.phi > phi_min)
            & (objects.phi < phi_max)
        )
        veto = veto | ak.any(hem_veto, axis=1)

    hem_filter = ~veto

    return hem_filter


def apply_hem_veto(events, *objects_list):
    hem_filter = get_hem_veto_filter(*objects_list)

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
    """Apply systematic uncertainty variations for cases in which recalculation
    of downstream quantities is needed e.g., JEC/JER variations.

    The following collections/branches are modified in place:
        * AK8 jets
        * MET
        * HT
    
    Args:
        events (ak.Array)
        variation (str): choose from "jec_up", "jec_down", "jer_up", "jer_down",
            "ue_up", "ue_down" for JEC, JER, unclustered energy up/down.

    Returns:
        events (ak.Array)
    """

    if variation is None: return events
    elif variation in ["jec_up", "jec_down", "jer_up", "jer_down", "ue_up", "ue_down"]:
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
            met_variation_name = "METUp" if "_up" in variation else "METDown"
            events = ak.with_field(
                events,
                events[met_variation_name][:, 1],
                "MET",
            )
            ht_variation_name = "HTJECup" if "_up" in variation else "HTJECdown"
            events = ak.with_field(
                events,
                events[ht_variation_name],
                "HT",
            )
        elif "jer" in variation:
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
            met_variation_name = "METUp" if "_up" in variation else "METDown"
            events = ak.with_field(
                events,
                events[met_variation_name][:, 0],
                "MET",
            )
            ht_variation_name = "HTJERup" if "_up" in variation else "HTJERdown"
            events = ak.with_field(
                events,
                events[ht_variation_name],
                "HT",
            )
        else:
            met_variation_name = "METUp" if "_up" in variation else "METDown"
            events = ak.with_field(
                events,
                events[met_variation_name][:, 5],
                "MET",
            )

        if "jec" in variation or "jer" in variation:
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


def __add_weight_variations(events, variation_up, variation_down, variation_name):
    """Add the up/down weight branches to the events.
    
    Args:
        events (ak.Array): the events erray to which to add the up/down weights variations
        variation_up (ak.Array): the up variations from the nominal weights
        variation_down (ak.Array): the down variations from the nominal weights
        variation_name (str): the name of the variation

    Returns:
        ak.Array, float, float: events, sumw up, and sumw down
    """

    weight_name = "Weight" if is_tree_maker(events) else "genWeight"
    nominal_weights = events[weight_name]

    weights_up = nominal_weights * variation_up
    weights_down = nominal_weights * variation_down

    # Create the new branches
    events[f"{weight_name}{variation_name}Up"] = weights_up
    events[f"{weight_name}{variation_name}Down"] = weights_down

    # Compute the sum of weights for the variations
    sumw_up = ak.sum(weights_up)
    sumw_down = ak.sum(weights_down)

    return events, sumw_up, sumw_down


def apply_scale_variations(events):
    """Calculate up/down renormalization and factorisation scale variation.
    
    Following definition here: https://github.com/TreeMaker/TreeMaker/blob/7a81115566ed1f2206eb4d447c9c7ba0870d88d0/Utils/src/PDFWeightProducer.cc#L167
    This should be done **before** any event selection is applied.

    The following collections/branches are added in place:
        * WeightScaleUp (TreeMaker) / genWeightScaleUp (PFNanoAOD)
        * WeightScaleDown (TreeMaker) / genWeightScaleDown (PFNanoAOD)

    The definition of the weights includes normalization factors, such that
    the normalization to unit luminosity is the same for the variations and
    the nominal weights.

    Args:
        events (ak.Array)

    Returns:
        ak.Array, float, float: events, sumw up, and sumw down
    """

    # Calculate up/down variations
    variation_up = ak.max(events.ScaleWeights[:, [i for i in range(9) if i not in (5, 7)]], axis=-1)
    variation_down = ak.min(events.ScaleWeights[:, [i for i in range(9) if i not in (5, 7)]], axis=-1)

    return __add_weight_variations(events, variation_up, variation_down, "Scale")


def apply_pdf_variations(events):
    """Calculate the PDF up/down variations.
    
    This should be done **before** any event selection is applied.

    The following collections/branches are added in place:
        * WeightPDFUp (TreeMaker) / genWeightPDFUp (PFNanoAOD)
        * WeightPDFDown (TreeMaker) / genWeightPDFDown (PFNanoAOD)

    The definition of the weights includes normalization factors, such that
    the normalization to unit luminosity is the same for the variations and
    the nominal weights.

    Args:
        events (ak.Array)

    Returns:
        ak.Array, float, float: events, sumw up, and sumw down
    """
    
    # Normalize the array of pdf weights by the first entry
    if is_tree_maker(events):
        pdf_variations = events.PDFweights.to_numpy()
        max_value = np.max(pdf_variations, where=pdf_variations<20, initial=1)
        pdf_variations = np.clip(pdf_variations, a_min=None, a_max=max_value)
        pdf_variations = pdf_variations / pdf_variations[:, :1]
    else:
        raise NotImplementedError()

    # Calculate the mean and standard deviation across replicas per event
    mean = np.mean(pdf_variations, axis=1)
    std = np.std(pdf_variations, axis=1)

    # Calculate the up/down variations across events
    variation_up = mean + std
    variation_down = mean - std

    return __add_weight_variations(events, variation_up, variation_down, "PDF")

