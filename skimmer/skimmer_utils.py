import awkward as ak
import numpy as np
import numba as nb
from coffea.nanoevents.methods import vector
import uproot

from utils.awkward_array_utilities import as_type
from utils.tree_maker.triggers import trigger_table as trigger_table_treemaker
from utils.systematics import calc_jec_variation, calc_jer_variation
from utils.Logger import *

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


def apply_phi_spike_filter(events, year, jet_eta_branch_name="Jet_eta", jet_phi_branch_name="Jet_phi", reverse=False):
    rad = 0.028816 # half the length of the diagonal of the eta-phi rectangular cell
    rad *= 0.35 # the factor of 0.35 was optimized from the signal vs. background sensitivity study

    eta_lead = None
    eta_sub = None
    phi_lead = None
    phi_sub = None
    if year == "2016":
        eta_lead = [0.048,0.24,1.488,1.584,-1.008]
        phi_lead = [-0.35,-0.35,-0.77,-0.77,-1.61]
        eta_sub = [-1.2,-0.912,-0.912,-0.816,-0.72,-0.72,-0.528,-0.432,-0.336,-0.24,-0.24,-0.144,-0.144,-0.048,0.144,0.912,0.912,1.008,1.296,-1.584,-0.816,-0.72,-0.144,-0.048,-0.048,0.048,1.104,1.488]
        phi_sub = [-1.19,2.03,3.01,-1.75,-2.17,-0.77,2.73,2.73,0.21,0.07,0.21,-2.59,0.77,0.91,1.75,1.75,2.87,0.63,-0.49,0.63,1.47,-2.31,0.07,-2.59,0.77,0.91,-3.15,2.73]
    elif year == "2017":
        eta_lead = [0.144,1.488,1.488,1.584,-0.624]
        phi_lead = [-0.35,-0.77,-0.63,-0.77,0.91]
        eta_sub = [-0.912,-0.912,-0.816,-0.72,-0.528,-0.336,-0.24,-0.24,-0.144,-0.144,-0.048,0.144,0.912,0.912,1.008,-1.2,-0.72,-0.72,-0.432,0.336,0.624,1.104,1.296]
        phi_sub = [2.03,3.01,-1.75,-0.77,2.73,0.21,0.07,0.21,-2.59,0.77,0.91,1.75,1.75,2.87,0.63,-1.19,-2.31,-2.17,2.73,-0.77,-0.77,-3.15,-0.49]
    elif year == "2018":
        eta_lead = [1.488,1.488,1.584]
        phi_lead = [-0.77,-0.63,-0.77]
        eta_sub = [-1.584,-1.2,-0.912,-0.912,-0.816,-0.816,-0.72,-0.72,-0.528,-0.432,-0.336,-0.24,-0.24,-0.144,-0.144,-0.144,-0.048,-0.048,0.144,0.912,0.912,1.008,1.296,-0.72,1.104,1.488,1.776]
        phi_sub = [0.63,-1.19,2.03,3.01,-1.75,-0.77,-2.17,-0.77,2.73,2.73,0.21,0.07,0.21,-2.59,0.07,0.77,0.77,0.91,1.75,1.75,2.87,0.63,-0.49,-2.31,-3.15,-0.21,0.77]
    else:
        raise ValueError("Invalid year")

    eta_lead = nb.typed.List(eta_lead)
    eta_sub = nb.typed.List(eta_sub)
    phi_lead = nb.typed.List(phi_lead)
    phi_sub = nb.typed.List(phi_sub)

    jets_eta = getattr(events, jet_eta_branch_name)
    jets_phi = getattr(events, jet_phi_branch_name)

    builder = ak.ArrayBuilder()
    phi_spike_filter = __get_phi_spike_filter(builder, eta_lead, phi_lead, eta_sub, phi_sub, rad, jets_eta, jets_phi, reverse=reverse).snapshot()

    events = events[phi_spike_filter]

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

    jet_hem_condition = (
        (jets.eta > -3.05)
        & (jets.eta < -1.35)
        & (jets.phi > -1.62)
        & (jets.phi < -0.82)
    )
    electron_hem_condition = (
        (electrons.eta > -3.05)
        & (electrons.eta < -1.35)
        & (electrons.phi > -1.62)
        & (electrons.phi < -0.82)
    )
    muon_hem_condition = (
        (muons.eta > -3.05)
        & (muons.eta < -1.35)
        & (muons.phi > -1.62)
        & (muons.phi < -0.82)
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


@nb.jit
def __get_phi_spike_filter(builder, eta_lead, phi_lead, eta_sub, phi_sub, rad, jets_eta, jets_phi, reverse):
    for jet_eta, jet_phi in zip(jets_eta, jets_phi):
        if len(jet_eta) < 2:
            builder.append(True)
        else:
            keep_event = True
            for iep in range(len(eta_lead)):
                if (eta_lead[iep] - jet_eta[0])**2 + (phi_lead[iep] - jet_phi[0])**2 < rad:
                    keep_event = False
                    break
            for iep in range(len(eta_sub)):
                if (eta_sub[iep] - jet_eta[1])**2 + (phi_sub[iep] - jet_phi[1])**2 < rad:
                    keep_event = False
                    break
            if reverse:
                builder.append(not keep_event)
            else:
                builder.append(keep_event)
    return builder


def apply_variation(events, variation):
    """Apply systematic uncertainty variations for cases in which recalculation
    of downstream quantities is needed e.g., JEC/JER variations.

    The following collections/branches are modified in place:
        * AK8 jets
        * MET
        * HT
    
    Args:
        events (ak.Array)
        variation (str): choose from "jec_up", "jec_down", "jer_up", "jer_down".
    """

    if variation is None: return events
    elif variation in ["jec_up", "jec_down", "jer_up", "jer_down"]:
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
    

def apply_scale_variations(events):
    '''
    Get up/down variations envelope, following definition
    here: https://github.com/TreeMaker/TreeMaker/blob/7a81115566ed1f2206eb4d447c9c7ba0870d88d0/Utils/src/PDFWeightProducer.cc#L167
    
    Adds "ScaleWeight_up" and "ScaleWeight_down" branches, plus stores the integrals
    for nominal and up/down variations in the cutflow
    '''

    envelope_up = ak.max(events.ScaleWeights[:,[i for i in range(9) if i not in (5, 7)]], axis=-1)
    envelope_down = ak.min(events.ScaleWeights[:, [i for i in range(9) if i not in (5, 7)]], axis=-1)

    # Calculate central norm and variations
    if is_tree_maker(events):
        sum_w_nominal = ak.sum(events.Weight)
        sum_w_up = ak.sum(events.Weight * envelope_up)
        sum_w_down = ak.sum(events.Weight * envelope_down)

        events["WeightScaleUp"] = events.Weight * envelope_up
        events["WeightScaleDown"] = events.Weight * envelope_down
    else:
        sum_w_nominal = ak.sum(events.genWeight)
        sum_w_up = ak.sum(events.genWeight * envelope_up)
        sum_w_down = ak.sum(events.genWeight * envelope_down)

        events["genWeightScaleUp"] = events.genWeight * envelope_up
        events["genWeightScalDown"] = events.genWeight * envelope_down

    return events, sum_w_nominal, sum_w_up, sum_w_down


def apply_pdf_variations(events):
    '''
    
    '''