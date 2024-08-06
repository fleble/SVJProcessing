import awkward as ak
from coffea.nanoevents.methods import vector

from utils.awkward_array_utilities import as_type
from utils.tree_maker.triggers import trigger_table
from utils.Logger import *
import uproot

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
        cut_flow[cut_name] += __get_number_of_events(events)
    else:
        cut_flow[cut_name] = __get_number_of_events(events)


def apply_trigger_cut(events, trigger_list):
    """Filter events using an or of all triggers.

    Args:
        events (ak.Array)
        trigger_list (list[str])

    Returns:
        ak.Array
    """

    trigger_filter = ak.zeros_like(events.EvtNum, dtype=bool)
    for trigger_name in trigger_list:
        trigger_index = trigger_table[trigger_name]
        trigger_branch = events.TriggerPass[:, trigger_index]
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
        met_filter = getattr(events, met_filter_name)
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


def is_mc(events):
    return "Weight" in events.fields


def is_data(events):
    return not is_mc(events)


def apply_hem_veto(events):
    """Apply HEM veto for the second part of 2018 data."""

    def hem_veto(events):

        jets = events.Jets
        electrons = events.Electrons
        muons = events.Muons

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
        return ~veto
    
    if is_mc(events):
        veto = hem_veto(events)
    else:
        veto = (
            ((events.RunNum >= 319077) & hem_veto(events))
            | (events.RunNum < 319077)
        )

    return events[veto]


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


def __get_number_of_events(events):
    if is_mc(events):
        return ak.sum(events.Weight)
    else:
        return len(events)


def get_cut_flow_from_skims(input_file, cut_flow_tree):
    f = uproot.open(input_file)
    cut_flow = f["CutFlow"].arrays(cut_flow_tree.keys(),  library="pd")
    return cut_flow.to_dict("list")

