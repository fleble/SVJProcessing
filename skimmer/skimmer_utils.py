import awkward as ak
from coffea.nanoevents.methods import vector

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


def __get_number_of_events(events):
    if is_mc(events):
        return ak.sum(events.Weight)
    else:
        return len(events)


def get_cut_flow_from_skims(input_file, cut_flow_tree):
    f = uproot.open(input_file)
    cut_flow = f["CutFlow"].arrays(cut_flow_tree.keys(),  library="pd")
    return cut_flow.to_dict("list")

