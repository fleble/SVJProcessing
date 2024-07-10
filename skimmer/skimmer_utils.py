import awkward as ak
import numpy as np
import numba as nb
from coffea.nanoevents.methods import vector

from utils.awkward_array_utilities import as_type
from utils.tree_maker.triggers import trigger_table as trigger_table_treemaker
from utils.nanoaod.triggers import trigger_table as trigger_table_nanoaod
from utils.Logger import *
import uproot

# Needed so that ak.zip({"pt": [...], "eta": [...], "phi": [...], "mass": [...]},
#                         with_name="PtEtaPhiMLorentzVector")
# is understood as a PtEtaPhiMLorentzVector from coffea.nanoevents.methods.vector_
ak.behavior.update(vector.behavior)


def calculate_number_of_objects(physics_objects):
    """Calculate number of a physics objects in all events.

    Args:
        physics_objects (awkward.Array): Jagged ak array where axis 0 is the
            event axis and axis 1 is the object axis.

    Returns:
        awkward.Array
    """

    return ak.num(physics_objects, axis=1)

def __get_pair_of_objects(physics_objects, indices=(0, 1), n_objects=None):
    """Returns a tuple of masked ak arrays corresponding to physics objects with required indices.

    Args:
        physics_objects (ak.Array)
        indices (tuple[int], optional): By default 2 leading objects
        n_objects (ak.Array, optional): Number of objects per event

    Returns:
        tuple[ak.Array]
    """

    if n_objects is None:
        n_objects = calculate_number_of_objects(physics_objects)

    max_index = max(indices)
    physics_objects = ak.mask(physics_objects, n_objects >= max_index + 1)

    obj0 = physics_objects[:, indices[0]]
    obj1 = physics_objects[:, indices[1]]

    return obj0, obj1


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

    #check if TriggerPass is inside events or not, and use it if it is (for treemaker skims)
    if "TriggerPass" in events.fields:
        trigger_filter = ak.zeros_like(events.EvtNum, dtype=bool)
        for trigger_name in trigger_list:
            trigger_index = trigger_table_treemaker[trigger_name]
            trigger_branch = events.TriggerPass[:, trigger_index]
            trigger_filter = trigger_filter | (trigger_branch == 1)
    else: #if TriggerPass is not inside events, nano trigger skim
        for idx,trigger_name in enumerate(trigger_list):
            trigger_index = trigger_table_nanoaod[trigger_name]
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
        events (EventsFromAkArray)
        cut_flow (dict[str, float])
        met_filter_names (list[str])

    Returns:
        None
    """

    for met_filter_name in met_filter_names:
        branch_name = "Flag_" + met_filter_name
        if branch_name in events.fields:
            met_filter = getattr(events, branch_name)
            events = events[met_filter]
            # Comment out the line below to have MET filter
            #update_cut_flow(cut_flow, utl.capitalize(met_filter_name), events)

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



def apply_hem_filter(events, year, jet_pt_branch_name="Jet_pt", jet_eta_branch_name="Jet_eta",
                     jet_phi_branch_name="Jet_phi"):

    hem_min_pt = 30
    hem_min_eta = -3.05
    hem_max_eta = -1.35
    hem_min_phi = -1.62
    hem_max_phi = -0.82

    if year != "2018":
        hem_filter = ak.ones_like(events.genWeight)
    else:
        jets_pt = getattr(events, jet_pt_branch_name)
        jets_eta = getattr(events, jet_eta_branch_name)
        jets_phi = getattr(events, jet_phi_branch_name)
        
        condition_per_jet = (jets_pt > hem_min_pt) \
            & (jets_eta > hem_min_eta) & (jets_eta < hem_max_eta) \
            & (jets_phi > hem_min_phi) & (jets_phi < hem_max_phi)
        hem_filter = ak.any(condition_per_jet, axis=1) == False

    if is_mc(events):
        events = events[hem_filter]
    else:
        hem_filter_data = (
            ((events.run >= 319077) & hem_filter)
            | (events.run < 319077)
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


def __get_number_of_events(events):
    if is_mc(events):
        return ak.sum(events.Weight)
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

