import awkward as ak
import numpy as np
import numba as nb
from coffea.nanoevents.methods import vector
from coffea.util import load, save
import uproot
import pickle  
import cachetools

from utils.awkward_array_utilities import as_type
from utils.tree_maker.triggers import trigger_table as trigger_table_treemaker
from utils.systematics import calc_jec_variation, calc_jer_variation, calc_jerc_variations_PFNano, calc_unclustered_met_variations_PFNano, calc_custom_svj_jes_variations_PFNano, apply_jercs_PFNano, apply_jecs_PFNano, propagate_jecs_to_MET_PFNano, propagate_jecs_to_METSig_PFNano
from utils.met_significance_factory_pfnano import MetSignificanceCalculator
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


#def __jet_var_i(var,i,pad_value=np.Inf):
#    padded_var = ak.fill_none(ak.pad_none(var,i+1),pad_value)
#    return padded_var[:,i]
#
#
#def __get_phi_spike_filter(hot_spots_dict,var_name,j_eta_i,j_phi_i,rad):
#    hot_etas_i, hot_phis_i = hot_spots_dict[var_name]
#    hot_etas_i_reshaped = np.reshape(hot_etas_i,(len(hot_etas_i),1))
#    hot_phis_i_reshaped = np.reshape(hot_phis_i,(len(hot_phis_i),1))
#    j_eta_i_reshaped = np.broadcast_to(list(j_eta_i),(len(hot_etas_i),len(j_eta_i)))
#    j_phi_i_reshaped = np.broadcast_to(list(j_phi_i),(len(hot_phis_i),len(j_phi_i)))
#    return np.prod((j_eta_i_reshaped - hot_etas_i_reshaped)**2 + (j_phi_i_reshaped - hot_phis_i_reshaped)**2 > rad, axis=0, dtype=bool)
#
#
#def apply_phi_spike_filter(
#        events,
#        year,
#        hot_spots_pkl,
#        n_jets,
#        jets_eta,
#        jets_phi,
#    ):
#
#    if year == "2016APV": year = "2016"
#    with open(hot_spots_pkl,"rb") as infile:
#        phi_spike_hot_spots = pickle.load(infile)
#    rad = 0.028816*0.35 # the factor of 0.35 was optimized from the signal vs. background sensitivity study for s-channel
#    hot_spots_dict = phi_spike_hot_spots[year]
#    conditions = np.ones(len(events), dtype=bool)
#    for i in range(n_jets):
#        conditions &= __get_phi_spike_filter(
#            hot_spots_dict,
#            f"j{i+1}Phivsj{i+1}Eta",
#            __jet_var_i(jets_eta, i),
#            __jet_var_i(jets_phi, i),
#            rad,
#        )
#    events = events[conditions]
#    return events

def apply_phi_spike_filter(events, year, jet_eta_branch_name="Jet_eta", jet_phi_branch_name="Jet_phi", reverse=False):
    rad = 0.028816 # half the length of the diagonal of the eta-phi rectangular cell
    rad *= 0.35 # the factor of 0.35 was optimized from the signal vs. background sensitivity study

    eta_lead = None
    eta_sub = None
    phi_lead = None
    phi_sub = None
    if year == "2016" or year == "2016APV":
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

def get_hem_veto_filter(*objects_list):

    eta_min = -3.05
    eta_max = -1.35
    phi_min = -1.62
    phi_max = -0.82

    veto = ak.ones_like(range(len(objects_list[0])), dtype=bool)
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
        variation (str): choose from "jec_up", "jec_down", "jer_up", "jer_down".

    Returns:
        events (ak.Array)
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


def __add_weight_variations(events, variation_up, variation_down, variation_name, computed_nominal_weights=None, multiply_by_pu_weights=False):
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

    sumw_nom = ak.sum(nominal_weights)

    if "PU" in variation_name:
        nominal_weights = nominal_weights*computed_nominal_weights
        events[f"{weight_name}{variation_name}"] = nominal_weights
    #    sumw_nom = ak.sum(nominal_weights)

    weights_up = nominal_weights * variation_up
    weights_down = nominal_weights * variation_down

    # Compute the sum of weights for the variations
    sumw_up = ak.sum(weights_up)
    sumw_down = ak.sum(weights_down)

    if multiply_by_pu_weights:
        nominal_weights = nominal_weights*events["genWeightPU"]

    # Create the new branches
    if multiply_by_pu_weights:
        events[f"{weight_name}{variation_name}UpPUNom"] = weights_up
        events[f"{weight_name}{variation_name}DownPUNom"] = weights_down
    else:
        events[f"{weight_name}{variation_name}Up"] = weights_up
        events[f"{weight_name}{variation_name}Down"] = weights_down

    return events, sumw_up, sumw_down, sumw_nom


def apply_scale_variations(events,is_nano=False, multiply_by_pu_weight=False):
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
    if is_tree_maker(events):
        variation_up = ak.max(events.ScaleWeights[:, [i for i in range(9) if i not in (5, 7)]], axis=-1)
        variation_down = ak.min(events.ScaleWeights[:, [i for i in range(9) if i not in (5, 7)]], axis=-1)
    elif is_nano:
        variation_up = ak.max(events.LHEScaleWeight[:, [i for i in range(9) if i not in (5, 7)]], axis=-1)
        variation_down = ak.min(events.LHEScaleWeight[:, [i for i in range(9) if i not in (5, 7)]], axis=-1)
    else:
        raise NotImplementedError()

    return __add_weight_variations(events, variation_up, variation_down, "Scale", computed_nominal_weights=None, multiply_by_pu_weights=multiply_by_pu_weight)


def apply_pdf_variations(events, is_nano=False, multiply_by_pu_weight=False):
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
        
    elif is_nano:
       pdf_variations = events.LHEPdfWeight.to_numpy()

    else:
        raise NotImplementedError()
    
    max_value = np.max(pdf_variations, where=pdf_variations<20, initial=1)
    pdf_variations = np.clip(pdf_variations, a_min=None, a_max=max_value)
    pdf_variations = pdf_variations / pdf_variations[:, :1]

    # Calculate the mean and standard deviation across replicas per event
    mean = np.mean(pdf_variations, axis=1)
    std = np.std(pdf_variations, axis=1)

    # Calculate the up/down variations across events
    variation_up = mean + std
    variation_down = mean - std

    return __add_weight_variations(events, variation_up, variation_down, "PDF", computed_nominal_weights=None, multiply_by_pu_weights=multiply_by_pu_weight)

def apply_ps_variations(events,is_nano=False,ps_type="ISR", multiply_by_pu_weight=False):
    """Calculate the PS up/down variations for ISR and FSR.
    """

    if is_nano:
        if ps_type == "ISR":
            variation_up = events.PSWeight[:, 0]
            variation_down = events.PSWeight[:, 2]
        if ps_type == "FSR":
            variation_up = events.PSWeight[:, 1]
            variation_down = events.PSWeight[:, 3]
       
    else:
        raise NotImplementedError()

    return __add_weight_variations(events, variation_up, variation_down, f"PS{ps_type}", computed_nominal_weights=None, multiply_by_pu_weights=multiply_by_pu_weight)

    
   


def apply_pu_variations(events, year, pfnano_sys_file=None , is_nano=False, multiply_by_pu_weight=False):

    # Normalize the array of pdf weights by the first entry
    if is_nano:
       pu_nTrueInt = events.Pileup_nTrueInt
       variations_factory = load(pfnano_sys_file)

    else:
        raise NotImplementedError()
    
    pu_nom, pu_up, pu_down   = variations_factory["get_pu_weight"](year, pu_nTrueInt)

    return __add_weight_variations(events,pu_up, pu_down, "PU", computed_nominal_weights=pu_nom, multiply_by_pu_weights=multiply_by_pu_weight)


###############################
####### PFNano section ########
###############################

def apply_variation_pfnano(events, variation, year, run, pfnano_sys_file):
    
    #if variation is None: return events
    #if variatyion is None set variation = "nominal"
    if variation is None:
        variation = "nominal"

        
    if variation in ["nominal","jec_up", "jec_down", "jer_up", "jer_down","SVJjec_up", "SVJjec_down"]:
        jerc_cache = cachetools.Cache(np.inf)
        #load the JEC/JER variations from the pfnano file
        if pfnano_sys_file is not None:
            jerc_variations = load(pfnano_sys_file)
        #CZZ: first compute the JEC variations for AK4 and AK8 jets
        for radius in [4, 8]:
            jet_coll = "Jet" if radius == 4 else "FatJet"

            if variation in ["nominal"]:  
                
                #here apply JECs and JERs (undo the original jec and apply the updated ones)
                corrected_JERC_jets  = apply_jercs_PFNano(
                        events,
                        year,
                        run,
                        jet_coll,
                        jerc_variations,
                        jerc_cache,
                        ) 
                
                
                if radius == 4:
                    #here unpack corrections
                    corrected_MET_updated_JECs = propagate_jecs_to_MET_PFNano(
                        events,
                        year,
                        run,
                        jet_coll,
                        jerc_variations,
                        jerc_cache,
                        ) 

                    met_ptcorr, met_phicorr = corrected_MET_updated_JECs
                    events["MET_pt"] = met_ptcorr
                    events["MET_phi"] = met_phicorr

                    #Here propagate the corrections to METSignificance
                    met_sig_corr_nom = propagate_jecs_to_METSig_PFNano(events,
                                                year,
                                                run,
                                                jet_coll,
                                                jerc_variations,
                                                jerc_cache,
                                                )
                    
                    events["MET_significance"] = met_sig_corr_nom


                #unpack corrected_JERC_jets
                jerc_corr_pt, jerc_corr_eta, jerc_corr_phi, jerc_corr_mass = corrected_JERC_jets

                #adding the JERC varied jets to the events
                jerc_permutation = ak.argsort(jerc_corr_pt, ascending=False)
                jerc_corr_pt = jerc_corr_pt[jerc_permutation]
                jerc_corr_eta = jerc_corr_eta[jerc_permutation]
                jerc_corr_phi = jerc_corr_phi[jerc_permutation]
                jerc_corr_mass = jerc_corr_mass[jerc_permutation]
                events[f"{jet_coll}_pt"] = jerc_corr_pt
                events[f"{jet_coll}_eta"] = jerc_corr_eta
                events[f"{jet_coll}_phi"] = jerc_corr_phi
                events[f"{jet_coll}_mass"] = jerc_corr_mass
                jerc_pf_cand_jet_idx = events[f"{jet_coll}PFCands_jetIdx"]
                jerc_sorted_pf_cand_jet_idx = ak.Array([p[idx] for idx, p in zip(jerc_pf_cand_jet_idx, jerc_permutation)])
                events[f"{jet_coll}PFCands_jetIdx"] = jerc_sorted_pf_cand_jet_idx


            if variation in ["jec_up", "jec_down", "SVJjec_up", "SVJjec_down","jer_up", "jer_down"]: 
                if variation in ["SVJjec_up", "SVJjec_down"]:

                    #need to apply to the nominal jets the updated corrections
                    corrected_JERC_jets  = apply_jercs_PFNano(
                        events,
                        year,
                        run,
                        jet_coll,
                        jerc_variations,
                        jerc_cache,
                        ) 
                    

                    #apply MET-T1 correction to MET
                    if radius == 4:
                        #here unpack corrections
                        corrected_MET_updated_JECs = propagate_jecs_to_MET_PFNano(
                            events,
                            year,
                            run,
                            jet_coll,
                            jerc_variations,
                            jerc_cache,
                            ) 

                        met_ptcorr, met_phicorr = corrected_MET_updated_JECs
                        events["MET_pt"] = met_ptcorr
                        events["MET_phi"] = met_phicorr

                    #unpack corrected_JERC_jets
                    jerc_corr_pt, jerc_corr_eta, jerc_corr_phi, jerc_corr_mass = corrected_JERC_jets

                    #adding the JEC varied jets to the events
                    jerc_permutation = ak.argsort(jerc_corr_pt, ascending=False)
                    jerc_corr_pt = jerc_corr_pt[jerc_permutation]
                    jerc_corr_eta = jerc_corr_eta[jerc_permutation]
                    jerc_corr_phi = jerc_corr_phi[jerc_permutation]
                    jerc_corr_mass = jerc_corr_mass[jerc_permutation]
                    events[f"{jet_coll}_pt"] = jerc_corr_pt
                    events[f"{jet_coll}_eta"] = jerc_corr_eta
                    events[f"{jet_coll}_phi"] = jerc_corr_phi
                    events[f"{jet_coll}_mass"] = jerc_corr_mass
                    jerc_pf_cand_jet_idx = events[f"{jet_coll}PFCands_jetIdx"]
                    jerc_sorted_pf_cand_jet_idx = ak.Array([p[idx] for idx, p in zip(jerc_pf_cand_jet_idx, jerc_permutation)])
                    events[f"{jet_coll}PFCands_jetIdx"] = jerc_sorted_pf_cand_jet_idx


                    #compute the custom jecs variations
                    corrected_jets_met_jercs = calc_custom_svj_jes_variations_PFNano(
                        events,
                        year,
                        run,
                        jet_coll,
                        variation,
                    )

                if variation in ["jer_up", "jer_down"]: 
                    corrected_jets_met_jercs = calc_jerc_variations_PFNano(
                        events,
                        year,
                        run,
                        jet_coll,
                        jerc_variations,
                        variation,
                        jerc_cache,
                    )


                    if radius == 4:
                        #here unpack corrections
                        corrected_MET_updated_JECs = propagate_jecs_to_MET_PFNano(
                            events,
                            year,
                            run,
                            jet_coll,
                            jerc_variations,
                            jerc_cache,
                            ) 

                        met_ptcorr, met_phicorr = corrected_MET_updated_JECs
                        events["MET_pt"] = met_ptcorr
                        events["MET_phi"] = met_phicorr

                        #Here propagate the corrections to METSignificance
                        met_sig_corr_nom = propagate_jecs_to_METSig_PFNano(events,
                                                    year,
                                                    run,
                                                    jet_coll,
                                                    jerc_variations,
                                                    jerc_cache,
                                                    )
                        
                        events["MET_significance"] = met_sig_corr_nom



                if variation in ["jec_up", "jec_down"]: 
                    
                    corrected_jets_met_jercs = calc_jerc_variations_PFNano(
                        events,
                        year,
                        run,
                        jet_coll,
                        jerc_variations,
                        variation,
                        jerc_cache,
                    )


                
                if variation in ["jec_up", "jec_down", "jer_up", "jer_down", "SVJjec_up", "SVJjec_down"]:
                    #here unpack corrections / here do not propagate the JER variations to MET
                    if (radius == 4) and (variation in ["jec_up", "jec_down", "SVJjec_up", "SVJjec_down"]):
                        corr_pt, corr_eta, corr_phi, corr_mass, met_ptcorr, met_phicorr, custom_corr_factor = corrected_jets_met_jercs
                    else:
                        corr_pt, corr_eta, corr_phi, corr_mass, _, _, custom_corr_factor = corrected_jets_met_jercs
                

                    #adding the JERC varied jets to the events
                    permutation = ak.argsort(corr_pt, ascending=False)
                    corr_pt = corr_pt[permutation]
                    corr_eta = corr_eta[permutation]
                    corr_phi = corr_phi[permutation]
                    corr_mass = corr_mass[permutation]
                    custom_corr_factor = custom_corr_factor[permutation]
                    events[f"{jet_coll}_pt"] = corr_pt
                    events[f"{jet_coll}_eta"] = corr_eta
                    events[f"{jet_coll}_phi"] = corr_phi
                    events[f"{jet_coll}_mass"] = corr_mass
                    #add corr factor for custom SVJ JEC
                    events[f"{jet_coll}_corrSVJJEC"] = custom_corr_factor
                    pf_cand_jet_idx = events[f"{jet_coll}PFCands_jetIdx"]
                    sorted_pf_cand_jet_idx = ak.Array([p[idx] for idx, p in zip(pf_cand_jet_idx, permutation)])
                    events[f"{jet_coll}PFCands_jetIdx"] = sorted_pf_cand_jet_idx
                
                    #add the MET corrections, here correct MET only if the variation is JEC (not JER)
                    if (radius == 4) and (variation in ["jec_up", "jec_down", "SVJjec_up", "SVJjec_down"]):
                        events["MET_pt"] = met_ptcorr
                        events["MET_phi"] = met_phicorr

                        #Here propagate the corrections to METSignificance
                        #CZZ: MET must come from jerc varied collections, otherwise nominal when running the unclustered energy variation
                        metSigCalc = MetSignificanceCalculator(events,
                                                year,
                                                run,
                                                ) 
                        events["MET_significance"] = metSigCalc.getSignificance()                 
        
    
    if variation in ["unclEn_up", "unclEn_down"]:

        #then do the same for MET with coffea libraries
        jerc_cache = cachetools.Cache(np.inf)
        #load the JEC/JER variations from the pfnano file
        if pfnano_sys_file is not None:
            jerc_variations = load(pfnano_sys_file)

        corrected_unclEn_met = calc_unclustered_met_variations_PFNano(
            events,
            year,
            run,
            jet_coll="Jet",
            jerc_variations=jerc_variations,
            variation=variation,
            jerc_cache=jerc_cache,
        )

        #unpack corrected_unclEn_met
        met_ptcorr, met_phicorr = corrected_unclEn_met

        #Here propagate the corrections to METSignificance
        met_sig_corr_nom = propagate_jecs_to_METSig_PFNano(events,
                                    year,
                                    run,
                                    jet_coll="Jet",
                                    jerc_variations=jerc_variations,
                                    jerc_cache=jerc_cache,
                                    make_unclustered_En_var = True,
                                    variation = variation,
                                    )
        
        events["MET_significance"] = met_sig_corr_nom

        #add the MET corrections
        events["MET_pt"] = met_ptcorr
        events["MET_phi"] = met_phicorr


        #CZZ: first compute the JEC variations for AK4 and AK8 jets
        for radius in [4, 8]:
            jet_coll = "Jet" if radius == 4 else "FatJet"
            
            #here apply JEC by default (undo JECs, and reapply them)
            corrected_JEC_jets  = apply_jercs_PFNano(
                    events,
                    year,
                    run,
                    jet_coll,
                    jerc_variations,
                    jerc_cache,
                    ) 
            
            #unpack corrected_JER_jets
            jec_corr_pt, jec_corr_eta, jec_corr_phi, jec_corr_mass = corrected_JEC_jets

            #adding the JER varied jets to the events
            jec_permutation = ak.argsort(jec_corr_pt, ascending=False)
            jec_corr_pt = jec_corr_pt[jec_permutation]
            jec_corr_eta = jec_corr_eta[jec_permutation]
            jec_corr_phi = jec_corr_phi[jec_permutation]
            jec_corr_mass = jec_corr_mass[jec_permutation]
            events[f"{jet_coll}_pt"] = jec_corr_pt
            events[f"{jet_coll}_eta"] = jec_corr_eta
            events[f"{jet_coll}_phi"] = jec_corr_phi
            events[f"{jet_coll}_mass"] = jec_corr_mass
            jec_pf_cand_jet_idx = events[f"{jet_coll}PFCands_jetIdx"]
            jec_sorted_pf_cand_jet_idx = ak.Array([p[idx] for idx, p in zip(jec_pf_cand_jet_idx, jec_permutation)])
            events[f"{jet_coll}PFCands_jetIdx"] = jec_sorted_pf_cand_jet_idx

         

    return events
