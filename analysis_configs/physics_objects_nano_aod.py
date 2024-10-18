import awkward as ak
from utils.Logger import *


import numpy as np
import awkward as ak
from coffea.nanoevents.methods import vector
# Needed so that ak.zip({"pt": [...], "eta": [...], "phi": [...], "mass": [...]},
#                         with_name="PtEtaPhiMLorentzVector")
# is understood as a PtEtaPhiMLorentzVector from coffea.nanoevents.methods.vector
ak.behavior.update(vector.behavior)


def rapidity(obj):
    return np.arcsinh(np.sinh(obj.eta)/np.sqrt(1+obj.mass**2/obj.pt**2))


def momentum(obj):
    return np.sqrt( obj.pt**2 + obj.pz**2 )


def get_collection(events, collection_name, variables=[], zero_variables=[]):
    """Get a collection.

    Args:
        events (awkward.Array): the Events TTree opened with uproot.
        collection_name (str): e.g. "MET", "FatJet", ...
        variables (list[str]): List of branches to read, e.g. pt, eta, ...
        zero_variables (list[str]): List of branches to create and fill with 0.

    Returns:
        awkward.Array: ak array with field pt, eta, phi, mass,
            that can be used as a coffea PtEtaPhiMLorentzVector.
    """

    mandatory_variables = ["pt", "eta", "phi", "mass"]
    for variable in mandatory_variables:
        if variable not in variables + zero_variables:
            log.critical(f"{variable} not in variables to read for collection {collection_name}.")
            log.critical(f"Mandatory varables are: {mandatory_variables}")
            exit(1)

    found_branch = False
    for field in events.fields:
        if collection_name in field:
            found_branch = True
            break

   
    if not found_branch:
        log.warning(f"{collection_name} not in events.fields!")
        collection = None

    else:
        #collection_branch = eval("events." + collection_name)
        collection_branch = collection_name + "_"
        zero = ak.zeros_like(events[collection_branch + variables[0]])
        collection = ak.zip(
            {
                **{variable: events[collection_branch + variable]
                   for variable in variables
                },
                **{zero_variable: zero
                   for zero_variable in zero_variables
                },
             },
             with_name="PtEtaPhiMLorentzVector",
        )
    collection = ak.with_field(collection, rapidity(collection), where="rapidity")     
    collection = ak.with_field(collection, momentum(collection), where="momentum")
    return collection


def get_collection_size(collection, attr="pt"):
    """Get the length of a collection per event.

    Args:
        collection (awkward.highlevel.Array): ak array
        attr (str, optional, default=pt): the name of a field of the ak array
    """

    return ak.count(getattr(collection, attr), axis=1)



def get_jets(events, collection_name, additional_variables=[]):
    variables = [
        "pt", "eta", "phi", "mass"
    ]


    variables += additional_variables
    jets = get_collection(events, collection_name, variables)
    return jets



def get_jet_pf_cands(events, jet_collection_name, pf_cands_collection_name, additional_variables=[]):
    """Get pf candidates clustered in jets for a given jet collection.

    The array returned has fields from the pf candidate collection and a
    jetIdx field, such that one knows to which jet belongs a pf candidate.
    There can be duplicate if a candidate belongs to several jets.
    If so, the reptitions of the candidate will have same variables but
    different jetIdx.

    Args:
        events (awkward.highlevel.Array): the Events TTree opened with uproot.
        jet_collection_name (str): usually "Jet" or "FatJet"
        pf_cands_collection_name (str): usually "PFCands"

    Returns:
        awkward.highlevel.Array: ak array with field jetIdx, pt, eta, phi, mass,
            charge and pdgId that can be used as a coffea
            PtEtaPhiMLorentzVector.
    """

    jet_pf_cands_collection_name = jet_collection_name + pf_cands_collection_name


    found_pf_cands_collection = False
    found_jet_pf_cands_collection = False

    for field in events.fields:
        if jet_pf_cands_collection_name in field:
            found_jet_pf_cands_collection = True
        if pf_cands_collection_name in field:
            found_pf_cands_collection = True
        
    if not found_pf_cands_collection and not found_jet_pf_cands_collection:
        pf_cands = None

    else:
        jet_pf_cands_idx_info = jet_pf_cands_collection_name + "_"
        jet_pfcands_idx = events[jet_pf_cands_idx_info + "pFCandsIdx"]
        variables = ["pt", "mass", "eta", "phi", "charge", "pdgId"]
        variables += additional_variables

        #check if first element of jet_pfcands_idx is an integer
        if not isinstance(jet_pfcands_idx[0][0], int):
            #convert to integer
            jet_pfcands_idx = ak.values_astype(jet_pfcands_idx, 'int64')

        pf_cands = ak.zip(
            {
                **{"jetIdx": events[jet_pf_cands_idx_info + "jetIdx"]},
                **{variable: events[pf_cands_collection_name + "_" + variable][jet_pfcands_idx]
                   for variable in variables
                  },
            },
            with_name="PtEtaPhiMLorentzVector",
        )
        pf_cands = ak.with_field(pf_cands, rapidity(pf_cands), where="rapidity")

    return pf_cands


def get_pf_cands(events, collection_name, additional_variables=[]):
    """Get pf candidates collection.

    The array returned does not contain any information about the jet the
    pf candidate belongs to. If this collection was defined such that it stores
    all pf candidates in the events, then this function returns all pf
    candidates in the event.

    Args:
        events (awkward.highlevel.Array): the Events TTree opened with uproot.
        collection_name (str): usually "PFCands"

    Returns:
        awkward.highlevel.Array: ak array with field pt, eta, phi, mass,
            charge and pdgId that can be used as a coffea
            PtEtaPhiMLorentzVector.
    """

    variables = ["pt", "mass", "eta", "phi", "charge", "pdgId", "d0", "dz", "d0Err", "dzErr"]
    variables += additional_variables
    pf_cands = get_collection(events, collection_name, variables)

    return pf_cands


