import numpy as np 
import awkward as ak

# convert phi values into spherical coordinate?
def phi_(x,y):
    phi = np.arctan2(y,x)
    return ak.where(phi < 0, phi + 2*np.pi, phi)

def delta_phi(jet_phi,met_phi):
    phi1 = phi_( np.cos(jet_phi), np.sin(jet_phi) )
    phi2 = phi_( np.cos(met_phi), np.sin(met_phi) )
    dphi = phi1 - phi2
    dphi_edited = ak.where(dphi < -np.pi, dphi + 2*np.pi, dphi)
    dphi_edited = ak.where(dphi_edited > np.pi, dphi_edited - 2*np.pi, dphi_edited)
    return dphi_edited

def delta_eta(eta0,eta1):
  return (eta0 - eta1)

def delta_R(eta0,eta1,phi0,phi1):
    dp = delta_phi(phi0,phi1)
    deta = delta_eta(eta0,eta1)
    deltaR2 = deta * deta + dp * dp
    return np.sqrt(deltaR2)

def run_jet_constituent_matching(events, orig_jets):
    jets = orig_jets[:]
    jet_constituents = events.JetsConstituents[:]
    jetsAK8_constituents_index = jets.constituentsIndex
    indices = ak.flatten(jetsAK8_constituents_index,axis=-1)
    jet_constituents_for_jets = jet_constituents[indices]
    jet_const_jet_array = ak.unflatten(jet_constituents_for_jets,ak.flatten(ak.num(jetsAK8_constituents_index,axis=-1)),axis=1)
    jets["Constituents"] = jet_const_jet_array
    return jets
