import numpy as np
import awkward as ak
import utils.awkward_array_utilities as akUtl
import numba as nb
from utils.Logger import *

NAN_VALUE = -1


def delta_phi(obj1, obj2):
    return obj1.delta_phi(obj2)

def delta_rapidity(obj1, obj2):
    return obj1.rapidity - obj2.rapidity

def delta_r_rapidity(obj1, obj2):
    
    dy = delta_rapidity(obj1, obj2)
    dphi = delta_phi(obj1, obj2)

    return np.sqrt(dy**2 + dphi**2)

def delta_r_pseudorapidity(obj1, obj2):
    
    return obj1.delta_r(obj2)


def delta_r(obj1, obj2, use_rapidity=False):
    if use_rapidity:
        return delta_r_rapidity(obj1, obj2)
    else:
        return delta_r_pseudorapidity(obj1, obj2)



def calculate_generalized_angularity(constituents, jets, jet_radius, beta, kappa, nan_value=NAN_VALUE):
    """Calculate generalized angularity.

    Args:
        constituents (awkward.Array or None):
            2D ak array of constituents per jet, where axis 0 is the jet axis,
            axis 1 is the constituents axis with fields pt, eta, phi, mass,
            and with name PtEtaPhiMLorentzVector.
            Can be None if kappa=0.
        jets (awkward.Array or None):
            1D ak array of jets, with fields pt, eta, phi and mass and with
            name PtEtaPhiMLorentzVector.
            Can be None if beta=0.
        beta (float): Angular distance exponent.
        kappa (float): Transverse momentum fraction exponent.
        nan_value (float, optional):
            Value to use when the angularity cannot be computed (should not
            happen!)

    Returns:
        awkward.Array

    Examples:
        >>> from tests.utilities.variablesComputation.awkwardArray.docTestHelper import get_data
        >>> from builtinUtilities import pretty_printer
        >>> constituents, jets = get_data()
        >>> pretty_printer(calculate_generalized_angularity(constituents, jets, 1, 0, 2))
        [0.500, -1.000, 0.297]
        >>> pretty_printer(calculate_generalized_angularity(constituents, jets, 1, 1, 1))
        [0.071, -1.000, 0.168]
    """

    if beta == 0 and kappa == 0:
        generalized_angularity = ak.num(constituents.pt, axis=1)

    else:
        if beta == 0:
            angular_term = 1.
        else:
            jet_broadcasted = akUtl.broadcast(jets, constituents)[0]
            delta_r = delta_r(constituents, jet_broadcasted)
            angular_term = delta_r**beta
    
        if kappa == 0:
            sum_constituents_pt = 1.
            momentum_term = 1.
        else:
            sum_constituents_pt = ak.sum(constituents.pt, axis=1)
            momentum_term = constituents.pt**kappa

        numerator = ak.sum(momentum_term * angular_term, axis=1)
        denominator = sum_constituents_pt**kappa * jet_radius**beta
        generalized_angularity = akUtl.divide_ak_arrays(numerator, denominator, division_by_zero_value=nan_value)
        generalized_angularity = ak.fill_none(generalized_angularity, nan_value)
    
    return generalized_angularity


def calculate_multiplicity(constituents, nan_value=NAN_VALUE):
    """Calculate number of constituents. See docstring of calculate_generalized_angularity."""

    return calculate_generalized_angularity(constituents, None, 1., 0, 0, nan_value=nan_value)


def calculate_delta_phi_with_met(jets, met):
    """Returns delta phi between jets and MET.

    Args:
        jets (ak.Array):
            Ak array where axis 0 is the event axis, axis 1 is the object axis
            with fields pt, eta, phi, mass and with name PtEtaPhiMLorentzVector.
        met (ak.Array):
            Ak array where axis 0 is the event axis with fields pt and phi and
            with name PtEtaPhiMLorentzVector.
    
    Returns:
        ak.Array
    """

    met = ak.broadcast_arrays(met, jets)[0]
    delta_phi = jets.delta_phi(met)

    return delta_phi


def calculate_lund_jet_plane_z_with_met(jets, met):
    """
    Args:
        jets (ak.Array):
            Ak array where axis 0 is the event axis, axis 1 is the object axis
            with fields pt, eta, phi, mass and with name PtEtaPhiMLorentzVector.
        met (ak.Array):
            Ak array where axis 0 is the event axis with fields pt and phi and
            with name PtEtaPhiMLorentzVector.
    
    Returns:
        ak.Array
    """

    met = ak.broadcast_arrays(met, jets)[0]
    z = ak.min([jets.pt, met.pt], axis=0) / (jets.pt + met.pt)

    return z


def calculate_invariant_mass_with_met(jets, met):
    """
    Args:
        jets (ak.Array):
            Ak array where axis 0 is the event axis, axis 1 is the object axis
            with fields pt, eta, phi, mass and with name PtEtaPhiMLorentzVector.
        met (ak.Array):
            Ak array where axis 0 is the event axis with fields pt and phi and
            with name PtEtaPhiMLorentzVector.
    
    Returns:
        ak.Array
    """

    met = ak.broadcast_arrays(met, jets)[0]
    term1 = np.sqrt(jets.mass**2 + jets.pt**2)
    term2 = np.cos(met.phi - jets.phi) * jets.pt
    mt = np.sqrt(jets.mass**2 + 2 * met.pt * (term1 - term2))

    return mt


def calculate_HEF(charged, constituents, nan_value=NAN_VALUE):
    """Calculate energy fraction for charged or neutral hadrons.

    Args:
        charged (bool)
        constituents (awkward.Array):
            2D ak array of constituents per jet, where axis 0 is the jet axis,
            axis 1 is the constituents axis with fields pdgId.
    """

    hadrons_constituents = constituents[akUtl.is_in_list(constituents.pdgId, pdg_id("hadrons"))]

    if charged:
        selected_hadrons_constituents = hadrons_constituents[(hadrons_constituents.charge != 0)]
    else:
        selected_hadrons_constituents = hadrons_constituents[(hadrons_constituents.charge == 0)]

    numerator = ak.sum(selected_hadrons_constituents.energy, axis=1)
    denominator = ak.sum(constituents.energy, axis=1)
    
    #energy_fraction = akUtl.divide_ak_arrays(numerator, denominator, division_by_zero_value=nan_value)
    #energy_fraction = ak.fill_none(energy_fraction, nan_value)
    
    energy_fraction = numerator/denominator
    return energy_fraction


def calculate_chHEF(constituents, nan_value=NAN_VALUE):
    """Return energy fraction for charged hadrons. See doctring of calculate_HEF."""

    return calculate_HEF(True, constituents, nan_value=nan_value)


def calculate_neHEF(constituents, nan_value=NAN_VALUE):
    """Return energy fraction for neutral hadrons. See doctring of calculate_HEF."""

    return calculate_HEF(False, constituents, nan_value=nan_value)


def calculate_chargedparticle_multiplicity(constituents, particle):
    """Calculate multiplicity for charged particles (electrons, muons, hadrons).

    Args:
        constituents (awkward.Array):
            2D ak array of constituents per jet, where axis 0 is the jet axis,
            axis 1 is the constituents axis with fields pdgId.
        particle (str): "electron", "muon"
    """

    pdg_ids = pdg_id(particle)
    if not isinstance(pdg_ids, list): pdg_ids = [pdg_ids]

    if pdg_ids == pdg_id("hadrons"):
        hadrons_constituents = constituents[akUtl.is_in_list(constituents.pdgId, pdg_ids)]
        selected_constituents = hadrons_constituents[(hadrons_constituents.charge != 0)]
    else:
        selected_constituents = constituents[akUtl.is_in_list(constituents.pdgId, pdg_ids)]

    return ak.num(selected_constituents.pt, axis=1)

def calculate_electron_multiplicity(constituents):
    return calculate_chargedparticle_multiplicity(constituents, particle = "electron")

def calculate_muon_multiplicity(constituents):
    return calculate_chargedparticle_multiplicity(constituents, particle = "muon")

def calculate_chargedhadron_multiplicity(constituents):
    return calculate_chargedparticle_multiplicity(constituents, particle = "hadrons")


def calculate_energy_fraction(constituents, particle, nan_value=NAN_VALUE):
    """Calculate energy fraction for charged or neutral hadrons.

    Args:
        constituents (awkward.Array):
            2D ak array of constituents per jet, where axis 0 is the jet axis,
            axis 1 is the constituents axis with fields pdgId.
        particle (str): "electron", "muon", "photon"
    """

    pdg_ids = pdg_id(particle)
    if not isinstance(pdg_ids, list): pdg_ids = [pdg_ids]
    selected_constituents = constituents[akUtl.is_in_list(constituents.pdgId, pdg_ids)]

    numerator = ak.sum(selected_constituents.energy, axis=1)
    denominator = ak.sum(constituents.energy, axis=1)
    energy_fraction = numerator/denominator

    return energy_fraction

def calculate_electron_energy_fraction(constituents, nan_value=NAN_VALUE):
    return calculate_energy_fraction(constituents, particle="electron", nan_value=nan_value)

def calculate_muon_energy_fraction(constituents, nan_value=NAN_VALUE):
    return calculate_energy_fraction(constituents, particle="muon", nan_value=nan_value)

def calculate_photon_energy_fraction(constituents, nan_value=NAN_VALUE):
    return calculate_energy_fraction(constituents, particle="photon", nan_value=nan_value)


def pdg_id(particle):

    pdg_id_dict = {
	## Quarks
	"d": 1,
	"u": 2,
	"s": 3,
	"c": 4,
	"b": 5,
	"t": 6,

	## Leptons
	"e-"    : 11,
	"e+"    : -11,
	"mu-"   : 13,
	"mu+"   : -13,
	"tau-"  : 15,
	"tau+"  : -15,
	"nu_e"  : 12,
	"nu_mu" : 14,
	"nu_tau": 16,

	## Gauge and Higgs bosons
	"g"     : 21,
	"photon": 22,
	"Z"     : 23,
	"W+"    : 24,
	"W-"    : -24,
	"H"     : 25,

	## Light I=1 mesons
	"pi0" : 111,
	"pi+" : 211,
	"pi-" : -211,
	"rho0": 113,
	"rho+": 213,
	"rho-": -213,

	## Strange mesons
	"KL0": 130,
    }

    d = pdg_id_dict   # shorthand

    pdg_id_dict["quarks"] = [d["d"], d["u"], d["s"], d["c"], d["b"], d["t"]]

    pdg_id_dict["electron"] = [d["e-"], d["e+"]]
    pdg_id_dict["muon"] = [d["mu-"], d["mu+"]]
    pdg_id_dict["charged_leptons"] = [d["e-"], d["e+"], d["mu-"], d["mu+"], d["tau-"], d["tau+"]]
    pdg_id_dict["neutral_leptons"] = [d["nu_e"], d["nu_mu"], d["nu_tau"]]
    pdg_id_dict["leptons"] =  pdg_id_dict["charged_leptons"] + pdg_id_dict["neutral_leptons"]

    pdg_id_dict["light_mesons"] = [d["pi0"], d["pi+"], d["pi-"], d["rho0"], d["rho+"], d["rho-"]]
    pdg_id_dict["strange_mesons"] = [d["KL0"]]
    pdg_id_dict["mesons"] = d["light_mesons"] + d["strange_mesons"]

    pdg_id_dict["hadrons"] = d["mesons"]


    return pdg_id_dict[particle]

def count_constituents(jet_indices, n_jets):
    """Count constituents per jet, with same jagged structure as the jet collection.

    Args:
        jet_indices (ak.Array): Jet constituents per event
        n_jets (ak.Array): Number of jets per event
    """

    #@nb.jit
    def __count_constituents(builder, jet_indices, n_jets):
        for indices_this_event, n_jet in zip(jet_indices, n_jets):
            builder.begin_list()
            for i_jet in range(n_jet):
                counter = 0
                for index in indices_this_event:
                    if index == i_jet:
                        counter += 1
                builder.append(counter)
            builder.end_list()
        return builder

    builder = ak.ArrayBuilder()
    return __count_constituents(builder, jet_indices, n_jets).snapshot()


def make_constituents_per_jet(constituents_per_event, n_jets, jet_idx_field_name="jetIdx"):
    """Make ak array of constituents per jet, from per event ak array.
    
    Args:
        constituents_per_event (ak.Array): 2D ak array with fields of jet
            constituents per event, with fields. The jet in which a given
            constituent belongs is indicated via a jet index field.
        n_jets (ak.Array): 1D ak array with number of jets per event.
        jet_idx_field_name (str, optional, default="jetIdx"):
            Jet index field name in the constituents_per_event ak array.
    """

    constituents_per_event = akUtl.sort_array_with_fields(constituents_per_event, jet_idx_field_name, ascending=True)
    counts = count_constituents(constituents_per_event[jet_idx_field_name], n_jets)

    flat_counts = ak.flatten(counts)
    flat_constituents = ak.flatten(constituents_per_event)
    constituents_per_jet = ak.unflatten(flat_constituents, flat_counts, axis=0)

    return constituents_per_jet, counts