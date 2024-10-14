import numpy as np
import awkward as ak
import utils.awkward_array_utilities as akUtl
import utils.physics_utilities as phUtl

NAN_VALUE = -1

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
            delta_r = vecUtl.delta_r(constituents, jet_broadcasted)
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

    hadrons_constituents = constituents[akUtl.is_in_list(constituents.pdgId, phUtl.pdg_id("hadrons"))]

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

    pdg_ids = phUtl.pdg_id(particle)
    if not isinstance(pdg_ids, list): pdg_ids = [pdg_ids]

    if pdg_ids == phUtl.pdg_id("hadrons"):
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

    pdg_ids = phUtl.pdg_id(particle)
    if not isinstance(pdg_ids, list): pdg_ids = [pdg_ids]
    selected_constituents = constituents[akUtl.is_in_list(constituents.pdgId, pdg_ids)]

    numerator = ak.sum(selected_constituents.energy, axis=1)
    denominator = ak.sum(constituents.energy, axis=1)
    #energy_fraction = akUtl.divide_ak_arrays(numerator, denominator, division_by_zero_value=nan_value)
    energy_fraction = numerator/denominator

    return energy_fraction

def calculate_electron_energy_fraction(constituents, nan_value=NAN_VALUE):
    return calculate_energy_fraction(constituents, particle="electron", nan_value=nan_value)

def calculate_muon_energy_fraction(constituents, nan_value=NAN_VALUE):
    return calculate_energy_fraction(constituents, particle="muon", nan_value=nan_value)

def calculate_photon_energy_fraction(constituents, nan_value=NAN_VALUE):
    return calculate_energy_fraction(constituents, particle="photon", nan_value=nan_value)