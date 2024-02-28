import numpy as np
import awkward as ak


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

