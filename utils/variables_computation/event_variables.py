import numpy as np
import awkward as ak


def calculate_number_of_objects(physics_objects):
    """Calculate number of a physics objects in all events.

    Args:
        physics_objects (awkward.Array): Jagged ak array where axis 0 is the
            event axis and axis 1 is the object axis.

    Returns:
        awkward.Array
    """

    return ak.num(physics_objects, axis=1)


def __get_pair_of_objects(physics_objects, indices=(0, 1)):
    """Returns a tuple of masked ak arrays corresponding to physics objects with required indices.

    Args:
        physics_objects (ak.Array)
        indices (tuple[int], optional): By default 2 leading objects

    Returns:
        tuple[ak.Array, ak.Array]
    """

    n_objects = calculate_number_of_objects(physics_objects)
    physics_objects_0 = ak.mask(physics_objects, n_objects >= indices[0] + 1)
    physics_objects_1 = ak.mask(physics_objects, n_objects >= indices[1] + 1)
    obj0 = physics_objects_0[:, indices[0]]
    obj1 = physics_objects_1[:, indices[1]]

    return obj0, obj1


def __get_objects_from_distance_with_met(physics_objects, met, mode="min"):
    """

    Args:
        physics_objects (ak.Array)
        n_objects (ak.Array, optional): Number of objects per event
        mode (str): "min", "max"

    Returns:
        ak.Array
    """

    assert mode == "min" or mode == "max"

    met = ak.broadcast_arrays(met, physics_objects)[0]
    distances = abs(physics_objects.delta_phi(met))
    indices = ak.argsort(distances, axis=1, ascending=mode=='min')
    ordered_physics_objects = physics_objects[indices]

    return ordered_physics_objects


def calculate_delta_eta(
        physics_objects,
        indices=(0, 1),
        nan_value=-9999,
        absolute_value=False,
    ):
    """Calculate delta eta between two physics objects.

    Args:
        physics_objects (awkward.Array): 
            Ak array where axis 0 is the event axis, axis 1 is the object axis
            with field eta and with name PtEtaPhiMLorentzVector.
            Physics objects can be jets, leptons etc...
        indices (tuple[int], optional):
            The indices of the object for which to compute delta eta.
            By default will do it the two leading objects.
        nan_value (float, optional, default=-9999):
            Value to use when the event has less objects than required.
            If None, Nones will not be replaced.
        absolute_value (bool, optional, default=False):
            Whether to return absolute or signed values

    Returns:
        awkward.Array
    """

    obj0, obj1 = __get_pair_of_objects(physics_objects, indices)
    delta_eta = obj0.eta - obj1.eta
    if absolute_value:
        delta_eta = abs(delta_eta)
    if nan_value is not None:
        delta_eta = ak.fill_none(delta_eta, nan_value)

    return delta_eta


def calculate_delta_phi(
        physics_objects,
        indices=(0, 1),
        nan_value=-9999,
        absolute_value=False,
    ):
    """Calculate delta phi between two physics objects.

    Args:
        physics_objects (awkward.Array): 
            Ak array where axis 0 is the event axis, axis 1 is the object axis
            with field phi and with name PtEtaPhiMLorentzVector.
            Physics objects can be jets, leptons etc...
        indices (tuple[int], optional):
            The indices of the object for which to compute delta eta.
            By default will do it the two leading objects.
        nan_value (float, optional, default=-9999):
            Value to use when the event has less objects than required.
            If None, Nones will not be replaced.
        absolute_value (bool, optional, default=False):
            Whether to return absolute or signed values

    Returns:
        awkward.Array
    """

    obj0, obj1 = __get_pair_of_objects(physics_objects, indices)
    delta_phi = obj0.delta_phi(obj1)
    if absolute_value:
        delta_phi = abs(delta_phi)
    if nan_value is not None:
        delta_phi = ak.fill_none(delta_phi, nan_value)

    return delta_phi


def delta_phi_dark_quark(phi1, phi2):
    """
    Calculate delta phi between two physics objects. 
    This function can handle both single values and arrays.
    Used to calculate delta phi between a dark quark and all jets. 
    (not sure if possible with existing delta phi function)
    
    Args:
        phi1, phi2 (array-like or float): 
            The two angles for which delta phi is to be calculated.
            Can be single values or arrays of the same shape.
    Returns:
        array-like or float
    """
    
    dphi = phi1 - phi2
    dphi = np.where(dphi > np.pi, dphi - 2 * np.pi, dphi)
    dphi = np.where(dphi < -np.pi, dphi + 2 * np.pi, dphi)
    return abs(dphi)


def delta_r_dark_quark(eta1, phi1, eta2, phi2):
    """
    Calculate delta R between two physics objects. 
    Used to calculate delta R between a dark quark and all jets.
    (not sure if possible with existing function)
    
    Args:
        eta, phi (array-like or float): 
            eta and phi values of the two objects for which delta R is to be calculated
            Can be single values or arrays of the same shape.
    Returns:
        array-like or float
    """
    
    delta_phi = delta_phi_dark_quark(phi1, phi2)
    delta_eta = abs(eta1 - eta2)
    delta_R = np.sqrt(delta_eta**2 + delta_phi**2)
    return delta_R


def calculate_delta_r(
        physics_objects,
        indices=(0, 1),
        nan_value=-9999,
    ):
    """Calculate delta R between two physics objects.

    Args:
        physics_objects (awkward.Array): 
            Ak array where axis 0 is the event axis, axis 1 is the object axis
            with field eta (or rapidity) and phi and with name PtEtaPhiMLorentzVector.
            Physics objects can be jets, leptons etc...
        indices (tuple[int], optional):
            The indices of the object for which to compute delta eta.
            By default will do it the two leading objects.
        nan_value (float, optional, default=-9999):
            Value to use when the event has less objects than required
            If None, Nones will not be replaced.

    Returns:
        awkward.Array
    """

    delta_phi = calculate_delta_phi(physics_objects, indices, nan_value=None)
    delta_eta = calculate_delta_eta(physics_objects, indices, nan_value=None)
    delta_r = np.sqrt(delta_phi ** 2 + delta_eta ** 2)
    if nan_value is not None:
        delta_r = ak.fill_none(delta_r, nan_value)

    return delta_r


def calculate_invariant_mass(
        physics_objects,
        indices=(0, 1),
        nan_value=-9999,
    ):
    """Calculate invariant mass between two physics objects.

    Args:
        physics_objects (awkward.Array): 
            Ak array where axis 0 is the event axis, axis 1 is the object axis
            with field eta (or rapidity) and phi and with name PtEtaPhiMLorentzVector.
            Physics objects can be jets, leptons etc...
        indices (tuple[int], optional):
            The indices of the object for which to compute delta eta.
            By default will do it the two leading objects.
        nan_value (float, optional, default=-9999):
            Value to use when the event has less objects than required
            If None, Nones will not be replaced.

    Returns:
        awkward.Array
    """

    obj0, obj1 = __get_pair_of_objects(physics_objects, indices)
    invariant_mass = (obj0 + obj1).mass
    if nan_value is not None:
        invariant_mass = ak.fill_none(invariant_mass, nan_value)

    return invariant_mass


def calculate_lund_jet_plane_z(
        physics_objects,
        indices=(0, 1),
        nan_value=-9999,
    ):
    """Calculate the Lund jet plane variable z between two physics objects.

    Args:
        physics_objects (awkward.Array): 
            Ak array where axis 0 is the event axis, axis 1 is the object axis
            with field phi and with name PtEtaPhiMLorentzVector.
            Physics objects can be jets, leptons etc...
        indices (tuple[int], optional):
            The indices of the object for which to compute delta eta.
            By default will do it the two leading objects.
        nan_value (float, optional, default=-9999):
            Value to use when the event has less objects than required.
            If None, Nones will not be replaced.

    Returns:
        awkward.Array
    """

    obj0, obj1 = __get_pair_of_objects(physics_objects, indices)
    min_pt = ak.min([obj0.pt, obj1.pt], axis=0)
    z = min_pt / (obj0.pt + obj1.pt)
    if nan_value is not None:
        z = ak.fill_none(z, nan_value)

    return z


def calculate_atlas_momentum_balance(jets, met, nan_value=-9999):
    """ATLAS momentum balance."""

    jets = __get_objects_from_distance_with_met(jets, met, mode="min")
    n_jets = calculate_number_of_objects(jets)
    jets = ak.mask(jets, n_jets > 0)
    jet_min = jets[:, 0]
    jet_max = jets[:, -1]
    dijet = jet_min + jet_max
    pt_balance = dijet.pt / (jet_min.pt + jet_max.pt)
    pt_balance = ak.mask(pt_balance, n_jets >= 2)
    pt_balance = ak.fill_none(pt_balance, nan_value)

    return pt_balance


def calculate_atlas_delta_phi_max_min(jets, met, nan_value=-9999):
    """ATLAS momentum balance."""

    jets = __get_objects_from_distance_with_met(jets, met, mode="min")
    n_jets = calculate_number_of_objects(jets)
    jets = ak.mask(jets, n_jets > 0)
    jet_min = jets[:, 0]
    jet_max = jets[:, -1]
    delta_phi_max_min = jet_max.delta_phi(jet_min)
    delta_phi_max_min = ak.mask(delta_phi_max_min, n_jets >= 2)
    delta_phi_max_min = ak.fill_none(delta_phi_max_min, nan_value)

    return delta_phi_max_min


def calculate_transverse_mass(jets, met, jet_indices=(0, 1), n_jets=None, nan_value=-9999):
    """Calculate transverse mass variable for all events using 2 leading jets.

    Args:
        jets (awkward.Array):
            Ak array where axis 0 is the event axis, axis 1 is the jet axis
            with fields pt, rapidity, phi and mass and with name
            PtEtaPhiMLorentzVector.
        met (awkward.Array):
            Missing transverse energy. 1D ak array with fields pt, eta, phi,
            mass and with name PtEtaPhiMLorentzVector.
        jet_indices (list[int], optional, default=[0, 1]):
            Indices of the two jets used to compute MT
        n_jets (awkward.Array, optional, default=None):
            Ak array with one axis with number of jets in each event.
            If None, will be computed from jets.
        nan_value (float, optional, default=0):
            Value to use when the event has less than 2 jets.

    Returns:
        awkward.Array
    """

    jet0, jet1 = __get_pair_of_objects(jets, jet_indices)
    dijet = jet0 + jet1
    mt = np.sqrt( dijet.mass**2 + 2 * ( np.sqrt(dijet.mass**2 + dijet.pt**2) * met.pt - met.dot(dijet) ) )
    mt = ak.fill_none(mt, nan_value)

    return mt


def calculate_delta_phi_min(jets, met, nan_value=-9999):
    """Calculate the minimum delta phi between the MET and the jets.

    Args:
        physics_objects (awkward.Array): 
            Ak array where axis 0 is the event axis, axis 1 is the object axis
            with field eta and with name PtEtaPhiMLorentzVector.
            Physics objects can be jets, leptons etc...
        indices (tuple[int], optional):
            The indices of the object for which to compute delta eta.
            By default will do it the tw leading objects.
        nan_value (float, optional, default=-9999):
            Value to use when the event has less objects than required

    Returns:
        awkward.Array
    """

    met = ak.broadcast_arrays(met, jets)[0]
    delta_phi_min = ak.min(abs(jets.delta_phi(met)), axis=1)

    return delta_phi_min

