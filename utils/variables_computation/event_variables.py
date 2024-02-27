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

