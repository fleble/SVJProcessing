import numba as nb
import awkward as ak
import utils.awkward_array_utilities as akUtl
from utils.Logger import *


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