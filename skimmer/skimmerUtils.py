import awkward as ak

from utils.Logger import *


def update_cut_flow(cut_flow, cut_name, events):
    """Update cut flow table in a coffea accumulator.

    Args:
        cut_flow (dict[str, float])
        cut_name (str): the name of the cut to appear in the cut flow tree
        events (EventsFromAkArray)
        use_raw_events (bool)
    """

    if cut_name in cut_flow.keys():
        cut_flow[cut_name][0] += get_number_of_events(events)
    else:
        cut_flow[cut_name] = [get_number_of_events(events)]


def get_number_of_events(events):
    return ak.sum(events.Weight)

