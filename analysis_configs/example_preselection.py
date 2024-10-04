import awkward as ak

from skimmer import skimmer_utils


def process(events, cut_flow, year, **kwargs):
    """Example pre-selection config."""

    # Example for object selection
    # This filters the jet collections, e.g. AK8 jets with pT < 50 GeV are removed from the events!
    # All cross reference are handled accordingly
    events["Jets"] = events.Jets[(events.Jets.pt > 30) & (abs(events.Jets.eta) < 2.4) & (events.Jets.ID == 1)]
    events["JetsAK8"] = events.JetsAK8[(events.JetsAK8.pt > 50) & (abs(events.JetsAK8.eta) < 2.4) & (events.JetsAK8.ID == 1)]

    # Example for event selection
    st = events.MET + ak.sum(events.Jets.pt, axis=1)
    events = events[st > 1300]
    skimmer_utils.update_cut_flow(cut_flow, "STGt1300GeV", events)

    return events, cut_flow

