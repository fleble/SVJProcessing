import awkward as ak

from skimmer import skimmer_utils


def process(events, cut_flow, year, primary_dataset="", pn_tagger=False):

    # Adding HT as event variable
    events = ak.with_field(
        events,
        ak.sum(events.Jet_pt, axis=-1),
        "HT",
    )

    # ST cut for triggers to be fully efficient
    events = events[events.HT > 800]
    skimmer_utils.update_cut_flow(cut_flow, "HTGt800GeV", events)

    return events, cut_flow
