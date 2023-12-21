import awkward as ak

from skimmer import skimmer_utils
from utils.awkward_array_utilities import as_type
import analysis_configs.triggers as trg
from analysis_configs.met_filters import met_filters


def process(events, cut_flow, year):
    """SVJ t-channel pre-selection."""

    # Trigger event selection
    triggers = eval(f"trg.t_channel_{year}")
    events = skimmer_utils.apply_trigger_cut(events, triggers)
    skimmer_utils.update_cut_flow(cut_flow, "Trigger", events)

    # Objects filter
    good_jet_filter = (
        (events.Jets.pt > 30)
        & (abs(events.Jets.eta) < 2.4)
        & (events.Jets.ID == 1)
    )
    good_ak8_jet_filter = (
        (events.JetsAK8.pt > 50)
        & (abs(events.JetsAK8.eta) < 2.4)
        & (events.JetsAK8.ID == 1)
    )
    good_electron_filter = (
        (events.Electrons.pt > 10)
        & (abs(events.Electrons.eta) < 2.4)
        & (abs(events.Electrons.iso) < 0.1)
    )
    good_muon_filter = (
        (events.Muons.pt > 10)
        & (abs(events.Muons.eta) < 2.4)
        & (abs(events.Muons.iso) < 0.4)
    )
    events["Jets"] = events.Jets[good_jet_filter]
    events["JetsAK8"] = events.JetsAK8[good_ak8_jet_filter]
    events["Electrons"] = events.Electrons[good_electron_filter]
    events["Muons"] = events.Muons[good_muon_filter]

    # ST cut for triggers to be fully efficient
    st = events.MET + ak.sum(events.Jets.pt, axis=1)
    events = events[st > 1300]
    skimmer_utils.update_cut_flow(cut_flow, "STGt1300GeV", events)

    # MET filter event selection
    events = skimmer_utils.apply_met_filters_cut(events, met_filters)
    skimmer_utils.update_cut_flow(cut_flow, "METFilters", events)

    # Veto events with mini-isolated leptons
    n_electrons = ak.count(events.Electrons.pt, axis=1)
    n_muons = ak.count(events.Muons.pt, axis=1)
    n_leptons = n_electrons + n_muons
    events = events[n_leptons == 0]
    skimmer_utils.update_cut_flow(cut_flow, "LeptonVeto", events)

    # # Delta phi min cut
    # met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
    #     pt=events.MET,
    #     phi=events.METPhi,
    # )
    # jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
    #     pt=events.JetsAK8.pt,
    #     eta=events.JetsAK8.eta,
    #     phi=events.JetsAK8.phi,
    #     mass=events.JetsAK8.mass,
    # )

    # met = ak.broadcast_arrays(met, jets)[0]
    # delta_phi_min = ak.min(abs(jets.delta_phi(met)), axis=1)
    # filter = delta_phi_min < 1.5
    # # Needed otherwise type is not defined and skim cannot be written
    # filter = as_type(filter, bool)
    # events = events[filter]
    # skimmer_utils.update_cut_flow(cut_flow, "DeltaPhiMinLt1p5", events)

    # Requiring at least 2 FatJets
    events = events[ak.count(events.JetsAK8.pt, axis=1) >= 2]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK8Gt2", events)

    return events, cut_flow

