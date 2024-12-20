import os
import awkward as ak

from skimmer import skimmer_utils
from utils.awkward_array_utilities import as_type
import analysis_configs.triggers as trg
import utils.variables_computation.event_variables as event_vars
from analysis_configs.met_filters import met_filters_nanoaod as met_filters
from analysis_configs import sequences_s_channel_leptons as sequences


def process(events, cut_flow, year, primary_dataset="", pn_tagger=False, **kwargs):
    """SVJ s-channel leptons pre-selection."""

    # Trigger event selection
    triggers = getattr(trg, f"s_channel_{year}")
    events = skimmer_utils.apply_trigger_cut(events, triggers)
    skimmer_utils.update_cut_flow(cut_flow, "Trigger", events)

    # Good jet filters
    events = sequences.apply_good_ak8_jet_filter(events)
    skimmer_utils.update_cut_flow(cut_flow, "GoodJetsAK8", events)


    # Adding JetsAK8_isGood branch already so that it can be used
    # in the rest of the pre-selection
    events = sequences.add_good_ak8_jet_branch(events)

    # Requiring at least 2 good FatJets
    filter = ak.count(events.FatJet_pt[events.FatJet_isGood], axis=1) >= 2
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK8Gt2", events)

    #apply RT filter (RT = MET over MT)
    if len(events) != 0:
        jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.FatJet_pt[events.FatJet_isGood],
            eta=events.FatJet_eta[events.FatJet_isGood],
            phi=events.FatJet_phi[events.FatJet_isGood],
            mass=events.FatJet_mass[events.FatJet_isGood],
        )
        met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET_pt,
            phi=events.MET_phi,
        )
        mt = event_vars.calculate_transverse_mass(jets, met)
        rt = events.MET_pt / mt
        filter_rt = rt > 0.15
        filter_rt = as_type(filter_rt, bool)   #not needed
        events = events[filter_rt]
    
    skimmer_utils.update_cut_flow(cut_flow, "RT selection", events)


    #apply DeltaEta filter
    if len(events) != 0:
        jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.FatJet_pt[events.FatJet_isGood],
            eta=events.FatJet_eta[events.FatJet_isGood],
            phi=events.FatJet_phi[events.FatJet_isGood],
            mass=events.FatJet_mass[events.FatJet_isGood],
        )
        delta_eta = abs(event_vars.calculate_delta_eta(jets))
        filter_deltaeta = delta_eta < 2.2
        filter_deltaeta = as_type(filter_deltaeta, bool)
        events = events[filter_deltaeta]
    
    skimmer_utils.update_cut_flow(cut_flow, "DeltaEtaj0j1 selection", events)

    #apply MT selection
    if len(events) != 0:
        jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.FatJet_pt[events.FatJet_isGood],
            eta=events.FatJet_eta[events.FatJet_isGood],
            phi=events.FatJet_phi[events.FatJet_isGood],
            mass=events.FatJet_mass[events.FatJet_isGood],
        )
        met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET_pt,
            phi=events.MET_phi,
        )
        mt = event_vars.calculate_transverse_mass(jets, met)
        filter_mt = mt > 1500
        filter_mt = as_type(filter_mt, bool)
        events = events[filter_mt]
    
    skimmer_utils.update_cut_flow(cut_flow, "MT selection", events)

    # MET filter event selection
    events = skimmer_utils.apply_met_filters_cut(events, met_filters)
    skimmer_utils.update_cut_flow(cut_flow, "METFilters", events)

    # Phi spike filter
    events = skimmer_utils.apply_phi_spike_filter(
        events,
        year,
        f"{os.environ['SVJ_PROCESSING_ROOT']}/analysis_configs/schannel_hot_spots.pkl",
        n_jets=2,
        jets_eta=events.Jet_eta,
        jets_phi=events.Jet_phi,
    )
    skimmer_utils.update_cut_flow(cut_flow, "PhiSpikeFilter", events)

    # apply HEM issue filter - to be applied only on 2018 data
    if year == "2018" and skimmer_utils.is_data(events):
        ak4_jets = ak.zip({
            "pt": events["Jet_pt"],
            "eta": events["Jet_eta"],
            "phi": events["Jet_phi"],
        })
        electrons = ak.zip({
            "pt":  events["Electron_pt"],
            "eta": events["Electron_eta"],
            "phi": events["Electron_phi"],
        })
        muons = ak.zip({
            "pt":  events["Muon_pt"],
            "eta": events["Muon_eta"],
            "phi": events["Muon_phi"],
        })
        events = skimmer_utils.apply_hem_veto(events, ak4_jets, electrons, muons)
        skimmer_utils.update_cut_flow(cut_flow, "HEMIssueFilter", events)

    # Delta phi min cut
    if len(events) != 0:
        # If needed because the selection crashes due to the special ak type
        met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET_pt,
            phi=events.MET_phi,
        )
        jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.FatJet_pt[events.FatJet_isGood],
            eta=events.FatJet_eta[events.FatJet_isGood],
            phi=events.FatJet_phi[events.FatJet_isGood],
            mass=events.FatJet_mass[events.FatJet_isGood],
        )

        met = ak.broadcast_arrays(met, jets)[0]
        delta_phi_min = ak.min(abs(jets.delta_phi(met)), axis=1)
        filter_deltaphi = delta_phi_min < 0.8
        # Needed otherwise type is not defined and skim cannot be written
        filter_deltaphi = as_type(filter_deltaphi, bool)
        events = events[filter_deltaphi]

    skimmer_utils.update_cut_flow(cut_flow, "DeltaPhiMin selection", events)
    
    events = sequences.add_analysis_branches(events)
    events = sequences.remove_collections(events)

    return events, cut_flow

