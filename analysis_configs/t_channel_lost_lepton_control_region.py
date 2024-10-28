import os
import awkward as ak

from skimmer import skimmer_utils
from utils.awkward_array_utilities import as_type
import analysis_configs.triggers as trg
from analysis_configs.met_filters import met_filters_treemaker as met_filters
from analysis_configs import objects_definition as obj
from analysis_configs import sequences


def process(events, cut_flow, year, primary_dataset="", pn_tagger=False, **kwargs):
    """SVJ t-channel lost lepton control region.
    
    Same selections as for preselection region, but requiring exactly 1 veto lepton.
    """

    if skimmer_utils.is_data(events):
        events = sequences.remove_primary_dataset_overlap(events, year, primary_dataset)
        skimmer_utils.update_cut_flow(cut_flow, "PrimaryDatasetOvelap", events)

    # Add new branches at the start so that the defined branches can be used thereafter
    events = sequences.add_good_ak8_jet_branch(events)
    events = sequences.add_good_ak4_jet_branch(events)
    events = sequences.add_is_veto_electron_branch(events)
    events = sequences.add_is_veto_muon_branch(events)
    events = sequences.add_analysis_branches(events)

    # Trigger event selection
    triggers = getattr(trg, f"t_channel_{year}")
    events = skimmer_utils.apply_trigger_cut(events, triggers)
    skimmer_utils.update_cut_flow(cut_flow, "Trigger", events)

    # ST cut for triggers to be fully efficient
    events = events[events.ST > 1300]
    skimmer_utils.update_cut_flow(cut_flow, "STGt1300GeV", events)

    # HEM veto
    if year == "2018" and skimmer_utils.is_data(events):
        good_ak4_jets = events.Jets[events.Jets.isGood]
        veto_electrons = events.Electrons[events.Electrons.isVeto]
        veto_muons = events.Muons[events.Muons.isVeto]
        events = skimmer_utils.apply_hem_veto(events, good_ak4_jets, veto_electrons, veto_muons)
        skimmer_utils.update_cut_flow(cut_flow, "HEMVeto", events)

    # Good jet filters
    events = sequences.apply_good_ak8_jet_filter(events)
    skimmer_utils.update_cut_flow(cut_flow, "GoodJetsAK8", events)

    # Requiring at least 2 good FatJets
    filter = ak.count(events.JetsAK8.pt[events.JetsAK8.isGood], axis=1) >= 2
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK8Gt2", events)

    # At least 1 veto lepton
    events = sequences.add_n_lepton_veto_branch(events)
    events = events[events["NVetoLeptons"] >= 1]
    skimmer_utils.update_cut_flow(cut_flow, "NVetoLeptonsGt1", events)

    # Delta phi min cut
    # as_type needed otherwise type is not defined and skim cannot be written
    filter = as_type(events.DeltaPhiMinGoodJetsAK8 < 1.5, bool)
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "DeltaPhiMinLt1p5", events)

    # MET cut
    events = events[events.MET > 200]
    skimmer_utils.update_cut_flow(cut_flow, "METGt200GeV", events)

    # MET filter event selection
    events = skimmer_utils.apply_met_filters_cut(events, met_filters)
    skimmer_utils.update_cut_flow(cut_flow, "METFilters", events)

    # Phi spike filter
    events = skimmer_utils.apply_phi_spike_filter(
        events,
        year,
        f"{os.environ['SVJ_PROCESSING_ROOT']}/analysis_configs/tchannel_phi_spike_hot_spots.pkl",
        n_jets=4,
        jets_eta=events.Jets[events.Jets.isGood].eta,
        jets_phi=events.Jets[events.Jets.isGood].phi,
    )
    skimmer_utils.update_cut_flow(cut_flow, "PhiSpikeFilter", events)

    if pn_tagger:
        events = sequences.add_particle_net_tagger(events)

    events = sequences.remove_collections(events)

    return events, cut_flow

