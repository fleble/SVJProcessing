import awkward as ak

from skimmer import skimmer_utils
from utils.awkward_array_utilities import as_type
import analysis_configs.triggers as trg
from analysis_configs.met_filters import met_filters_treemaker as met_filters
from analysis_configs import sequences


def process(events, cut_flow, year, primary_dataset="", pn_tagger=False):
    """SVJ t-channel lost lepton control region.
    
    Same selections as for preselection region, but requiring exactly 1 veto lepton.
    """

    if not skimmer_utils.is_mc(events):
        events = sequences.remove_primary_dataset_overlap(events, year, primary_dataset)
        skimmer_utils.update_cut_flow(cut_flow, "PrimaryDatasetOvelap", events)

    # Trigger event selection
    triggers = getattr(trg, f"t_channel_{year}")
    events = skimmer_utils.apply_trigger_cut(events, triggers)
    skimmer_utils.update_cut_flow(cut_flow, "Trigger", events)

    # ST cut for triggers to be fully efficient
    st = events.MET + events.HT
    events = events[st > 1300]
    skimmer_utils.update_cut_flow(cut_flow, "STGt1300GeV", events)

    # MET filter event selection
    events = skimmer_utils.apply_met_filters_cut(events, met_filters)
    skimmer_utils.update_cut_flow(cut_flow, "METFilters", events)

    # Good jet filters
    events = sequences.apply_good_ak8_jet_filter(events)
    skimmer_utils.update_cut_flow(cut_flow, "GoodJetsAK8", events)

    # Adding JetsAK8_isGood branch already so that it can be used
    # in the rest of the pre-selection
    events = sequences.add_good_ak8_jet_branch(events)

    # Requiring at least 2 good FatJets
    filter = ak.count(events.JetsAK8.pt[events.JetsAK8.isGood], axis=1) >= 2
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK8Gt2", events)

    # Exactly 1 veto lepton
    events = sequences.require_n_veto_leptons(events, n=1)
    skimmer_utils.update_cut_flow(cut_flow, "OneVetoLepton", events)

    # Delta phi min cut
    if len(events) != 0:
        # If needed because the selection crashes due to the special ak type
        met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET,
            phi=events.METPhi,
        )
        good_jets_ak8 = events.JetsAK8[events.JetsAK8.isGood]
        jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=good_jets_ak8.pt,
            eta=good_jets_ak8.eta,
            phi=good_jets_ak8.phi,
            mass=good_jets_ak8.mass,
        )

        met = ak.broadcast_arrays(met, jets)[0]
        delta_phi_min = ak.min(abs(jets.delta_phi(met)), axis=1)
        filter = delta_phi_min < 1.5
        # Needed otherwise type is not defined and skim cannot be written
        filter = as_type(filter, bool)
        events = events[filter]
        
        delta_phi_min = delta_phi_min[filter]

    skimmer_utils.update_cut_flow(cut_flow, "DeltaPhiMinLt1p5", events)

    # MET cut
    events = events[events.MET > 200]
    skimmer_utils.update_cut_flow(cut_flow, "METGt200GeV", events)


    events = sequences.add_analysis_branches(events)

    if pn_tagger:
        events = sequences.add_particle_net_tagger(events)

    events = sequences.remove_collections(events)

    return events, cut_flow

