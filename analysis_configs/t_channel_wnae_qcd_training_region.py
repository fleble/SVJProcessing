import awkward as ak

from skimmer import skimmer_utils
from utils.awkward_array_utilities import as_type
import analysis_configs.triggers as trg
import analysis_configs.met_filters as met_filters
from analysis_configs import sequences


def process(events, cut_flow, year, **kwargs):
    """SVJ t-channel WNAE QCD training region targetting gamma + jets events."""

    # Trigger event selection
    triggers = getattr(trg, f"photon_{year}")
    events = skimmer_utils.apply_trigger_cut(events, triggers)
    skimmer_utils.update_cut_flow(cut_flow, "Trigger", events)

    # Adding branches already so that it can be used in the rest of the selection
    events = sequences.add_st(events)
    events = sequences.add_good_ak8_jet_branch(events)
    events = sequences.add_good_ak4_jet_branch(events)
    events = sequences.add_is_veto_electron_branch(events)
    events = sequences.add_is_veto_muon_branch(events)
    events = sequences.add_good_photon_branch(events)

    # Require at least 1 photon
    filter = ak.count(events.Photons.pt[events.Photons.isGood], axis=1) >= 1
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, f"nPhotonPtGt1", events)

    # Cuts for the triggers to be fully efficient
    if year in ["2017", "2018"]:
        min_photon_pt = 250
    else:
        min_photon_pt = 200
    filter = events.Photons.pt[:, 0] > min_photon_pt
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, f"LeadingPhotonPtGt{min_photon_pt}GeV", events)

    events = events[events.ST > 800]
    skimmer_utils.update_cut_flow(cut_flow, "STGt800GeV", events)

    # MET filter event selection
    met_filters_names = getattr(met_filters, f"met_filters_treemaker_{year}")
    events = skimmer_utils.apply_met_filters_cut(events, met_filters_names)
    skimmer_utils.update_cut_flow(cut_flow, "METFilters", events)

    # HEM veto
    good_ak4_jets = events.Jets[events.Jets.isGood]
    good_photons = events.Photons[events.Photons.isGood]
    if year == "2018" and skimmer_utils.is_data(events):
        events = skimmer_utils.apply_hem_veto(events, good_ak4_jets, good_photons)
        skimmer_utils.update_cut_flow(cut_flow, "HEMVeto", events)
    if year == "2018" and skimmer_utils.is_mc(events):
        filter = skimmer_utils.get_hem_veto_filter(good_ak4_jets, good_photons)
        events["HEMVeto"] = filter

    # Define as training jets the AK8 jets with angular separation
    # with the tag photon greater than 0.8
    jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=events.JetsAK8.pt,
        eta=events.JetsAK8.eta,
        phi=events.JetsAK8.phi,
        mass=events.JetsAK8.mass,
    )
    photons = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=events.Photons.pt,
        eta=events.Photons.eta,
        phi=events.Photons.phi,
        mass=events.Photons.mass,
    )
    tag_photon = photons[:, 0:1]

    jets_, tag_photon_ = ak.unzip(ak.cartesian([jets, tag_photon], axis=-1, nested=True))
    min_delta_r = ak.min(jets_.delta_r(tag_photon_), axis=-1)
    is_training_jet = (min_delta_r > 0.8) & events.JetsAK8.isGood
    is_training_jet = as_type(is_training_jet, bool)
    events["JetsAK8"] = ak.with_field(
        events["JetsAK8"],
        is_training_jet,
        "isTrainingJet",
    )

    # Veto events with mini-isolated leptons
    events = sequences.apply_lepton_veto(events)
    skimmer_utils.update_cut_flow(cut_flow, "LeptonVeto", events)

    # Requiring at least 1 training AK8 jet
    filter = ak.count(events.JetsAK8.pt[events.JetsAK8.isTrainingJet], axis=1) >= 1
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nTrainingJetsAK8Gt1", events)

    events = sequences.remove_collections(events)

    return events, cut_flow

