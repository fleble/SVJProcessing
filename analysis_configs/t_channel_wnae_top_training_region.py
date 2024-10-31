import awkward as ak

from skimmer import skimmer_utils
import analysis_configs.triggers as trg
from analysis_configs.met_filters import met_filters_treemaker as met_filters
from analysis_configs import objects_definition as obj
from analysis_configs import sequences


def process(events, cut_flow, year, primary_dataset="", **kwargs):
    """SVJ t-channel WNAE top training region targetting semi-leptonic ttbar events."""

    if not skimmer_utils.is_mc(events):
        events = sequences.remove_single_lepton_primary_dataset_overlap(events, year, primary_dataset)
        skimmer_utils.update_cut_flow(cut_flow, "PrimaryDatasetOvelap", events)

    # Adding branches already so that it can be used in the rest of the selection
    events = sequences.add_good_ak8_jet_branch(events)
    events = sequences.add_good_ak4_jet_branch(events)
    events = sequences.add_is_veto_electron_branch(events)
    events = sequences.add_is_veto_muon_branch(events)

    # Trigger event selection
    triggers = getattr(trg, f"single_lepton_{year}")
    events = skimmer_utils.apply_trigger_cut(events, triggers)
    skimmer_utils.update_cut_flow(cut_flow, "Trigger", events)
    
    # ST cut for the training phase space to closer to the preselection phase space
    events = sequences.add_st(events)
    events = events[events.ST > 600]
    skimmer_utils.update_cut_flow(cut_flow, "STGt600GeV", events)


    # MET filter event selection
    events = skimmer_utils.apply_met_filters_cut(events, met_filters)
    skimmer_utils.update_cut_flow(cut_flow, "METFilters", events)


    # HEM veto
    good_ak4_jets = events.Jets[events.Jets.isGood]
    veto_electrons = events.Electrons[events.Electrons.isVeto]
    veto_muons = events.Muons[events.Muons.isVeto]
    if year == "2018" and skimmer_utils.is_data(events):
        events = skimmer_utils.apply_hem_veto(events, good_ak4_jets, veto_electrons, veto_muons)
        skimmer_utils.update_cut_flow(cut_flow, "HEMVeto", events)
    if year == "2018" and skimmer_utils.is_mc(events):
        filter = skimmer_utils.get_hem_veto_filter(good_ak4_jets, veto_electrons, veto_muons)
        events["HEMVeto"] = filter


    # Require exactly 1 "tag lepton" (medium ID, tight mini-iso electron or muon)
    tag_electrons = events.Electrons[obj.is_tag_electron(events.Electrons)]
    tag_muons = events.Muons[obj.is_tag_muon(events.Muons)]
    n_tag_electrons = ak.count(tag_electrons.pt, axis=1)
    n_tag_muons = ak.count(tag_muons.pt, axis=1)
    n_tag_leptons = n_tag_electrons + n_tag_muons
    events = events[n_tag_leptons == 1]
    skimmer_utils.update_cut_flow(cut_flow, f"nTagLeptonEq1", events)
    

    # Veto events with additional veto leptons
    events = sequences.apply_lepton_veto(
        events,
        electron_extra_condition=~obj.is_tag_electron(events.Electrons),
        muon_extra_condition=~obj.is_tag_muon(events.Muons),
    )
    skimmer_utils.update_cut_flow(cut_flow, "LeptonVeto", events)


    # Define clean AK4 and AK8 jets
    cleaning_electrons = events.Electrons[obj.is_cleaning_electron(events.Electrons)]
    cleaning_muons = events.Muons[obj.is_cleaning_muon(events.Muons)]
    cleaning_electrons = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=cleaning_electrons.pt,
        eta=cleaning_electrons.eta,
        phi=cleaning_electrons.phi,
        mass=cleaning_electrons.mass,
    )
    cleaning_muons = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=cleaning_muons.pt,
        eta=cleaning_muons.eta,
        phi=cleaning_muons.phi,
        mass=cleaning_muons.mass,
    )
    cleaning_leptons = ak.concatenate((cleaning_electrons, cleaning_muons), axis=1)

    is_clean_ak4_jet = skimmer_utils.is_clean(
        events.Jets,
        cleaning_leptons,
        radius=0.4,
    )
    events["Jets"] = ak.with_field(
        events["Jets"],
        is_clean_ak4_jet,
        "isClean",
    )

    is_clean_ak8_jet = skimmer_utils.is_clean(
        events.JetsAK8,
        cleaning_leptons,
        radius=0.8,
    )
    events["JetsAK8"] = ak.with_field(
        events["JetsAK8"],
        is_clean_ak8_jet,
        "isClean",
    )


    # Require at least 4 AK4 jets
    is_good_clean_ak4_jets = events.Jets.isGood & events.Jets.isClean
    n_ak4_jets = ak.sum(is_good_clean_ak4_jets, axis=1)
    events = events[n_ak4_jets >= 4]
    skimmer_utils.update_cut_flow(cut_flow, "nAK4JetsGt4", events)


    # Require exactly 2 b-tagged jets among 4 leading AK4 jets
    b_tagging_score = (
        events.Jets.bJetTagDeepFlavourprobb
        + events.Jets.bJetTagDeepFlavourprobbb
        + events.Jets.bJetTagDeepFlavourproblepb
    )
    is_b_tagged = b_tagging_score >= skimmer_utils.get_b_tagging_wp(year)
    is_b_tagged = is_b_tagged & (ak.local_index(is_b_tagged) < 4)
    events["Jets"] = ak.with_field(
        events["Jets"],
        is_b_tagged,
        "isBTagged",
    )
    n_b_tag = ak.sum(events.Jets.isBTagged & events.Jets.isGood & events.Jets.isClean, axis=1)
    events = events[n_b_tag == 2]
    skimmer_utils.update_cut_flow(cut_flow, "nBTagEq2", events)


    # Do not consider as training jets the AK8 jets with angular separation
    # less than 0.8 from the b-tagged AK4 jet closest to the tag lepton
    tag_electrons = events.Electrons[obj.is_tag_electron(events.Electrons)]
    tag_muons = events.Muons[obj.is_tag_muon(events.Muons)]
    tag_leptons = ak.concatenate((tag_electrons, tag_muons), axis=1)

    is_b_tagged = events.Jets.isBTagged & events.Jets.isGood & events.Jets.isClean
    b_tagged_jets = events.Jets[is_b_tagged]

    mapping = skimmer_utils.collections_matching(
        collection1=tag_leptons,
        collection2=b_tagged_jets,
    )
    matched_b_tagged_jet = b_tagged_jets[mapping]
    is_not_tag_side = skimmer_utils.is_clean(
        events.JetsAK8,
        matched_b_tagged_jet,
        radius=0.8,
    )

    is_training_jet = (
        events.JetsAK8.isGood
        & events.JetsAK8.isClean
        & is_not_tag_side
    )

    events["JetsAK8"] = ak.with_field(
        events["JetsAK8"],
        is_training_jet,
        "isTrainingJet",
    )


    # Require at least 1 training AK8 jet
    filter = ak.count(events.JetsAK8.pt[events.JetsAK8.isTrainingJet], axis=1) >= 1
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nTrainingJetsAK8Gt1", events)


    # Remove branches added for convenience
    branch_names = ["isGood", "isClean", "isBTagged"]
    for branch_name in branch_names:
        events["Jets"] = events["Jets"][[x for x in events["Jets"].fields if x != branch_name]]
    branch_names = ["isClean"]
    for branch_name in branch_names:
        events["JetsAK8"] = events["JetsAK8"][[x for x in events["JetsAK8"].fields if x != branch_name]]

    
    events = sequences.remove_collections(events)


    return events, cut_flow

