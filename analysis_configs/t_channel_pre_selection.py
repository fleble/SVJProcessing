import awkward as ak

from skimmer import skimmer_utils
from utils.awkward_array_utilities import as_type
from utils.variables_computation import event_variables
import analysis_configs.triggers as trg
from analysis_configs.met_filters import met_filters


def __is_good_jet(jets_ak8):
    return jets_ak8.ID == 1


def __is_analysis_jet(jets_ak8):
    filter = (
        (jets_ak8.pt > 50)
        & (abs(jets_ak8.eta) < 2.4)
    )
    return filter


def process(events, cut_flow, year):
    """SVJ t-channel pre-selection."""

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
    analysis_jets = events.JetsAK8[__is_analysis_jet(events.JetsAK8)]
    good_jets_filter = ak.all(__is_good_jet(analysis_jets), axis=1)
    events = events[good_jets_filter]
    skimmer_utils.update_cut_flow(cut_flow, "GoodJetsAK8", events)

    # Adding JetsAK8_idGood branch already so that it can be used
    # in the rest of the pre-selection
    is_good_analysis_jet = (
        __is_analysis_jet(events.JetsAK8)
        & __is_good_jet(events.JetsAK8)
    )
    events["JetsAK8"] = ak.with_field(
        events["JetsAK8"],
        is_good_analysis_jet,
        "isGood",
    )

    # Requiring at least 2 good FatJets
    filter = ak.count(events.JetsAK8.pt[events.JetsAK8.isGood], axis=1) >= 2
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK8Gt2", events)

    # Veto events with mini-isolated leptons
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
    good_electrons = events.Electrons[good_electron_filter]
    good_muons = events.Muons[good_muon_filter]
    n_electrons = ak.count(good_electrons.pt, axis=1)
    n_muons = ak.count(good_muons.pt, axis=1)
    n_leptons = n_electrons + n_muons
    events = events[n_leptons == 0]
    skimmer_utils.update_cut_flow(cut_flow, "LeptonVeto", events)

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


    # Adding new branches
    
    # Jets AK8 variables
    new_branches = {}
    jets_ak8 = events.JetsAK8
    gen_jets_ak8 = events.GenJetsAK8
    jets_ak8_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=jets_ak8.pt,
        eta=jets_ak8.eta,
        phi=jets_ak8.phi,
        mass=jets_ak8.mass,
    )
    met_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=events.MET,
        phi=events.METPhi,
    )

    met_bc = ak.broadcast_arrays(met_lv, jets_ak8_lv)[0]
    new_branches["deltaPhiMET"] = jets_ak8_lv.delta_phi(met_bc)

    if "hvCategory" in gen_jets_ak8.fields:  # If signal samples
        gen_jet_index = ak.mask(jets_ak8.genIndex, jets_ak8.genIndex >= 0)
        hv_category = gen_jets_ak8.hvCategory[gen_jet_index]
        hv_category = ak.fill_none(hv_category, -9999)
        new_branches["hvCategory"] = hv_category

        new_branches["isDarkJetTightNoMix"] = (
            (hv_category == 1)
            | (hv_category == 3)
            | (hv_category == 5)
            | (hv_category == 9)
        )

        new_branches["isDarkJetTight"] = (
            (hv_category != -9999)
            & (hv_category != 0)
            & (hv_category < 16)
        )

        new_branches["isDarkJetMedium"] = (
            (hv_category != -9999)
            & (hv_category != 0)
            & (hv_category != 16)
            & (hv_category != 17)
        )

        new_branches["isDarkJetLoose"] = (
            (hv_category != -9999)
            & (hv_category != 0)
            & (hv_category != 16)
        )

    for branch_name, branch in new_branches.items():
        events["JetsAK8"] = ak.with_field(
            events["JetsAK8"],
            branch,
            branch_name,
        )
    
    # Event variables
    nan_value = -9999
    good_jets_ak8 = events.JetsAK8[events.JetsAK8.isGood]
    good_jets_ak8_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=good_jets_ak8.pt,
        eta=good_jets_ak8.eta,
        phi=good_jets_ak8.phi,
        mass=good_jets_ak8.mass,
    )

    for index_0 in [0, 1, 2, 3]:
        for index_1 in [0, 1, 2, 3]:
            if index_1 == index_0: continue
            delta_eta = event_variables.calculate_delta_eta(
                physics_objects=good_jets_ak8_lv,
                indices=(index_0, index_1),
                absolute_value=False,
                nan_value=None,
            )
            delta_phi = event_variables.calculate_delta_phi(
                physics_objects=good_jets_ak8_lv,
                indices=(index_0, index_1),
                absolute_value=False,
                nan_value=None,
            )
            delta_r = event_variables.calculate_delta_r(
                physics_objects=good_jets_ak8_lv,
                indices=(index_0, index_1),
                nan_value=nan_value,
            )
            delta_eta_abs = abs(delta_eta)
            delta_phi_abs = abs(delta_phi)
            delta_eta = ak.fill_none(delta_eta, nan_value)
            delta_phi = ak.fill_none(delta_eta, nan_value)
            delta_eta_abs = ak.fill_none(delta_eta, nan_value)
            delta_phi_abs = ak.fill_none(delta_eta, nan_value)

            events[f"DeltaEta{index_0}{index_1}GoodJetsAK8"] = delta_eta
            events[f"DeltaPhi{index_0}{index_1}GoodJetsAK8"] = delta_phi
            events[f"DeltaR{index_0}{index_1}GoodJetsAK8"] = delta_r
            events[f"DeltaEtaAbs{index_0}{index_1}GoodJetsAK8"] = delta_eta_abs
            events[f"DeltaPhiAbs{index_0}{index_1}GoodJetsAK8"] = delta_phi_abs
   
    events["DeltaPhiMinGoodJetsAK8"] = ak.min(abs(good_jets_ak8.deltaPhiMET), axis=1)
    events["ST"] = events.MET + events.HT
    events["ATLASDeltaPhiMinMax"] = event_variables.calculate_atlas_delta_phi_max_min(
        jets=good_jets_ak8_lv,
        met=met_lv,
        nan_value=nan_value,
    )
    events["ATLASPtBalance"] = event_variables.calculate_atlas_momentum_balance(
        jets=good_jets_ak8_lv,
        met=met_lv,
        nan_value=nan_value,
    )

    # Removing un-necessary collections
    events = events[[x for x in events.fields if x != "JetsAK15"]]


    return events, cut_flow

