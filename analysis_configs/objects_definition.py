def is_good_ak8_jet(jets_ak8):
    return jets_ak8.ID == 1


def is_analysis_ak8_jet(jets_ak8):
    return (
        (jets_ak8.pt > 50)
        & (abs(jets_ak8.eta) < 2.4)
    )


def is_good_ak4_jet(jets):
    return (
        (jets.pt > 30)
        & (abs(jets.eta) < 2.4)
        & (jets.ID == 1)
    )


def is_analysis_electron(electrons):
    return (
        (electrons.pt > 10)
        & (abs(electrons.eta) < 2.4)
    )


def is_analysis_muon(muons):
    return (
        (muons.pt > 10)
        & (abs(muons.eta) < 2.4)
    )


def is_veto_electron(electrons):
    return (
        is_analysis_electron(electrons)
        & (abs(electrons.iso) < 0.1)
    )


def is_veto_muon(muons):
    return (
        is_analysis_muon(muons)
        & (abs(muons.iso) < 0.4)
    )


def is_tag_electron(electrons):
    return (
        is_analysis_electron(electrons)
        & (electrons.mediumID == 1)
        & (electrons.iso < 0.1)  # mini-isolation, tight WP
    )


def is_tag_muon(muons):
    return (
        is_analysis_muon(muons)
        & (muons.mediumID == 1)
        & (muons.iso < 0.1)  # mini-isolation, tight WP
    )


def __is_isolated_electron(electrons):
    return electrons.pfRelIso < 0.15


def __is_isolated_muon(muons):
    return muons.pfRelIso < 0.15


def is_cleaning_electron(electrons):
    return (
        is_tag_electron(electrons)
        & __is_isolated_electron(electrons)
    )


def is_cleaning_muon(muons):
    return (
        is_tag_muon(muons)
        & __is_isolated_muon(muons)
    )


def is_good_photon(photons):
    filter = (
        (abs(photons.eta) < 2.5)  # Photon in ECAL
        & ((abs(photons.eta) < 1.4442) | (abs(photons.eta) > 1.566))  # Exclude barrel-endcap transition 
        & (photons.cutBasedID >= 2)  # Medium ID
        & (photons.hasPixelSeed == False)  # No pixel seed -> conversion veto
    )
    return filter

