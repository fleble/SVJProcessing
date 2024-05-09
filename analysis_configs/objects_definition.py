def is_good_jet(jets_ak8):
    return jets_ak8.ID == 1


def is_analysis_jet(jets_ak8):
    filter = (
        (jets_ak8.pt > 50)
        & (abs(jets_ak8.eta) < 2.4)
    )
    return filter


def is_veto_electron(electrons):
    filter = (
        (electrons.pt > 10)
        & (abs(electrons.eta) < 2.4)
        & (abs(electrons.iso) < 0.1)
    )
    return filter


def is_veto_muon(muons):
    filter = (
        (muons.pt > 10)
        & (abs(muons.eta) < 2.4)
        & (abs(muons.iso) < 0.4)
    )
    return filter


def is_good_photon(photons):
    filter = (
        (abs(photons.eta) < 2.5)  # Photon in ECAL
        & ((abs(photons.eta) < 1.4442) | (abs(photons.eta) > 1.566))  # Exclude barrel-endcap transition 
        & (photons.cutBasedID >= 2)  # Medium ID
        & (photons.hasPixelSeed == False)  # No pixel seed -> conversion veto
    )
    return filter

