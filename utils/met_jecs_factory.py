import awkward
import numpy
from copy import copy

#here function from coffea METCorrector corrector to be used to update T1 MET with current JECs applied
def update_met_t1_corr(met_pt, met_phi, jet_pt, jet_phi, jet_pt_orig):
    sj, cj = numpy.sin(jet_phi), numpy.cos(jet_phi)

    #invert signs
    x = met_pt * numpy.cos(met_phi) - awkward.sum(
        jet_pt * cj - jet_pt_orig * cj, axis=1
    )
    y = met_pt * numpy.sin(met_phi) - awkward.sum(
        jet_pt * sj - jet_pt_orig * sj, axis=1
    )    
    return awkward.zip({"pt": numpy.hypot(x, y), "phi": numpy.arctan2(y, x)})

def apply_uncl_variation_to_met_t1(met_pt, met_phi,positive=None,dx=None,dy=None):

    x = met_pt * numpy.cos(met_phi) 
    y = met_pt * numpy.sin(met_phi) 
    if positive is not None and dx is not None and dy is not None:
        x = x + dx if positive else x - dx
        y = y + dy if positive else y - dy

    return awkward.zip(
        {"pt": numpy.hypot(x, y), "phi": numpy.arctan2(y, x)}, depth_limit=1
    )