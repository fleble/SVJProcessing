import awkward as ak
import numpy as np
import utils.gen_matching_tools as gen_matching_tools
from utils.jet_energy_scale_svj_factory import SVJCustomJESCalculator 
from coffea.lookup_tools import extractor
from coffea.jetmet_tools import CorrectedMETFactory

def calc_jec_variation(
        pt, eta, phi, energy,
        jer_factor, jec_unc, orig_idx,
        variation_orig_idx, variation_jer_factor
        ):
    """
    Applies a JEC variation (up or down) on a 4-vector.

    Note there are 3 'ordering levels':
    - "Final": the final ordering of jets after centrally applied corrections
    - "Original": the ordering of jets _before_ any corrections
    - "Variation": the final ordering of jets after the applying the correction of the
        _variation_

    The algorithm below first creates a map to reorder "Final" to "Variation", then
    applies the correction after ordering everything in "Variation" ordering.

    Args:
        pt (ak.Array): jet pt
        eta (ak.Array): jet eta
        phi (ak.Array): jet phi
        energy (ak.Array): jet energy
        jer_factor (ak.Array): the JER factor that was applied centrally to obtain the
            final jet
        jec_unc (ak.Array): the JEC uncertainty
        orig_idx (ak.Array): mapping of final corrected jet ordering back to 'original'
            ordering
        variation_orig_idx (ak.Array): mapping of variation ordering back to 'original'
            ordering
        variation_jer_factor (ak.Array): the variation's JER factor

    Returns:
        (ak.Array, ak.Array, ak.Array, ak.Array, ak.Array) : pt, eta, phi, energy after
            applying the variation and reordered by the pT after variation; permutation
            to propagate the reordering to the whole collection
    """

    # Create a map to reorder final corrected jets to the ordering of the variation
    map_orig_idx_to_var_idx = ak.argsort(variation_orig_idx, axis=-1)
    map_final_idx_to_orig_idx = orig_idx
    reorder_final_to_var = map_final_idx_to_orig_idx[map_orig_idx_to_var_idx]

    # Reorder everything that is in "Final" order to "Variation" order
    pt = pt[reorder_final_to_var]
    eta = eta[reorder_final_to_var]
    phi = phi[reorder_final_to_var]
    energy = energy[reorder_final_to_var]
    jer_factor = jer_factor[reorder_final_to_var]
    jec_unc = jec_unc[reorder_final_to_var]

    corr = 1. / jer_factor * (1.+jec_unc) * variation_jer_factor
    return pt*corr, eta, phi, energy*corr, reorder_final_to_var


def calc_jer_variation(
        pt, eta, phi, energy,
        jer_factor, orig_idx,
        variation_orig_idx, variation_jer_factor
        ):
    """
    Applies a JER variation (up or down) on a 4-vector.

    Note there are 3 'ordering levels':
    - "Final": the final ordering of jets after centrally applied corrections
    - "Original": the ordering of jets _before_ any corrections
    - "Variation": the final ordering of jets after the applying the correction of the
        _variation_

    The algorithm below first creates a map to reorder "Final" to "Variation", then
    applies the correction after ordering everything in "Variation" ordering.

    Args:
        pt (ak.Array): jet pt
        eta (ak.Array): jet eta
        phi (ak.Array): jet phi
        energy (ak.Array): jet energy
        jer_factor (ak.Array): the JER factor that was applied centrally to obtain the
            final jet
        orig_idx (ak.Array): mapping of final corrected jet ordering back to 'original'
            ordering
        variation_orig_idx (ak.Array): mapping of variation ordering back to 'original'
            ordering
        variation_jer_factor (ak.Array): the variation's JER factor

    Returns:
        (ak.Array, ak.Array, ak.Array, ak.Array) : pt, eta, phi, and energy after
            applying the variation and reordered by the pT after variation.
    """

    # Create a map to reorder final corrected jets to the ordering of the variation
    map_orig_idx_to_var_idx = ak.argsort(variation_orig_idx, axis=-1)
    map_final_idx_to_orig_idx = orig_idx
    reorder_final_to_var = map_final_idx_to_orig_idx[map_orig_idx_to_var_idx]

    # Reorder everything that is in "Final" order to "Variation" order
    pt = pt[reorder_final_to_var]
    eta = eta[reorder_final_to_var]
    phi = phi[reorder_final_to_var]
    energy = energy[reorder_final_to_var]
    jer_factor = jer_factor[reorder_final_to_var]

    corr = 1. / jer_factor * variation_jer_factor
    return pt*corr, eta, phi, energy*corr, reorder_final_to_var



###############################
####### PFNano section ########
###############################


def make_jets_for_jerc(events,jet_coll, event_rho, correction_key):
        
        jets = {}
        jets["pt_raw"] = (1 - events[f"{jet_coll}_rawFactor"])*events[f"{jet_coll}_pt"]
        jets["mass_raw"] = (1 - events[f"{jet_coll}_rawFactor"])*events[f"{jet_coll}_mass"]
        jets["event_rho"] = ak.broadcast_arrays(event_rho, events[f"{jet_coll}_pt"])[0]
        if "NOJER" not in correction_key:
            if "FatJet" in jet_coll:
                m_genMatch_dR2max = 0.4*0.4
                jets["pt_gen"] = gen_matching_tools.get_matched_gen_jets(events[f"{jet_coll}_pt"], 
                                                                         events[f"{jet_coll}_eta"], 
                                                                         events[f"{jet_coll}_phi"], 
                                                                         events[f"GenJetAK8_pt"], 
                                                                         events[f"GenJetAK8_eta"], 
                                                                         events[f"GenJetAK8_phi"], 
                                                                         events[f"FatJet_genJetAK8Idx"],
                                                                         m_genMatch_dR2max,
                                                                         correction_key.replace("NOJEC",""),
                                                                         jet_coll,
                                                                         jets["event_rho"]
                                                                         )
            else:
                m_genMatch_dR2max = 0.2*0.2
                jets["pt_gen"] = gen_matching_tools.get_matched_gen_jets(events[f"{jet_coll}_pt"], 
                                                                         events[f"{jet_coll}_eta"], 
                                                                         events[f"{jet_coll}_phi"], 
                                                                         events[f"Gen{jet_coll}_pt"], 
                                                                         events[f"Gen{jet_coll}_eta"], 
                                                                         events[f"Gen{jet_coll}_phi"], 
                                                                         events[f"{jet_coll}_genJetIdx"],
                                                                         m_genMatch_dR2max,
                                                                         correction_key.replace("NOJEC",""),
                                                                         jet_coll,
                                                                         jets["event_rho"]
                                                                         )


        #make jets an ak array
        if "NOJER" not in correction_key:
            jets = ak.zip({
                "pt": events[f"{jet_coll}_pt"],
                "eta": events[f"{jet_coll}_eta"],
                "phi": events[f"{jet_coll}_phi"],
                "mass": events[f"{jet_coll}_mass"],
                "area": events[f"{jet_coll}_area"],
                "pt_raw": jets["pt_raw"],
                "mass_raw": jets["mass_raw"],
                "event_rho": jets["event_rho"],
                "pt_gen": ak.values_astype(ak.fill_none(jets["pt_gen"], 0), np.float32),
            })
        else:
            jets = ak.zip({
                "pt": events[f"{jet_coll}_pt"],
                "eta": events[f"{jet_coll}_eta"],
                "phi": events[f"{jet_coll}_phi"],
                "mass": events[f"{jet_coll}_mass"],
                "area": events[f"{jet_coll}_area"],
                "pt_raw": jets["pt_raw"],
                "mass_raw": jets["mass_raw"],
                "event_rho": jets["event_rho"],
            })
             
        return jets

#CZZ: hardcoded PFMET for now
def make_met_for_jerc(events):
        
        met = {}

        met["pt"] = events["MET_pt"]
        met["phi"] = events["MET_phi"]
        met["MetUnclustEnUpDeltaX"] = events["MET_MetUnclustEnUpDeltaX"]
        met["MetUnclustEnUpDeltaY"] = events["MET_MetUnclustEnUpDeltaY"]

        #make met an ak array
        met = ak.zip({
            "pt": met["pt"],
            "phi": met["phi"],
            "MetUnclustEnUpDeltaX": met["MetUnclustEnUpDeltaX"],
            "MetUnclustEnUpDeltaY": met["MetUnclustEnUpDeltaY"],
        })
        
             
        return met



def calc_jerc_variations_PFNano(
    events: ak.Array,
    year: str,
    run: str,
    jet_coll: str,
    jerc_variations: dict,
    variation: str,
    jerc_cache: dict,
    ) -> ak.Array:
    
    jet_factory = jerc_variations[f'{jet_coll.lower()}_factory']
    met_factory = jerc_variations['met_factory']


    # calculate all variables needed as inputs
    jerc_key_label = ""
    access_jerc_corr_jets =  ""
    if "jec" in variation:
        jerc_key_label = "NOJER"
        access_jerc_corr_jets = "JES_jes"
    if "jer" in variation:
        jerc_key_label = "NOJEC"
        access_jerc_corr_jets = "JER"


    #build jet corrections
    correction_key = None
    if "2016" in year:
        if "APV" in year:
            correction_key = f"{year.replace('APV','preVFP')}{run.lower()}{jerc_key_label}" 
        else:
            correction_key = f"{year.replace(year,'2016postVFP')}{run.lower()}{jerc_key_label}"
    else:
        correction_key = f"{year}{run.lower()}{jerc_key_label}" 

    rhos = events.fixedGridRhoFastjetAll
    jets_corrected = jet_factory[correction_key].build(make_jets_for_jerc(events,jet_coll, rhos, correction_key), jerc_cache)

    met =  None
    if jet_coll == "Jet":
        met = met_factory.build(make_met_for_jerc(events), jets_corrected, {})
    
        
    #extract the direction of the variation
    direction = "up" if "up" in variation else "down"

    #make dummy x_ratio
    x_ratio = ak.ones_like(jets_corrected.pt)

    #extract corrections to jets
    pt_corr = eval(f"jets_corrected.{access_jerc_corr_jets}.{direction}.pt")
    eta_corr = eval(f"jets_corrected.{access_jerc_corr_jets}.{direction}.eta")
    phi_corr = eval(f"jets_corrected.{access_jerc_corr_jets}.{direction}.phi")
    mass_corr = eval(f"jets_corrected.{access_jerc_corr_jets}.{direction}.mass")

    #extract corrections to met
    met_ptcorr = None
    met_phicorr = None
    if jet_coll == "Jet":
        met_ptcorr = eval(f"met.{access_jerc_corr_jets}.{direction}.pt")
        met_phicorr = eval(f"met.{access_jerc_corr_jets}.{direction}.phi")

    return pt_corr, eta_corr, phi_corr, mass_corr, met_ptcorr, met_phicorr, x_ratio


def calc_unclustered_met_variations_PFNano(
    events: ak.Array,
    year: str,
    run: str,
    jet_coll: str,
    jerc_variations: dict,
    variation: str,
    jerc_cache: dict,
    ) -> ak.Array:
    
    jet_factory = jerc_variations[f'{jet_coll.lower()}_factory']
    met_factory = jerc_variations['met_factory']

    #build jet corrections
    correction_key = None
    if "2016" in year:
        if "APV" in year:
            correction_key = f"{year.replace('APV','preVFP')}{run.lower()}" 
        else:
            correction_key = f"{year.replace(year,'2016postVFP')}{run.lower()}"
    else:
        correction_key = f"{year}{run.lower()}"

    rhos = events.fixedGridRhoFastjetAll
    jets_corrected = jet_factory[correction_key].build(make_jets_for_jerc(events,jet_coll, rhos, correction_key), jerc_cache)
    met = met_factory.build(make_met_for_jerc(events), jets_corrected, {})
    
        
    #extract the direction of the variation
    direction = "up" if "up" in variation else "down"

    #extract corrections to met
    met_ptcorr = None
    met_phicorr = None

    met_ptcorr = eval(f"met.MET_UnclusteredEnergy.{direction}.pt")
    met_phicorr = eval(f"met.MET_UnclusteredEnergy.{direction}.phi")

    return met_ptcorr, met_phicorr


def calc_custom_svj_jes_variations_PFNano(
                events: ak.Array,
                year: str,
                run: str,
                jet_coll: str,
                jerc_variations: dict,
                variation: str,
                ) -> ak.Array:
    

    #build jet corrections
    correction_key = None
    if "2016" in year:
        if "APV" in year:
            correction_key = f"{year.replace('APV','preVFP')}{run.lower()}" 
        else:
            correction_key = f"{year.replace(year,'2016postVFP')}{run.lower()}"
    else:
        correction_key = f"{year}{run.lower()}"

    #extract the direction of the variation
    direction = "up" if "up" in variation else "down"

    #build custom jec calculator
    svjJESCalc = SVJCustomJESCalculator(
        events,
        jet_coll,
        correction_key
    )

    #fetch correction
    svj_jecs, x_ratio = svjJESCalc.getVariation(direction)
    jets_corrected = ak.zip({
        "pt": events[f"{jet_coll}_pt"]*svj_jecs,
        "phi": events[f"{jet_coll}_phi"],
        "pt_raw": events[f"{jet_coll}_pt"],
    })


    #now propagate custom jes to met (T1-like correction)
    met_factory = jerc_variations['met_factory']
    if jet_coll == "Jet":
        met = met_factory.build(make_met_for_jerc(events), jets_corrected, {})


    #extract corrections to met
    met_ptcorr = None
    met_phicorr = None
    if jet_coll == "Jet":
        met_ptcorr = eval(f"met.MET_UnclusteredEnergy.{direction}.pt")
        met_phicorr = eval(f"met.MET_UnclusteredEnergy.{direction}.phi")


    return jets_corrected.pt, events[f"{jet_coll}_eta"], jets_corrected.phi, events[f"{jet_coll}_mass"], met_ptcorr, met_phicorr, x_ratio
    



    