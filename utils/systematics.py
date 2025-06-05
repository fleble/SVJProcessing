import awkward as ak
import numpy as np
import copy as cp
import utils.gen_matching_tools as gen_matching_tools
from utils.jet_energy_scale_svj_factory import SVJCustomJESCalculator 
from utils.met_jecs_factory import update_met_t1_corr, apply_uncl_variation_to_met_t1
from coffea.lookup_tools import extractor
from coffea.jetmet_tools import CorrectedMETFactory
from utils.met_significance_factory_pfnano import MetSignificanceCalculator

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

        met["pt"] = events["RawMET_pt"]    #RawMET_pt
        met["phi"] = events["RawMET_phi"] #RawMET_phi
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




def apply_jers_PFNano(
    events: ak.Array,
    year: str,
    run: str,
    jet_coll: str,
    jerc_variations: dict,
    jerc_cache: dict,
    ) -> ak.Array:
    
    jet_factory = jerc_variations[f'{jet_coll.lower()}_factory']


    # calculate all variables needed as inputs
    jerc_key_label = ""

    jerc_key_label = "NOJEC"

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

    #extract corrections to jets
    pt_corr = jets_corrected.pt_jer
    eta_corr = events[f"{jet_coll}_eta"]
    phi_corr = events[f"{jet_coll}_phi"]
    mass_corr = jets_corrected.mass_jer


    return pt_corr, eta_corr, phi_corr, mass_corr


def apply_jecs_PFNano(
    events: ak.Array,
    year: str,
    run: str,
    jet_coll: str,
    jerc_variations: dict,
    jerc_cache: dict,
    ) -> ak.Array:
    
    jet_factory = jerc_variations[f'{jet_coll.lower()}_factory']


    # calculate all variables needed as inputs
    jerc_key_label = "NOJER"


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


    #extract corrections to jets
    pt_corr = jets_corrected.pt_jec
    eta_corr = events[f"{jet_coll}_eta"]
    phi_corr = events[f"{jet_coll}_phi"]
    mass_corr = jets_corrected.mass_jec

    #define new raw factor
    pt_raw = (1. - events[f"{jet_coll}_rawFactor"])*events[f"{jet_coll}_pt"]
    raw_factor = 1. - pt_raw/pt_corr


    return pt_corr, eta_corr, phi_corr, mass_corr,raw_factor

#Apply JECs and JERs to nominal jets 
def apply_jercs_PFNano(
    events: ak.Array,
    year: str,
    run: str,
    jet_coll: str,
    jerc_variations: dict,
    jerc_cache: dict,
    ) -> ak.Array:
    
    jet_factory = jerc_variations[f'{jet_coll.lower()}_factory']

    # calculate all variables needed as inputs, apply both JEC and JER
    jerc_key_label = ""

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


    #extract corrections to jets
    pt_corr = jets_corrected.pt
    eta_corr = jets_corrected.eta
    phi_corr = jets_corrected.phi
    mass_corr = jets_corrected.mass


    return pt_corr, eta_corr, phi_corr, mass_corr


#Propage JECs to MET
def propagate_jecs_to_MET_PFNano(
                events: ak.Array,
                year: str,
                run: str,
                jet_coll: str,
                jerc_variations: dict,
                jerc_cache: dict,
                ) -> ak.Array:
    
    jet_factory = jerc_variations[f'{jet_coll.lower()}_factory']


    # calculate all variables needed as inputs
    jerc_key_label = "NOJER"


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
    jets_corrected_nom = jet_factory[correction_key].build(make_jets_for_jerc(events,jet_coll, rhos, correction_key), jerc_cache)

    #extract corrected jets
    jet_pt_corr_nom  = jets_corrected_nom.pt
    jet_phi_corr_nom  = events[f"{jet_coll}_phi"]

    #extract raw jet pt
    jet_pt_raw = jets_corrected_nom.pt_raw

    met_pt_raw = events["RawMET_pt"]
    met_phi_raw = events["RawMET_phi"]

    corr_t1_met =  update_met_t1_corr(met_pt_raw, met_phi_raw, jet_pt_corr_nom, jet_phi_corr_nom, jet_pt_raw)
    
    return corr_t1_met.pt, corr_t1_met.phi



#Propage JECs to MET
def propagate_jecs_to_METSig_PFNano(
                events: ak.Array,
                year: str,
                run: str,
                jet_coll: str,
                jerc_variations: dict,
                jerc_cache: dict,
                make_unclustered_En_var: bool = False,
                variation: str = None,
                ) -> ak.Array:
    
    jet_factory = jerc_variations[f'{jet_coll.lower()}_factory']

    # calculate all variables needed as inputs
    jerc_key_label = "NOJER"

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
    jets_corrected_nom = jet_factory[correction_key].build(make_jets_for_jerc(events,jet_coll, rhos, correction_key), jerc_cache)

    #extract corrected jets
    jet_pt_corr_nom  = jets_corrected_nom.pt

    #copy events to avoid modifying the original array
    events_copy = cp.deepcopy(events)

    #correct the copy of events with corrected jets
    #adding the JEC varied jets to the events
    jes_permutation = ak.argsort(jet_pt_corr_nom, ascending=False)
    jes_corr_pt = jet_pt_corr_nom[jes_permutation]
    jes_corr_eta = events_copy[f"{jet_coll}_eta"][jes_permutation]
    jes_corr_phi = events_copy[f"{jet_coll}_phi"][jes_permutation]
    jes_corr_mass = events_copy[f"{jet_coll}_mass"][jes_permutation]
    events_copy[f"{jet_coll}_pt"] = jes_corr_pt
    events_copy[f"{jet_coll}_eta"] = jes_corr_eta
    events_copy[f"{jet_coll}_phi"] = jes_corr_phi
    events_copy[f"{jet_coll}_mass"] = jes_corr_mass
    jes_pf_cand_jet_idx = events_copy[f"{jet_coll}PFCands_jetIdx"]
    jes_sorted_pf_cand_jet_idx = ak.Array([p[idx] for idx, p in zip(jes_pf_cand_jet_idx, jes_permutation)])
    events_copy[f"{jet_coll}PFCands_jetIdx"] = jes_sorted_pf_cand_jet_idx

    variation_direction = ""
    #extract variation direction
    if make_unclustered_En_var:
        #extract from variation
        variation_direction = "up" if "up" in variation else "down"

    met_sig_recalculator = MetSignificanceCalculator(events_copy,
                            year,
                            run,
                            make_unclustered_En_var,
                            variation_direction,
                        ) 
 
     
    met_sig_corr = met_sig_recalculator.getSignificance() 

    #delete the copy of events
    del events_copy

    return met_sig_corr




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

    # calculate all variables needed as inputs
    jerc_key_label = ""
    access_jerc_corr_jets =  ""

    if "jec" in variation:
        access_jerc_corr_jets = "JES_jes"
    if "jer" in variation:
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

        
    #extract the direction of the variation
    direction = "up" if "up" in variation else "down"

    #make dummy x_ratio
    x_ratio = ak.ones_like(jets_corrected.pt)

    #extract corrections to jets
    pt_corr = eval(f"jets_corrected.{access_jerc_corr_jets}.{direction}.pt")
    eta_corr = eval(f"jets_corrected.{access_jerc_corr_jets}.{direction}.eta")
    phi_corr = eval(f"jets_corrected.{access_jerc_corr_jets}.{direction}.phi")
    mass_corr = eval(f"jets_corrected.{access_jerc_corr_jets}.{direction}.mass")

    pt_raw_var = eval(f"jets_corrected.{access_jerc_corr_jets}.{direction}.pt_raw")

    met =  None
    if jet_coll == "Jet":
        met = update_met_t1_corr(events["RawMET_pt"], events["RawMET_phi"], pt_corr, phi_corr, pt_raw_var)

    #extract corrections to met
    met_ptcorr = None
    met_phicorr = None
    if jet_coll == "Jet":
        met_ptcorr = met.pt
        met_phicorr = met.phi


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

    # calculate all variables needed as inputs
    jerc_key_label = "NOJER"


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
    jets_corrected_nom = jet_factory[correction_key].build(make_jets_for_jerc(events,jet_coll, rhos, correction_key), jerc_cache)

    #extract corrected jets
    jet_pt_corr_nom  = jets_corrected_nom.pt
    jet_phi_corr_nom  = events[f"{jet_coll}_phi"]

    #extract raw jet pt
    jet_pt_raw = jets_corrected_nom.pt_raw

    met_pt_raw = events["RawMET_pt"]
    met_phi_raw = events["RawMET_phi"]

    #get T1 corrected MET (nominal)
    corr_t1_met = update_met_t1_corr(met_pt_raw, met_phi_raw, jet_pt_corr_nom, jet_phi_corr_nom, jet_pt_raw)

    #extract the direction of the variation
    direction = "up" if "up" in variation else "down"

    positive = None
    if (direction == "up"):
        positive = True
    if (direction == "down"):
        positive = False

    #unpack corr_t1_met
    met_pt_t1_nom = corr_t1_met.pt
    met_phi_t1_nom = corr_t1_met.phi

    #extracting correction due to unclustered energy variation
    corr_t1_met_uncl_var = apply_uncl_variation_to_met_t1(met_pt_t1_nom,met_phi_t1_nom, positive=positive, dx=events["MET_MetUnclustEnUpDeltaX"], dy=events["MET_MetUnclustEnUpDeltaY"])

    return corr_t1_met_uncl_var.pt, corr_t1_met_uncl_var.phi



def calc_custom_svj_jes_variations_PFNano(
                events: ak.Array,
                year: str,
                run: str,
                jet_coll: str,
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
    #need to retrieve pT raw pT to propagate to MET 
    svj_jecs, x_ratio = svjJESCalc.getVariation(direction)
    jets_corrected = ak.zip({
        "pt": events[f"{jet_coll}_pt"]*svj_jecs,
        "phi": events[f"{jet_coll}_phi"],
        "pt_raw": events[f"{jet_coll}_pt"],
    })


    #now propagate custom jes to met (T1-like correction) if jet_coll is Jet
    met_ptcorr = None
    met_phicorr = None
    if jet_coll == "Jet":
        corr_met_custom_jecs = update_met_t1_corr(events["MET_pt"], events["MET_phi"], jets_corrected.pt, jets_corrected.phi, jets_corrected.pt_raw)
        met_ptcorr = corr_met_custom_jecs.pt
        met_phicorr = corr_met_custom_jecs.phi

    return jets_corrected.pt, events[f"{jet_coll}_eta"], jets_corrected.phi, events[f"{jet_coll}_mass"]*svj_jecs, met_ptcorr, met_phicorr, x_ratio
    



    