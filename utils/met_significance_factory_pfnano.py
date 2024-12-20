import awkward as ak
import numba as nb
import numpy as np
import utils.jet_resolution_utils_pfnano as jet_res_utils
import skimmer.skimmer_utils as skimmer_utils


#Class adapted from CMSSW implementation: https://github.com/cms-sw/cmssw/blob/48600daa76c0dd3925eb34d216bc756ed556016d/RecoMET/METAlgorithms/src/METSignificance.cc
#Variations procedure for unclustered energy taken from: https://github.com/Lvigilante/nanoAOD-tools/blob/MetSig/python/postprocessing/modules/jme/METSigProducer.py 
#Metsignificance tuning parameters taken from:  https://github.com/cms-sw/cmssw/blob/92333e3acb2db4776c5772fbd71e3ef325a9b46e/RecoMET/METProducers/python/METSignificanceParams_cfi.py#L19 

nb.jit(nopython=True)
class MetSignificanceCalculator:
    
    def __init__(self, 
                 events,
                 year,
                 run,
                 make_unclustered_En_var = False,
                 variation_direction = "",
                 jet_coll="Jet",
                 met_coll="MET",
                 electron_coll="Electron",
                 muon_coll="Muon",
                 dr_match=0.4,
                 jet_threshold=15.,
                 jet_params=[1.39,1.26,1.21,1.23,1.28],
                 pjet_params=[-0.2586,0.6173],
                 ):

                correction_key = None
                if "2016" in year:
                    if "APV" in year:
                        correction_key = f"{year.replace('APV','preVFP')}{run.lower()}" 
                    else:
                        correction_key = f"{year.replace(year,'2016postVFP')}{run.lower()}"
                else:
                    correction_key = f"{year}{run.lower()}" 


                #get objects for calculation from from events
                
                #CZZ: MET must come from jerc varied collections, otherwise nominal when running the unclustered energy variation
                met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
                    pt=eval(f"events.{met_coll}_pt"),
                    phi=eval(f"events.{met_coll}_phi"),
                )

                #get jets and leptons
                jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
                    pt=eval(f"events.{jet_coll}_pt"),
                    eta=eval(f"events.{jet_coll}_eta"),
                    phi=eval(f"events.{jet_coll}_phi"),
                    mass=eval(f"events.{jet_coll}_mass"),
                )

                #access resolution objects
                rho = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, events[f"{jet_coll}_pt"])[0]
                resobj_pt = jet_res_utils.fetch_resolution(correction_key,jet_coll,var="pt")
                resobj_phi = jet_res_utils.fetch_resolution(correction_key,jet_coll,var="phi")

                #get pt resolution and phi resolution
                pt_resolutions = resobj_pt.getResolution(JetEta=jets.eta, Rho=rho, JetPt=jets.pt)
                phi_resolutions = resobj_phi.getResolution(JetEta=jets.eta, Rho=rho, JetPt=jets.pt)

                #add pt and phi resolutions to the jets
                jets = ak.zip(
                    {
                        "pt": jets.pt,
                        "eta": jets.eta,
                        "phi": jets.phi,
                        "mass": jets.mass,
                        "dpt": pt_resolutions,
                        "dphi": phi_resolutions,
                    }
                )

                electrons = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
                    pt=eval(f"events.{electron_coll}_pt"),
                    eta=eval(f"events.{electron_coll}_eta"),
                    phi=eval(f"events.{electron_coll}_phi"),
                    mass=eval(f"events.{electron_coll}_mass"),
                )

                muons = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
                    pt=eval(f"events.{muon_coll}_pt"),
                    eta=eval(f"events.{muon_coll}_eta"),
                    phi=eval(f"events.{muon_coll}_phi"),
                    mass=eval(f"events.{muon_coll}_mass"),
                )

                self.events = events
                self.jets  = jets
                self.electrons = electrons
                self.muons = muons
                self.met = met
                self.sum_pt_unclustered = eval(f"events.{met_coll}_sumPtUnclustered")
                self.met_unclust_EnUp_DeltaX = eval(f"events.{met_coll}_MetUnclustEnUpDeltaX")
                self.met_unclust_EnUp_DeltaY = eval(f"events.{met_coll}_MetUnclustEnUpDeltaY")
                self.dR2match = dr_match ** 2
                self.jetThreshold = jet_threshold
                self.jetParams = jet_params
                self.pjetParams = pjet_params
                self.make_unclustered_En_var = make_unclustered_En_var
                self.variation_direction = variation_direction


    
    def _getBin(self, abseta):
        etabins = [0.8,1.3,1.9,2.5,100]
        for i, a in enumerate(etabins):
            if abseta < a:
                return int(i)
                break
                

    def _deltaPhi(self, phi1, phi2):
        dphi = phi2-phi1
        if  dphi > np.pi:
            dphi -= 2.0*np.pi
        if dphi <= -np.pi:
            dphi += 2.0*np.pi
        return abs(dphi)

    def _deltaR2(self, l1, l2):
        return self._deltaPhi(l1.phi, l2.phi)**2 + (l1.eta - l2.eta)**2

    def _deltaR(self, l1, l2):
        return np.sqrt(self.deltaR2(l1,l2))
    
    def _cleanJet(self, jet, leptons):

        for lep in leptons:
            if self._deltaR2(lep, jet) < self.dR2match and lep.pt > 10:
                return False
        
        return True
    
    
    def _getCovarianceEv(self,ev_idx):
        
        #get jets and leptons
        jets_ev = self.jets[ev_idx]
        all_leptons_ev = [self.electrons[ev_idx],self.muons[ev_idx]]

        #get unclustered energy
        sumPtUnclustered_ev = self.sum_pt_unclustered[ev_idx]
        unclust_EnUp_DeltaX_ev = self.met_unclust_EnUp_DeltaX[ev_idx]
        unclust_EnUp_DeltaY_ev = self.met_unclust_EnUp_DeltaY[ev_idx]

        #define total sum pt from unclustered
        totalSumPt_ev = sumPtUnclustered_ev

        #get cleaned jets
        cleanedJets_ev = []
        for j in jets_ev:
                for leptons in all_leptons_ev:
                    is_cleaned_jet = self._cleanJet(j, leptons)
                if is_cleaned_jet:
                    cleanedJets_ev += [j]

        # metsig covariance
        cov_xx = 0
        cov_xy = 0
        cov_yy = 0
        
        for j in cleanedJets_ev:
                eta_bin_idx = self._getBin(abs(j.eta))
                # split into high/low pt jets
                if j.pt > self.jetThreshold:
                    scale = 0
                    if (eta_bin_idx == 0):
                        scale = self.jetParams[0]
                    elif (eta_bin_idx == 1):
                        scale = self.jetParams[1]
                    elif (eta_bin_idx == 2):
                        scale = self.jetParams[2]
                    elif (eta_bin_idx == 3):
                        scale = self.jetParams[3]
                    else:
                        scale = self.jetParams[4]

                    
                    cj = np.cos(j.phi)
                    sj = np.sin(j.phi)
                    dpt = scale * j.pt * j.dpt
                    dph = j.pt * j.dphi

                    cov_xx += dpt * dpt * cj * cj + dph * dph * sj * sj
                    cov_xy += (dpt * dpt - dph * dph) * cj * sj
                    cov_yy += dph * dph * cj * cj + dpt * dpt * sj * sj

                else:
                    #add cleaned jets not passing the threshold to the unclustered energy component
                    totalSumPt_ev += j.pt
        
        if self.make_unclustered_En_var:
            sumPtUnclustEnDeltaUp = np.sqrt(unclust_EnUp_DeltaX_ev*unclust_EnUp_DeltaX_ev  + unclust_EnUp_DeltaY_ev*unclust_EnUp_DeltaY_ev )

            if self.variation_direction == "up":
                totalSumPt_ev += sumPtUnclustEnDeltaUp 

            if self.variation_direction == "down":
                totalSumPt_ev -= sumPtUnclustEnDeltaUp


        #protection against unphysical events
        if (sumPtUnclustered_ev < 0):
            sumPtUnclustered_ev = 0


        #add unclustered energy to covariance matrix
        cov_tt = self.pjetParams[0] * self.pjetParams[0] + self.pjetParams[1] * self.pjetParams[1] * totalSumPt_ev

        cov_xx += cov_tt
        cov_yy += cov_tt

        cov = [[cov_xx, cov_xy], [cov_xy, cov_yy]]

        return cov
    
    
    def _getSignificanceEv(self,cov, ev_idx):

        # Covariance matrix determinant
        det = cov[0][0] * cov[1][1] - cov[0][1] * cov[1][0]

        # Inverse matrix
        ncov_xx = cov[1][1] / det
        ncov_xy = -cov[0][1] / det
        ncov_yy = cov[0][0] / det

        # Product of met and inverse of covariance
        sig = self.met.px[ev_idx] * self.met.px[ev_idx] * ncov_xx + 2 * self.met.px[ev_idx] * self.met.py[ev_idx] * ncov_xy + self.met.py[ev_idx] * self.met.py[ev_idx] * ncov_yy

        return sig
    
    
    def _build(self):
        builder = ak.ArrayBuilder()
        for ev_idx,_ in enumerate(self.events):
            builder.begin_list()
            cov_matrix_ev = self._getCovarianceEv(ev_idx)
            significance_ev = self._getSignificanceEv(cov_matrix_ev, ev_idx)
            builder.append(significance_ev)
            builder.end_list()
        return builder



    def getSignificance(self):
        
        return self._build()
    
    
 
            


