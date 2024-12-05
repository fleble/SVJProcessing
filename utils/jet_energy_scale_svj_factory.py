import awkward as ak
import numba as nb
import numpy as np
import utils.gen_matching_tools as gen_matching_tools

nb.jit(nopython=True)
class SVJCustomJESCalculator:

    def __init__(self, 
                 events,
                 jet_coll,
                 std_jer_key
                 ):
        


        self.events = events
        self.jet_coll = jet_coll
        self.std_jer_key = std_jer_key

    

    def _build(self):
        
        jets = {}
        jets["event_rho"] = ak.broadcast_arrays(self.events.fixedGridRhoFastjetAll, self.events[f"{self.jet_coll}_pt"])[0]
        jets["pt_reco"] = self.events[f"{self.jet_coll}_pt"]

        if "FatJet" in self.jet_coll:
            m_genMatch_dR2max = 0.4*0.4
            jets["pt_gen"] = gen_matching_tools.get_matched_gen_jets(jets["pt_reco"], 
                                                                    self.events[f"{self.jet_coll}_eta"], 
                                                                    self.events[f"{self.jet_coll}_phi"], 
                                                                    self.events[f"GenJetAK8_pt"], 
                                                                    self.events[f"GenJetAK8_eta"], 
                                                                    self.events[f"GenJetAK8_phi"], 
                                                                    self.events[f"FatJet_genJetAK8Idx"],
                                                                    m_genMatch_dR2max,
                                                                    self.std_jer_key.replace("NOJEC",""),
                                                                    self.jet_coll,
                                                                    jets["event_rho"]
                                                                    )
            
        else:
            m_genMatch_dR2max = 0.2*0.2
            jets["pt_gen"] = gen_matching_tools.get_matched_gen_jets(jets["pt_reco"], 
                                                                    self.events[f"{self.jet_coll}_eta"], 
                                                                    self.events[f"{self.jet_coll}_phi"], 
                                                                    self.events[f"Gen{self.jet_coll}_pt"], 
                                                                    self.events[f"Gen{self.jet_coll}_eta"], 
                                                                    self.events[f"Gen{self.jet_coll}_phi"], 
                                                                    self.events[f"{self.jet_coll}_genJetIdx"],
                                                                    m_genMatch_dR2max,
                                                                    self.std_jer_key.replace("NOJEC",""),
                                                                    self.jet_coll,
                                                                    jets["event_rho"]
                                                                    )

        
        jets = ak.zip({
                "pt": self.events[f"{self.jet_coll}_pt"],
                "pt_gen": ak.values_astype(ak.fill_none(jets["pt_gen"], 0), np.float32),
            })

        #compute ratio of gen pt to reco pt
        x_jes = np.abs(jets.pt_gen / jets.pt - ak.ones_like(jets.pt))
        svj_jecsup = ak.ones_like(jets.pt) + x_jes
        svj_jecdown = ak.ones_like(jets.pt) - x_jes
        #if svj_jecdown has negative values, set them to 0
        svj_jecdown = ak.where(svj_jecdown < 0, ak.zeros_like(svj_jecdown), svj_jecdown)
        

        return svj_jecsup , svj_jecdown


    def getVariation(self, direction):
        svj_jecsup , svj_jecdown =  self._build()
        if direction == "up":
            return svj_jecsup
        elif direction == "down":
            return svj_jecdown

        
        


    