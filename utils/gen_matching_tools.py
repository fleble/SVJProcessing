import os
import numpy as np
import awkward as ak
import numba as nb
import utils.jet_resolution_utils_pfnano as jet_res_utils



def phi_mpi_pi(angle):
        
    # Initial check
    if -np.pi <= angle <= np.pi:
        return angle
    
    # Normalization process
    if angle > 0:
        n = int(0.5 * (angle / np.pi + 1))
        angle -= 2 * n * np.pi
    else:
        n = int(-0.5 * (angle / np.pi - 1))
        angle += 2 * n * np.pi
    
    return angle



def get_dr2(phi1, eta1, phi2, eta2):

    dphi = phi_mpi_pi(phi1 - phi2)
    deta = eta1 - eta2
    return dphi*dphi + deta*deta


nb.jit(nopython=True)
def get_matched_gen_jets_(builder, jet_pt_reco , jet_eta_reco, jet_phi_reco,  jet_pt_gen, jet_eta_gen, jet_phi_gen, genJetIdx, m_genMatch_dR2max, pt_resolution):
    for j_pts_reco, j_etas_reco, j_phis_reco, j_pts_gen, j_etas_gen, j_phis_gen, j_reco_genIdxs, j_resos_reco in zip(jet_pt_reco , jet_eta_reco, jet_phi_reco, jet_pt_gen, jet_eta_gen, jet_phi_gen, genJetIdx, pt_resolution):
        builder.begin_list()
        for reco_j_idx,j_genIdx in enumerate(j_reco_genIdxs):
            if (j_genIdx >= 0) and (j_genIdx < len(j_pts_gen)):
                dr2 = get_dr2(j_phis_reco[reco_j_idx], j_etas_reco[reco_j_idx], j_phis_gen[j_genIdx], j_etas_gen[j_genIdx])
                if ((dr2 < m_genMatch_dR2max) and jet_res_utils.check_pt_resolution(j_pts_reco[reco_j_idx], j_pts_gen[j_genIdx],j_resos_reco[reco_j_idx])):
                    builder.append(j_pts_gen[j_genIdx])
                else:
                    builder.append(None)
            else:
                builder.append(None)
        builder.end_list()
    return builder


def get_matched_gen_jets(jet_pt_reco , jet_eta_reco, jet_phi_reco, jet_pt_gen, jet_eta_gen, jet_phi_gen, genJetIdx, m_genMatch_dR2max, correction_key,jet_coll,rho):
    pt_resolution = jet_res_utils.get_jets_resolution(jet_pt_reco, jet_eta_reco, rho, correction_key,jet_coll,var="pt")
    builder = ak.ArrayBuilder()
    return get_matched_gen_jets_(builder, jet_pt_reco , jet_eta_reco, jet_phi_reco,  jet_pt_gen, jet_eta_gen, jet_phi_gen, genJetIdx, m_genMatch_dR2max, pt_resolution).snapshot()