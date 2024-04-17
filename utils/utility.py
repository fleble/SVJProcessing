import numpy as np
import awkward as ak
# from mt2 import mt2
from . import objects as ob
from . import samples as s
from itertools import combinations
import time

def awkwardReshape(akArray,npArray):
    if len(akArray) == 0:
        return ak.Array([])
    else:
        return ak.broadcast_arrays(akArray.pt,1.0)[1] * npArray

def arrayConcatenate(array1,array2):
    if ak.any(array1) != True:
        return [array2]
    else:
        return ak.concatenate((array1, [array2]),axis=0)

# the two functions below return infinity when the event doesn't have the required
# number of jets or return the correct value for the jet variable
def jetVar_i(var,i,padValue=np.Inf):
    paddedVar = ak.fill_none(ak.pad_none(var,i+1),padValue)
    return paddedVar[:,i]

# convert phi values into spherical coordinate?
def phi_(x,y):
    phi = np.arctan2(y,x)
    return ak.where(phi < 0, phi + 2*np.pi, phi)

def deltaPhi(jetphiL,metPhiL):
    phi1 = phi_( np.cos(jetphiL), np.sin(jetphiL) )
    phi2 = phi_( np.cos(metPhiL), np.sin(metPhiL) )
    dphi = phi1 - phi2
    dphi_edited = ak.where(dphi < -np.pi, dphi + 2*np.pi, dphi)
    dphi_edited = ak.where(dphi_edited > np.pi, dphi_edited - 2*np.pi, dphi_edited)
    return abs(dphi_edited)

def deltaEta(eta0,eta1):
    return abs(eta0 - eta1)

def delta_R(eta0,eta1,phi0,phi1):
    dp = deltaPhi(phi0,phi1)
    deta = deltaEta(eta0,eta1)
    deltaR2 = deta * deta + dp * dp
    return np.sqrt(deltaR2)

def tauRatio(tau_a,tau_b,i):
    Ji_tau_a = jetVar_i(tau_a,i)
    Ji_tau_b = jetVar_i(tau_b,i)
    Ji_tau_ab = Ji_tau_a/Ji_tau_b
    return Ji_tau_ab

def mtMETLepton(met, metPhi, lepton):
    """Calculate the transverse mass of the MET and the leading Lepton"""
    leptonPhi = jetVar_i(lepton.phi,0)
    leptonPt = jetVar_i(lepton.pt,0)
    dphi = deltaPhi(leptonPhi,metPhi)
    mt = np.sqrt(2 * leptonPhi * met * (1 - np.cos(dphi)))
    return mt

## Need to update mt2 code.
# def lorentzVector(pt,eta,phi,mass,i):
#     px = pt[i]*np.cos(phi[i])
#     py = pt[i]*np.sin(phi[i])
#     pz = pt[i]*np.sinh(eta[i])
#     p2 = px**2 + py**2 + pz**2
#     m2 = mass[i]**2
#     energy = np.sqrt(m2+p2)
#     return np.array([energy,px,py,pz])
#
# def mass4vec(m4vec):
#     E = m4vec[0]
#     px = m4vec[1]
#     py = m4vec[2]
#     pz = m4vec[3]
#     return np.sqrt(E**2 - (px**2 + py**2 + pz**2))
#
# def eta4vec(m4vec):
#     px = m4vec[1]
#     py = m4vec[2]
#     pz = m4vec[3]
#     pt = np.sqrt(px**2 + py**2)
#     return np.arcsinh(pz/pt)
#
# def phi4vec(m4vec):
#     px = m4vec[1]
#     py = m4vec[2]
#     pt = np.sqrt(px**2 + py**2)
#     return np.arcsin(py/pt)
#
# def M_2J(j1,j2):
#     totJets = j1+j2
#     return mass4vec(totJets)
#
# def MT2Cal(FDjet0,FSMjet0,FDjet1,FSMjet1,met,metPhi):
#     Fjet0 = FDjet0 + FSMjet0
#     Fjet1 = FDjet1 + FSMjet1
#     METx = met*np.cos(metPhi)
#     METy = met*np.sin(metPhi)
#     MT2v = mt2(
#     mass4vec(Fjet0), Fjet0[1], Fjet0[2],
#     mass4vec(Fjet1), Fjet1[1], Fjet1[2],
#     METx, METy, 0.0, 0.0, 0
#     )
#     return MT2v
#
# def f4msmCom(pt,eta,phi,mass,met,metPhi,cut):
#     MT2 = 0
#     if len(pt) >= 4:
#         List4jets_3Com = [[0,1,2,3],[0,2,1,3],[0,3,1,2]]
#         diffList = []
#         jetList = []
#         for c in List4jets_3Com:
#             jc1 = lorentzVector(pt,eta,phi,mass,c[0])
#             jc2 = lorentzVector(pt,eta,phi,mass,c[1])
#             jc3 = lorentzVector(pt,eta,phi,mass,c[2])
#             jc4 = lorentzVector(pt,eta,phi,mass,c[3])
#             jetList.append([jc1,jc2,jc3,jc4])
#             diffList.append(abs(M_2J(jc1,jc2) - M_2J(jc3,jc4)))
#         comIndex = np.argmin(diffList)
#         msmJet = jetList[comIndex]
#         # applying angular cut
#         dEtaCut = 1.8
#         dPhiCut = 0.9
#         dRCutLow = 1.5
#         dRCutHigh = 4.0
#         dPhiMETCut = 1.5
#
#         if cut == "dEta":
#             if deltaEta(eta4vec(msmJet[0]),eta4vec(msmJet[1])) < dEtaCut and deltaEta(eta4vec(msmJet[2]),eta4vec(msmJet[3])) < dEtaCut:
#                 MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
#         elif cut == "dPhi":
#             if deltaPhiji(phi4vec(msmJet[0]),phi4vec(msmJet[1])) < dPhiCut and deltaPhiji(phi4vec(msmJet[2]),phi4vec(msmJet[3])) > dPhiCut:
#                 MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
#         elif cut == "dR":
#             if (dRCutLow < delta_R(eta4vec(msmJet[0]),eta4vec(msmJet[1]),phi4vec(msmJet[0]),phi4vec(msmJet[1])) < dRCutHigh) and (dRCutLow < delta_R(eta4vec(msmJet[2]),eta4vec(msmJet[3]),phi4vec(msmJet[2]),phi4vec(msmJet[3])) < dRCutHigh):
#                 MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
#         elif cut == "":
#             MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
#     return MT2
#
# def f4msmCom_vec(pt,eta,phi,mass,met,metPhi,cut):
#     vfunc = np.vectorize(lambda pt,eta,phi,mass,met,metPhi,cut: f4msmCom(pt,eta,phi,mass,met,metPhi,cut))
#     return vfunc(pt,eta,phi,mass,met,metPhi,cut)

def decode(hvCat,nlabel,clabel,catList):
    if hvCat >= nlabel:
        hvCat -= nlabel
        catList.append(clabel)
    return hvCat

def tch_hvCat_decode(hvCat):
    catList = []
    hvCat = decode(hvCat,16,"QsM",catList)
    hvCat = decode(hvCat,8,"QdM",catList)
    hvCat = decode(hvCat,4,"Gd",catList)
    hvCat = decode(hvCat,2,"Qd",catList)
    hvCat = decode(hvCat,1,"stableD",catList)
    return catList

def baselineVar(dataset,events,hemPeriod,sFactor):
    varVal = {}
    dataKeys = ["HTMHTData","JetHTData","METData","SingleElectronData","SingleMuonData","SinglePhotonData","EGammaData"]
    scaleFactor = s.sfGetter(dataset,scaleOn=sFactor)
    print("dataset = {} scaleon = {} scaleFactor = {}".format(dataset, sFactor, scaleFactor))
    isData = False
    for dKey in dataKeys:
        if dKey in dataset:
            isData = True
            break
    evtw = np.ones(len(events))
    if not isData:
         # 2018 lumi
        if "2016" in dataset:
            luminosity = 35921.036
        elif "2017" in dataset:
            luminosity = 41521.331
        elif "2018" in dataset:
            if hemPeriod == "PreHEM":
                luminosity = 21071.460
            elif hemPeriod == "PostHEM":
                luminosity = 38621.232
            else:
                luminosity = 59692.692
        evtw = luminosity*events.Weight*scaleFactor
    eCounter = np.where(evtw >= 0, 1, -1)
    obj = ob.Objects(events)
    jets = obj.goodJets()
    bjets = obj.goodBJets(dataset,jets)
    fjets = obj.goodFatJets()
    electrons = obj.goodElectrons()
    muons = obj.goodMuons()
    nonIsoMuons = obj.nonIsoMuons()
    met = events.MET
    metPhi = events.METPhi
    mtAK8 = events.MT_AK8
    ht = ak.sum(jets.pt,axis=1)
    st = ht + met
    # AK4 Jet Variables
    jetPhi = jets.phi
    dPhij = deltaPhi(jetPhi,metPhi)
    dPhiMinj = ak.min(dPhij,axis=1,mask_identity=False)
    # AK8 Jet Variables
    jetAK8Phi = fjets.phi
    dPhijAK8 = deltaPhi(jetAK8Phi,metPhi)
    dPhiMinjAK8 = ak.min(dPhijAK8,axis=1,mask_identity=False)
    # cr region 
    crElectrons = obj.crElectrons()
    crMuons     = obj.crMuons()

    if len(bjets) > 0:
        nBJets = ak.num(bjets)
    else:
        nBJets = np.zeros(len(evtw))

    isSignal = 0
    if "mMed" in dataset:
        if "s-channel" in dataset:
            isSignal = 1
        else:
            isSignal = 2
    varVal["isSignal"] = isSignal
    ## getting jet category for each jet
    jetCats = []
    if isSignal == 1:
        jetCats = ak.where(fjets.isHV,9,0) # treating SVJ from s-channel as QM jets and non-SVJ as SM jets

    elif isSignal == 2:
        GenJetsAK8 = events.GenJetsAK8
        jetsAK8GenInd = fjets.genIndex
        fjets = fjets[jetsAK8GenInd != -1]
        jetsAK8GenInd = jetsAK8GenInd[jetsAK8GenInd != -1]
        genHVCategory = GenJetsAK8.hvCategory
        jetCats = genHVCategory[jetsAK8GenInd]
    else:
        jetCats = awkwardReshape(fjets,np.ones(len(events))*-1)

    ############## keeping only certain jets for training ##############
    hvCond = ak.zeros_like(jetCats,dtype=bool)
    jetCatsUsed = [-1]
    if isSignal == 1:
        jetCatsUsed = [9] # s-channel SVJs
    elif isSignal == 2:
        jetCatsUsed = [3,9,11,5,7,13] # t-channel, all dark jets
        # jetCatsUsed = [3,9]
    for jetCat in jetCatsUsed:
        hvCond = hvCond | (jetCats == jetCat)
    djets = fjets[hvCond]
    ########################################################################

    ##### HLT muon matching for trigger study ####
    # hltMuons = obj.hltMuons
    # triggerOfflineMuons = obj.triggerOfflineMuons()
    # offMuonPhi = triggerOfflineMuons.phi
    # offMuonEta = triggerOfflineMuons.eta
    # hltMuonPhi = hltMuons.phi
    # hltMuonEta = hltMuons.eta

    # anum = ak.num(offMuonEta,axis=1)
    # bnum = ak.num(hltMuonEta,axis=1)
    # maxNum = ak.max(ak.concatenate([anum,bnum]))
    # numEvent = len(hltMuonEta)
    # if maxNum == 0:
    #     nHLTMatchedMuons = np.zeros(numEvent)
    # else:
    #     offMuonPhi = ak.fill_none(ak.pad_none(offMuonPhi,maxNum),np.Inf)
    #     offMuonEta = ak.fill_none(ak.pad_none(offMuonEta,maxNum),np.Inf)
    #     hltMuonPhi = ak.fill_none(ak.pad_none(hltMuonPhi,maxNum),np.Inf)
    #     hltMuonEta = ak.fill_none(ak.pad_none(hltMuonEta,maxNum),np.Inf)

    #     matchedList = []
    #     # for the single muon dataset, maxNum = 1, so this loop is not going to take that long to run
    #     for i in range(maxNum):
    #         matched = np.zeros(numEvent,dtype=bool)
    #         for j in range(maxNum):
    #             deltaRMatch = delta_R(hltMuonEta[:,i],offMuonEta[:,j],hltMuonPhi[:,i],offMuonPhi[:,j])
    #             cond = (deltaRMatch<0.2) & np.isfinite(deltaRMatch)
    #             matched = matched | cond
    #         matchedList.append(matched)
    #     matchedList = np.transpose(matchedList)
    #     nHLTMatchedMuons = np.sum(matchedList,axis=1)
    # varVal['nOffMuons'] = ak.num(triggerOfflineMuons)
    # varVal['nHLTMatchedMuons'] = nHLTMatchedMuons
    ############################################

    varVal['djets'] = djets
    varVal['fjets'] = fjets
    varVal['jets'] = jets
    varVal['bjets'] = bjets
    varVal['electrons'] = electrons
    varVal['muons'] = muons
    varVal['nonIsoMuons'] = nonIsoMuons
    varVal['eCounter'] = eCounter
    varVal['evtw'] = evtw
    varVal['nl'] = (ak.num(electrons) + ak.num(muons))
    varVal['nnim'] = ak.num(nonIsoMuons)
    varVal['njets'] = ak.num(jets)
    varVal['njetsAK8'] = ak.num(fjets)
    varVal['nb'] = nBJets
    varVal['met'] = met
    varVal['metPhi'] = metPhi
    ##### TTStitch ######
    # varVal['madHT'] = events.madHT
    #####################
    varVal['mT'] = mtAK8
    varVal['ht'] = ht
    varVal['st'] = st
    varVal['METrHT_pt30'] = met/ht
    varVal['METrST_pt30'] = met/st
    varVal['dPhiMinjMET'] = dPhiMinj
    varVal['dPhiMinjMETAK8'] = dPhiMinjAK8
    varVal['J1AK8Pt'] = jetVar_i(fjets.pt,0)
    if isData == False:
        varVal['GenJetsAK8'] = events.GenJetsAK8
        varVal['GenParticles'] = events.GenParticles
    # cr variables to be stored
    varVal['crElectrons']   = crElectrons
    varVal['crMuons']       = crMuons

    return varVal

def jConstVarGetter(dataset,events,varVal,cut):
    evtw = varVal["evtw"][cut]
    fjets = varVal["fjets"][cut]
    jetCats = []
    fjw = awkwardReshape(fjets,evtw)
    evtNum = events.EvtNum
    fjEvtNum = awkwardReshape(fjets,evtNum)
    bkgKeys = ["QCD","TTJets","WJets","ZJets"]
    isSignal = 0
    if "mMed" in dataset:
        if "s-channel" in dataset:
            isSignal = 1
        else:
            isSignal = 2
    # getting jet category for each jet
    if isSignal == 1:
        jetCats = ak.where(fjets.isHV,9,0) # treating SVJ from s-channel as QM jets and non-SVJ as SM jets

    elif isSignal == 2:
        GenJetsAK8 = events.GenJetsAK8
        jetsAK8GenInd = fjets.genIndex
        fjets = fjets[jetsAK8GenInd != -1]
        jetsAK8GenInd = jetsAK8GenInd[jetsAK8GenInd != -1]
        genHVCategory = GenJetsAK8.hvCategory
        jetCats = genHVCategory[jetsAK8GenInd]
    else:
        jetCats = awkwardReshape(fjets,np.ones(len(events))*-1)

    jetConstituents = events.JetsConstituents
    JetsAK8_constituentsIndex = fjets.constituentsIndex
    jCstPt = jetConstituents.pt
    indices = ak.flatten(JetsAK8_constituentsIndex,axis=-1)
    jetConstituentsForJets = jetConstituents[indices]
    goodJetConst = ak.unflatten(jetConstituentsForJets,ak.flatten(ak.num(JetsAK8_constituentsIndex,axis=-1)),axis=1)

    jCst4vec = {}
    jCstVar = {}
    jCst4vec["jCstPt"] = goodJetConst.pt
    jCst4vec["jCstEta"] = goodJetConst.eta
    jCst4vec["jCstPhi"] = goodJetConst.phi
    jCst4vec["jCstEnergy"] = goodJetConst.energy
    jCst4vec["jCstPdgId"] = goodJetConst.PdgId
    jCst4vec["jCstdxy"] = goodJetConst.dxy
    jCst4vec["jCstdxysig"] = goodJetConst.dxysig
    jCst4vec["jCstdz"] = goodJetConst.dz
    jCst4vec["jCstdzsig"] = goodJetConst.dzsig
    jCst4vec["jCstPuppiWeight"] = goodJetConst.PuppiWeight

    fjets["jetCats"] = jetCats
    fjets["fjw"] = fjw
    fjets["fjEvtNum"] = fjEvtNum
    fjets["fJNum"] = ak.local_index(fjw)

    #fjets_const = ak.broadcast_arrays(fjets,goodJetConst)[0]
    goodJetConstPt = goodJetConst.pt
    jCstVar["jCstPtAK8"] = ak.broadcast_arrays(fjets.pt,goodJetConstPt)[0]
    jCstVar["jCstEtaAK8"] = ak.broadcast_arrays(fjets.eta,goodJetConstPt)[0]
    jCstVar["jCstPhiAK8"] = ak.broadcast_arrays(fjets.phi,goodJetConstPt)[0]
    jCstVar["jCstEnergyAK8"] = ak.broadcast_arrays(fjets.energy,goodJetConstPt)[0]
    jCstVar["jCstAxismajorAK8"] = ak.broadcast_arrays(fjets.axismajor,goodJetConstPt)[0]
    jCstVar["jCstAxisminorAK8"] = ak.broadcast_arrays(fjets.axisminor,goodJetConstPt)[0]
    jCstVar["jCstTau1AK8"] = ak.broadcast_arrays(fjets.NsubjettinessTau1,goodJetConstPt)[0]
    jCstVar["jCstTau2AK8"] = ak.broadcast_arrays(fjets.NsubjettinessTau2,goodJetConstPt)[0]
    jCstVar["jCstTau3AK8"] = ak.broadcast_arrays(fjets.NsubjettinessTau3,goodJetConstPt)[0]
    jCstVar["jCstPtDAK8"] = ak.broadcast_arrays(fjets.ptD,goodJetConstPt)[0]
    jCstVar["jCstSoftDropMassAK8"] = ak.broadcast_arrays(fjets.softDropMass,goodJetConstPt)[0]
    jCstVar["jCsthvCategory"] = ak.broadcast_arrays(fjets.jetCats,goodJetConstPt)[0]
    jCstVar["jCstWeightAK8"] = ak.broadcast_arrays(fjets.fjw,goodJetConstPt)[0]
    jCstVar["jCstEvtNum"] = ak.broadcast_arrays(fjets.fjEvtNum,goodJetConstPt)[0]
    jCstVar["jCstJNum"] = ak.broadcast_arrays(fjets.fJNum,goodJetConstPt)[0]

    return jCst4vec,jCstVar

def varGetter(dataset,events,varVal,cut,jNVar=False):
    jets = varVal['jets'][cut] 
    bjets = varVal['bjets'][cut] 
    fjets = varVal['fjets'][cut] 
    djets = varVal['djets'][cut] 
    electrons = varVal['electrons'][cut] 
    muons = varVal['muons'][cut] 
    nonIsoMuons = varVal['nonIsoMuons'][cut] 
    evtw = varVal['evtw'][cut] 
    eCounter = varVal['eCounter'][cut] 
    nBJets = varVal['nb'][cut]
    ##### TTStitch ######
    # madHT= varVal['madHT'][cut]
    #####################
    met = varVal['met'][cut]
    metPhi = varVal['metPhi'][cut]
    mtAK8 = varVal['mT'][cut]
    ht = varVal['ht'][cut]
    st = varVal['st'][cut]
    dPhiMinj = varVal['dPhiMinjMET'][cut]
    dPhiMinjAK8 = varVal['dPhiMinjMETAK8'][cut]
    # cr leptons 
    crElectrons = varVal['crElectrons'][cut]
    crMuons     = varVal['crMuons'][cut]

    jetAK8Eta = fjets.eta
    jetAK8Phi = fjets.phi
    j1_etaAK8 = jetVar_i(jetAK8Eta,0)
    j2_etaAK8 = jetVar_i(jetAK8Eta,1)
    j1_phiAK8 = jetVar_i(jetAK8Phi,0)
    j2_phiAK8 = jetVar_i(jetAK8Phi,1)

    ## GenJetsAK8_hvCategory is only present in the signal samples, not the background
    jetCats = []
    jetDarkPtFracs = []
    bkgKeys = ["QCD","TTJets","WJets","ZJets"]
    isSignal = False
    if "mMed" in dataset:
        isSignal = True        

    if isSignal:
        ## Calculating the number of N-med events
        GenParticlesPdgId = abs(varVal['GenParticles'][cut].PdgId)
        hvCond = ak.zeros_like(GenParticlesPdgId,dtype=bool)
        medIDs = [4900001,4900002,4900003,4900004,4900005,4900006]
        for medID in medIDs:
            hvCond = hvCond | (GenParticlesPdgId == medID)
        num_of_med = ak.sum(hvCond,axis=1)
        GenJetsAK8 = events.GenJetsAK8
        jetsAK8GenInd = fjets.genIndex
        fjets = fjets[jetsAK8GenInd != -1]
        genHVCategory = GenJetsAK8.hvCategory
        genDarkPtFrac = GenJetsAK8.darkPtFrac
        jetsAK8GenInd = jetsAK8GenInd[jetsAK8GenInd != -1]
        jetCats = genHVCategory[jetsAK8GenInd]
        jetDarkPtFracs = genDarkPtFrac[jetsAK8GenInd]
    else:
        num_of_med = np.zeros(len(events)) 
        jetCats = awkwardReshape(fjets,np.ones(len(evtw))*-1)
        jetDarkPtFracs = awkwardReshape(fjets,np.ones(len(evtw))*-1)
        # GenJetsAK8 = np.ones(len(events))

    varVal['JetsAK8_hvCategory'] = jetCats
    varVal['JetsAK8_darkPtFrac'] = jetDarkPtFracs
    # varVal['GenJetsAK8'] = ak.num(GenJetsAK8)

    for i in range(4):
        varVal['J{}_hvCategory'.format(i+1)] = jetVar_i(varVal['JetsAK8_hvCategory'],i)
        varVal['J{}_darkPtFrac'.format(i+1)] = jetVar_i(varVal['JetsAK8_darkPtFrac'],i)
    ew = awkwardReshape(electrons,evtw)
    mw = awkwardReshape(muons,evtw)
    nimw = awkwardReshape(nonIsoMuons,evtw)
    jw = awkwardReshape(jets,evtw)
    fjw = awkwardReshape(fjets,evtw)
    crew = awkwardReshape(crElectrons,evtw)
    crmw = awkwardReshape(crMuons,evtw)

    # AK4 Jet Variables
    jetPhi = jets.phi
    jetEta = jets.eta
    j1_eta = jetVar_i(jetEta,0)
    j2_eta = jetVar_i(jetEta,1)
    j1_phi = jetVar_i(jetPhi,0)
    j2_phi = jetVar_i(jetPhi,1)
    dPhij1 = deltaPhi(j1_phi,metPhi)
    dPhij2 = deltaPhi(j2_phi,metPhi)
    dPhij1rdPhij2 = dPhij1/dPhij2
    dPhiMinj = ak.min(deltaPhi(jetPhi,metPhi),axis=1,mask_identity=False)
    dEtaj12 = deltaEta(j1_eta,j2_eta)
    deltaR12j = delta_R(j1_eta,j2_eta,j1_phi,j2_phi)

    # AK8 Jet Variables
    jetAK8pT = fjets.pt
    jetAK8Phi = fjets.phi
    jetAK8Eta = fjets.eta
    jetAK8M = fjets.mass
    j1_etaAK8 = jetVar_i(jetAK8Eta,0)
    j2_etaAK8 = jetVar_i(jetAK8Eta,1)
    j1_phiAK8 = jetVar_i(jetAK8Phi,0)
    j2_phiAK8 = jetVar_i(jetAK8Phi,1)
    dPhij1AK8 = deltaPhi(j1_phiAK8,metPhi)
    dPhij2AK8 = deltaPhi(j2_phiAK8,metPhi)
    dPhij1rdPhij2AK8 = dPhij1AK8/dPhij2AK8
    dPhijAK8 = deltaPhi(jetAK8Phi,metPhi)
    dPhiMinjAK8 = ak.min(dPhijAK8,axis=1,mask_identity=False)
    dEtaj12AK8 = deltaEta(j1_etaAK8,j2_etaAK8)
    deltaR12jAK8 = delta_R(j1_etaAK8,j2_etaAK8,j1_phiAK8,j2_phiAK8)
    tau1 = fjets.NsubjettinessTau1
    tau2 = fjets.NsubjettinessTau2
    tau3 = fjets.NsubjettinessTau3
    J_tau21 = tau2/tau1
    J_tau32 = tau3/tau2
    J1_tau21 = tauRatio(tau2,tau1,0)
    J1_tau32 = tauRatio(tau3,tau2,0)
    J2_tau21 = tauRatio(tau2,tau1,1)
    J2_tau32 = tauRatio(tau3,tau2,1)

    varVal['eCounter'] = eCounter
    varVal['evtw'] = evtw
    varVal['jw'] = jw
    varVal['fjw'] = fjw
    varVal['ew'] = ew
    varVal['mw'] = mw
    varVal['nimw'] = nimw
    varVal['njets'] = ak.num(jets)
    varVal['njetsAK8'] = ak.num(fjets)
    varVal['nTruthSVJ'] = ak.num(djets)
    varVal['nb'] = ak.num(bjets)
    varVal['nl'] = varVal['nl'][cut]
    varVal['nnim'] = varVal['nnim'][cut]
    varVal['met'] = met
    varVal['logMET'] = np.log(met)
    varVal['metPhi'] = metPhi
    varVal['metSig'] = events.METSignificance
    varVal['logMETSig'] = np.log(varVal['metSig'])
    varVal['mT'] = mtAK8
    varVal['ht'] = ht
    varVal['st'] = st
    varVal['METrHT_pt30'] = varVal['METrHT_pt30'][cut]
    varVal['METrST_pt30'] = varVal['METrST_pt30'][cut]
    varVal['dPhiMinjMET'] = dPhiMinj
    varVal['dPhiMinjMETAK8'] = dPhiMinjAK8
    ##### TTStitch ######
    # varVal['madHT'] = madHT
    #####################
    varVal['jPt'] = jets.pt
    varVal['jEta'] = jetEta
    varVal['jPhi'] = jetPhi
    varVal['jE'] = jets.energy
    varVal['jAxismajor'] = jets.axismajor
    varVal['jAxisminor'] = jets.axisminor
    varVal['jPtD'] = jets.ptD
    varVal['dPhiMinjMET'] = dPhiMinj
    varVal['jPtAK8'] = fjets.pt
    varVal['jEtaAK8'] = jetAK8Eta
    varVal['jPhiAK8'] = jetAK8Phi
    varVal['jEAK8'] = fjets.energy
    varVal['jAxismajorAK8'] = fjets.axismajor
    varVal['jAxisminorAK8'] = fjets.axisminor
    varVal['jChEMEFractAK8'] = fjets.chargedEmEnergyFraction
    varVal['jChHadEFractAK8'] = fjets.chargedHadronEnergyFraction
    varVal['jChHadMultAK8'] = fjets.chargedHadronMultiplicity
    varVal['jChMultAK8'] = fjets.chargedMultiplicity
    varVal['jecfN2b1AK8'] = fjets.ecfN2b1
    varVal['jecfN2b2AK8'] = fjets.ecfN2b2
    varVal['jecfN3b1AK8'] = fjets.ecfN3b1
    varVal['jecfN3b2AK8'] = fjets.ecfN3b2
    varVal['jEleEFractAK8'] = fjets.electronEnergyFraction
    varVal['jEleMultAK8'] = fjets.electronMultiplicity
    varVal['jGirthAK8'] = fjets.girth
    varVal['jHfEMEFractAK8'] = fjets.hfEMEnergyFraction
    varVal['jHfHadEFractAK8'] = fjets.hfHadronEnergyFraction
    varVal['jMultAK8'] = fjets.multiplicity
    varVal['jMuEFractAK8'] = fjets.muonEnergyFraction
    varVal['jMuMultAK8'] = fjets.muonMultiplicity
    varVal['jNeuEmEFractAK8'] = fjets.neutralEmEnergyFraction
    varVal['jNeuHadEFractAK8'] = fjets.neutralHadronEnergyFraction
    varVal['jNeuHadMultAK8'] = fjets.neutralHadronMultiplicity
    varVal['jNeuMultAK8'] = fjets.neutralMultiplicity
    varVal['jTau1AK8'] = tau1
    varVal['jTau2AK8'] = tau2
    varVal['jTau3AK8'] = tau3
    varVal['jTau21AK8'] = J_tau21
    varVal['jTau32AK8'] = J_tau32
    varVal['jPhoEFractAK8'] = fjets.photonEnergyFraction
    varVal['jPhoMultAK8'] = fjets.photonMultiplicity
    varVal['jPtDAK8'] = fjets.ptD
    varVal['jSoftDropMassAK8'] = fjets.softDropMass
    varVal['dPhijMETAK8'] = dPhijAK8
    varVal['dPhiMinjMETAK8'] = dPhiMinjAK8
    varVal['dEtaj12AK8'] = dEtaj12AK8
    varVal['dRJ12AK8'] = deltaR12jAK8
    varVal['dPhij1rdPhij2AK8'] = dPhij1rdPhij2AK8
    varVal['electronsIso'] = electrons.iso


    # Variable for the Control region cut
    varVal['electronPT'] = electrons.pt
    varVal['electronPhi'] = electrons.phi
    varVal['electronEta'] = electrons.eta
    varVal['muonPT'] = muons.pt
    varVal['muonPhi'] = muons.phi
    varVal['muonEta'] = muons.eta

    varVal['crew']          = crew
    varVal['crmw']          = crmw
    varVal['crElectronPT']  = crElectrons.pt
    varVal['crElectronPhi'] = crElectrons.phi
    varVal['crElectronEta'] = crElectrons.eta
    varVal['crMuonPT']     = crMuons.pt
    varVal['crMuonPhi']     = crMuons.phi
    varVal['crMuonEta']     = crMuons.eta
    varVal['dPhiMinJAK8crElectron1'] = ak.min(deltaPhi(jetAK8Phi,jetVar_i(crElectrons.phi,0)),axis=1,mask_identity=False)
    varVal['dPhiMinJAK8crMuon1']     = ak.min(deltaPhi(jetAK8Phi,jetVar_i(crMuons.phi,0)),axis=1,mask_identity=False)
    mtMETCRMuon = mtMETLepton(met,metPhi,crMuons)
    mtMETCRElectron = mtMETLepton(met,metPhi,crElectrons)
    # print(mtMETCRMuon)
    # print(mtMETCRElectron)
    varVal['mtMETCRMuon'] = mtMETCRMuon
    varVal['mtMETCRElectron'] = mtMETCRElectron




    #varVal['muonsIso'] = muons.iso
    #varVal['nonIsoMuonsPt'] = nonIsoMuons.pt
    #varVal['nonIsoMuonsIso'] = nonIsoMuons.iso
    ## Save first AK8 jet info 
    


    if jNVar:
        # preparing histograms for jN variables
        # print("jetVar_i(fjets.pt,0) ",jetVar_i(fjets.pt,1))
        # maxN = len(fjets) if len(fjets) < 3 else 3
        maxN = 4
        # print("len(fjets) = ",len(fjets))
        for i in range(maxN):
            # i_th_mask = fjets
            varVal['j{}Pt'.format(i+1)] = jetVar_i(jets.pt,i)
            varVal['j{}Eta'.format(i+1)] = jetVar_i(jetEta,i)
            varVal['j{}Phi'.format(i+1)] = jetVar_i(jetPhi,i)
            varVal['j{}E'.format(i+1)] = jetVar_i(jets.energy,i)
            varVal['j{}Axismajor'.format(i+1)] = jetVar_i(jets.axismajor,i)
            varVal['j{}Axisminor'.format(i+1)] = jetVar_i(jets.axisminor,i)
            varVal['j{}PtD'.format(i+1)] = jetVar_i(jets.ptD,i)
            varVal['dPhij{}MET'.format(i+1)] = deltaPhi(jetVar_i(jetPhi,i),metPhi)
            varVal['j{}PtAK8'.format(i+1)] = jetVar_i(fjets.pt,i)
            varVal['j{}EtaAK8'.format(i+1)] = jetVar_i(jetAK8Eta,i)
            varVal['j{}PhiAK8'.format(i+1)] = jetVar_i(jetAK8Phi,i)
            varVal['j{}EAK8'.format(i+1)] = jetVar_i(fjets.energy,i)
            varVal['j{}AxismajorAK8'.format(i+1)] = jetVar_i(fjets.axismajor,i)
            varVal['j{}AxisminorAK8'.format(i+1)] = jetVar_i(fjets.axisminor,i)
            varVal['j{}GirthAK8'.format(i+1)] = jetVar_i(fjets.girth,i)
            varVal['j{}PtDAK8'.format(i+1)] = jetVar_i(fjets.ptD,i)
            varVal['j{}Tau1AK8'.format(i+1)] = jetVar_i(tau1,i)
            varVal['j{}Tau2AK8'.format(i+1)] = jetVar_i(tau2,i)
            varVal['j{}Tau3AK8'.format(i+1)] = jetVar_i(tau3,i)
            varVal['j{}Tau21AK8'.format(i+1)] = tauRatio(tau2,tau1,i)
            varVal['j{}Tau32AK8'.format(i+1)] = tauRatio(tau3,tau2,i)
            varVal['j{}SoftDropMassAK8'.format(i+1)] = jetVar_i(fjets.softDropMass,i)
            varVal['dPhij{}METAK8'.format(i+1)] = deltaPhi(jetVar_i(jetAK8Phi,i),metPhi)
            varVal['dRj{}AK8crMuon1'.format(i+1)]= delta_R(jetVar_i(jetAK8Eta,i),jetVar_i(crMuons.eta,0),jetVar_i(jetAK8Phi,i),jetVar_i(crMuons.phi,0))
            varVal['dRj{}AK8crElectron1'.format(i+1)]= delta_R(jetVar_i(jetAK8Eta,i),jetVar_i(crElectrons.eta,0),jetVar_i(jetAK8Phi,i),jetVar_i(crElectrons.phi,0))
        
        allComs = list(combinations(range(maxN),2))
        for com in allComs:
            j1 = com[0]
            j2 = com[1]
            j1_eta = jetVar_i(jetEta,j1)
            j2_eta = jetVar_i(jetEta,j2)
            j1_phi = jetVar_i(jetPhi,j1)
            j2_phi = jetVar_i(jetPhi,j2)
            dPhij1 = deltaPhi(j1_phi,metPhi)
            dPhij2 = deltaPhi(j2_phi,metPhi)
            j1_etaAK8 = jetVar_i(jetAK8Eta,j1)
            j2_etaAK8 = jetVar_i(jetAK8Eta,j2)
            j1_phiAK8 = jetVar_i(jetAK8Phi,j1)
            j2_phiAK8 = jetVar_i(jetAK8Phi,j2)
            dPhij1AK8 = deltaPhi(j1_phiAK8,metPhi)
            dPhij2AK8 = deltaPhi(j2_phiAK8,metPhi)
            varVal['dEtaj{}{}'.format(j1+1,j2+1)] = deltaEta(j1_eta,j2_eta)
            varVal['dPhij{}{}'.format(j1+1,j2+1)] = deltaPhi(j1_phi,j2_phi)
            varVal['dRj{}{}'.format(j1+1,j2+1)] = delta_R(j1_eta,j2_eta,j1_phi,j2_phi)
            varVal['dPhij{}rdPhij{}'.format(j1+1,j2+1)] = dPhij1/dPhij2
            varVal['dEtaj{}{}AK8'.format(j1+1,j2+1)] = deltaEta(j1_etaAK8,j2_etaAK8)
            varVal['dPhij{}{}AK8'.format(j1+1,j2+1)] = deltaPhi(j1_phiAK8,j2_phiAK8)
            varVal['dRj{}{}AK8'.format(j1+1,j2+1)] = delta_R(j1_etaAK8,j2_etaAK8,j1_phiAK8,j2_phiAK8)
            varVal['dPhij{}rdPhij{}AK8'.format(j1+1,j2+1)] = dPhij1AK8/dPhij2AK8
    varVal['nNMedEvent'] = np.array(num_of_med)
    # varVal['mT2_f4_msm'] = f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"")
    # varVal['mT2_f4_msm_dEta'] = f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dEta")
    # varVal['mT2_f4_msm_dPhi'] = f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dPhi")
    # varVal['mT2_f4_msm_dR'] = f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dR")
    # varVal['GenMT2_AK8'] = GenMT2_AK8
