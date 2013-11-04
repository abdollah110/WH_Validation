/* 
 * File:   zh_Tree.h
 * Author: abdollah
 *
 * Created on February 6, 2013, 1:56 PM
 */

#ifndef ZH_TREE_H
#define	ZH_TREE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "myevent.h"
#include "myobject.h"
#include "Leptons_IdIso.h"
#include "zh_Auxiliary.h"
#include "Leptons_PreSelection.h"
#include "zh_Functions.h"





unsigned int Channel = 0;
unsigned int Run = 0;
unsigned int Lumi = 0;
unsigned int Event = 0;
float IMass = 0;
float ZMass = 0;
float HMass = 0;
float met, metPhi, pfmet, pfmetPhi;
float l1M, l1Px, l1Py, l1Pz, l1E, l1Pt, l1Phi, l1Eta, l1Charge, l1_muIso, l1_eleIso, l1_eleMVANonTrg, l1_eleNumHit = -10;
float l2M, l2Px, l2Py, l2Pz, l2E, l2Pt, l2Phi, l2Eta, l2Charge, l2_muIso, l2_eleIso, l2_eleMVANonTrg, l2_eleNumHit = -10;
float l3M, l3Px, l3Py, l3Pz, l3E, l3Pt, l3Phi, l3Eta, l3Charge, l3_muIso, l3_eleIso, l3_eleMVANonTrg, l3_eleNumHit, l3_tauIsoMVA2raw = -10;
float l3_RefJetPt, l3_RefJetEta, l3_RefJetPhi = -10;
bool l1_muId, l1_eleId, l2_muId, l2_eleId, l3_muId_Loose, l3_muId_Tight, l3_eleId_Loose, l3_eleId_Tight;
bool l3_tauIsoL, l3_tauIsoM, l3_tauIsoT, l3_tauRejMuL, l3_tauRejMuM, l3_tauRejMuT, l3_tauRejEleL, l3_tauRejEleM, l3_tauRejEleMVA;
bool l3_tauIso3HitL, l3_tauIso3HitM, l3_tauIso3HitT, l3_tauRejMu2L, l3_tauRejMu2M, l3_tauRejMu2T, l3_tauRejEleMVA3L, l3_tauRejEleMVA3M, l3_tauRejEleMVA3T;
bool l3_tauIsoVL, l3_tauIsoMVA2L, l3_tauIsoMVA2M, l3_tauIsoMVA2T;
bool l1_muTrgObjMatch, l2_muTrgObjMatch, l1_eleTrgObjMatch, l2_eleTrgObjMatch;

float covMet11, covMet12, covMet21, covMet22;
float pfcovMet11, pfcovMet12, pfcovMet21, pfcovMet22;
float eff_Correction, pu_Weight;
int num_PV, num_bjet, num_goodjet;
int mu_Size, BareMuon_Size, electron_Size, BareElectron_Size, tau_Size, BareTau_Size;
int mu_partTight_Size, ele_partTight_Size;
float l3_CloseJetPt;
float l3_CloseJetEta;
float l3_CloseJetPhi;
float l1_CloseJetPt;
float l1_CloseJetEta;
float l1_CloseJetPhi;
float l2_CloseJetPt;
float l2_CloseJetEta;
float l2_CloseJetPhi;
bool l1_passConversionVeto, l2_passConversionVeto;
bool l1_isGsfCtfScPixChargeConsistent;
bool l1_isGsfScPixChargeConsistent;
bool l1_isGsfCtfChargeConsistent;
bool l2_isGsfCtfScPixChargeConsistent;
bool l2_isGsfScPixChargeConsistent;
bool l2_isGsfCtfChargeConsistent;

void fillTree(TTree * Run_Tree, myevent *m, float cor_eff, float PU_Weight, int channel, myobject obj1, myobject obj2, myobject obj3) {



    vector<myobject> Met = m->RecMVAMet;
    vector<myobject> PFMet = m->RecPFMet;
    Channel = channel;
    Run = m->runNumber;
    Lumi = m->lumiNumber;
    Event = m->eventNumber;
    IMass = InvarMass_4(obj1, obj2, obj3, obj3);
    ZMass = InvarMass_2(obj1, obj2);
    HMass = InvarMass_2(obj3, obj3);
    covMet11 = m->MVAMet_sigMatrix_00;
    covMet12 = m->MVAMet_sigMatrix_01;
    covMet21 = m->MVAMet_sigMatrix_10;
    covMet22 = m->MVAMet_sigMatrix_11;
    pfcovMet11 = m->MET_sigMatrix_00;
    pfcovMet12 = m->MET_sigMatrix_01;
    pfcovMet21 = m->MET_sigMatrix_10;
    pfcovMet22 = m->MET_sigMatrix_11;
    met = Met.front().pt;
    pfmet = PFMet.front().pt;
    metPhi = Met.front().phi;
    pfmetPhi = PFMet.front().phi;
    eff_Correction = cor_eff;
    mu_Size = myCleanLepton(m, "mu").size();
    BareMuon_Size = myCleanBareLepton(m, "mu").size();
    electron_Size = myCleanLepton(m, "ele").size();
    BareElectron_Size = myCleanBareLepton(m, "ele").size();
    tau_Size = myCleanLepton(m, "tau").size();
    BareTau_Size = myCleanBareLepton(m, "tau").size();
    vector<myobject> Vertex = m->Vertex;
    num_PV = Vertex.size();
    num_bjet = bjet_Multiplicity(m);
    pu_Weight = PU_Weight;
    num_goodjet = GoodJet(m).size();

    l1M = obj1.mass;
    l1Px = obj1.px;
    l1Py = obj1.py;
    l1Pz = obj1.pz;
    l1E = obj1.E;
    l1Pt = obj1.pt;
    l1Phi = obj1.phi;
    l1Eta = obj1.eta;
    l1_muId = Id_Mu_Loose(obj1);
    l1_muIso = Iso_Mu_dBeta(obj1);
    //    l1_eleId = EleMVANonTrigId_Loose(obj1);
    l1_eleId = EleMVANonTrigId_Tight(obj1);
    l1_eleIso = Iso_Ele_dBeta(obj1);
    l1_eleMVANonTrg = obj1.Id_mvaNonTrg;
    l1_eleNumHit = obj1.numHitEleInner;
    l1Charge = obj1.charge;
    l1_muTrgObjMatch = obj1.hasTrgObject_loose;
    l1_eleTrgObjMatch = obj1.hasTrgObject_loose;
    l1_passConversionVeto = obj1.passConversionVeto;
    l1_isGsfCtfScPixChargeConsistent = obj1.isGsfCtfScPixChargeConsistent;
    l1_isGsfScPixChargeConsistent = obj1.isGsfScPixChargeConsistent;
    l1_isGsfCtfChargeConsistent = obj1.isGsfCtfChargeConsistent;
    l1_CloseJetPt = Find_Closet_Jet(obj1, m)[0];
    l1_CloseJetEta = Find_Closet_Jet(obj1, m)[1];
    l1_CloseJetPhi = Find_Closet_Jet(obj1, m)[2];

    l2M = obj2.mass;
    l2Px = obj2.px;
    l2Py = obj2.py;
    l2E = obj2.E;
    l2Pt = obj2.pt;
    l2Phi = obj2.phi;
    l2Pz = obj2.pz;
    l2Eta = obj2.eta;
    l2_muId = Id_Mu_Loose(obj2);
    l2_muIso = Iso_Mu_dBeta(obj2);
    //    l2_eleId = EleMVANonTrigId_Loose(obj2);
    l2_eleId = EleMVANonTrigId_Tight(obj2);
    l2_eleIso = Iso_Ele_dBeta(obj2);
    l2_eleMVANonTrg = obj2.Id_mvaNonTrg;
    l2_eleNumHit = obj2.numHitEleInner;
    l2Charge = obj2.charge;
    l2_muTrgObjMatch = obj2.hasTrgObject_loose;
    l2_eleTrgObjMatch = obj2.hasTrgObject_loose;
    l2_passConversionVeto = obj2.passConversionVeto;
    l2_isGsfCtfScPixChargeConsistent = obj2.isGsfCtfScPixChargeConsistent;
    l2_isGsfScPixChargeConsistent = obj2.isGsfScPixChargeConsistent;
    l2_isGsfCtfChargeConsistent = obj2.isGsfCtfChargeConsistent;
    l2_CloseJetPt = Find_Closet_Jet(obj2, m)[0];
    l2_CloseJetEta = Find_Closet_Jet(obj2, m)[1];
    l2_CloseJetPhi = Find_Closet_Jet(obj2, m)[2];

    l3M = obj3.mass;
    l3Px = obj3.px;
    l3Py = obj3.py;
    l3Pz = obj3.pz;
    l3E = obj3.E;
    l3Pt = obj3.pt;
    l3Phi = obj3.phi;
    l3Eta = obj3.eta;
    l3_muId_Loose = Id_Mu_Loose(obj3);
    l3_muId_Tight = Id_Mu_Tight(obj3);
    l3_muIso = Iso_Mu_dBeta(obj3);
    l3_eleId_Loose = EleMVANonTrigId_Loose(obj3);
    l3_eleId_Tight = EleMVANonTrigId_Tight(obj3);
    l3_eleIso = Iso_Ele_dBeta(obj3);
    l3_eleMVANonTrg = obj3.Id_mvaNonTrg;
    l3_eleNumHit = obj3.numHitEleInner;
    l3_tauIsoMVA2raw = obj3.byIsolationMVA2raw;
    l3_tauIsoVL = obj3.byVLooseCombinedIsolationDeltaBetaCorr;
    l3_tauIsoL = obj3.byLooseCombinedIsolationDeltaBetaCorr;
    l3_tauIso3HitL = obj3.byLooseCombinedIsolationDeltaBetaCorr3Hits;
    l3_tauIsoMVA2L = obj3.byLooseIsolationMVA2;
    l3_tauIsoM = obj3.byMediumCombinedIsolationDeltaBetaCorr;
    l3_tauIso3HitM = obj3.byMediumCombinedIsolationDeltaBetaCorr3Hits;
    l3_tauIsoMVA2M = obj3.byMediumIsolationMVA2;
    l3_tauIsoT = obj3.byTightCombinedIsolationDeltaBetaCorr;
    l3_tauIso3HitT = obj3.byTightCombinedIsolationDeltaBetaCorr3Hits;
    l3_tauIsoMVA2T = obj3.byTightIsolationMVA2;
    l3_tauRejMuL = obj3.discriminationByMuonLoose;
    l3_tauRejMu2L = obj3.discriminationByMuonLoose2;
    l3_tauRejMuM = obj3.discriminationByMuonMedium;
    l3_tauRejMu2M = obj3.discriminationByMuonMedium2;
    l3_tauRejMuT = obj3.discriminationByMuonTight;
    l3_tauRejMu2T = obj3.discriminationByMuonTight2;
    l3_tauRejEleL = obj3.discriminationByElectronLoose;
    l3_tauRejEleM = obj3.discriminationByElectronMedium;
    l3_tauRejEleMVA = obj3.discriminationByElectronMVA;
    l3_tauRejEleMVA3L = obj3.discriminationByElectronMVA3Loose;
    l3_tauRejEleMVA3M = obj3.discriminationByElectronMVA3Medium;
    l3_tauRejEleMVA3T = obj3.discriminationByElectronMVA3Tight;
    l3Charge = obj3.charge;
    l3_CloseJetPt = Find_Closet_Jet(obj3, m)[0];
    l3_CloseJetEta = Find_Closet_Jet(obj3, m)[1];
    l3_CloseJetPhi = Find_Closet_Jet(obj3, m)[2];
    l3_RefJetPt = obj3.jetPt;
    l3_RefJetEta = obj3.jetEta;
    l3_RefJetPhi = obj3.jetPhi;




    Run_Tree->Fill();
}



#endif	/* ZH_TREE_H */

