// The code to do teh ZH totautau Analysis
// to make it excutable run: ./Make.sh ZH_Analysis.cc

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
//needed to make the executable
#include "interface/myevent.h"
#include "interface/LinkDef.h"
#include "interface/myobject.h"
// needed header files
#include "interface/makeHisto.h"
#include "interface/Leptons_PreSelection.h"
#include "interface/zh_Auxiliary.h"
#include "interface/zh_Corrector.h"
#include "interface/zh_Trigger.h"
#include "interface/zh_Tree.h"
#include "interface/Leptons_IdIso.h"
#include "interface/zh_Functions.h"
#include "interface/LumiReweightingStandAlone.h"

int main(int argc, char** argv) {

    using namespace std;
    //define myevent class
    myevent *m = new myevent;
    //define 1D and 2D histogram
    myMap1 = new map<string, TH1F*>();
    myMap2 = new map<string, TH2F*>();

    cout << "\n######################### Analysis is initializing ####################################### " << endl;

    //#################################################################################################
    //################# First Argument, Data or MC, which type of data or MC    #######################
    //#################################################################################################

    string status_sample = *(argv + 1);
    cout << "*** First Argument, Data or MC, which type of data or MC ***" << endl;
    cout << status_sample.c_str() << endl;

    bool mc12 = (status_sample.compare("mc12") == 0 ? true : false);
    bool mc11 = (status_sample.compare("mc11") == 0 ? true : false);
    bool data12 = (status_sample.compare("data12") == 0 ? true : false);
    bool data11 = (status_sample.compare("data11") == 0 ? true : false);
    if (!(mc12 || mc11 || data12 || data11))
        cout << "xxxxx Error, please slecet between: mc12 || mc11 || data12 || data11 " << endl;

    //#################################################################################################
    //############## Second Argument, Run over just di-ele, just di-mu (for data) or total (for MC) ###
    //#################################################################################################

    string status_type = *(argv + 2);
    cout << "**** Second Argument, Run over just di-ele, just di-mu (for data) or total (for MC) ***" << endl;
    cout << status_type.c_str() << endl;
    bool is_tot = (status_type.compare("Tot") == 0 ? true : false);
    bool is_ele = (status_type.compare("Ele") == 0 ? true : false);
    bool is_mueg = (status_type.compare("MuEG") == 0 ? true : false);
    bool is_mu = (status_type.compare("Mu") == 0 ? true : false);
    if (!(is_tot || is_ele || is_mu || is_mueg))
        cout << "xxxxx Error, please slecet between: Tot || Ele || MuEG || Mu " << endl;

    //#################################################################################################
    //############## Third anad Forth Argument,   OutPut Name/ Input Files                         ########################
    //#################################################################################################

    string out = *(argv + 3);

    std::vector<string> fileNames;
    for (int f = 4; f < argc; f++) {
        fileNames.push_back(*(argv + f));
        // printing the input NAME
        cout << "\n INPUT NAME IS:   " << fileNames[f - 4] << "\t";
    }
    //#################################################################################################
    //############## defining an out_file name need on the given argument  information  ###############
    //#################################################################################################

    string outname = status_sample + "_" + status_type + "_" + out;
    //PRINTING THE OUTPUT name
    cout << "\n\n\n OUTPUT NAME IS:    " << outname << endl;
    TFile *fout = TFile::Open(outname.c_str(), "RECREATE");

    //#################################################################################################
    //############## initializing the PU correction                                    ###############
    //#################################################################################################

    reweight::LumiReWeighting* LumiWeights_12;
    LumiWeights_12 = new reweight::LumiReWeighting("interface/Summer12_PU.root", "interface/dataPileUpHistogram_True_2012.root", "mcPU", "pileup");
    reweight::LumiReWeighting* LumiWeights_11;
    LumiWeights_11 = new reweight::LumiReWeighting("interface/Fall11_PU_observed.root", "interface/dataPileUpHistogram_Observed_2011.root", "mcPU", "pileup");

    //#################################################################################################
    //############## defining Tree Branches Filled via fillTree function                ###############
    //#################################################################################################
    TTree *Run_Tree = new TTree("RLE_tree", "RLE_tree");
    //    To force a memory-resident Tree
    Run_Tree->SetDirectory(0);

    Run_Tree->Branch("Channel", &Channel, "Channel/I");
    Run_Tree->Branch("Run", &Run, "Run/I");
    Run_Tree->Branch("Lumi", &Lumi, "Lumi/I");
    Run_Tree->Branch("Event", &myEvent, "Event/I");
    Run_Tree->Branch("IMass", &IMass, "IMass/F");
    Run_Tree->Branch("ZMass", &ZMass, "ZMass/F");
    Run_Tree->Branch("HMass", &HMass, "HMass/F");
    Run_Tree->Branch("met", &met, "met/F");
    Run_Tree->Branch("pfmet", &pfmet, "pfmet/F");
    Run_Tree->Branch("metPhi", &metPhi, "metPhi/F");
    Run_Tree->Branch("pfmetPhi", &pfmetPhi, "pfmetPhi/F");
    Run_Tree->Branch("covMet11", &covMet11, "covMet11/F");
    Run_Tree->Branch("covMet12", &covMet12, "covMet12/F");
    Run_Tree->Branch("covMet21", &covMet21, "covMet21/F");
    Run_Tree->Branch("covMet22", &covMet22, "covMet22/F");
    Run_Tree->Branch("pfcovMet11", &pfcovMet11, "pfcovMet11/F");
    Run_Tree->Branch("pfcovMet12", &pfcovMet12, "pfcovMet12/F");
    Run_Tree->Branch("pfcovMet21", &pfcovMet21, "pfcovMet21/F");
    Run_Tree->Branch("pfcovMet22", &pfcovMet22, "pfcovMet22/F");
    Run_Tree->Branch("num_PV", &num_PV, "num_PV/I");
    Run_Tree->Branch("num_bjet", &num_bjet, "num_bjet/I");
    Run_Tree->Branch("num_goodjet", &num_goodjet, "num_goodjet/I");
    Run_Tree->Branch("eff_Correction", &eff_Correction, "eff_Correction/F");
    Run_Tree->Branch("pu_Weight", &pu_Weight, "pu_Weight/F");

    Run_Tree->Branch("mu_Size", &mu_Size, "mu_Size/I");
    Run_Tree->Branch("BareMuon_Size", &BareMuon_Size, "BareMuon_Size/I");
    Run_Tree->Branch("electron_Size", &electron_Size, "electron_Size/I");
    Run_Tree->Branch("BareElectron_Size", &BareElectron_Size, "BareElectron_Size/I");
    Run_Tree->Branch("tau_Size", &tau_Size, "tau_Size/I");
    Run_Tree->Branch("BareTau_Size", &BareTau_Size, "BareTau_Size/I");
    Run_Tree->Branch("mu_partTight_Size", &mu_partTight_Size, "mu_partTight_Size/I");
    Run_Tree->Branch("ele_partTight_Size", &ele_partTight_Size, "ele_partTight_Size/I");

    Run_Tree->Branch("l1M", &l1M, "l1M/F");
    Run_Tree->Branch("l1E", &l1E, "l1E/F");
    Run_Tree->Branch("l1Px", &l1Px, "l1Px/F");
    Run_Tree->Branch("l1Py", &l1Py, "l1Py/F");
    Run_Tree->Branch("l1Pz", &l1Pz, "l1Pz/F");
    Run_Tree->Branch("l1Pt", &l1Pt, "l1Pt/F");
    Run_Tree->Branch("l1Eta", &l1Eta, "l1Eta/F");
    Run_Tree->Branch("l1Phi", &l1Phi, "l1Phi/F");
    Run_Tree->Branch("l1Charge", &l1Charge, "l1Charge/F");
    Run_Tree->Branch("l1_muId", &l1_muId, "l1_muId/O");
    Run_Tree->Branch("l1_muIso", &l1_muIso, "l1_muIso/F");
    Run_Tree->Branch("l1_eleId", &l1_eleId, "l1_eleId/O");
    Run_Tree->Branch("l1_eleIso", &l1_eleIso, "l1_eleIso/F");
    Run_Tree->Branch("l1_eleMVANonTrg", &l1_eleMVANonTrg, "l1_eleMVANonTrg/F");
    Run_Tree->Branch("l1_eleNumHit", &l1_eleNumHit, "l1_eleNumHit/F");
    Run_Tree->Branch("l1_muTrgObjMatch", &l1_muTrgObjMatch, "l1_muTrgObjMatch/O");
    Run_Tree->Branch("l1_eleTrgObjMatch", &l1_eleTrgObjMatch, "l1_eleTrgObjMatch/O");
    Run_Tree->Branch("l1_muTrgObjMatchMed", &l1_muTrgObjMatchMed, "l1_muTrgObjMatchMed/O");
    Run_Tree->Branch("l1_eleTrgObjMatchMed", &l1_eleTrgObjMatchMed, "l1_eleTrgObjMatchMed/O");
    Run_Tree->Branch("l1_passConversionVeto", &l1_passConversionVeto, "l1_passConversionVeto/O");
    Run_Tree->Branch("l1_isGsfCtfScPixChargeConsistent", &l1_isGsfCtfScPixChargeConsistent, "l1_isGsfCtfScPixChargeConsistent/O");
    Run_Tree->Branch("l1_isGsfScPixChargeConsistent", &l1_isGsfScPixChargeConsistent, "l1_isGsfScPixChargeConsistent/O");
    Run_Tree->Branch("l1_isGsfCtfChargeConsistent", &l1_isGsfCtfChargeConsistent, "l1_isGsfCtfChargeConsistent/O");
    Run_Tree->Branch("l1_CloseJetPt", &l1_CloseJetPt, "l1_CloseJetPt/F");
    Run_Tree->Branch("l1_CloseJetEta", &l1_CloseJetEta, "l1_CloseJetEta/F");
    Run_Tree->Branch("l1_CloseJetPhi", &l1_CloseJetPhi, "l1_CloseJetPhi/F");


    Run_Tree->Branch("l2M", &l2M, "l2M/F");
    Run_Tree->Branch("l2E", &l2E, "l2E/F");
    Run_Tree->Branch("l2Px", &l2Px, "l2Px/F");
    Run_Tree->Branch("l2Py", &l2Py, "l2Py/F");
    Run_Tree->Branch("l2Pz", &l2Pz, "l2Pz/F");
    Run_Tree->Branch("l2Pt", &l2Pt, "l2Pt/F");
    Run_Tree->Branch("l2Eta", &l2Eta, "l2Eta/F");
    Run_Tree->Branch("l2Phi", &l2Phi, "l2Phi/F");
    Run_Tree->Branch("l2Charge", &l2Charge, "l2Charge/F");
    Run_Tree->Branch("l2_muId", &l2_muId, "l2_muId/O");
    Run_Tree->Branch("l2_muIso", &l2_muIso, "l2_muIso/F");
    Run_Tree->Branch("l2_eleId", &l2_eleId, "l2_eleId/O");
    Run_Tree->Branch("l2_eleIso", &l2_eleIso, "l2_eleIso/F");
    Run_Tree->Branch("l2_eleMVANonTrg", &l2_eleMVANonTrg, "l2_eleMVANonTrg/F");
    Run_Tree->Branch("l2_eleNumHit", &l2_eleNumHit, "l2_eleNumHit/F");
    Run_Tree->Branch("l2_muTrgObjMatch", &l2_muTrgObjMatch, "l2_muTrgObjMatch/O");
    Run_Tree->Branch("l2_eleTrgObjMatch", &l2_eleTrgObjMatch, "l2_eleTrgObjMatch/O");
    Run_Tree->Branch("l2_muTrgObjMatchMed", &l2_muTrgObjMatchMed, "l2_muTrgObjMatchMed/O");
    Run_Tree->Branch("l2_eleTrgObjMatchMed", &l2_eleTrgObjMatchMed, "l2_eleTrgObjMatchMed/O");
    Run_Tree->Branch("l2_passConversionVeto", &l2_passConversionVeto, "l2_passConversionVeto/O");
    Run_Tree->Branch("l2_isGsfCtfScPixChargeConsistent", &l2_isGsfCtfScPixChargeConsistent, "l2_isGsfCtfScPixChargeConsistent/O");
    Run_Tree->Branch("l2_isGsfScPixChargeConsistent", &l2_isGsfScPixChargeConsistent, "l2_isGsfScPixChargeConsistent/O");
    Run_Tree->Branch("l2_isGsfCtfChargeConsistent", &l2_isGsfCtfChargeConsistent, "l2_isGsfCtfChargeConsistent/O");
    Run_Tree->Branch("l2_CloseJetPt", &l2_CloseJetPt, "l2_CloseJetPt/F");
    Run_Tree->Branch("l2_CloseJetEta", &l2_CloseJetEta, "l2_CloseJetEta/F");
    Run_Tree->Branch("l2_CloseJetPhi", &l2_CloseJetPhi, "l2_CloseJetPhi/F");

    Run_Tree->Branch("l3M", &l3M, "l3M/F");
    Run_Tree->Branch("l3E", &l3E, "l3E/F");
    Run_Tree->Branch("l3Px", &l3Px, "l3Px/F");
    Run_Tree->Branch("l3Py", &l3Py, "l3Py/F");
    Run_Tree->Branch("l3Pz", &l3Pz, "l3Pz/F");
    Run_Tree->Branch("l3Pt", &l3Pt, "l3Pt/F");
    Run_Tree->Branch("l3Eta", &l3Eta, "l3Eta/F");
    Run_Tree->Branch("l3Phi", &l3Phi, "l3Phi/F");
    Run_Tree->Branch("l3Charge", &l3Charge, "l3Charge/F");
    Run_Tree->Branch("l3_CloseJetPt", &l3_CloseJetPt, "l3_CloseJetPt/F");
    Run_Tree->Branch("l3_CloseJetEta", &l3_CloseJetEta, "l3_CloseJetEta/F");
    Run_Tree->Branch("l3_CloseJetPhi", &l3_CloseJetPhi, "l3_CloseJetPhi/F");
    Run_Tree->Branch("l3_muId_Loose", &l3_muId_Loose, "l3_muId_Loose/O");
    Run_Tree->Branch("l3_muId_Tight", &l3_muId_Tight, "l3_muId_Tight/O");
    Run_Tree->Branch("l3_eleId_Loose", &l3_eleId_Loose, "l3_eleId_Loose/O");
    Run_Tree->Branch("l3_eleId_Tight", &l3_eleId_Tight, "l3_eleId_Tight/O");
    Run_Tree->Branch("l3_muIso", &l3_muIso, "l3_muIso/F");
    Run_Tree->Branch("l3_eleIso", &l3_eleIso, "l3_eleIso/F");
    Run_Tree->Branch("l3_eleMVANonTrg", &l3_eleMVANonTrg, "l3_eleMVANonTrg/F");
    Run_Tree->Branch("l3_eleNumHit", &l3_eleNumHit, "l3_eleNumHit/F");
    Run_Tree->Branch("l3_tauIsoVL", &l3_tauIsoVL, "l3_tauIsoVL/O");
    Run_Tree->Branch("l3_tauIso3HitL", &l3_tauIso3HitL, "l3_tauIso3HitL/O");
    Run_Tree->Branch("l3_tauIso3HitM", &l3_tauIso3HitM, "l3_tauIso3HitM/O");
    Run_Tree->Branch("l3_tauIso3HitT", &l3_tauIso3HitT, "l3_tauIso3HitT/O");
    Run_Tree->Branch("l3_tauIsoL", &l3_tauIsoL, "l3_tauIsoL/O");
    Run_Tree->Branch("l3_tauIsoM", &l3_tauIsoM, "l3_tauIsoM/O");
    Run_Tree->Branch("l3_tauIsoT", &l3_tauIsoT, "l3_tauIsoT/O");
    Run_Tree->Branch("l3_tauIsoMVA2L", &l3_tauIsoMVA2L, "l3_tauIsoMVA2L/O");
    Run_Tree->Branch("l3_tauIsoMVA2M", &l3_tauIsoMVA2M, "l3_tauIsoMVA2M/O");
    Run_Tree->Branch("l3_tauIsoMVA2T", &l3_tauIsoMVA2T, "l3_tauIsoMVA2T/O");
    Run_Tree->Branch("l3_tauIsoMVA2raw", &l3_tauIsoMVA2raw, "l3_tauIsoMVA2raw/F");
    Run_Tree->Branch("l3_tauRejMuL", &l3_tauRejMuL, "l3_tauRejMuL/O");
    Run_Tree->Branch("l3_tauRejMu2L", &l3_tauRejMu2L, "l3_tauRejMu2L/O");
    Run_Tree->Branch("l3_tauRejMuM", &l3_tauRejMuM, "l3_tauRejMuM/O");
    Run_Tree->Branch("l3_tauRejMu2M", &l3_tauRejMu2M, "l3_tauRejMu2M/O");
    Run_Tree->Branch("l3_tauRejMuT", &l3_tauRejMuT, "l3_tauRejMuT/O");
    Run_Tree->Branch("l3_tauRejMu2T", &l3_tauRejMu2T, "l3_tauRejMu2T/O");
    Run_Tree->Branch("l3_tauRejEleL", &l3_tauRejEleL, "l3_tauRejEleL/O");
    Run_Tree->Branch("l3_tauRejEleM", &l3_tauRejEleM, "l3_tauRejEleM/O");
    Run_Tree->Branch("l3_tauRejEleMVA", &l3_tauRejEleMVA, "l3_tauRejEleMVA/O");
    Run_Tree->Branch("l3_tauRejEleMVA3L", &l3_tauRejEleMVA3L, "l3_tauRejEleMVA3L/O");
    Run_Tree->Branch("l3_tauRejEleMVA3M", &l3_tauRejEleMVA3M, "l3_tauRejEleMVA3M/O");
    Run_Tree->Branch("l3_tauRejEleMVA3T", &l3_tauRejEleMVA3T, "l3_tauRejEleMVA3T/O");
    Run_Tree->Branch("l3_RefJetPt", &l3_RefJetPt, "l3_RefJetPt/F");
    Run_Tree->Branch("l3_RefJetEta", &l3_RefJetEta, "l3_RefJetEta/F");
    Run_Tree->Branch("l3_RefJetPhi", &l3_RefJetPhi, "l3_RefJetPhi/F");




    //#################################################################################################
    //###################      Starting the analysis, making loop over files    #######################
    //#################################################################################################
    //running over the
    for (int k = 0; k < fileNames.size(); k++) {

        TChain *rootTree = new TChain("t");
        rootTree->Add(fileNames[k].c_str());
        int nev = int(rootTree->GetEntries());
        TBranch* branch = rootTree->GetBranch("myevent");
        branch->SetAddress(&m);
        cout << "number of entries is : " << nev << endl;


        // running over the root files
        for (int i = 0; i < nev; i++) {
            rootTree->GetEvent(i);
            if (i % 1000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nev);
            fflush(stdout);

            //*********************************************************************************************
            //****************************    Object definitions    ***************************************
            //*********************************************************************************************
            vector<myobject> BareMuon = myCleanBareLepton(m, "mu");
            vector<myobject> BareElectron = myCleanBareLepton(m, "ele");
            vector<myobject> BareTau = myCleanBareLepton(m, "tau");


            //Number of B-jets
            int num_Bjet = bjet_Multiplicity(m);
            //*********************************************************************************************
            //****************************    PileUp re weighting    ***************************************
            //*********************************************************************************************
            int num_PU = 1;
            float PU_Weight = 1;

            if (mc12) {
                num_PU = m->PUInfo_true;
                PU_Weight = LumiWeights_12->weight(num_PU);
            }
            if (mc11) {
                num_PU = m->PUInfo;
                PU_Weight = LumiWeights_11->weight(num_PU);
            }
            //*********************************************************************************************
            //****************************    Trigger      ************************************************
            //*********************************************************************************************
            bool Trigger;
            if (mc12) Trigger = Trg_MC_12(m);
            if (mc11) Trigger = Trg_MC_11(m);
            if (data12) Trigger = Trg_Data_12(m);
            if (data11) Trigger = Trg_Data_11(m);
            //*********************************************************************************************
            //*********************************************************************************************
            //*******************    Default Values     ***************************************************
            //*********************************************************************************************
            float HighPt_Lep = 20;
            float Cor_eff = 1;
            //#################################################################################################
            //#################################################################################################
            //###############    2l2tau Analysis       #########################################################
            //#################################################################################################
            //#################################################################################################
            //#################################################################################################
            //#################################################################################################
            //#########################        Fake Rate Estimation        ##########################################
            //#################################################################################################
            //#################################################################################################
            //#################################################################################################
            if (is_mu || is_tot) {
                //##############################################################################
                // MMT
                //##############################################################################

                if (BareTau.size() > 0 && BareMuon.size() > 1 && Trigger) {
                    for (int i = 0; i < BareMuon.size(); i++) {
                        for (int j = i + 1; j < BareMuon.size(); j++) {
                            for (int k = 0; k < BareTau.size(); k++) {



                                bool first_l_HighPt = BareMuon[i].pt > HighPt_Lep;
                                bool Overlap_Dz = OverLap(BareMuon[i], BareMuon[j], BareTau[k]);
                                bool bjet_num = num_Bjet < 1;
                                Cor_eff = getCorrFactor("mmt", status_sample.c_str(), BareMuon[i], BareMuon[j], BareTau[k]);

                                bool preSelection = bjet_num && Overlap_Dz && first_l_HighPt;

                                if (preSelection) {
                                    fillTree(Run_Tree, m, PU_Weight, Cor_eff, 1, BareMuon[i], BareMuon[j], BareTau[k]);
                                }

                            }
                        }
                    }
                }
            }
            //##############################################################################
            // MET
            //##############################################################################
            if (is_mueg || is_tot) {
                if (BareTau.size() > 0 && BareElectron.size() > 0 && BareMuon.size() > 0 && Trigger) {
                    for (int i = 0; i < BareMuon.size(); i++) {
                        for (int k = 0; k < BareElectron.size(); k++) {
                            for (int l = 0; l < BareTau.size(); l++) {


                                bool first_l_HighPt = BareMuon[i].pt > HighPt_Lep || BareElectron[k].pt > HighPt_Lep;
                                bool Overlap_Dz = OverLap(BareMuon[i], BareElectron[k], BareTau[l]);
                                bool bjet_num = num_Bjet < 1;
                                Cor_eff = getCorrFactor("emt", status_sample.c_str(), BareMuon[i], BareElectron[k], BareTau[l]);

                                bool preSelection = Overlap_Dz && bjet_num && first_l_HighPt;

                                if (preSelection) {
                                    fillTree(Run_Tree, m, PU_Weight, Cor_eff, 2, BareMuon[i], BareElectron[k], BareTau[l]);
                                }
                            }
                        }
                    }
                }
            }





            //##############################################################################
            //##############################################################################
            //##############################################################################
            //##############################################################################
            //##############################################################################
            //##############################################################################
            //##############################################################################
            //##############################################################################
            //##############################################################################
            //##############################################################################

            //##############################################################################
            // BEST ZEE SELECTION


        }//loop over events


        delete rootTree;
    }
    delete m;


    fout->cd();

    Run_Tree->Write();

    map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
    map<string, TH1F*>::const_iterator jMap1 = myMap1->end();

    for (; iMap1 != jMap1; ++iMap1)
        nplot1(iMap1->first)->Write();

    map<string, TH2F*>::const_iterator iMap2 = myMap2->begin();
    map<string, TH2F*>::const_iterator jMap2 = myMap2->end();

    for (; iMap2 != jMap2; ++iMap2)
        nplot2(iMap2->first)->Write();



    fout->Close();
    return 0;
}
