#include "tr_Tree.h"

int main(int argc, char** argv) {


    std::string out = *(argv + 1);
    std::string input = *(argv + 2);

    //PRINTING THE OUTPUT name
    std::cout << "\n\n\n OUTPUT NAME IS:    " << out << endl;
    TFile *fout = TFile::Open(out.c_str(), "RECREATE");

    using namespace std;

    myMap1 = new std::map<std::string, TH1F*>();
    myMap2 = new map<string, TH2F*>();
    //
    TFile *f_Double = new TFile(input.c_str());
    TTree *Run_Tree = (TTree*) f_Double->Get("RLE_tree");
    Run_Tree->AddFriend("Mass_tree");



    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(1);

    Run_Tree->SetBranchAddress("Channel", &Channel);
    Run_Tree->SetBranchAddress("Run", &Run);
    Run_Tree->SetBranchAddress("Lumi", &Lumi);
    Run_Tree->SetBranchAddress("Event", &Event);
    Run_Tree->SetBranchAddress("IMass", &IMass);
    Run_Tree->SetBranchAddress("ZMass", &ZMass);
    Run_Tree->SetBranchAddress("HMass", &HMass);
    Run_Tree->SetBranchAddress("met", &met);
    Run_Tree->SetBranchAddress("pfmet", &pfmet);
    Run_Tree->SetBranchAddress("metPhi", &metPhi);
    Run_Tree->SetBranchAddress("pfmetPhi", &pfmetPhi);
    Run_Tree->SetBranchAddress("covMet11", &covMet11);
    Run_Tree->SetBranchAddress("covMet12", &covMet12);
    Run_Tree->SetBranchAddress("covMet21", &covMet21);
    Run_Tree->SetBranchAddress("covMet22", &covMet22);
    Run_Tree->SetBranchAddress("pfcovMet11", &pfcovMet11);
    Run_Tree->SetBranchAddress("pfcovMet12", &pfcovMet12);
    Run_Tree->SetBranchAddress("pfcovMet21", &pfcovMet21);
    Run_Tree->SetBranchAddress("pfcovMet22", &pfcovMet22);
    Run_Tree->SetBranchAddress("num_PV", &num_PV);
    Run_Tree->SetBranchAddress("num_bjet", &num_bjet);
    Run_Tree->SetBranchAddress("num_goodjet", &num_goodjet);
    Run_Tree->SetBranchAddress("eff_Correction", &eff_Correction);
    Run_Tree->SetBranchAddress("pu_Weight", &pu_Weight);

    Run_Tree->SetBranchAddress("mu_Size", &mu_Size);
    Run_Tree->SetBranchAddress("BareMuon_Size", &BareMuon_Size);
    Run_Tree->SetBranchAddress("electron_Size", &electron_Size);
    Run_Tree->SetBranchAddress("BareElectron_Size", &BareElectron_Size);
    Run_Tree->SetBranchAddress("tau_Size", &tau_Size);
    Run_Tree->SetBranchAddress("BareTau_Size", &BareTau_Size);
    Run_Tree->SetBranchAddress("mu_partTight_Size", &mu_partTight_Size);
    Run_Tree->SetBranchAddress("ele_partTight_Size", &ele_partTight_Size);

    Run_Tree->SetBranchAddress("l1M", &l1M);
    Run_Tree->SetBranchAddress("l1E", &l1E);
    Run_Tree->SetBranchAddress("l1Px", &l1Px);
    Run_Tree->SetBranchAddress("l1Py", &l1Py);
    Run_Tree->SetBranchAddress("l1Pz", &l1Pz);
    Run_Tree->SetBranchAddress("l1Pt", &l1Pt);
    Run_Tree->SetBranchAddress("l1Eta", &l1Eta);
    Run_Tree->SetBranchAddress("l1Phi", &l1Phi);
    Run_Tree->SetBranchAddress("l1Charge", &l1Charge);
    Run_Tree->SetBranchAddress("l1_muId", &l1_muId);
    Run_Tree->SetBranchAddress("l1_muIso", &l1_muIso);
    Run_Tree->SetBranchAddress("l1_eleId", &l1_eleId);
    Run_Tree->SetBranchAddress("l1_eleIso", &l1_eleIso);
    Run_Tree->SetBranchAddress("l1_eleMVANonTrg", &l1_eleMVANonTrg);
    Run_Tree->SetBranchAddress("l1_eleNumHit", &l1_eleNumHit);

    Run_Tree->SetBranchAddress("l1_muTrgObjMatch", &l1_muTrgObjMatch);
    Run_Tree->SetBranchAddress("l1_eleTrgObjMatch", &l1_eleTrgObjMatch);
    Run_Tree->SetBranchAddress("l1_passConversionVeto", &l1_passConversionVeto);
    Run_Tree->SetBranchAddress("l1_isGsfCtfScPixChargeConsistent", &l1_isGsfCtfScPixChargeConsistent);
    Run_Tree->SetBranchAddress("l1_isGsfScPixChargeConsistent", &l1_isGsfScPixChargeConsistent);
    Run_Tree->SetBranchAddress("l1_isGsfCtfChargeConsistent", &l1_isGsfCtfChargeConsistent);




    Run_Tree->SetBranchAddress("l2M", &l2M);
    Run_Tree->SetBranchAddress("l2E", &l2E);
    Run_Tree->SetBranchAddress("l2Px", &l2Px);
    Run_Tree->SetBranchAddress("l2Py", &l2Py);
    Run_Tree->SetBranchAddress("l2Pz", &l2Pz);
    Run_Tree->SetBranchAddress("l2Pt", &l2Pt);
    Run_Tree->SetBranchAddress("l2Eta", &l2Eta);
    Run_Tree->SetBranchAddress("l2Phi", &l2Phi);
    Run_Tree->SetBranchAddress("l2Charge", &l2Charge);
    Run_Tree->SetBranchAddress("l2_muId", &l2_muId);
    Run_Tree->SetBranchAddress("l2_muIso", &l2_muIso);
    Run_Tree->SetBranchAddress("l2_eleId", &l2_eleId);
    Run_Tree->SetBranchAddress("l2_eleIso", &l2_eleIso);
    Run_Tree->SetBranchAddress("l2_eleMVANonTrg", &l2_eleMVANonTrg);
    Run_Tree->SetBranchAddress("l2_eleNumHit", &l2_eleNumHit);

    Run_Tree->SetBranchAddress("l2_muTrgObjMatch", &l2_muTrgObjMatch);
    Run_Tree->SetBranchAddress("l2_eleTrgObjMatch", &l2_eleTrgObjMatch);
    Run_Tree->SetBranchAddress("l2_passConversionVeto", &l2_passConversionVeto);
    Run_Tree->SetBranchAddress("l2_isGsfCtfScPixChargeConsistent", &l2_isGsfCtfScPixChargeConsistent);
    Run_Tree->SetBranchAddress("l2_isGsfScPixChargeConsistent", &l2_isGsfScPixChargeConsistent);
    Run_Tree->SetBranchAddress("l2_isGsfCtfChargeConsistent", &l2_isGsfCtfChargeConsistent);

    Run_Tree->SetBranchAddress("l3M", &l3M);
    Run_Tree->SetBranchAddress("l3E", &l3E);
    Run_Tree->SetBranchAddress("l3Px", &l3Px);
    Run_Tree->SetBranchAddress("l3Py", &l3Py);
    Run_Tree->SetBranchAddress("l3Pz", &l3Pz);
    Run_Tree->SetBranchAddress("l3Pt", &l3Pt);
    Run_Tree->SetBranchAddress("l3Eta", &l3Eta);
    Run_Tree->SetBranchAddress("l3Phi", &l3Phi);
    Run_Tree->SetBranchAddress("l3Charge", &l3Charge);
    Run_Tree->SetBranchAddress("l3_CloseJetPt", &l3_CloseJetPt);
    Run_Tree->SetBranchAddress("l3_CloseJetEta", &l3_CloseJetEta);
    Run_Tree->SetBranchAddress("l3_CloseJetPhi", &l3_CloseJetPhi);
    Run_Tree->SetBranchAddress("l3_muId_Loose", &l3_muId_Loose);
    Run_Tree->SetBranchAddress("l3_muId_Tight", &l3_muId_Tight);
    Run_Tree->SetBranchAddress("l3_muIso", &l3_muIso);
    Run_Tree->SetBranchAddress("l3_eleId_Loose", &l3_eleId_Loose);
    Run_Tree->SetBranchAddress("l3_eleId_Tight", &l3_eleId_Tight);
    Run_Tree->SetBranchAddress("l3_eleIso", &l3_eleIso);
    Run_Tree->SetBranchAddress("l3_eleMVANonTrg", &l3_eleMVANonTrg);
    Run_Tree->SetBranchAddress("l3_eleNumHit", &l3_eleNumHit);
    Run_Tree->SetBranchAddress("l3_tauIsoVL", &l3_tauIsoVL);
    Run_Tree->SetBranchAddress("l3_tauIso3HitL", &l3_tauIso3HitL);
    Run_Tree->SetBranchAddress("l3_tauIso3HitM", &l3_tauIso3HitM);
    Run_Tree->SetBranchAddress("l3_tauIso3HitT", &l3_tauIso3HitT);
    Run_Tree->SetBranchAddress("l3_tauIsoL", &l3_tauIsoL);
    Run_Tree->SetBranchAddress("l3_tauIsoM", &l3_tauIsoM);
    Run_Tree->SetBranchAddress("l3_tauIsoT", &l3_tauIsoT);
    Run_Tree->SetBranchAddress("l3_tauIsoMVA2L", &l3_tauIsoMVA2L);
    Run_Tree->SetBranchAddress("l3_tauIsoMVA2M", &l3_tauIsoMVA2M);
    Run_Tree->SetBranchAddress("l3_tauIsoMVA2T", &l3_tauIsoMVA2T);
    Run_Tree->SetBranchAddress("l3_tauIsoMVA2raw", &l3_tauIsoMVA2raw);
    Run_Tree->SetBranchAddress("l3_tauRejMuL", &l3_tauRejMuL);
    Run_Tree->SetBranchAddress("l3_tauRejMu2L", &l3_tauRejMu2L);
    Run_Tree->SetBranchAddress("l3_tauRejMuM", &l3_tauRejMuM);
    Run_Tree->SetBranchAddress("l3_tauRejMu2M", &l3_tauRejMu2M);
    Run_Tree->SetBranchAddress("l3_tauRejMuT", &l3_tauRejMuT);
    Run_Tree->SetBranchAddress("l3_tauRejMu2T", &l3_tauRejMu2T);
    Run_Tree->SetBranchAddress("l3_tauRejEleL", &l3_tauRejEleL);
    Run_Tree->SetBranchAddress("l3_tauRejEleM", &l3_tauRejEleM);
    Run_Tree->SetBranchAddress("l3_tauRejEleMVA", &l3_tauRejEleMVA);
    Run_Tree->SetBranchAddress("l3_tauRejEleMVA3L", &l3_tauRejEleMVA3L);
    Run_Tree->SetBranchAddress("l3_tauRejEleMVA3M", &l3_tauRejEleMVA3M);
    Run_Tree->SetBranchAddress("l3_tauRejEleMVA3T", &l3_tauRejEleMVA3T);
    Run_Tree->SetBranchAddress("l3_RefJetPt", &l3_RefJetPt);
    Run_Tree->SetBranchAddress("l3_RefJetEta", &l3_RefJetEta);
    Run_Tree->SetBranchAddress("l3_RefJetPhi", &l3_RefJetPhi);


    //New BG_Tree
    TTree * BG_Tree = new TTree("BG_Tree", "BG_Tree");
    //    To force a memory-resident Tree
    BG_Tree->SetDirectory(0);

    BG_Tree->Branch("Channel_", &Channel_, "Channel_/I");
    BG_Tree->Branch("subChannel_", &subChannel_, "subChannel_/I");
    BG_Tree->Branch("Run_", &Run_, "Run_/I");
    BG_Tree->Branch("Lumi_", &Lumi_, "Lumi_/I");
    BG_Tree->Branch("Event_", &Event_, "Event_/I");

    BG_Tree->Branch("HMass_", &HMass_, "HMass_/F");
    BG_Tree->Branch("SVMass_", &SVMass_, "SVMass_/F");

    BG_Tree->Branch("l3Pt_", &l3Pt_, "l3Pt_/F");
    BG_Tree->Branch("l3Eta_", &l3Eta_, "l3Eta_/F");
    BG_Tree->Branch("l3_CloseJetPt_", &l3_CloseJetPt_, "l3_CloseJetPt_/F");
    BG_Tree->Branch("l3_CloseJetEta_", &l3_CloseJetEta_, "l3_CloseJetEta_/F");


    //###############################################################################################
    //######################################################
    //Booking MVA just once
    //######################################################
    TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");


    reader->AddVariable("N_LT", &N_LT);
    reader->AddVariable("N_MET", &N_MET);
    reader->AddVariable("N_l3_tauIsoMVA2raw", &N_IsoTau1);
    reader->AddVariable("N_l4_tauIsoMVA2raw", &N_IsoTau2);


    TString dir = "Weight_emu/";
    TString prefix = "TMVAClassification";

    // Book method(s)
    TString methodName = TString("KNN") + TString(" method");
    TString weightfile = dir + prefix + TString("_") + TString("KNN") + TString(".weights.xml");
    reader->BookMVA(methodName, weightfile);
    //######################################################
    //TREE FOR MVA
    //######################################################

    float N_LT, N_MET;
    float N_l3_tauIsoMVA2raw, N_l4_tauIsoMVA2raw, N_IsoTot;
    float N_l3_LepId, N_l3_LepIso, N_TMass;

    // MVA for LLTT ######################################################
    TTree * MVATree = new TTree("MVATree", "MVATree");
    // To force a memory-resident Tree
    MVATree->SetDirectory(0);
    MVATree->Branch("N_LT", &N_LT, "N_LT/F");
    MVATree->Branch("N_MET", &N_MET, "N_MET/F");
    MVATree->Branch("N_l3_tauIsoMVA2raw", &N_l3_tauIsoMVA2raw, "N_l3_tauIsoMVA2raw/F");
    MVATree->Branch("N_l4_tauIsoMVA2raw", &N_l4_tauIsoMVA2raw, "N_l4_tauIsoMVA2raw/F");
    MVATree->Branch("N_IsoTot", &N_IsoTot, "N_IsoTot/F");

    //###############################################################################################
    //Just each categori should be filled once
    int Event_Double[8][9];
    memset(Event_Double, 0, sizeof (Event_Double[0][0]) * 8 * 9);
    int Ev_double_tau = 0; //need to be checked gloat instead of integer!!!!!!!!
    int Ev_double_ele = 0;
    int Ev_double_mu = 0;
    int Ev_double_Ltau = 0;
    //###############################################################################################

    Int_t nentries_wtn = (Int_t) Run_Tree->GetEntries();
    int y = 0;
    for (Int_t i = 0; i < nentries_wtn; i++) {
        Run_Tree->GetEntry(i);
        if (i % 10000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
        fflush(stdout);

        //###############################################################################################
        // Default values
        //###############################################################################################
        float VisibleMass = HMass;
        //###############################################################################################
        //###########################        Optimum    #################################################
        //###############################################################################################

        //####################################################
        // MMT FakeRateation
        //####################################################

        bool SS_dimuon = l1Charge * l2Charge > 0;
        bool Mu1Id = l1_muId > 0;
        bool Mu1Iso = (abs(l1Eta) < 1.479 && l1_muIso < 0.15) || (fabs(l1Eta) > 1.479 && l1_muIso < 0.10);
        bool Mu1MatchedTrg = l1_muTrgObjMatch;
        bool Mu2Id = l2_muId > 0;
        bool Mu2Iso = (fabs(l2Eta) < 1.479 && l2_muIso < 0.20) || (fabs(l2Eta) > 1.479 && l2_muIso < 0.15);
        bool Mu2MatchedTrg = l2_muTrgObjMatch;
        bool tauIso = l3_tauIso3HitL && l3_tauRejEleMVA3L && l3_tauRejMu2T;
        bool Mass_m1m2 = InvarMass_F(l1E, l2E, l1Px, l2Px, l1Py, l2Py, l1Pz, l2Pz) > 20;
        bool Mass_m2tau = InvarMass_F(l3E, l2E, l3Px, l2Px, l3Py, l2Py, l3Pz, l2Pz) > 20;
        bool Jest1Lep = mu_Size == 2 && electron_Size == 0 && tau_Size == 1;
        bool OS_Charge_lt = l1Charge * l3Charge < 0;
        bool selectMMT = OS_Charge_lt && Jest1Lep && SS_dimuon && Mu1Id && Mu1Iso && Mu1MatchedTrg && Mu2Id && Mu2Iso && Mu2MatchedTrg && tauIso && Mass_m1m2 && Mass_m2tau;


        if (Channel == 1 && selectMMT && (Event != Event_Double[1][1])) {

            fillTreeN(BG_Tree, Channel, 0, VisibleMass, SVMass, Run, Lumi, Event, l3Pt, l3Eta, l3_CloseJetPt, l3_CloseJetEta);
            Event_Double[1][1] = Event;


        }
        //####################################################
        // EMT FakeRateation
        //####################################################

        bool SS_dilep = l1Charge * l2Charge > 0;
        Mu1Id = l1_muId > 0;
        Mu1Iso = (abs(l1Eta) < 1.479 && l1_muIso < 0.15) || (fabs(l1Eta) > 1.479 && l1_muIso < 0.10);
        //        bool Mu1MatchedTrg = l1_muTrgObjMatch;
        bool Ele2Id = l2_eleId > 0;
        bool Ele2Iso = (fabs(l2Eta) < 1.479 && l2_muIso < 0.15) || (fabs(l2Eta) > 1.479 && l2_eleIso < 0.10);
        bool Ele2Conversion = l2_passConversionVeto;
        bool Ele2MissingHit = l2_eleNumHit < 1;
        bool Ele2_chargeConsistancy1 = l2_isGsfCtfScPixChargeConsistent == l2_isGsfCtfChargeConsistent;
        bool Ele2_chargeConsistancy2 = l2_isGsfCtfScPixChargeConsistent == l2_isGsfScPixChargeConsistent;
        //        bool Mu2MatchedTrg = l2_muTrgObjMatch;
        tauIso = l3_tauIso3HitL && l3_tauRejEleMVA3L && l3_tauRejMu2T;
        Mass_m1m2 = InvarMass_F(l1E, l2E, l1Px, l2Px, l1Py, l2Py, l1Pz, l2Pz) > 20;
        bool Mass_ele2tau = InvarMass_F(l3E, l2E, l3Px, l2Px, l3Py, l2Py, l3Pz, l2Pz) > 20;
        float Mass_ele2tauDiff = fabs(InvarMass_F(l3E, l2E, l3Px, l2Px, l3Py, l2Py, l3Pz, l2Pz) - 90);
        bool MassReq = (Mass_ele2tauDiff < 20 && l3_tauRejEleMVA3M) || (Mass_ele2tauDiff > 20 && l3_tauRejEleMVA3L);
        Jest1Lep = mu_Size == 1 && electron_Size == 1 && tau_Size == 1;
        OS_Charge_lt = l1Charge * l3Charge < 0;
        bool selectEMT = OS_Charge_lt && MassReq && Ele2Conversion && Ele2MissingHit && Ele2_chargeConsistancy1 && Ele2_chargeConsistancy2 && Mass_ele2tau && Jest1Lep && SS_dilep && Mu1Id && Mu1Iso && Ele2Id && Ele2Iso && tauIso && Mass_m1m2 && Mass_ele2tau;


        if (Channel == 2 && selectEMT && (Event != Event_Double[1][2])) {

            fillTreeN(BG_Tree, Channel, 0, VisibleMass, SVMass, Run, Lumi, Event, l3Pt, l3Eta, l3_CloseJetPt, l3_CloseJetEta);
            Event_Double[1][2] = Event;


        }






    }


    fout->cd();
    BG_Tree->Write();

    map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
    map<string, TH1F*>::const_iterator jMap1 = myMap1->end();

    for (; iMap1 != jMap1; ++iMap1)
        nplot1(iMap1->first)->Write();

    map<string, TH2F*>::const_iterator iMap2 = myMap2->begin();
    map<string, TH2F*>::const_iterator jMap2 = myMap2->end();

    for (; iMap2 != jMap2; ++iMap2)
        nplot2(iMap2->first)->Write();

    fout->Close();


}

