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



    //    cout.setf(ios::fixed, ios::floatfield);
    //    cout.precision(4);

    Run_Tree->SetBranchAddress("Channel", &Channel);
    Run_Tree->SetBranchAddress("Run", &Run);
    Run_Tree->SetBranchAddress("Lumi", &Lumi);
    Run_Tree->SetBranchAddress("Event", &myEvent);
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

    Run_Tree->SetBranchAddress("l1_CloseJetPt", &l1_CloseJetPt);
    Run_Tree->SetBranchAddress("l1_CloseJetEta", &l1_CloseJetEta);
    Run_Tree->SetBranchAddress("l1_CloseJetPhi", &l1_CloseJetPhi);

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
    Run_Tree->SetBranchAddress("l2_CloseJetPt", &l2_CloseJetPt);
    Run_Tree->SetBranchAddress("l2_CloseJetEta", &l2_CloseJetEta);
    Run_Tree->SetBranchAddress("l2_CloseJetPhi", &l2_CloseJetPhi);


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
    BG_Tree->Branch("myEvent_", &Event_, "myEvent_/I");

    BG_Tree->Branch("HMass_", &HMass_, "HMass_/F");
    BG_Tree->Branch("SVMass_", &SVMass_, "SVMass_/F");

    BG_Tree->Branch("l3Pt_", &l3Pt_, "l3Pt_/F");
    BG_Tree->Branch("l3Eta_", &l3Eta_, "l3Eta_/F");
    BG_Tree->Branch("l3_CloseJetPt_", &l3_CloseJetPt_, "l3_CloseJetPt_/F");
    BG_Tree->Branch("l3_CloseJetEta_", &l3_CloseJetEta_, "l3_CloseJetEta_/F");


    //###############################################################################################
    //    //######################################################
    //    //Booking MVA just once eleEM
    //    //######################################################
    TMVA::Reader *reader_eleEM = new TMVA::Reader("!Color:!Silent");

    reader_eleEM->AddVariable("electronJetPt", &electronJetPt);
    reader_eleEM->AddVariable("electronPt", &electronPt);
    reader_eleEM->AddVariable("numJets20", &numJets20);


    TString dir_eleEM = "Weight_emu/";
    TString prefix_eleEM = "em_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k100.";

    // Book method(s)
    TString methodName_eleEM = TString("KNN") + TString(" method");
    TString weightfile_eleEM = dir_eleEM + prefix_eleEM + TString("KNN") + TString(".weights.xml");
    reader_eleEM->BookMVA(methodName_eleEM, weightfile_eleEM);
    //        //######################################################
    //    //######################################################
    //    //Booking MVA just once muEM
    //    //######################################################
    TMVA::Reader *reader_muEM = new TMVA::Reader("!Color:!Silent");

    reader_muEM->AddVariable("muonJetPt", &muonJetPt);
    reader_muEM->AddVariable("muonPt", &muonPt);
    reader_muEM->AddVariable("numJets20", &numJets20);


    TString dir_muEM = "Weight_emu/";
    TString prefix_muEM = "em_Mwjets_pt10_h2taucuts_muonInfo_k100.";

    // Book method(s)
    TString methodName_muEM = TString("KNN") + TString(" method");
    TString weightfile_muEM = dir_muEM + prefix_muEM + TString("KNN") + TString(".weights.xml");
    reader_muEM->BookMVA(methodName_muEM, weightfile_muEM);
    //        //######################################################
    //    //######################################################
    //    //Booking MVA just once muMM
    //    //######################################################
    TMVA::Reader *reader_muMM = new TMVA::Reader("!Color:!Silent");

    reader_muMM->AddVariable("muonJetPt", &muonJetPt);
    reader_muMM->AddVariable("muonPt", &muonPt);
    reader_muMM->AddVariable("numJets20", &numJets20);


    TString dir_muMM = "Weight_mumu/";
    TString prefix_muMM = "mm_wjets_pt10_h2taucuts_muonInfo_k100.";

    // Book method(s)
    TString methodName_muMM = TString("KNN") + TString(" method");
    TString weightfile_muMM = dir_muMM + prefix_muMM + TString("KNN") + TString(".weights.xml");
    reader_muMM->BookMVA(methodName_muMM, weightfile_muMM);
    //        //######################################################
    //        //TREE FOR MVA
    //        //######################################################
    //
    //        float N_LT, N_MET;
    //        float N_l3_tauIsoMVA2raw, N_l4_tauIsoMVA2raw, N_IsoTot;
    //        float N_l3_LepId, N_l3_LepIso, N_TMass;
    //
    //        // MVA for LLTT ######################################################
    //        TTree * MVATree = new TTree("MVATree", "MVATree");
    //        // To force a memory-resident Tree
    //        MVATree->SetDirectory(0);
    //        MVATree->Branch("N_LT", &N_LT, "N_LT/F");
    //        MVATree->Branch("N_MET", &N_MET, "N_MET/F");
    //        MVATree->Branch("N_l3_tauIsoMVA2raw", &N_l3_tauIsoMVA2raw, "N_l3_tauIsoMVA2raw/F");
    //        MVATree->Branch("N_l4_tauIsoMVA2raw", &N_l4_tauIsoMVA2raw, "N_l4_tauIsoMVA2raw/F");
    //        MVATree->Branch("N_IsoTot", &N_IsoTot, "N_IsoTot/F");

    //###############################################################################################
    //Just each categori should be filled once
    int myEvent_Double[8][9];
    memset(myEvent_Double, 0, sizeof (myEvent_Double[0][0]) * 8 * 9);
    int Ev_double_tau = 0; //need to be checked gloat instead of integer!!!!!!!!
    int Ev_double_ele = 0;
    int Ev_double_mu = 0;
    int Ev_double_Ltau = 0;
    //###############################################################################################

    Int_t nentries_wtn = (Int_t) Run_Tree->GetEntries();
    int y = 0;
    for (Int_t i = 0; i < nentries_wtn; i++) {
        Run_Tree->GetEntry(i);
        if (i % 10000 == 0) fprintf(stdout, "\r  Processed myEvents: %8d of %8d ", i, nentries_wtn);
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
        bool Mass_m1m2 = (InvarMass_F(l1E, l2E, l1Px, l2Px, l1Py, l2Py, l1Pz, l2Pz) > 20);
        bool Mass_m2tau = (InvarMass_F(l3E, l2E, l3Px, l2Px, l3Py, l2Py, l3Pz, l2Pz) > 20);
        bool Jest2mu = (mu_Size == 2);
        bool Jest0ele = (electron_Size == 0);
        bool Jest1tau = (tau_Size == 1);
        bool OS_Charge_lt = (l1Charge * l3Charge < 0);
        //        bool selectMMT = OS_Charge_lt && Jest2mu && Jest0ele && Jest1tau && SS_dimuon && Mu1Id && Mu1Iso && Mu1MatchedTrg && Mu2Id && Mu2Iso && Mu2MatchedTrg && tauIso && Mass_m1m2 && Mass_m2tau;


        bool commonCuts = OS_Charge_lt && Jest0ele && Jest1tau && SS_dimuon && Mu1MatchedTrg && Mu2MatchedTrg && tauIso && Mass_m1m2 && Mass_m2tau;
        bool selectMMT = commonCuts && Mu1Id && Mu1Iso && Mu2Id && Mu2Iso && Jest2mu;
        bool selectMMT_Cat0 = commonCuts && (!Mu1Id || !Mu1Iso) && (!Mu2Id || !Mu2Iso);
        bool selectMMT_Cat1 = commonCuts && (Mu1Id && Mu1Iso) && (!Mu2Id || !Mu2Iso);
        bool selectMMT_Cat2 = commonCuts && (!Mu1Id || !Mu1Iso) && (Mu2Id && Mu2Iso);


        // Event Selection
        if (selectMMT && (myEvent != myEvent_Double[1][1])) {
            fillTreeN(BG_Tree, Channel, 1, VisibleMass, SVMass, Run, Lumi, myEvent, l3Pt, l3Eta, l3_CloseJetPt, l3_CloseJetEta);
            plotFill("MMT", 1, 5, 0, 5, 1);
            myEvent_Double[1][1] = myEvent;
        }
        //Category 0
        if (selectMMT_Cat0 && (myEvent != myEvent_Double[2][1])) {
            fillTreeN(BG_Tree, Channel, 2, VisibleMass, SVMass, Run, Lumi, myEvent, l3Pt, l3Eta, l3_CloseJetPt, l3_CloseJetEta);
            float weight_mu1 = KNNforMuon(reader_muMM, l1_CloseJetPt, l1Pt, num_goodjet);
            float weight_mu2 = KNNforMuon(reader_muMM, l2_CloseJetPt, l2Pt, num_goodjet);
            //            float weight_mu1 = KNNforMuon(reader_muMM, l1Pt + l1_muIso*l1Pt, l1Pt, BareTau_Size + 2);
            //            float weight_mu2 = KNNforMuon(reader_muMM, l2Pt + l2_muIso*l2Pt, l2Pt, BareTau_Size + 2);
            plotFill("MMT", 2, 5, 0, 5, weight_mu2 * weight_mu1);
            myEvent_Double[2][1] = myEvent;
        }
        //Category 1
        if (selectMMT_Cat1 && (myEvent != myEvent_Double[3][1])) {
            fillTreeN(BG_Tree, Channel, 3, VisibleMass, SVMass, Run, Lumi, myEvent, l3Pt, l3Eta, l3_CloseJetPt, l3_CloseJetEta);
            float weight_mu2 = KNNforMuon(reader_muMM, l2_CloseJetPt, l2Pt, num_goodjet);
            //            float weight_mu2 = KNNforMuon(reader_muMM, l2Pt + l2_muIso*l2Pt, l2Pt, BareTau_Size + 2);
            plotFill("MMT", 3, 5, 0, 5, weight_mu2);
            myEvent_Double[3][1] = myEvent;
        }
        //Category 2
        if (selectMMT_Cat2 && (myEvent != myEvent_Double[4][1])) {
            fillTreeN(BG_Tree, Channel, 4, VisibleMass, SVMass, Run, Lumi, myEvent, l3Pt, l3Eta, l3_CloseJetPt, l3_CloseJetEta);
            float weight_mu1 = KNNforMuon(reader_muMM, l1_CloseJetPt, l1Pt, num_goodjet);
            //            float weight_mu1 = KNNforMuon(reader_muMM, l1Pt + l1_muIso*l1Pt, l1Pt, BareTau_Size + 2);
            plotFill("MMT", 4, 5, 0, 5, weight_mu1);
            myEvent_Double[4][1] = myEvent;
        }



        //####################################################
        // EMT FakeRateation
        //####################################################
        if (Channel == 2) {

            bool SS_dilep = l1Charge * l2Charge > 0;
            bool Mu1Id = l1_muId > 0;
            bool Mu1Iso = (abs(l1Eta) < 1.479 && l1_muIso < 0.15) || (fabs(l1Eta) > 1.479 && l1_muIso < 0.10);
            bool Mu1MatchedTrg = l1_muTrgObjMatch || l1_eleTrgObjMatch;
            bool Ele2Id = l2_eleId > 0;
            bool Ele2Iso = (fabs(l2Eta) < 1.479 && l2_muIso < 0.15) || (fabs(l2Eta) > 1.479 && l2_eleIso < 0.10);
            bool Ele2Conversion = l2_passConversionVeto;
            bool Ele2MissingHit = l2_eleNumHit < 1;
            bool Ele2_chargeConsistancy1 = l2_isGsfCtfScPixChargeConsistent == l2_isGsfCtfChargeConsistent;
            bool Ele2_chargeConsistancy2 = l2_isGsfCtfScPixChargeConsistent == l2_isGsfScPixChargeConsistent;
            bool Ele2MatchedTrg = l2_eleTrgObjMatch || l2_muTrgObjMatch;
            bool tauIso_1 = l3_tauIso3HitL;
//                    bool tauIso_2 = l3_tauRejEleMVA3L;
//            bool tauIso_3 = 1;
            bool tauIso_3 = l3_tauRejMuL;
            bool Mass_m1m2 = InvarMass_F(l1E, l2E, l1Px, l2Px, l1Py, l2Py, l1Pz, l2Pz) > 20;
            float Mass_mu2tau = InvarMass_F(l3E, l1E, l3Px, l1Px, l3Py, l1Py, l3Pz, l1Pz);
            float Mass_ele2tau = InvarMass_F(l3E, l2E, l3Px, l2Px, l3Py, l2Py, l3Pz, l2Pz);
            bool l2TauMass = (l2Pt < l1Pt ? Mass_ele2tau > 20 : Mass_mu2tau > 20);
            float Mass_ele2tauDiff = fabs(InvarMass_F(l3E, l2E, l3Px, l2Px, l3Py, l2Py, l3Pz, l2Pz) - 91);
            bool MassReq = (Mass_ele2tauDiff < 20 && l3_tauRejEleMVA3M) || (Mass_ele2tauDiff > 20 && l3_tauRejEleMVA3L);
            bool Jest1mu = mu_Size == 1;
            bool Jest1ele = electron_Size == 1;
            bool Jest1tau = tau_Size == 1;
            bool OS_Charge_lt = l1Charge * l3Charge < 0;

            bool commonCuts = Jest1tau && OS_Charge_lt && MassReq && Ele2Conversion && Ele2MissingHit && Ele2_chargeConsistancy1 && Ele2_chargeConsistancy2 && SS_dilep && tauIso_1 && tauIso_3 && Mass_m1m2 && l2TauMass ;
            bool selectEMT = commonCuts && Mu1Id && Mu1Iso && Ele2Id && Ele2Iso && Jest1mu && Jest1ele;
            bool selectEMT_Cat0 = commonCuts && (!Mu1Id || !Mu1Iso) && (!Ele2Id || !Ele2Iso);
            bool selectEMT_Cat1 = commonCuts && (Mu1Id && Mu1Iso) && (!Ele2Id || !Ele2Iso) && Jest1mu;
            bool selectEMT_Cat2 = commonCuts && (!Mu1Id || !Mu1Iso) && (Ele2Id && Ele2Iso) && Jest1ele;


            // Event Selection
            if (selectEMT && (myEvent != myEvent_Double[1][2])) {
                fillTreeN(BG_Tree, Channel, 1, VisibleMass, SVMass, Run, Lumi, myEvent, l3Pt, l3Eta, l3_CloseJetPt, l3_CloseJetEta);
                plotFill("EMT", 1, 5, 0, 5, 1);
                myEvent_Double[1][2] = myEvent;
            }
            //Category 0 
            if (selectEMT_Cat0 && (myEvent != myEvent_Double[2][2])) {
                fillTreeN(BG_Tree, Channel, 2, VisibleMass, SVMass, Run, Lumi, myEvent, l3Pt, l3Eta, l3_CloseJetPt, l3_CloseJetEta);
                float weight_mu = KNNforMuon(reader_muEM, l1_CloseJetPt, l1Pt, num_goodjet);
                float weight_ele = KNNforElectron(reader_eleEM, l2_CloseJetPt, l2Pt, num_goodjet);
//                float weight_mu = KNNforMuon(reader_muEM, l1Pt + l1_muIso*l1Pt, l1Pt, BareTau_Size + 2);
//                float weight_ele = KNNforElectron(reader_eleEM, l2Pt + l2_eleIso*l2Pt, l2Pt, BareTau_Size + 2);
                plotFill("EMT", 2, 5, 0, 5, weight_mu * weight_ele);
                myEvent_Double[2][2] = myEvent;
            }
            //Category 1
            if (selectEMT_Cat1 && (myEvent != myEvent_Double[3][2])) {
                fillTreeN(BG_Tree, Channel, 3, VisibleMass, SVMass, Run, Lumi, myEvent, l3Pt, l3Eta, l3_CloseJetPt, l3_CloseJetEta);
                float weight_ele = KNNforElectron(reader_eleEM, l2_CloseJetPt, l2Pt, num_goodjet);
//                float weight_ele = KNNforElectron(reader_eleEM, l2_CloseJetPt, l2Pt, num_goodjet);
                plotFill("EMT", 3, 5, 0, 5, weight_ele);
                myEvent_Double[3][2] = myEvent;
            }
            //Category 2
            if (selectEMT_Cat2 && (myEvent != myEvent_Double[4][2])) {
                fillTreeN(BG_Tree, Channel, 4, VisibleMass, SVMass, Run, Lumi, myEvent, l3Pt, l3Eta, l3_CloseJetPt, l3_CloseJetEta);
                float weight_mu = KNNforMuon(reader_muEM, l1_CloseJetPt, l1Pt, num_goodjet);
//                float weight_mu = KNNforMuon(reader_muEM, l1_CloseJetPt, l1Pt, num_goodjet);
                plotFill("EMT", 4, 5, 0, 5, weight_mu);
                myEvent_Double[4][2] = myEvent;
            }




            //        if (Run == 1 && Lumi == 6414 && myEvent == 1923808) { //1111111111111111011
            //        if (Run == 1 && Lumi == 4722 && myEvent == 1416276) { // ?
            //        if (Run == 1 && Lumi == 11624 && myEvent == 3486306) { //  1111111111111111011
            //        if (Run == 1 && Lumi == 3931 && myEvent == 1179002) { //
            //        if (Run == 1 && Lumi == 13211 && myEvent == 3962445) { // ?
//            1:6414:1923808
//            if (Run == 1 && Lumi == 8227 && myEvent == 2467446) { // 111111111111111110
//            if (Run == 1 && Lumi == 3952 && myEvent == 1185308) { //
//                if (Run == 1 && Lumi == 4722 && myEvent == 1416276) {
//if (Run == 1 && Lumi == 13211 && myEvent == 3962445) {
//if (Run == 1 && Lumi == 3734 && myEvent == 1119716) {
//if (Run == 1 && Lumi == 8227 && myEvent == 2467446) {
//if (Run == 1 && Lumi == 86 && myEvent == 25557) {
//if (Run == 1 && Lumi == 3952 && myEvent == 1185308) {
//if (Run == 1 && Lumi == 15248 && myEvent == 4573251) {?
//if (Run == 1 && Lumi == 528 && myEvent == 158151) {?
//if (Run == 1 && Lumi == 15188 && myEvent == 4555454) {  101111101111101111++++++Mass= 93.2521
//if (Run == 1 && Lumi == 6204 && myEvent == 1860747) {101111100111001111++++++Mass= 54.8918
//if (Run == 1 && Lumi == 16461 && myEvent == 4937349) {?
//if (Run == 1 && Lumi == 7003 && myEvent == 2100294) {?
//if (Run == 1 && Lumi == 2864 && myEvent == 858957) {?
//if (Run == 1 && Lumi == 8514 && myEvent == 2553549) {?
//if (Run == 1 && Lumi == 13174 && myEvent == 3951356) {?
//if (Run == 1 && Lumi == 849 && myEvent == 254662) {?
//if (Run == 1 && Lumi == 14640 && myEvent == 4391012) {?
//if (Run == 1 && Lumi == 13678 && myEvent == 4102416) {?
//if (Run == 1 && Lumi == 5869 && myEvent == 1760171) {?
//if (Run == 1 && Lumi == 6127 && myEvent == 1837553) {?
if (Run == 1 && Lumi == 9632 && myEvent == 2888996) {

                cout << Run << " " << Lumi << " " << myEvent << "\n";
                cout << Jest1tau << OS_Charge_lt << MassReq << Ele2Conversion << Ele2MissingHit << Ele2_chargeConsistancy1 << Ele2_chargeConsistancy2 << SS_dilep << tauIso_1 << tauIso_3 << Mass_m1m2 << l2TauMass << Mu1Id << Mu1Iso << Ele2Id << Ele2Iso << Jest1mu << Jest1ele;
                cout << "++++++Mass= " << InvarMass_F(l3E, l2E, l3Px, l2Px, l3Py, l2Py, l3Pz, l2Pz) << "\n";

            }

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

