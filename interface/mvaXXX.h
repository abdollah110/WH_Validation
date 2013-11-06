/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "/Users/abdollahmohammadi/Downloads/root/tmva/test/TMVAGui.C"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

//#if not defined(__CINT__) || defined(__MAKECINT__)
//#endif

using namespace TMVA;

float muonJetPt, muonPt, numJets20 ;
float KNNforMuon(TMVA::Reader *reader, float muonJetPt_ = 20, float muonPt_ = 20, float numJets20_ = 20) {
    gROOT->ProcessLine(".O0"); // turn off optimization in CINT
    muonJetPt = muonJetPt_;
    muonPt = muonPt_;
    numJets20 = numJets20_;
    return reader->EvaluateMVA("KNN method");
}

float electronJetPt, electronPt;
float KNNforElectron(TMVA::Reader *reader, float electronJetPt_ = 20, float electronPt_ = 20, float numJets20_ = 20) {
    gROOT->ProcessLine(".O0"); // turn off optimization in CINT
    electronJetPt = electronJetPt_;
    electronPt = electronPt_;
    numJets20 = numJets20_;
    return reader->EvaluateMVA("KNN method");
}
