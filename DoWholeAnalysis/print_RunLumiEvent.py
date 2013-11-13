#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="abdollahmohammadi"
__date__ ="$Nov 3, 2013 11:06:39 AM$"
import ROOT
from ROOT import Double
from ROOT import TCanvas
from ROOT import TFile
from ROOT import TH1F
from ROOT import TH2F
from ROOT import TLatex
from ROOT import TLegend
from ROOT import TNtuple
from ROOT import TProfile
from ROOT import TTree
from ROOT import gBenchmark
from ROOT import gROOT
from ROOT import gRandom
from ROOT import gStyle
from ROOT import gSystem
import numpy as n
import os


def _Print_RLE(file_,tree_,channel):
    myFile = ROOT.TFile.Open(file_, "READ")
    myTree = myFile.Get(tree_)
    
    

    
    for entry in xrange(myTree.GetEntries()):
                myTree.GetEntry(entry)
                if (myTree.Channel_==channel and myTree.subChannel_==1):
                    print int(myTree.Run_),":",int(myTree.Lumi_),":",int(myTree.myEvent_)
#                    sys.stdout.write('.')
    return 0
    



if __name__ == "__main__":
#    _Print_RLE("test_wz.root", "BG_Tree", 2)
    _Print_RLE("outMuEG_Saturday.root", "BG_Tree", 2)

#    _Print_RLE("test.root", "BG_Tree", 1,"rle_mmt.root")
