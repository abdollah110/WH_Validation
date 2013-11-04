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


def _Print_RLE(file_,tree_,channel,Out):
    myFile = ROOT.TFile.Open(file_, "READ")
    myTree = myFile.Get(tree_)
    
    
    Run = n.zeros(1, dtype=float)
    Lumi = n.zeros(1, dtype=float)
    Event = n.zeros(1, dtype=float)
    
    for entry in xrange(myTree.GetEntries()):
#                if (entry % 10000 == 0): print "Entry is : ", entry
                myTree.GetEntry(entry)
                if (myTree.Channel_==2):
                    Run[0] = myTree.Run_
                    Lumi[0] = myTree.Lumi_
                    Event[0] = myTree.Event_
                    print int(Run[0]),":",int(Lumi[0]),":",int(Event[0])
    return 0
    



if __name__ == "__main__":
    _Print_RLE("test.root", "BG_Tree", 1,"rle_mmt.root")
#    _print_RLE("file.root", "BG_Tree", 1,"rle_mmt.root")
