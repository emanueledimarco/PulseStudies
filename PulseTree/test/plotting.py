#!/usr/bin/env python
from math import *
import re
import os, os.path
from array import array

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] tree.root")
    (options, args) = parser.parse_args()

    formats = ["pdf","png"]

    ROOT.gROOT.cd()
    tfile = ROOT.TFile(args[0])
    tree = tfile.Get("PulseTree/pulses")
    
    resolte = ROOT.TH1F("resoltemp","",40,-2,2)
    resolutions = {"EB":ROOT.TH1F("resol_EB","",40,-0.2,0.2),
                   "EE":ROOT.TH1F("resol_EE","",40,-0.2,1.)
                   }
    cuts = {"EB":"abs(eta)<1.479",
            "EE":"abs(eta)>1.479"}
    drawopt = {"EB":"",
               "EE":"same"}

    common_cuts = "gain[5]>1"

    MfVsW = ROOT.TProfile("MfVsW","",40,2000,9000)
    
    canv = ROOT.TCanvas("canv","",600,600)
    for k,cut in cuts.iteritems():
        print "Running cut = ",cut
        ROOT.gStyle.SetOptStat(1111)
        tree.Draw( "(ampl_multifit2[5]-ampl_multifit[5])/ampl_multifit2[5]*100. >> resol_%s" % k, "%s && %s" % (cut,common_cuts))
        resolutions[k].SetLineColor(ROOT.kBlack)
        resolutions[k].SetLineWidth(2)
        resolutions[k].GetXaxis().SetTitle("(mf_{fix}-mf_{std})/mf_{fix} [%]")
        resolutions[k].GetYaxis().SetTitle("entries")
        resolutions[k].Draw()
        [canv.SaveAs("plots/resol_%s.%s" % (k,f)) for f in formats]

        ROOT.gStyle.SetOptStat(0)
        MfVsWx = MfVsW.Clone("MfVsW_%s" % k)
        tree.Draw("ampl_multifit2[5]/ampl_weights:ampl_weights >> MfVsW_%s" % k, cut)
        MfVsWx.GetXaxis().SetTitle("ampl_weights")
        MfVsWx.GetYaxis().SetTitle("mf_{fix}/weights")
        MfVsWx.GetYaxis().SetRangeUser(0.99,1.01)
        MfVsWx.SetLineColor(ROOT.kBlack)
        MfVsWx.SetLineWidth(2)
        MfVsWx.SetMarkerStyle(ROOT.kOpenCircle)
        MfVsWx.Draw(drawopt[k])
        if k=="EE": [canv.SaveAs("plots/techtest.%s" % f) for f in formats]

    profEE = ROOT.TProfile("profEE","",10,1.479,3.0)
    tree.Draw( "(ampl_multifit2[5]-ampl_multifit[5])/ampl_multifit2[5]*100:abs(eta) >> profEE", common_cuts + " && " + cuts["EE"]) 
    ROOT.gStyle.SetOptStat(0)
    profEE.SetLineColor(ROOT.kBlack)
    profEE.SetLineWidth(2)
    profEE.SetMarkerStyle(ROOT.kOpenCircle)
    profEE.GetXaxis().SetTitle("|#eta|")
    profEE.GetYaxis().SetTitle("(mf_{fix}-mf_{std})/mf_{fix} [%]")
    profEE.Draw()
    [canv.SaveAs("plots/resolEEvsEta.%s" % f) for f in formats]


