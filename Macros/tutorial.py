"""
\author Luca Micheletti <luca.micheletti@cern.ch>, CERN
"""
import argparse
import json
from array import array
import os
from os import path

from DQFitter import DQFitter
from ROOT import TF1, TH1F, TFile, TTree, gRandom, TCanvas, TLegend, kRed

import re

def GenerateTutorialSample():
    """
    This method create the sample for the tutorial
    """
    nEvents = 1000000
    SigOverBkg1 = 0.03
    SigOverBkg2 = SigOverBkg1 / 10.
    fOut = TFile("tutorial.root", "RECREATE")

    funcMassBkg = TF1("funcMassBkg", "expo", 2.0, 5.0)
    funcMassBkg.SetParameter(0, 0.00)
    funcMassBkg.SetParameter(1, -0.5)

    funcMassSig1 = TF1("funcMassSig1", "gaus(0) + gaus(3)", 2.0, 5.0)
    funcMassSig1.SetParameter(0, 1.0)
    funcMassSig1.SetParameter(1, 3.1)
    funcMassSig1.SetParameter(2, 0.07)
    funcMassSig1.SetParameter(3, 1.0)
    funcMassSig1.SetParameter(4, 3.1)
    funcMassSig1.SetParameter(5, 0.09)

    funcMassSig2 = TF1("funcMassSig2", "gaus(0) + gaus(3)", 2.0, 5.0)
    funcMassSig2.SetParameter(0, 1.0)
    funcMassSig2.SetParameter(1, 3.686)
    funcMassSig2.SetParameter(2, 1.05 * 0.07)
    funcMassSig2.SetParameter(3, 1.0)
    funcMassSig2.SetParameter(4, 3.686)
    funcMassSig2.SetParameter(5, 1.05 * 0.09)

    histMass = TH1F("histMass", "histMass", 100, 2.0, 5.0)
    histMass.FillRandom("funcMassBkg", int(nEvents - (nEvents * SigOverBkg1)))
    histMass.FillRandom("funcMassSig1", int(nEvents * SigOverBkg1))
    histMass.FillRandom("funcMassSig2", int(nEvents * SigOverBkg2))
    histMass.Write()

    m = array("f", [0.0])
    tree = TTree("data", "data")
    tree.Branch("m", m, "m/F")

    for iEvent in range(0, nEvents):
        seed = gRandom.Rndm()
        if seed > SigOverBkg1:
            m[0] = funcMassBkg.GetRandom()
        else:
            if seed > SigOverBkg2:
                m[0] = funcMassSig1.GetRandom()
            else:
                m[0] = funcMassSig2.GetRandom()
        tree.Fill()
    tree.Write()

    fOut.Close()


def main():
    parser = argparse.ArgumentParser(description="Arguments to pass")
    parser.add_argument(
        "cfgFileName", metavar="text", default="configFit.json", help="config file name"
    )
    parser.add_argument(
        "--gen_tutorial", help="generate tutorial sample", action="store_true"
    )
    parser.add_argument("--run_fit", help="run the multi trial", action="store_true")
    parser.add_argument("--run_fit_projections", help="run the fitter for multiple x projections", action="store_true")
    args = parser.parse_args()

    print("Loading task configuration: ...", end="\r")
    with open(args.cfgFileName, "r") as jsonCfgFile:
        inputCfg = json.load(jsonCfgFile)
    print("Loading task configuration: Done!")

    if args.gen_tutorial:
        GenerateTutorialSample()

    if args.run_fit:
        if not path.isdir(inputCfg["output"]["output_file_name"]):
            os.system("mkdir -p %s" % (inputCfg["output"]["output_file_name"]))
        dqFitter = DQFitter(
            inputCfg["input"]["input_file_name"], inputCfg["input"]["input_name"], inputCfg["output"]["output_file_name"], 2, 5
        )
        dqFitter.SetFitConfig(inputCfg["input"]["pdf_dictionary"])
        dqFitter.MultiTrial()
    if args.run_fit_projections:
        print("Running fit for multiple x projections")
        if not path.isdir(inputCfg["output"]["output_file_name"]):
            os.system("mkdir -p %s" % (inputCfg["output"]["output_file_name"]))
        InputFileName = inputCfg["input"]["input_file_name"]
        FileIn = TFile(InputFileName, "UPDATE")
        InputHist = FileIn.Get(inputCfg["input"]["input_name"])
        Nbinsz=InputHist.GetNbinsY()
        Test = False
        if Test:
            Nbinsz = 3
        print("Nbinsz: ", Nbinsz)
        BinResultsDicts = []
        for binZ in range(1, Nbinsz + 1):
            print("Running fit for bin index %d" % binZ)
            MinDataRange = inputCfg["input"]["pdf_dictionary"]["fitRangeMin"][0] #the range of the data is the greater range of the fit (the greater has to be the first?)
            MaxDataRange = inputCfg["input"]["pdf_dictionary"]["fitRangeMax"][0]
            print("MinDataRange: ", MinDataRange)
            print("MaxDataRange: ", MaxDataRange)
            dqFitter = DQFitter(
                inputCfg["input"]["input_file_name"], inputCfg["input"]["input_name"], inputCfg["output"]["output_file_name"], MinDataRange, MaxDataRange, inputCfg["input"]["desired_bins"], True, binZ
            )
            dqFitter.SetFitConfig(inputCfg["input"]["pdf_dictionary"])
            BinResultsDicts.append(dqFitter.MultiTrial())
        print(f"BinResultsDicts:{BinResultsDicts} \n\n\n\n")
        FileIn.Close()
        FileOut = TFile("{}{}.root".format(inputCfg["output"]["output_file_name"], inputCfg["input"]["input_name"].rsplit("/", 1)[-1]), "UPDATE")
        FileOut.cd()
        InitialBin = 2 #=0 to get the whole x range
        for RangeName in BinResultsDicts[0]: #Iterate through dictionary keys ("fit" ranges). Assuming all z bins use same ranges for mass fit!
            print("RangeName: ", RangeName)
            # pTJetMin = inputCfg["input"]["input_name"][-15] # string
            # pTJetMax = inputCfg["input"]["input_name"][-7]
            pTBounds = re.findall(r'\d+', inputCfg["input"]["input_name"])
            pTRangeTitle = pTBounds[0][:-3] + " to " + pTBounds[1][:-3] + " GeV/c"
            hist_z = TH1F(f"h_z_Range{RangeName}", f"J/#Psi z fit range {RangeName[10:]} & pT {pTRangeTitle}", 10-InitialBin, (InitialBin*0.1), 1)
            zBinHist = 1
            for zBinDict in range(InitialBin, Nbinsz):
                print("zBinDict: ", zBinDict)
                print("zBin: ", zBinHist)
                NJPsi = BinResultsDicts[zBinDict][RangeName]["sig_Jpsi"]
                NJPsiErr = BinResultsDicts[zBinDict][RangeName]["sig_Jpsi_err"] 
                print("NJpsi: ", NJPsi)
                hist_z.SetBinContent(zBinHist, NJPsi)
                hist_z.SetBinError(zBinHist, NJPsiErr)
                zBinHist += 1
            hist_z.SetTitleSize(0.01, "t")
            hist_z.SetMarkerStyle(24) 
            hist_z.SetMarkerSize(1.0)  
            hist_z.SetMarkerColor(kRed-3)
            hist_z.SetLineColor(kRed-3)
            hist_z.SetLineWidth(2)
            hist_z.SetStats(0)
            hist_z.GetYaxis().SetTitle("#frac{1}{N^{J/#Psi}_{jets}} #frac{dN}{dz^{ch.}}")
            hist_z.GetXaxis().SetTitle("z^{ch.}_{J/#Psi}")
            # hist_z.Scale(1. / NJPsi) #Normalize to 1
            hist_z.Write()
            cz = TCanvas(f"c{RangeName}", f"Canvas{RangeName}", 1400, 1080)
            hist_z.Draw()
            NJPsiTot = hist_z.Integral()
            legend = TLegend(0.7, 0.7, 0.9, 0.9)
            legend.AddEntry(hist_z, "pp #sqrt{s} = 13.6 TeV, inclusive J/#Psi", "")
            legend.AddEntry(hist_z, "p_{T, J/#Psi} > 1 GeV/c", "")
            legend.AddEntry(hist_z, pTBounds[0][:-3] + " < p_{T, Jet} < " + pTBounds[1][:-3] + " GeV/c, R = 0.4", "")
            legend.AddEntry(hist_z, "data, |y| < 0.9", "pe")
            legend.AddEntry(hist_z, "N_{{J/#Psi}} = {:.2f}".format(NJPsiTot), "")
            # legend.SetEntrySeparation(0.1)
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            legend.Draw()
            cz.Update()
            cz.Write()
            czNormalized = TCanvas(f"c{RangeName}norm", f"Canvas{RangeName} Normalized", 1400, 1080)
            hist_z_norm = hist_z.Clone()
            binWidth = hist_z.GetBinWidth(1) #assuming all z bins have the same width!
            print("binWidth: ", binWidth)
            hist_z_norm.Scale(1. / (NJPsiTot*binWidth)) #Normalize to 1
            hist_z_norm.Draw("P")
            legend.Draw()
            czNormalized.Update()
            czNormalized.Write()


            # czNormalizedSliced = TCanvas(f"c{RangeName}normSliced", f"Canvas{RangeName} Normalized", 1920, 1080)
            # hist_z_normlSliced  = TH1F(hist_z, f"JPsi z With Mass Fit From Range {RangeName[10:]} Normalized", 3, hist_z.GetNbinsX())
            # NJPsiTotSliced = hist_z_normlSliced.Integral()
            # hist_z_normlSliced.Scale(1. / (NJPsiTotSliced*binWidth)) #Normalize to 1
            # hist_z_normlSliced.Draw()
            # legend.Draw()
            # czNormalizedSliced.Update()
            # czNormalizedSliced.Write()            
        FileOut.Close()
        print("Fit for multiple x projections: Done!")
        print("Output created: {}{}.root".format(inputCfg["output"]["output_file_name"], inputCfg["input"]["input_name"].rsplit("/", 1)[-1]))

main()
