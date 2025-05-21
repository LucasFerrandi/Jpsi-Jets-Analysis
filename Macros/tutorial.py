"""
\author Luca Micheletti <luca.micheletti@cern.ch>, CERN
"""
import argparse
import json
from array import array
import os
from os import path

from DQFitter import DQFitter
from ROOT import TF1, TH1F, TFile, TTree, gRandom, TCanvas, TLegend, kRed, kBlue, kGreen, kOrange, kGray, TF1, TGraph, TLine, TMath
import re

import numpy as np

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

def crystal_ball(x, params):
    A, B, C, D, E, F = params
    t = (x - A) / B
    if C < 0:
        t = -t
    
    absAlpha = abs(C)
    absAlpha2 = abs(E)
    
    if -absAlpha <= t < absAlpha2:
        return np.exp(-0.5 * t * t)
    
    if t < -absAlpha:
        a = (D/absAlpha)**D * np.exp(-0.5 * absAlpha**2)
        b = D/absAlpha - absAlpha
        return a / (b - t)**D
    
    if t >= absAlpha2:
        c = (F/absAlpha2)**F * np.exp(-0.5 * absAlpha2**2)
        d = F/absAlpha2 - absAlpha2
        return c / (d + t)**F
    
    return 0.0

# ROOT function wrapper
class CB2Pdf:
    def __call__(self, x, par):
        return crystal_ball(x[0], [par[0], par[1], par[2], par[3], par[4], par[5]])

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
        FileIn = TFile(InputFileName, "READ")
        JPsiDirName = "j-psi-fragmentation-function-task"
        JPsiDir = FileIn.GetDirectory(JPsiDirName)
        for hist in JPsiDir.GetListOfKeys():
            if "h_diel_z" in hist.GetName():
                print("Found h_diel_z histogram: ", hist.GetName())
                input_name = JPsiDirName + "/" + hist.GetName()
                InputHist = FileIn.Get(input_name)
                Nbinsz=InputHist.GetNbinsY()
                Test = False
                if Test:
                    Nbinsz = 5
                skipEnergyAnalysis = True
                if skipEnergyAnalysis:
                    if "Energy" in input_name:
                        print("skipEnergyAnalysis = True, then skipping Energy analysis")
                        continue
                BinResultsDicts = [] # List of z_size dictionaries, each one with 3 dicts (1 for each fit range)
                '''
                if "Energy" in input_name:
                    print("Found Energy in histogram name")
                    if "SemiInclusive" in input_name:
                        print("Found SemiInclusive in histogram name")
                        continue #temporary. Later include the same analysis for semiinclusive
                    # EDir.cd()
                    continue #temporary. Later include the same analysis for energy for analysis of Xi(z, E)
                '''
                # FileIn.Close() #Maybe this should not be commented. previously this was after "for binZ..."
                for binZ in range(1, Nbinsz + 1):
                    print("Running fit for bin index %d" % binZ)
                    MinDataRange = inputCfg["input"]["pdf_dictionary"]["fitRangeMin"][0] #the range of the data is the greater range of the fit (the larger has to be the first?)
                    MaxDataRange = inputCfg["input"]["pdf_dictionary"]["fitRangeMax"][0]
                    print("MinDataRange: ", MinDataRange)
                    print("MaxDataRange: ", MaxDataRange)
                    dqFitter = DQFitter(
                        inputCfg["input"]["input_file_name"], input_name, inputCfg["output"]["output_file_name"], MinDataRange, MaxDataRange, inputCfg["input"]["desired_bins"], True, binZ
                    )
                    dqFitter.SetFitConfig(inputCfg["input"]["pdf_dictionary"])
                    BinResultsDicts.append(dqFitter.MultiTrial())
                # print(f"BinResultsDicts:{BinResultsDicts} \n\n\n\n")
                FileOut = TFile("{}FitterResults_{}".format(inputCfg["output"]["output_file_name"], inputCfg["input"]["input_file_name"].rsplit("/", 1)[-1]), "UPDATE") #Has to be opened again due to the fact that it was closed in dqFitter.MultiTrial()
                FileOut.cd()
                if "pT" in input_name:
                    print("Found pT in histogram name")
                    DirName = "Analysis_in_pT_Ranges"
                if "Energy" in input_name:
                    print("Found Energy in histogram name")
                    DirName = "Analysis_in_Energy_Ranges"
                Dir = FileOut.GetDirectory(DirName) or FileOut.mkdir(DirName)
                # if "SemiInclusive" in input_name:
                #     print("Found SemiInclusive in histogram name")
                #     continue #temporary. Later include the same analysis for semiinclusive
                Dir.cd()
                thisRangeDir = Dir.GetDirectory(input_name.split("/")[-1]) or Dir.mkdir(input_name.split("/")[-1])
                thisRangeDir.cd() #maybe exit the directory later
                '''
                if "Energy" in input_name:
                    DirName = "Analysis_in_Energy_Ranges"
                    Dir = FileOut.GetDirectory(DirName) or FileOut.mkdir(DirName)
                    print("Found Energy in histogram name")
                    # if "SemiInclusive" in input_name:
                    #     print("Found SemiInclusive in histogram name")
                    #     continue #temporary. Later include the same analysis for semiinclusive
                    Dir.cd()
                    thisRangeDir = Dir.GetDirectory(input_name.split("/")[-1]) or Dir.mkdir(input_name.split("/")[-1])
                    thisRangeDir.cd()
                '''
                InitialBin = 2 #"=0" to get the whole x range
                for RangeName in BinResultsDicts[0]: #Iterate through "fit" ranges (0 is arbitrary). Assuming all z bins use same ranges for mass fit!
                    print("RangeName: ", RangeName)
                    # pTJetMin = input_name[-15] # string
                    # pTJetMax = input_name[-7]
                    print("input_name: ", input_name)
                    pTBounds = re.findall(r'\d+', input_name) # list of sequences of digits
                    print("pTBounds: ", pTBounds)
                    pTRangeTitle = pTBounds[0][:-3] + " to " + pTBounds[1][:-3] + " GeV"
                    hist_z_title = "J/#Psi z, p_{T} " + pTRangeTitle + " (fit range " + RangeName[10:] + ")"
                    hist_z = TH1F(f"h_z_{RangeName}", hist_z_title, Nbinsz-InitialBin, (InitialBin/Nbinsz), 1)
                    zBinHist = 1
                    for zBinDict in range(InitialBin, Nbinsz):
                        # print("zBinDict: ", zBinDict)
                        # print("zBin: ", zBinHist)
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
                
                # Canvas for each parameter. One curve for each range 
                parCanvases = []
                parDir = thisRangeDir.GetDirectory("Fit Parameters") or thisRangeDir.mkdir("Fit Parameters")
                parDir.cd()
                for i in range(len(inputCfg["input"]["pdf_dictionary"]["parName"])):
                    for j, parName in enumerate(inputCfg["input"]["pdf_dictionary"]["parName"][i]):
                    # if True: #only for testing
                        # parName = "mean_Jpsi" #only for testing
                        print("parName: ", parName)
                        c = TCanvas(f"canvas_{parName}", f"{parName}", 1000, 800)
                        c.cd()
                        # c.SetTitle(f"Fit Parameters for {parName};Bin;{parName}")
                        graphs = []
                        Colors =[kRed-3, kOrange-2, kBlue-2] #kBlue+3 #kGreen-1
                        legend = TLegend()
                        legend.SetHeader("Fit Mass Ranges")
                        legend.SetFillStyle(0)
                        legend.SetBorderSize(0)
                        legend.SetTextSize(0.02)
                        graphsY = [] #Y values for all (3) ranges. Important for calculation of mean
                        for r, RangeName in enumerate(BinResultsDicts[0]): # "DictRange_1.8 to 4.2 GeV", etc.
                            graphY = [] #Y values for this range
                            for zBinDict in range(InitialBin, Nbinsz):
                                graphY.append(BinResultsDicts[zBinDict][RangeName][parName]) #Value of the parameter for each z bin
                            parGraph = TGraph(Nbinsz-InitialBin, array("f", np.arange(InitialBin*binWidth, Nbinsz*binWidth, binWidth)), array("f", graphY))
                            # parGraph.SetTitle(f"{RangeName[10:]}")
                            parGraph.SetMarkerStyle(21)  # Different marker styles (20, 21, 22...)
                            parGraph.SetMarkerColor(Colors[r])
                            parGraph.SetMarkerSize(1.0)
                            parGraph.SetLineColor(Colors[r])
                            # parGraph.SetLineWidth(2)
                            legend.AddEntry(parGraph, f"{RangeName[10:]}", "p")
                            graphs.append(parGraph)
                            graphsY.append(graphY)
                        graphs[0].SetTitle(f"{parName} values from fit;z Bin (lower bound);{parName}")
                        graphs[0].Draw("APL")
                        graphsY = np.array(graphsY)
                        meanGraphY = [] #Graph for the mean of the parameter for different ranges
                        for k, zBinDict in enumerate(range(InitialBin, Nbinsz)):
                            graphsY_column = graphsY[:, k]
                            mean = np.mean(graphsY_column)
                            meanGraphY.append(mean)
                            # std_dev = TMath.StdDev(len(graphsY_column), graphsY_column) #PENSAR EM COMO FAZER!
                            # mean_error = std_dev / TMath.Sqrt(len(graphsY_column)) #PENSAR EM COMO FAZER!
                        meanGraph = TGraph(Nbinsz-InitialBin, array("f", np.arange(InitialBin*binWidth, Nbinsz*binWidth, binWidth)), array("f", meanGraphY))
                        meanGraph.SetMarkerStyle(34)
                        legend.AddEntry(meanGraph, "Mean", "p")
                        meanGraph.Draw("PL SAME")
                        
                        # Horizontal lines representing min, max and initial values for parameters fit
                        y_Min = inputCfg["input"]["pdf_dictionary"]["parLimMin"][i][j]
                        lineMin = TLine(graphs[0].GetXaxis().GetXmin(), y_Min,
                        graphs[0].GetXaxis().GetXmax(), y_Min)
                        lineMin.SetLineColor(kGray+3)
                        lineMin.SetLineStyle(9)
                        lineMin.SetLineWidth(2)
                        lineMin.Draw("same")
                        y_Max = inputCfg["input"]["pdf_dictionary"]["parLimMax"][i][j]
                        lineMax = TLine(graphs[0].GetXaxis().GetXmin(), y_Max,
                        graphs[0].GetXaxis().GetXmax(), y_Max)
                        lineMax.SetLineColor(kGray+3)
                        lineMax.SetLineStyle(5)
                        lineMax.SetLineWidth(2)
                        lineMax.Draw("same")
                        y_Ini = inputCfg["input"]["pdf_dictionary"]["parVal"][i][j]
                        lineIni = TLine(graphs[0].GetXaxis().GetXmin(), y_Ini,
                        graphs[0].GetXaxis().GetXmax(), y_Ini)
                        lineIni.SetLineColor(kGray+3)
                        lineIni.SetLineStyle(2)
                        lineIni.SetLineWidth(2)
                        lineIni.Draw("same")
                        legend.AddEntry(lineMax, f"Max Value", "l")
                        legend.AddEntry(lineIni, f"Initial Value", "l")
                        legend.AddEntry(lineMin, f"Min Value", "l")
                        legend.Draw()
                        for graph in graphs[1:]:
                        # for graph in graphs:
                            graph.Draw("PL SAME")
                        # c.Update()
                        parCanvases.append(c)
                        c.Write()
                # for canva in parCanvases:
                #     canva.Write()

                # Plot function with parameters from json (without fitting)
                c = TCanvas("c", "Double Crystal Ball", 1000, 800)
                f = TF1("cb2", CB2Pdf(), 0, 5, 6)
                parValIni = inputCfg["input"]["pdf_dictionary"]["parVal"][0]
                f.SetParameters(*parValIni)
                f.SetParNames("A", "B", "C", "D", "E", "F")
                f.SetLineColor(kBlue+3)
                f.SetLineWidth(2)
                f.SetTitle("CB2 Function With Initial-Values Parameters;x;Probability Density")  
                f.SetMaximum(1.2)     
                f.Draw()
                c.Write()  # Save canvas to ROOT file
                f.Write()  # Save function to ROOT file
                c.SaveAs("CB2_3.png")

                FileOut.Close() #Maybe this shouldnt be commented
                print("Fit for multiple x projections: Done!")
        print("Output created: {}FitterResults_{}".format(inputCfg["output"]["output_file_name"], inputCfg["input"]["input_file_name"].rsplit("/", 1)[-1]))
            # dqFitter.PlotInitialParameters("JpsiPdf", "initial_signal.png")

main()
