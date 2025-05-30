"""
\author Luca Micheletti <luca.micheletti@cern.ch>, CERN
"""
import argparse
import json
from array import array
import os
from os import path

from DQFitter import DQFitter
from ROOT import (
    TF1,
    TH1F,
    TFile,
    TTree,
    gRandom,
    TCanvas,
    TLegend,
    kRed,
    kBlue,
    kGreen,
    kBlack,
    kOrange,
    kGray,
    TF1,
    TGraph,
    TLine,
    TPaveText,
    THStack
)
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

def vwg(x, params):
    A, B, C = params
    sigma = B + C * ((x - A) / A)
    return np.exp(-(x - A) * (x - A) / (2. * sigma * sigma))

def funcsSum(x, paramscb2, paramsvwg):
    return crystal_ball(x, paramscb2) + vwg(x, paramsvwg)

# ROOT function wrapper
class PDF:
    def __init__(self, pdf_name):
        self.pdf_name = pdf_name

    def __call__(self, x, par):
        if self.pdf_name == "CB2":
            return crystal_ball(x[0], [par[0], par[1], par[2], par[3], par[4], par[5]])
        if self.pdf_name == "VWG":
            return vwg(x[0], [par[0], par[1], par[2]])
        
LEGEND_SIZE_DEFAULT = 0.02
LEGEND_COORDS_DEFAULT = [0.2, 0.7, 0.4, 0.9]
def createLegend(legLines, legCoords=LEGEND_COORDS_DEFAULT, textSize=LEGEND_SIZE_DEFAULT, isMC=False):
    legend = TLegend(*legCoords)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(textSize)
    for legLine in legLines:
        legend.AddEntry(legend, legLine, "")
    legend.AddEntry(legend, "Data, |y| < 0.9", "")
    if isMC:
        legend.AddEntry(legend, "MC, |y| < 0.9", "")
    return legend

TITLE_SIZE_DEFAULT = 0.033
def createPaveText(text):
    pave_text = TPaveText(0.1, 0.96, 0.9, 0.99, "NDC")
    pave_text.AddText(text)
    pave_text.SetFillColor(0)
    pave_text.SetFillStyle(0)  # Transparent background
    pave_text.SetBorderSize(0)
    pave_text.SetTextFont(42)
    pave_text.SetTextAlign(22)
    pave_text.SetTextSize(TITLE_SIZE_DEFAULT)

    return pave_text

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

    # if args.run_fit_projections:
    print("Running fit for multiple x projections")
    legLine1 = "pp #sqrt{s} = 13.6 TeV"
    legLine2 = "p_{T}^{J/#Psi} > 1 GeV/c"
    legLines = [legLine1, legLine2, ""]
    markerStyle = 21 # Square. For circle, = 20
    Colors =[kRed-2, kOrange-2, kBlue-2] # kBlue+3 kGreen-1 kRed-3
    ColorsExt = [kGreen-1, kRed-2, kOrange-2, kBlue-2] # [kRed-2, kOrange+2, kOrange-2, kGreen-1, kBlue-2]
    if not path.isdir(inputCfg["output"]["output_file_name"]):
        os.system("mkdir -p %s" % (inputCfg["output"]["output_file_name"]))
    InputFileName = inputCfg["input"]["input_file_name"]
    FileIn = TFile(InputFileName, "READ")
    JPsiDirName = "j-psi-fragmentation-function-task"
    JPsiDir = FileIn.GetDirectory(JPsiDirName)
    hists_z_mean_norm = []
    histspTInput = [hist for hist in JPsiDir.GetListOfKeys() if "h_diel_z" in hist.GetName() and "pT" in hist.GetName()]
    skipEnergyAnalysis = True
    if skipEnergyAnalysis:
        print("skipEnergyAnalysis = True, then skipping Energy analysis")
        histsInput = histspTInput
    else:
        histsEInput = [hist for hist in JPsiDir.GetListOfKeys() if "h_diel_z" in hist.GetName() and "Energy" in hist.GetName()]
        for histEInput in histsEInput:
            print("Found Energy histogram: ", histEInput.GetName())
        # TODO: Implement the analysis for energy histograms aiming Xi analysis
    onlySemiInclusivepT = False # For testing
    if onlySemiInclusivepT:
        print("Only SemiInclusive pT histograms will be analyzed")
        histsInput = [hist for hist in JPsiDir.GetListOfKeys() if "h_diel_z" in hist.GetName() and "SemiInclusive" in hist.GetName() and "pT" in hist.GetName()]
    for histpT in histsInput:
        print("Found h_diel_z histogram: ", histpT.GetName())
        input_name = JPsiDirName + "/" + histpT.GetName()
        InputHist = FileIn.Get(input_name)
        lastBinZ=InputHist.GetNbinsY()
        Nbinsz = lastBinZ
        Test = False
        if Test:
            Nbinsz = 5
            if "Energy" in input_name:
                print("Test = True, then skipping Energy analysis")
                continue
        BinResultsDicts_inv = [] # List of z_size dictionaries, each one with 3 dicts (1 for each fit range)
        # FileIn.Close() #Maybe this should not be commented. previously this was after "for binZ..."
        #for binZ in range(1, Nbinsz + 1):
        for binZ in range(lastBinZ, lastBinZ - Nbinsz, -1): # from higher bins to lower
            print("Running fit for bin index %d" % binZ)
            MinDataRange = inputCfg["input"]["pdf_dictionary"]["fitRangeMin"][0] #the range of the data is the greater range of the fit (the larger has to be the first?)
            MaxDataRange = inputCfg["input"]["pdf_dictionary"]["fitRangeMax"][0]
            dqFitter = DQFitter(
                inputCfg["input"]["input_file_name"], input_name, inputCfg["output"]["output_file_name"], MinDataRange, MaxDataRange, inputCfg["input"]["desired_bins"], True, binZ
            )
            if binZ != lastBinZ: # Fit-parameters initial values given by the result of previous fit (except for the first bin)
                print("binZ :", binZ, " != lastBinZ-1:", lastBinZ-1)
                NparPDFs = [len(parNamePDF) for parNamePDF in inputCfg["input"]["pdf_dictionary"]["parName"]]
                print("NparPDF:", NparPDFs)
                for PDF_i, parNamesPDF in enumerate(inputCfg["input"]["pdf_dictionary"]["parName"]):
                    for PDFpar_i, parNamePdf in enumerate(parNamesPDF):
                            inputCfg["input"]["pdf_dictionary"]["parVal"][PDF_i][PDFpar_i] = BinResultsDicts_inv[-1]["DictRange_1.8 to 4.2 GeV"][parNamePdf]
            dqFitter.SetFitConfig(inputCfg["input"]["pdf_dictionary"])
            BinResultsDicts_inv.append(dqFitter.MultiTrial()) # Fitting happens here
        # print(f"BinResultsDicts_inv:{BinResultsDicts_inv} \n\n\n\n")
        BinResultsDicts = BinResultsDicts_inv[::-1] # Reverse the list to have the bins in ascending order
        FileOutName = "{}{}FitterResults.root".format(inputCfg["output"]["output_file_name"], inputCfg["input"]["input_file_name"].rsplit("/", 1)[-1][9:-5])
        FileOut = TFile(FileOutName, "UPDATE") # Has to be openned again due to the fact that it was closed in dqFitter.
        FileOut.cd()
        if "pT" in input_name:
            print("Found pT in histogram name")
            DirName = "Analysis_in_pT_Ranges"
        if "Energy" in input_name:
            print("Found Energy in histogram name")
            DirName = "Analysis_in_Energy_Ranges"
        Dir = FileOut.GetDirectory(DirName) or FileOut.mkdir(DirName)
        Dir.cd()
        thisRangeDir = Dir.GetDirectory(input_name.split("/")[-1]) or Dir.mkdir(input_name.split("/")[-1])
        thisRangeDir.cd() #maybe exit the directory later
        ''' Energy analysis
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
        InitialBin = 2 # Exclude low-stat bins. "=0" to get the whole x range
        NJPsi_AllFitRanges = []
        NJPsiErr_AllFitRanges = []
        hists_z = []
        hists_z_titles = []
        pTBounds = re.findall(r'\d+', input_name) # list of sequences of digits
        legLine3 = pTBounds[0][:-3] + " < p_{T}^{Jet} < " + pTBounds[1][:-3] + " GeV/c, R = 0.4"
        legLines[2] = legLine3
        pTRangeTitle = pTBounds[0][:-3] + " <  p_{T}^{jet} < " + pTBounds[1][:-3] + " GeV"
        hist_z_titlepT = "J/#Psi z, " + pTRangeTitle
        NJPsiTots = []
        for RangeName in BinResultsDicts[0]: #Iterate through "fit" ranges (0 is arbitrary). Assuming all z bins use same ranges for mass fit!
            hist_z_titleRange =  " (fit range " + RangeName[10:] + ")"
            hist_z_title = hist_z_titlepT + hist_z_titleRange
            hist_z = TH1F(f"h_z_{RangeName}", "", Nbinsz-InitialBin, (InitialBin/Nbinsz), 1)
            zBinHist = 1
            NJPsi_ThisFitRange = []
            NJPsiErr_ThisFitRange = []
            for zBinDict in range(InitialBin, Nbinsz):
                NJPsi = BinResultsDicts[zBinDict][RangeName]["sig_Jpsi"]
                NJPsiErr = BinResultsDicts[zBinDict][RangeName]["sig_Jpsi_err"]
                NJPsi_ThisFitRange.append(NJPsi)
                NJPsiErr_ThisFitRange.append(NJPsiErr)
                hist_z.SetBinContent(zBinHist, NJPsi)
                hist_z.SetBinError(zBinHist, NJPsiErr)
                zBinHist += 1
            NJPsi_AllFitRanges.append(NJPsi_ThisFitRange)
            NJPsiErr_AllFitRanges.append(NJPsiErr_ThisFitRange)

            # Histogram style and propaganda canvas
            hist_z.SetTitleSize(0.01, "t")
            hist_z.SetMarkerStyle(markerStyle) 
            hist_z.SetMarkerSize(1.0)
            hist_z.SetMarkerColor(Colors[0])
            hist_z.SetLineColor(Colors[0])
            hist_z.SetLineWidth(2)
            hist_z.SetStats(0)
            hist_z.GetYaxis().SetTitle("#frac{1}{N^{J/#Psi}_{jets}} #frac{dN}{dz^{ch.}}")
            hist_z.GetXaxis().SetTitle("z^{ch.}_{J/#Psi}")
            # hist_z.GetYaxis().SetLabelSize(0.03)
            # hist_z.GetXaxis().SetLabelSize(0.03)
            # hist_z.GetYaxis().SetTitleSize(0.03)
            # hist_z.GetXaxis().SetTitleSize(0.03)
            hist_z.Write()
            cz = TCanvas(f"c{RangeName}", f"Canvas{RangeName}", 1400, 1080)
            hist_z.Draw()
            NJPsiTot = hist_z.Integral()
            NJPsiTots.append(NJPsiTot)
            legend_z = createLegend(legLines)
            legend_z.AddEntry(hist_z, "N_{{J/#Psi}} = {:.2f}".format(NJPsiTot), "")
            legend_z.Draw()
            cztitle = TPaveText(0.1, 0.96, 0.9, 0.99, "NDC")
            cztitle = createPaveText(hist_z_title)
            cztitle.Draw()
            cz.Update()
            cz.Write()
            czNormalized = TCanvas(f"c{RangeName}norm", f"Canvas{RangeName} Normalized", 1400, 1080)
            hist_z_norm = hist_z.Clone()
            binWidth = hist_z.GetBinWidth(1) #assuming all z bins have the same width!
            hist_z_norm.Scale(1. / (NJPsiTot*binWidth))
            hist_z_norm.Draw("P")
            legend_z.Draw()
            cztitle.Draw()
            czNormalized.Update()
            czNormalized.Write()
            hists_z.append(hist_z)
            hists_z_titles.append(hist_z_titleRange)
        czsNormalized = TCanvas(hist_z_titlepT[7:] + ", Fit Ranges", hist_z_titlepT[7:] + ", Fit Ranges", 1400, 1080)
        legend_zs = createLegend(legLines)
        for range_i, hist_z in enumerate(hists_z):
            hist_z.SetMarkerStyle(markerStyle)
            hist_z.SetMarkerColor(Colors[range_i])
            hist_z.SetLineColor(Colors[range_i])
            if range_i == 0:
                hist_z.Draw("P")
            else:
                hist_z.Draw("P same")
            legend_zs.AddEntry(hist_z, hists_z_titles[range_i], "pe")
        legend_zs.Draw()
        czstitle = TPaveText(0.1, 0.96, 0.9, 0.99, "NDC")
        czstitle = createPaveText(hist_z_titlepT + ", Fit Ranges")
        czstitle.Draw()
        czsNormalized.Update()
        czsNormalized.Write()

        # Final z distribution
        hist_z_mean = TH1F(f"h_z_mean", "", Nbinsz-InitialBin, (InitialBin/Nbinsz), 1)
        # Uncertainties
        hist_z_stat = hist_z_mean.Clone("hist_z_stat")
        hist_z_syst = hist_z_mean.Clone("hist_z_syst")
        hist_z_tot = hist_z_mean.Clone("hist_z_tot")
        for h in [hist_z_stat, hist_z_syst, hist_z_tot]:
            h.SetTitle("")
        NJPsiMeans = np.mean(NJPsi_AllFitRanges, axis = 0)
        NJPsiErrMeans = np.mean(NJPsiErr_AllFitRanges, axis = 0) # Errors from the fitting (stat unc)
        NJPsiMeans_StdDev = np.std(NJPsi_AllFitRanges, ddof = 1, axis=0)
        NJPsiMeans_StdErr = NJPsiMeans_StdDev/np.sqrt(len(NJPsi_AllFitRanges)) # Errors from varying fit ranges (sys unc)
        NJPsiMeans_TotUnc = np.sqrt(NJPsiErrMeans**2 + NJPsiMeans_StdErr**2) # Total Uncertainty
        for iBin in range(len(NJPsiMeans)):
            hist_z_mean.SetBinContent(iBin + 1, NJPsiMeans[iBin])
            hist_z_mean.SetBinError(iBin + 1, NJPsiMeans_TotUnc[iBin])
            hist_z_stat.SetBinContent(iBin + 1, NJPsiErrMeans[iBin]*100/NJPsiMeans[iBin])
            hist_z_syst.SetBinContent(iBin + 1, NJPsiMeans_StdErr[iBin]*100/NJPsiMeans[iBin])
            hist_z_tot.SetBinContent(iBin + 1, NJPsiMeans_TotUnc[iBin]*100/NJPsiMeans[iBin])
        NJPsiTotmean = hist_z_mean.Integral()

        # Histogram style and propaganda canvas
        hist_z_mean.SetTitleSize(0.01, "t")
        hist_z_mean.SetMarkerStyle(markerStyle) 
        hist_z_mean.SetMarkerSize(1.0)  
        hist_z_mean.SetMarkerColor(Colors[0])
        hist_z_mean.SetLineColor(Colors[0])
        hist_z_mean.SetLineWidth(2)
        hist_z_mean.SetStats(0)
        hist_z_mean.GetYaxis().SetTitle("#frac{1}{N^{J/#Psi}_{jets}} #frac{dN}{dz^{ch.}}")
        hist_z_mean.GetXaxis().SetTitle("z^{ch.}_{J/#Psi}")
        cz_mean = TCanvas(hist_z_titlepT[7:] + ", Mean", hist_z_titlepT[7:] + ", Mean", 1400, 1080)
        hist_z_mean.Draw()

        hist_z_stat.SetTitleSize(0.01, "t")
        hist_z_stat.SetLineColor(Colors[0])
        hist_z_stat.SetMarkerStyle(markerStyle)
        hist_z_stat.SetLineWidth(2)
        hist_z_stat.SetStats(0)
        hist_z_stat.GetYaxis().SetTitle("Relative Unc. (%)")
        hist_z_stat.GetXaxis().SetTitle("z^{ch.}_{J/#Psi}")
        hist_z_syst.SetLineColor(Colors[1])
        hist_z_syst.SetMarkerStyle(markerStyle)
        hist_z_syst.SetLineWidth(2)
        hist_z_syst.SetStats(0)
        hist_z_tot.SetLineColor(Colors[2])
        hist_z_tot.SetMarkerStyle(markerStyle)
        hist_z_tot.SetLineWidth(2)
        hist_z_tot.SetLineStyle(2)
        hist_z_tot.SetStats(0)
        legendUnc = createLegend(legLines)
        legendUnc.AddEntry(hist_z_stat, "Stat Unc. (from the fitting)", "l")
        legendUnc.AddEntry(hist_z_syst, "Syst Unc. (from varying mass-fit range)", "l")
        legendUnc.AddEntry(hist_z_tot, "Tot Unc.", "l")
        legendMean = createLegend(legLines)
        legendMean.AddEntry(hist_z_mean, "N_{{J/#Psi}} = {:.2f}".format(NJPsiTotmean), "")
        legendMean.Draw()
        cz_meantitle = TPaveText(0.1, 0.96, 0.9, 0.99, "NDC")
        cz_meantitle = createPaveText(hist_z_titlepT)
        cz_meantitle.Draw()
        cz_mean.Update()
        cz_mean.Write()
        cz_meanNormalized = TCanvas(hist_z_titlepT[7:] + ", Mean Norm", hist_z_titlepT[7:] + ", Mean Normalized", 1400, 1080)
        hist_z_mean_norm = hist_z_mean.Clone()
        binWidth = hist_z_mean.GetBinWidth(1) #assuming all z bins have the same width!
        hist_z_mean_norm.Scale(1. / (NJPsiTotmean*binWidth)) #Normalize to 1
        hist_z_mean_norm.Draw("P")
        hist_z_mean_norm.SetDirectory(0)
        hists_z_mean_norm.append([hist_z_mean_norm, hist_z_titlepT[7:], NJPsiTotmean])
        legendMean.Draw()
        cz_meantitle.Draw()
        cz_meanNormalized.Update()
        cz_meanNormalized.Write()

        cz_meanUncs = TCanvas(hist_z_titlepT[7:] + ", Uncs", hist_z_titlepT[7:] + ", Uncertainties", 1400, 1080)
        cz_meanUncstitle = TPaveText(0.1, 0.96, 0.9, 0.99, "NDC")
        cz_meanUncstitle = createPaveText(f"Relative Uncertainties for z, p_{{T}}^{{jet}} {pTRangeTitle}")

        hist_z_stat.GetYaxis().SetRangeUser(0, 100.0)
        hist_z_stat.Draw("HIST")
        hist_z_syst.Draw("HIST SAME")
        hist_z_tot.Draw("HIST SAME")
        legendUnc.Draw()

        cz_meanUncstitle.Draw()
        cz_meanUncs.Update()
        cz_meanUncs.Write()

        # Canvas for each parameter. One curve for each range 
        parCanvases = []
        parDir = thisRangeDir.GetDirectory("Fit Parameters") or thisRangeDir.mkdir("Fit Parameters")
        parDir.cd()
        for i in range(len(inputCfg["input"]["pdf_dictionary"]["parName"])):
            for j, parName in enumerate(inputCfg["input"]["pdf_dictionary"]["parName"][i]):
                c = TCanvas(f"canvas_{parName}", f"{parName}", 1000, 800)
                c.cd()
                graphs = []
                legendPar = createLegend(legLines)
                legendPar.SetHeader("Fit Mass Ranges")
                graphsY = [] #Y values for all (3) ranges. Important for calculation of mean
                for r, RangeName in enumerate(BinResultsDicts[0]): # "DictRange_1.8 to 4.2 GeV", etc.
                    graphY = [] #Y values for this range
                    for zBinDict in range(InitialBin, Nbinsz):
                        graphY.append(BinResultsDicts[zBinDict][RangeName][parName]) #Value of the parameter for each z bin
                    parGraph = TGraph(Nbinsz-InitialBin, array("f", np.arange(InitialBin*binWidth, Nbinsz*binWidth, binWidth)), array("f", graphY))
                    # parGraph.SetTitle(f"{RangeName[10:]}")
                    parGraph.SetMarkerStyle(markerStyle)  # Different marker styles (20, 21, 22...)
                    parGraph.SetMarkerColor(Colors[r])
                    parGraph.SetMarkerSize(1.0)
                    parGraph.SetLineColor(Colors[r])
                    # parGraph.SetLineWidth(2)
                    legendPar.AddEntry(parGraph, f"{RangeName[10:]}", "p")
                    graphs.append(parGraph)
                    graphsY.append(graphY)
                graphs[0].SetTitle(f";z Bin (lower bound);{parName}") # Axis titles
                graphs[0].Draw("APL")
                parTitle = createPaveText(f"{parName} values from fit")
                parTitle.Draw()
                graphsY = np.array(graphsY)
                meanGraphY = [] #Graph for the mean of the parameter for different ranges
                for k, zBinDict in enumerate(range(InitialBin, Nbinsz)):
                    graphsY_column = graphsY[:, k]
                    mean = np.mean(graphsY_column)
                    meanGraphY.append(mean)
                meanGraph = TGraph(Nbinsz-InitialBin, array("f", np.arange(InitialBin*binWidth, Nbinsz*binWidth, binWidth)), array("f", meanGraphY))
                meanGraph.SetMarkerStyle(34)
                legendPar.AddEntry(meanGraph, "Mean", "p")
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
                legendPar.AddEntry(lineMax, f"Max Value", "l")
                legendPar.AddEntry(lineIni, f"Initial Value", "l")
                legendPar.AddEntry(lineMin, f"Min Value", "l")
                legendPar.Draw()
                for graph in graphs[1:]:
                    graph.Draw("PL SAME")
                # c.Update()
                parCanvases.append(c)
                c.Write()

        # Plot function with parameters from json (without fitting)
        canvasRawPDFs = TCanvas("canvasRawPDFs", "canvasRawPDFs", 1400, 1080)
        fitSumLeg = createLegend(legLines)
        funcs = []
        for pdf_i, pdf_name in enumerate(inputCfg["input"]["pdf_dictionary"]["pdf"][:-1]): # excludes the "SUM"
            fitRangeMin = inputCfg["input"]["pdf_dictionary"]["fitRangeMin"][0] # The first range has to be the largest, as in DQFitter
            fitRangeMax = inputCfg["input"]["pdf_dictionary"]["fitRangeMax"][0]
            n_pars = len(inputCfg["input"]["pdf_dictionary"]["parName"][pdf_i]) 
            f = TF1(pdf_name, PDF(pdf_name), fitRangeMin, fitRangeMax, n_pars)
            parValIni = inputCfg["input"]["pdf_dictionary"]["parVal"][pdf_i]
            f.SetParameters(*parValIni)
            f.SetParNames(*inputCfg["input"]["pdf_dictionary"]["parName"][pdf_i])
            f.SetLineColor(Colors[pdf_i])
            f.SetLineWidth(2)
            f.SetTitle(f"PDFs With Initial-Values Parameters;x;Probability Density")  
            f.SetMaximum(2.0)
            print("pdf_i: ", pdf_i)
            if pdf_i == 0:
                f.Draw()
            else:
                f.Draw("same")
            # f.Write()
            fitSumLeg.AddEntry(f, pdf_name, "l")
            funcs.append(f)
        # Plot the sum of all pdfs
        fitSumMin = inputCfg["input"]["pdf_dictionary"]["fitRangeMin"][0]
        fitSumMax = inputCfg["input"]["pdf_dictionary"]["fitRangeMax"][0]
        npoints = 500
        xs = [fitSumMin + i*(fitSumMax-fitSumMin)/(npoints-1) for i in range(npoints)]
        ys = []
        for x in xs:
            total = sum(f.Eval(x) for f in funcs)   # funcs contém [CB2, VWG, ...]
            ys.append(total)
        xa = np.array(xs, dtype=np.float64)
        ya = np.array(ys, dtype=np.float64)
        fitSum = TGraph(npoints, xa, ya)
        fitSum.SetLineColor(Colors[len(funcs)])
        fitSum.SetLineWidth(2)
        fitSum.SetTitle("PDFs Sum;x;PDF")
        fitSum.Draw("L SAME")
        fitSumLeg.AddEntry(fitSum, "PDFs Sum", "l")
        fitSumLeg.Draw()
        canvasRawPDFs.Update()
        canvasRawPDFs.Write()
        FileOut.Close()
        print("Fit for multiple x projections: Done!")
    # All z distributions for different pT ranges
    FileOut = TFile(FileOutName, "UPDATE")
    FileOut.cd()
    Dir = FileOut.GetDirectory(DirName)
    Dir.cd()
    czspTsNormalized = TCanvas("JPsi z, All pT Ranges", "JPsi z, All pT Ranges", 1400, 1080)
    legend_zspTs = createLegend(legLines[:2])
    hists_z_mean_norm = np.array(hists_z_mean_norm, dtype=object)
    hists_z_mean_norm_hs = THStack("hists_z_mean_norm", "")
    skipSemiInlclusive = True
    hist_z_mean_norm_ini = 0
    if skipSemiInlclusive == True:
        hist_z_mean_norm_ini = 1
    for range_pT_i, hist_z in enumerate(hists_z_mean_norm[hist_z_mean_norm_ini:, 0], hist_z_mean_norm_ini): # starts from hist_z_mean_norm_ini
        hist_z.SetMarkerStyle(markerStyle)
        hist_z.SetMarkerColor(ColorsExt[range_pT_i])
        hist_z.SetLineColor(ColorsExt[range_pT_i])
        hists_z_mean_norm_hs.Add(hist_z)
        hist_z_leg = f"{hists_z_mean_norm[range_pT_i, 1]}, N_{{J/#Psi}} = {hists_z_mean_norm[range_pT_i, 2]:.2f}"
        legend_zspTs.AddEntry(hist_z, hist_z_leg, "pe")
    hists_z_mean_norm_hs.Draw("nostack")
    hists_z_mean_norm_hs.GetXaxis().SetTitle(hists_z_mean_norm[0,0].GetXaxis().GetTitle())
    hists_z_mean_norm_hs.GetYaxis().SetTitle(hists_z_mean_norm[0,0].GetYaxis().GetTitle())
    hists_z_mean_norm_hs.SetMinimum(-0.5)
    legend_zspTs.Draw()
    czspTstitle = TPaveText(0.1, 0.96, 0.9, 0.99, "NDC")
    czspTstitle = createPaveText("JPsi z, All pT Ranges")
    czspTstitle.Draw()
    czspTsNormalized.Update()
    czspTsNormalized.Write()
    print("Output created: ", FileOutName)

main()
