"""
\author Luca Micheletti <luca.micheletti@cern.ch>, CERN
\author Victor Valencia <valencia@subatech.in2p3.fr>, subatech
"""
from os.path import exists
from plot_library import DoResidualPlot, DoPullPlot, DoCorrMatPlot, DoPropagandaPlot

import ROOT
from ROOT import (
    TH1F,
    RooArgSet,
    RooDataHist,
    RooDataSet,
    RooRealVar,
    RooWorkspace,
    TCanvas,
    TFile,
    gPad,
    gROOT,
)


class DQFitter:
    def __init__(self, fInName, fInputName, fOutPath, fMinDataRange, fMaxDataRange, desired_Nbins=700, run_projections=False, binY=None):
        if not exists(fInName):
            print("The input file does not exist, exit...")
            exit()
        self.fPdfDict = {}
        self.fOutPath = fOutPath
        # self.fFileOut = TFile("{}{}.root".format(fOutPath, fInputName.rsplit("/", 1)[-1]), "UPDATE")
        self.fInputName = fInputName
        self.fFileIn = TFile.Open(fInName)
        self.fInput = self.fFileIn.Get(fInputName)
        print("self.fInput.ClassName() = ", self.fInput.ClassName())
        print("self.fInput = ", self.fInput)
        if binY is None or binY == 1: # if None, not analysis in bins. If 1st bin, the file has to be created
            self.fFileOut = TFile("{}{}.root".format(fOutPath, fInputName.rsplit("/", 1)[-1]), "RECREATE")
            self.fFileOut.cd()
            self.fInput.Write()
            h_inclusive_mass = self.fFileIn.Get("j-psi-fragmentation-function-task/h_diel_mass")
            h_inclusive_mass.Write()
        else:
            self.fFileOut = TFile("{}{}.root".format(fOutPath, fInputName.rsplit("/", 1)[-1]), "UPDATE") 
        self.fFileOut.ls()
        self.run_projections = run_projections
        if run_projections:
            print("Projecting 2D histogram to 1D histograms")
            print(f"Projecting bin {binY}")
            y_max = self.fInput.GetYaxis().GetXmax()
            y_min = self.fInput.GetYaxis().GetXmin()
            NBinsY = self.fInput.GetNbinsY()
            binYValueMin = y_min + ((binY-1)*(y_max - y_min) / NBinsY)
            binYValueMax = y_min + ((binY)*(y_max - y_min) / NBinsY)
            self.subDirName = "bin_" + str(binYValueMin) + "_to_" + str(binYValueMax)
            print(f"Subdirectory name: {self.subDirName}")
            projected_hist = self.fInput.ProjectionX("proj_x", binY, binY)
            original_bins = projected_hist.GetNbinsX()
            rebin_factor = original_bins / desired_Nbins
            rebinned_hist = projected_hist.Rebin(int(rebin_factor), f"rebinned_proj_x{binY}")
            print(f"Rebinned bins: {rebinned_hist.GetNbinsX()}")
            # self.fFileOut.cd()
            # rebinned_hist.Write()
            self.fInput = rebinned_hist
        self.fRooWorkspace = RooWorkspace("w", "workspace")
        self.fParNames = []
        self.fFitRangeMin = []
        self.fFitRangeMax = []
        self.fTrialName = ""
        self.fRooMass = RooRealVar("m", "#it{M} (GeV/#it{c}^{2})", fMinDataRange, fMaxDataRange)

    def SetFitConfig(self, pdfDict):
        """
        Method to set the fit PDFs
        """
        self.fPdfDict = pdfDict
        self.fFitRangeMin = pdfDict["fitRangeMin"]
        self.fFitRangeMax = pdfDict["fitRangeMax"]

        pdfList = []
        for pdf in self.fPdfDict["pdf"][:-1]:
            self.fTrialName = self.fTrialName + pdf + "_"

        for i in range(0, len(self.fPdfDict["pdf"])):
            if not self.fPdfDict["pdf"][i] == "SUM":
                gROOT.ProcessLineSync(
                    ".x fit_library/{}Pdf.cxx+".format(self.fPdfDict["pdf"][i])
                )

        for i in range(0, len(self.fPdfDict["pdf"])):
            parVal = self.fPdfDict["parVal"][i]
            parLimMin = self.fPdfDict["parLimMin"][i]
            parLimMax = self.fPdfDict["parLimMax"][i]
            parName = self.fPdfDict["parName"][i]

            if not self.fPdfDict["pdf"][i] == "SUM":
                for j in range(0, len(parVal)):
                    if ("sum" in parName[j]) or ("prod" in parName[j]):
                        self.fRooWorkspace.factory("{}".format(parName[j]))
                        r1 = parName[j].find("::") + 2
                        r2 = parName[j].find("(", r1)
                        parName[j] = parName[j][r1:r2]
                        if (parLimMin == parLimMax):
                            self.fRooWorkspace.factory("{}[{}]".format(parName[j], parVal[j]))
                    else:
                        if (parLimMin == parLimMax):
                            self.fRooWorkspace.factory("{}[{}]".format(parName[j], parVal[j]))
                        else:
                            self.fRooWorkspace.factory("{}[{},{},{}]".format(parName[j], parVal[j], parLimMin[j], parLimMax[j]))
                        self.fParNames.append(parName[j])
                nameFunc = self.fPdfDict["pdf"][i]
                nameFunc += "Pdf::{}Pdf(m[{},{}]".format(self.fPdfDict["pdfName"][i],self.fPdfDict["fitRangeMin"][0],self.fPdfDict["fitRangeMax"][0])
                pdfList.append(self.fPdfDict["pdfName"][i])
                for j in range(0, len(parVal)):
                    nameFunc += ",{}".format(parName[j])
                nameFunc += ")"
                self.fRooWorkspace.factory(nameFunc)
            else:
                nameFunc = self.fPdfDict["pdf"][i]
                nameFunc += "::sum("
                for j in range(0, len(pdfList)):
                    nameFunc += "{}[{},{},{}]*{}Pdf".format(parName[j], parVal[j], parLimMin[j], parLimMax[j], pdfList[j])
                    self.fParNames.append(parName[j])
                    if not j == len(pdfList) - 1:
                        nameFunc += ","
                nameFunc += ")"
                self.fRooWorkspace.factory(nameFunc)



    def FitInvMassSpectrum(self, fitRangeMin, fitRangeMax, run_projections=False, subDirName=None, DictRange=None):
        print("TEST1")
        """
        Method to perform binned / unbinned fit to a ROOT histogram / tree
        """
        trialName = self.fTrialName + "_" + str(fitRangeMin) + "_" + str(fitRangeMax)
        binPlotsName = self.subDirName + "_" + trialName
        self.fRooWorkspace.Print()
        pdf = self.fRooWorkspace.pdf("sum")
        self.fRooMass.setRange("range", fitRangeMin, fitRangeMax)
        print("fitRangeMin = ", fitRangeMin)
        print("fitRangeMax = ", fitRangeMax)
        if "TTree" in self.fInput.ClassName():
            print("Perform unbinned fit")
            rooDs = RooDataSet(
                "data",
                "data",
                RooArgSet(self.fRooMass),
                ROOT.RooFit.Import(self.fInput),
            )
        else:
            print("Perform binned fit")
            # trimmed_hist = self.fInput.Clone("trimmed_hist")
            # trimmed_hist.GetXaxis().SetRangeUser(fitRangeMin, fitRangeMax)
            # trimmed_hist = trimmed_hist.Rebin(
            #     trimmed_hist.GetNbinsX(),
            #     "rebinned_hist",
            #     trimmed_hist.GetXaxis().GetXbins().GetArray())
            rooDs = RooDataHist(
                "data",
                "data",
                RooArgSet(self.fRooMass),
                ROOT.RooFit.Import(self.fInput),
                # ROOT.RooFit.Import(trimmed_hist),
                # ROOT.RooFit.Range("range")
            )
        rooFitRes = ROOT.RooFitResult(pdf.fitTo(rooDs, ROOT.RooFit.Save(True)))
        fRooPlot = self.fRooMass.frame(ROOT.RooFit.Title(binPlotsName))
        fRooPlotCopy = self.fRooMass.frame(ROOT.RooFit.Title(binPlotsName))
        print("fitRangeMinPLOT = ", fitRangeMin)
        print("fitRangeMaxPLOT = ", fitRangeMax)
        rooDs.plotOn(fRooPlot, ROOT.RooFit.MarkerStyle(20), ROOT.RooFit.MarkerSize(0.6), ROOT.RooFit.Range(fitRangeMin, fitRangeMax))
        pdf.plotOn(fRooPlot, ROOT.RooFit.LineColor(ROOT.kRed+1), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Range(fitRangeMin, fitRangeMax))

        index = 1
        histResults = TH1F(
            "fit_results_{}".format(binPlotsName),
            "fit_results_{}".format(binPlotsName),
            len(self.fParNames),
            0.0,
            len(self.fParNames),
        )
        for parName in self.fParNames:
            histResults.GetXaxis().SetBinLabel(index, parName)
            histResults.SetBinContent(index, self.fRooWorkspace.var(parName).getVal())
            print(f"self.fRooWorkspace.var({parName}).getVal(), subdirname {subDirName} = ", self.fRooWorkspace.var(parName).getVal())
            DictRange[parName] = (self.fRooWorkspace.var(parName).getVal())
            DictRange[parName+"_err"] = (self.fRooWorkspace.var(parName).getError())
            histResults.SetBinContent(index, self.fRooWorkspace.var(parName).getError())
            index += 1

        for i in range(0, len(self.fPdfDict["pdf"])):
            print("TEST2")
            print('self.fPdfDict["pdfName"][i] = ', self.fPdfDict["pdfName"][i])
            if not self.fPdfDict["pdfName"][i] == "SUM":
                pdf.plotOn(fRooPlot, ROOT.RooFit.Components("{}Pdf".format(self.fPdfDict["pdfName"][i])), ROOT.RooFit.LineColor(self.fPdfDict["pdfColor"][i]), ROOT.RooFit.LineStyle(self.fPdfDict["pdfStyle"][i]), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Range(fitRangeMin, fitRangeMax))

        extraText = []
        paveText = ROOT.TPaveText(0.85, 0.45, 0.99, 0.94, "brNDC")
        paveText.SetTextFont(42)
        paveText.SetTextSize(0.015)
        paveText.SetFillColor(ROOT.kWhite)
        for parName in self.fParNames:
            print("TEST3")
            paveText.AddText("{} = {:.4f} #pm {:.4f}".format(parName, self.fRooWorkspace.var(parName).getVal(), self.fRooWorkspace.var(parName).getError()))
            if self.fPdfDict["parForPropagandaPlot"].count(parName) > 0:
                text = self.fPdfDict["parNameForPropagandaPlot"][self.fPdfDict["parForPropagandaPlot"].index(parName)]
                if "sig" in parName:
                    extraText.append("{} = {:.0f} #pm {:.0f}".format(text, self.fRooWorkspace.var(parName).getVal(), self.fRooWorkspace.var(parName).getError()))
                else:
                    extraText.append("{} = {:.3f} #pm {:.3f}".format(text, self.fRooWorkspace.var(parName).getVal(), self.fRooWorkspace.var(parName).getError()))
            for i in range(0, len(self.fPdfDict["pdfName"])):
                if self.fPdfDict["pdfName"][i] in parName:
                    (paveText.GetListOfLines().Last()).SetTextColor(self.fPdfDict["pdfColor"][i])

        nPars = rooFitRes.floatParsFinal().getSize()
        if "TTree" in self.fInput.ClassName():
            print("TEST4")
            # Convert RooDataSet into RooDataHist to extract the Chi2 value
            rooDh = RooDataHist("rooDh", "binned version of rooDs", RooArgSet(self.fRooMass), rooDs)
            chi2 = ROOT.RooChi2Var("chi2", "chi2", pdf, rooDh)
            #chi2 = ROOT.RooChi2Var("chi2", "chi2", pdf, rooDh, False, 1)
            nbinsperGev = rooDh.numEntries() / (self.fPdfDict["fitRangeMax"][0] - self.fPdfDict["fitRangeMin"][0])
            nBins = (fitRangeMax - fitRangeMin) * nbinsperGev
            ndof = nBins - nPars
            reduced_chi2 = chi2.getVal() / ndof
            paveText.AddText("#bf{#chi^{2}/dof = %3.2f}" % (reduced_chi2))
            extraText.append("#chi^{2}/dof = %3.2f" % reduced_chi2)
        else:
            print("TEST5")
            # To Do : Find a way to get the number of bins differently. The following is a temparary solution.
            # WARNING : The largest fit range has to come first in the config file otherwise it does not work
            chi2 = ROOT.RooChi2Var("chi2", "chi2", pdf, rooDs, False, 1)
            # chi2 = ROOT.RooChi2Var("chi2", "chi2", pdf, rooDs)
            nbinsperGev = rooDs.numEntries() / (self.fPdfDict["fitRangeMax"][0] - self.fPdfDict["fitRangeMin"][0])
            nBins = (fitRangeMax - fitRangeMin) * nbinsperGev
            ndof = nBins - nPars
            reduced_chi2 = chi2.getVal() / ndof
            #reduced_chi2 = 3 #chi2.getVal() / ndof
            paveText.AddText("n Par = %3.2f" % (nPars))
            paveText.AddText("n Bins = %3.2f" % (nBins))
            paveText.AddText("#bf{#chi^{2}/dof = %3.2f}" % reduced_chi2)
            fRooPlot.addObject(paveText)
            extraText.append("#chi^{2}/dof = %3.2f" % reduced_chi2)
            print("TEST6")
        
        # Fit plot
        canvasFit = TCanvas(
            "fit_plot_{}".format(binPlotsName), "fit_plot_{}".format(binPlotsName), 1400, 1080
        )
        canvasFit.SetLeftMargin(0.15)
        gPad.SetLeftMargin(0.15)
        fRooPlot.GetYaxis().SetTitleOffset(1.4)
        fRooPlot.Draw()

        # Residual plot
        print("TEST7")
        canvasResidual = DoResidualPlot(fRooPlot, self.fRooMass, binPlotsName)

        # Pull plot
        print("TEST8")
        canvasPull = DoPullPlot(fRooPlot, self.fRooMass, binPlotsName)

        # Correlation matrix plot
        print("TEST9")
        canvasCorrMat = DoCorrMatPlot(rooFitRes, binPlotsName)

        # Propaganda plot
        print("TEST10")
        if self.fPdfDict["doPropagandaPlot"]:
            DoPropagandaPlot(rooDs, pdf, fRooPlotCopy, self.fPdfDict, self.fInputName, binPlotsName, self.fOutPath, extraText)

        # Save results
        self.fFileOut.cd()
        if run_projections:
            if not self.fFileOut.GetDirectory(self.subDirName):
                self.fFileOut.mkdir(self.subDirName) # create 1 subdir for each z bin
            self.fFileOut.cd(self.subDirName)
        print("self.fFileOut = ", self.fFileOut)
        self.fFileOut.ls()
        canvasFit.Write()
        canvasResidual.Write()
        canvasPull.Write()
        canvasCorrMat.Write()
        histResults.Write()

    def MultiTrial(self):
        """
        Method to run multiple fits of the same invariant mass distribution
        """
        # ResultsDict=dict.fromkeys(self.fParNames) 
        # Resultsbin = {}
        BinResultsDict = {} # Dictionary that will contain the dictionaries for each pT range
        for iRange in range(0, len(self.fFitRangeMin)):
            print("iRangemin = ",self.fFitRangeMin[iRange])
            DictRange = {} # Dictionary for each pT range
            self.FitInvMassSpectrum(
                self.fFitRangeMin[iRange], self.fFitRangeMax[iRange], self.run_projections, self.subDirName, DictRange
            )
            RangeName = f"DictRange_{self.fFitRangeMin[iRange]}_to_{self.fFitRangeMax[iRange]}GeV"
            BinResultsDict[RangeName] = DictRange
            print("DictRange = ", DictRange)
        for i in BinResultsDict:
            print(f"BinResultsDict[{i}] = \n", BinResultsDict[i])
        self.fFileOut.Close()
        return BinResultsDict
