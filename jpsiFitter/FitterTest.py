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

File = TFile.Open("/home/ferrandi/alice/O2Physics/PWGJE/Tasks/JPsiWorkDir/ResultsLHC24am_pass1_DIELC_R4_4/AnalysisCompleteSubJobsLHC24am_apass1_25-03-21.root")
InputHist = File.Get("j-psi-fragmentation-function-task/h_diel_mass")

# fMinDataRange = 2
# fMaxDataRange = 5
fMinDataRange = 2
fMaxDataRange = 5
fFitRangeMin = 2.8
fFitRangeMax = 3.2
fRooMass = RooRealVar("m", "#it{M} (GeV/#it{c}^{2})", fMinDataRange, fMaxDataRange)
fRooMass.setRange("fit_range", fFitRangeMin, fFitRangeMax)
trimmed_hist = InputHist.Clone("trimmed_hist")
trimmed_hist.GetXaxis().SetRangeUser(fFitRangeMin, fFitRangeMax)
rooDs = RooDataHist(
                "data",
                "data",
                RooArgSet(fRooMass),
                ROOT.RooFit.Import(trimmed_hist),
                # ROOT.RooFit.Import(fRooMass),
                ROOT.RooFit.Range("fit_range")
            )


a0 = ROOT.RooRealVar("a0", "a0", 0.0, -10, 10)  # Initial value = 0, range [-10, 10]
a1 = ROOT.RooRealVar("a1", "a1", 0.0, -10, 10)
px = ROOT.RooPolynomial("px", "px", fRooMass, ROOT.RooArgList(a0, a1))
mx = ROOT.RooRealVar("mx", "mean", 3.1, 3, 3.2)
sigma = ROOT.RooRealVar("sigma", "sigma", 0.1, 0.01, 1.0)
gx = ROOT.RooGaussian("gx", "gx", fRooMass, mx, sigma)
f = ROOT.RooRealVar("f", "f", 0.0, 1.0)
model = ROOT.RooAddPdf("model", "model", [gx, px], [f])


rooFitRes = model.fitTo(rooDs, ROOT.RooFit.Range("fit_range"), ROOT.RooFit.Save(True))
fRooPlot = fRooMass.frame(Title="Fitting a sub range")
rooDs.plotOn(fRooPlot)
rooFitRes.Print()
model.plotOn(fRooPlot, ROOT.RooFit.LineColor(ROOT.kRed))
c = ROOT.TCanvas("rf203_ranges", "rf203_ranges", 1920, 1080)
fRooPlot.Draw()

c.SaveAs("ResultFitterTest.png")