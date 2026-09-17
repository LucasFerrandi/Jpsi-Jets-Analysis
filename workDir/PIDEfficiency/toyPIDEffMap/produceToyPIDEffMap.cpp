#include <cmath>
#include <iostream>
#include "TFile.h"
#include "TH2D.h"

// mesma função de eficiência
double elecPidEff(double pt, double eta)
{
    double ptTerm  = 2.0 / (1.0 + std::exp(0.15*(pt - 1.0)));
    double etaTerm = std::exp(0.5 * std::pow(0.5*eta, 2))-0.5;
    return 0.9 * ptTerm * etaTerm;
}

int main()
{
    const int nPtBins  = 20;
    const int nEtaBins = 16;

    double ptMin = 0.0, ptMax = 20.0;
    double etaMin = -1.6, etaMax = 1.6;

    // cria histo
    TH2D* hEff = new TH2D("hElecPidEff",
                         "Toy electron PID efficiency; p_{T} (GeV/c); #eta",
                         nPtBins, ptMin, ptMax,
                         nEtaBins, etaMin, etaMax);

    // preenche
    for (int i = 1; i <= nPtBins; ++i) {
        double pt = hEff->GetXaxis()->GetBinCenter(i);
        for (int j = 1; j <= nEtaBins; ++j) {
            double eta = hEff->GetYaxis()->GetBinCenter(j);
            double eff = elecPidEff(pt, eta);
            hEff->SetBinContent(i, j, eff);
        }
    }

    // salva
    TFile f("toyElectronPidEff.root", "RECREATE");
    hEff->Write();
    f.Close();

    std::cout << "Saved toyElectronPidEff.root\n";
}
