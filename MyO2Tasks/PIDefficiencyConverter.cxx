#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

// #include "Common/Core/TableHelper.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
// #include "Framework/ASoAHelpers.h"
// #include "Framework/AnalysisHelpers.h"

#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

using MySimpleEvents = aod::ReducedEvents;
using MySimpleBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelLabels>;

struct AnalysisPIDEfficiency{
    HistogramRegistry registry{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
    Configurable<bool> test{"test", false, "Limits number of events"};
    Configurable<int> Nevts{"NevtsForTest", 100, "Number of events per DF (for testing)"};
    Configurable<std::string> effMapPath{"effMapPath", "$HOME/ferrandi/alice/Jpsi-Jets-Analysis/JpsiWorkDir/PIDEfficiency/PIDEfficiencyConverter/toyPIDEffMap.root", "Path to root file with PID-efficiency map"};
    TH2D* hPMElecEffMap = nullptr;
    void init(o2::framework::InitContext const&)
    {
        registry.add("hEvtCounter", "Event Counter", kTH1F, {{2, -0.5, 1.5}});
        registry.add("ptResolution", "PM Electrons ptResolution", kTH2F, {{100, 0.f, 20.f, "p_{T} (GeV/#it{c})"}, {200, -1.f, 1.f, "#Delta(p_{T}) (GeV/#it{c})"}});
        registry.add("hPtPMElec", "PM Electrons pT", kTH1F, {{200, 0, 20, "p_{T} (GeV/#it{c})"}});
        registry.add("hEtaPMElec", "PM Electrons Eta", kTH1F, {{30, -1.5, +1.5, "#eta"}});
        registry.add("hPtJpsi", "Gen Jpsi pT", kTH1F, {{100, 0, 20, "p_{T} (GeV/#it{c})"}});
        registry.add("hEtaJpsi", "Gen Jpsi Eta", kTH1F, {{30, -1.5, +1.5, "#eta"}});
        registry.add("hPhiJpsi", "Gen Jpsi Phi", kTH1F, {{30, 0.0, 4.0, "#varphi"}});
        registry.add("hDaughtersPt", "Gen Daughters pT",kTH1F, {{200, 0.f, 20.f, "p_{T} (GeV/#it{c})"}});
        registry.add("hDaughtersEta", "Gen Daughters Eta",kTH1F, {{30, -1.5, +1.5, "#eta"}});
        registry.add("hDaughtersPhi", "Gen Daughters Phi",kTH1F, {{30, 0.0, +3.15,"#varphi"}});
        registry.add("hElecDaughtersPt", "Gen ElecDaughters pT",kTH1F, {{200, 0.f, 20.f, "p_{T} (GeV/#it{c})"}});
        registry.add("hElecDaughtersEta", "Gen ElecDaughters Eta",kTH1F, {{30, -1.5, +1.5, "#eta"}});
        registry.add("hElecDaughtersPhi", "Gen ElecDaughters Phi",kTH1F, {{30, 0.0, +3.15,"#varphi"}});
        registry.add("hPositDaughtersPt", "Gen PositDaughters pT",kTH1F, {{200, 0.f, 20.f, "p_{T} (GeV/#it{c})"}});
        registry.add("hPositDaughtersEta", "Gen PositDaughters Eta",kTH1F, {{30, -1.5, +1.5, "#eta"}});
        registry.add("hPositDaughtersPhi", "Gen PositDaughters Phi",kTH1F, {{30, 0.0, +3.15,"#varphi"}});
        registry.add("hDeltaDaughterTrackId", "Delta Track ID", kTH1F, {{5000, -0.5, 5000.5, "(DaughterId - GenElecId)"}});
        registry.add("hDaughtersSize", "Number of Daughters",kTH1F, {{10, -0.5, 9.5}});
        registry.add("hPtRecoElecDaughter", "Electrons Daughter pT (All with J/psi mother)", kTH1F, {{100, 0, 10, "p_{T} (GeV/#it{c})"}});
        registry.add("hPINRecoElecDaughter", "Electrons Daughter pIN (All with J/psi mother)", kTH1F, {{100, 0, 10, "p_{IN} (GeV/#it{c})"}});
        registry.add("hEtaRecoElecDaughter", "Electrons Daughter Eta (All with J/psi mother)", kTH1F, {{30, -1.5, +1.5, "#eta"}});
        registry.add("hPtRecoPositDaughter", "Positrons Daughter pT", kTH1F, {{100, 0, 20, "p_{T} (GeV/#it{c})"}});
        registry.add("hPINRecoPositDaughter", "Positrons Daughter pIN", kTH1F, {{100, 0, 20, "p_{IN} (GeV/#it{c})"}});
        registry.add("hEtaRecoPositDaughter", "Positrons Daughter Eta", kTH1F, {{30, -1.5, +1.5, "#eta"}});
        registry.add("hPtRecoMatchedElecDaughter", "Electrons Daughter pT (Only With Matched Positron)", kTH1F, {{100, 0, 20, "p_{T} (GeV/#it{c})"}});
        registry.add("hPINRecoMatchedElecDaughter", "Electrons Daughter pIN (Only With Matched Positron)", kTH1F, {{100, 0, 20, "p_{IN} (GeV/#it{c})"}});
        registry.add("hEtaRecoMatchedElecDaughter", "Electrons Daughter Eta (Only With Matched Positron)", kTH1F, {{30, -1.5, +1.5, "#eta"}});
        registry.add("hElecEffDist", "Electron PID Efficiency", kTH1F, {{100, 0.0, +1.0, "PID Efficiency"}});
        registry.add("hPositEffDist", "Positron PID Efficiency", kTH1F, {{100, 0.0, +1.0, "PID Efficiency"}});
        registry.add("hJpsiEffDist", "Positron PID Efficiency", kTH1F, {{100, 0.0, +1.0, "PID Efficiency"}});
        registry.add("hMatchedJpsiPt", "Matched-J/#psi p_{T}", kTH1F, {{400, 0.f, 20.f, "p_{T} (GeV/#it{c})"}});
        registry.add("hMatchedJpsiPtEta", "Matched-J/#psi p_{T} and #eta", kTH2F, {{400, 0.f, 20.f, "p_{T} (GeV/#it{c})"}, {400, -1.5f, 1.5f, "#eta"}});
        registry.add("hWeightedJpsiEffPt", "Weighted J/#psi PID Efficiency", kTH1F, {{400, 0.f, 20.f, "p_{T} (GeV/#it{c})"}});
        registry.add("hWeightedJpsiEffPtEta", "Weighted J/#psi PID Efficiency", kTH2F, {{400, 0.f, 20.f, "p_{T} (GeV/#it{c})"}, {400, -1.5f, 1.5f, "#eta"}});
        // registry.add("hNormJpsiEffpT", "Normalized J/#psi PID Efficiency", kTH1F, {{400, 0.f, 20.f, "p_{T} (GeV/#it{c})"}});
        // registry.add("hNormJpsiEffpTEta", "Normalized J/#psi PID Efficiency", kTH2F, {{400, 0.f, 20.f, "p_{T} (GeV/#it{c})"}, {400, -1.5f, 1.5f, "#eta"}});
        // registry.add("hJpsiEff2D", "J/#psi PID Efficiency", kTH2F, {{400, 0.f, 20.f, "p_{T} (GeV/#it{c})"}, {400, -1.5f, 1.5f, "#eta"}});
        TFile* fEffMap = TFile::Open(effMapPath.value.c_str());
        hPMElecEffMap = (TH2D*) fEffMap->Get("hElecPidEff");
        if (!hPMElecEffMap) {
            LOG(fatal) << "Electron PID map not found";
        }
    }

    void processDummy(MySimpleEvents const&)
    {
    }
    PROCESS_SWITCH(AnalysisPIDEfficiency, processDummy, "dummy task", true);

    void processJpsiPIDEfficiency(MySimpleEvents const& events, MySimpleBarrelTracks const& tracks, ReducedMCEvents const&, ReducedMCTracks const&)
    {
        int i = 0;
        for (auto& event : events) {
            registry.fill(HIST("hEvtCounter"), 0);
            if (test.value){
                if (i % 10 == 0) {
                    std::cout << "Event ID: " << event.globalIndex() << "\n";
                    std::cout << "Event i: " << i << "\n";
                }
                i++;
                if (i>Nevts.value) {
                    std::cout << "Reached maximum number of events for testing: " << Nevts.value << ". Stopping.\n";
                    return;
                };
            }
            for (auto& track : tracks) {
                if (track.has_reducedMCTrack()){
                    auto mcPart = track.reducedMCTrack();
                    if (mcPart.pdgCode() == 11) { // electron. for electron or positron, use std::abs. \\ TODO: think if this way, some reconstruction bias interfere in the efficiency
                        registry.fill(HIST("ptResolution"), track.pt(), track.pt() - mcPart.pt());
                        registry.fill(HIST("hEtaPMElec"), track.eta());
                        registry.fill(HIST("hPtPMElec"), track.pt());
                        if (mcPart.has_mothers()){
                            auto mother = mcPart.mothers_first_as<aod::ReducedMCTracks>();
                            if (mother.pdgCode() == 443) {
                                registry.fill(HIST("hPtJpsi"), mother.pt());
                                registry.fill(HIST("hEtaJpsi"), mother.eta());
                                registry.fill(HIST("hPhiJpsi"), mother.phi());
                                registry.fill(HIST("hPtRecoElecDaughter"), track.pt());
                                registry.fill(HIST("hPINRecoElecDaughter"), track.tpcInnerParam());
                                registry.fill(HIST("hEtaRecoElecDaughter"), track.eta());
                                auto genDaughters = mother.daughters_as<aod::ReducedMCTracks>();
                                registry.fill(HIST("hDaughtersSize"), genDaughters.size());
                                for (const auto& daughter : genDaughters) {
                                    if (daughter.pdgCode() == 11) {
                                        registry.fill(HIST("hElecDaughtersPt"),  daughter.pt());
                                        registry.fill(HIST("hElecDaughtersEta"), daughter.eta());
                                        registry.fill(HIST("hElecDaughtersPhi"), daughter.phi());
                                    }
                                    if (daughter.pdgCode() == -11) { //Positron
                                        registry.fill(HIST("hPositDaughtersPt"),  daughter.pt());
                                        registry.fill(HIST("hPositDaughtersEta"), daughter.eta());
                                        registry.fill(HIST("hPositDaughtersPhi"), daughter.phi());
                                        for (auto& positRecoDaughterCand : tracks) { // Looking for reco-level positron daughter. Needed because O2 stores links from reco to gen, not viceversa
                                            if (positRecoDaughterCand.has_reducedMCTrack()){
                                                auto genDaughterCand = positRecoDaughterCand.reducedMCTrack();
                                                if (genDaughterCand.pdgCode() == -11) { //Positron
                                                    if (genDaughterCand.globalIndex() == daughter.globalIndex()){
                                                        registry.fill(HIST("hPtRecoPositDaughter"), positRecoDaughterCand.pt());
                                                        registry.fill(HIST("hPINRecoPositDaughter"), positRecoDaughterCand.tpcInnerParam());
                                                        registry.fill(HIST("hEtaRecoPositDaughter"), positRecoDaughterCand.eta());
                                                        registry.fill(HIST("hPtRecoMatchedElecDaughter"), track.pt());
                                                        registry.fill(HIST("hPINRecoMatchedElecDaughter"), track.tpcInnerParam());
                                                        registry.fill(HIST("hEtaRecoMatchedElecDaughter"), track.eta());
                                                        int positMapBin = hPMElecEffMap->FindBin(positRecoDaughterCand.tpcInnerParam(), positRecoDaughterCand.eta());
                                                        double positEff = hPMElecEffMap->GetBinContent(positMapBin);
                                                        int elecMapBin = hPMElecEffMap->FindBin(track.tpcInnerParam(), track.eta());
                                                        double elecEff = hPMElecEffMap->GetBinContent(elecMapBin);
                                                        double jpsiEff = positEff * elecEff;
                                                        registry.fill(HIST("hElecEffDist"), elecEff);
                                                        registry.fill(HIST("hPositEffDist"), positEff);
                                                        registry.fill(HIST("hJpsiEffDist"), jpsiEff);
                                                        registry.fill(HIST("hMatchedJpsiPt"), mother.pt());
                                                        registry.fill(HIST("hMatchedJpsiPtEta"), mother.pt(), mother.eta());
                                                        registry.fill(HIST("hWeightedJpsiEffPt"), mother.pt(), jpsiEff);
                                                        registry.fill(HIST("hWeightedJpsiEffPtEta"), mother.pt(), mother.eta(), jpsiEff);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    registry.fill(HIST("hDaughtersPt"), daughter.pt());
                                    registry.fill(HIST("hDaughtersEta"), daughter.eta());
                                    registry.fill(HIST("hDaughtersPhi"), daughter.phi());
                                    int deltaDaughterTrack = std::abs(daughter.globalIndex() - mcPart.globalIndex());
                                    registry.fill(HIST("hDeltaDaughterTrackId"), deltaDaughterTrack);
                                }
                            }
                        }
                    }   
                }
            }
        }
    }
    PROCESS_SWITCH(AnalysisPIDEfficiency, processJpsiPIDEfficiency, "Calculation of Jpsi PID efficiency in terms of electron PID efficiency", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisPIDEfficiency>(cfgc)};
}