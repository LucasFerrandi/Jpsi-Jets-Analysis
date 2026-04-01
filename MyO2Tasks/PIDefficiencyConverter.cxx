// This code loops over reco-level electrons, find their gen-level J/psi mother, find their associated reco-level positron and, given an electron PID efficiency map, applies it to the electrons p_IN and eta in order to calculate the weighted J/psi PID efficiency, which must be normalized afterwards.

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

    ConfigurableAxis binPt{"binPt", {160, 0.f, 40.f}, "Binning of the pT axis of particles"};
    ConfigurableAxis binP_IN{"binP_IN", {160, 0.f, 40.f}, "Binning of the p_IN axis of particles"};
    ConfigurableAxis binEta{"binEta", {120, -1.5f, 1.5f}, "Binning of the eta axis of particles"};
    ConfigurableAxis binPhi{"binPhi", {80, 0.f, 4.f}, "Binning of the phi axis of particles"};
    ConfigurableAxis binPIDEff{"binPIDEff", {100, 0.f, 1.f}, "Binning of the PID efficiency axis"};

    Configurable<bool> test{"test", false, "Limits number of events"};
    Configurable<int> Nevts{"NevtsForTest", 100, "Number of events per DF (for testing)"};
    // Configurable<std::string> effMapPath{"effMapPath", "$HOME/alice/Jpsi-Jets-Analysis/JpsiWorkDir/PIDEfficiency/PIDEfficiencyConverter/toyPIDEffMap/toyElectronPidEff.root", "Path to root file with PID-efficiency map"};
    Configurable<std::string> effMapPath{"effMapPath", "$HOME/alice/Jpsi-Jets-Analysis/JpsiWorkDir/PIDEfficiency/PIDEfficiencyConverter/IdasElectronMaps_DQ_LHC24_pass1_skimmed_V0candidates/TrackBarrel_Conversions_withPID_nSigmaEl-4-4/effMap_TrackBarrel_Conversions_withPID_nSigmaEl-4-4.root", "Path to root file with PID-efficiency map"};
    TH2D* hPMElecEffMap = nullptr;
    void init(o2::framework::InitContext const&)
    {
        const AxisSpec axisPt{binPt, "p_{T} (GeV/c)"};
        const AxisSpec axisP_IN{binP_IN, "p_{IN} (GeV/c)"};
        const AxisSpec axisEta{binEta, "#eta"};
        const AxisSpec axisPhi{binPhi, "#varphi"};
        const AxisSpec axisPIDEff{binPIDEff, "PID Efficiency"};

        registry.add("hEvtCounter", "Event Counter", kTH1F, {{2, -0.5, 1.5}});
        registry.add("hEvtCounter2", "Event Counter", kTH1F, {{2, -0.5, 1.5}});
        registry.add("hPtEtaTracks", "Tracks p_{T} and #eta", kTH2F, {{axisPt}, {axisEta}});
        registry.add("hPtEtaTracksWithMC", "p_{T} and #eta of Tracks with MC", kTH2F, {{axisPt}, {axisEta}});
        registry.add("ptResolution", "PM Electrons ptResolution", kTH2F, {{axisPt}, {2000, -1.f, 1.f, "#Delta(p_{T}) (GeV/c)"}});
        registry.add("hPtPMElec", "PM Electrons p_{T}",kTH1F, {{axisPt}});
        registry.add("hEtaPMElec", "PM Electrons Eta", kTH1F, {{axisEta}});
        registry.add("hPhiPMElec", "PM Electrons Phi", kTH1F, {{axisPhi}});
        registry.add("hPtEtaPMElec", "PM Electrons p_{T} and #eta", kTH2F, {{axisPt}, {axisEta}});

        registry.add("hPtJpsi", "Gen Jpsi p_{T}", kTH1F, {{axisPt}});
        registry.add("hEtaJpsi", "Gen Jpsi Eta", kTH1F, {{axisEta}});
        registry.add("hPhiJpsi", "Gen Jpsi Phi", kTH1F, {{axisPhi}});
        registry.add("hPtEtaJpsi", "Gen Jpsi p_{T} and #eta", kTH2F, {{axisPt}, {axisEta}});

        registry.add("hPtPromptJpsi", "Gen Jpsi p_{T}", kTH1F, {{axisPt}});
        registry.add("hEtaPromptJpsi", "Gen Jpsi Eta", kTH1F, {{axisEta}});
        registry.add("hPhiPromptJpsi", "Gen Jpsi Phi", kTH1F, {{axisPhi}});

        registry.add("hPtNonPromptJpsi", "Gen Jpsi p_{T}", kTH1F, {{axisPt}});
        registry.add("hEtaNonPromptJpsi", "Gen Jpsi Eta", kTH1F, {{axisEta}});
        registry.add("hPhiNonPromptJpsi", "Gen Jpsi Phi", kTH1F, {{axisPhi}});

        registry.add("hDaughtersPt", "Gen Daughters p_{T}",kTH1F, {{axisPt}});
        registry.add("hDaughtersEta", "Gen Daughters Eta",kTH1F, {{axisEta}});
        registry.add("hDaughtersPhi", "Gen Daughters Phi",kTH1F, {{axisPhi}});
        registry.add("hElecDaughtersPt", "Gen ElecDaughters p_{T}",kTH1F, {{axisPt}});
        registry.add("hElecDaughtersEta", "Gen ElecDaughters Eta",kTH1F, {{axisEta}});
        registry.add("hElecDaughtersPhi", "Gen ElecDaughters Phi",kTH1F, {{axisPhi}});
        registry.add("hPositDaughtersPt", "Gen PositDaughters p_{T}",kTH1F, {{axisPt}});
        registry.add("hPositDaughtersEta", "Gen PositDaughters Eta",kTH1F, {{axisEta}});
        registry.add("hPositDaughtersPhi", "Gen PositDaughters Phi",kTH1F, {{axisPhi}});

        registry.add("hGenElecPt", "Gen Elec p_{T}",kTH1F, {{axisPt}});
        registry.add("hGenPositPt", "Gen Posit p_{T}",kTH1F, {{axisPt}});
        registry.add("hGenPionPlusPt", "Gen #pi^{+} p_{T}",kTH1F, {{axisPt}});
        registry.add("hGenPionMinusPt", "Gen #pi^{-} p_{T}",kTH1F, {{axisPt}});
        registry.add("hGenPionPlusPtEta", "Gen #pi^{+} p_{T} and #eta", kTH2F, {{axisPt}, {axisEta}});
        registry.add("hGenPionPlusPtEtaLowStat", "Gen #pi^{+} p_{T} and #eta (only 10% of pions)", kTH2F, {{axisPt}, {axisEta}});
        registry.add("hGenPionPlusPtEtaLowerStat", "Gen #pi^{+} p_{T} and #eta (only 1% of pions)", kTH2F, {{axisPt}, {axisEta}});
        registry.add("hGenProtonPt", "Gen Proton p_{T}",kTH1F, {{axisPt}});
        registry.add("hGenProtonPlusPtEta", "Gen proton p_{T} and #eta", kTH2F, {{axisPt}, {axisEta}});
        registry.add("hGenAntiprotonPt", "Gen Antiproton p_{T}",kTH1F, {{axisPt}});
        
        registry.add("hDeltaDaughterTrackId", "Delta Track ID", kTH1F, {{5000, -0.5, 5000.5, "(DaughterId - GenElecId)"}});
        registry.add("hDaughtersSize", "Number of Daughters",kTH1F, {{10, -0.5, 9.5}});
        registry.add("hPtRecoElecDaughter", "Electrons Daughter pT (All with J/psi mother)", kTH1F, {{axisPt}});
        registry.add("hPINRecoElecDaughter", "Electrons Daughter pIN (All with J/psi mother)", kTH1F, {{axisP_IN}});
        registry.add("hEtaRecoElecDaughter", "Electrons Daughter Eta (All with J/psi mother)", kTH1F, {{axisEta}});
        registry.add("hPtRecoPositDaughter", "Positrons Daughter pT", kTH1F, {{axisPt}});
        registry.add("hPINRecoPositDaughter", "Positrons Daughter pIN", kTH1F, {{axisP_IN}});
        registry.add("hEtaRecoPositDaughter", "Positrons Daughter Eta", kTH1F, {{axisEta}});
        registry.add("hPtRecoMatchedElecDaughter", "Electrons Daughter pT (Only With Matched Positron)", kTH1F, {{axisPt}});
        registry.add("hEtaRecoMatchedElecDaughter", "Electrons Daughter Eta (Only With Matched Positron)", kTH1F, {{axisEta}});
        registry.add("hPINRecoMatchedElecDaughter", "Electrons Daughter pIN (Only With Matched Positron)", kTH1F, {{axisP_IN}});
        registry.add("hPtRecoMatchedPositDaughter", "Electrons Daughter pT (Only With Matched Positron)", kTH1F, {{axisPt}});
        registry.add("hPINRecoMatchedPositDaughter", "Electrons Daughter pIN (Only With Matched Positron)", kTH1F, {{axisP_IN}});
        registry.add("hEtaRecoMatchedPositDaughter", "Electrons Daughter Eta (Only With Matched Positron)", kTH1F, {{axisEta}});
        registry.add("hElecEffDist", "Electron PID Efficiency", kTH1F, {{axisPIDEff}});
        registry.add("hPositEffDist", "Positron PID Efficiency", kTH1F, {{axisPIDEff}});
        registry.add("hJpsiEffDist", "Jpsi PID Efficiency", kTH1F, {{axisPIDEff}});

        registry.add("hMatchedJpsiPt", "Gen Matched J/#psi p_{T}", kTH1F, {{axisPt}});
        registry.add("hMatchedJpsiEta", "Gen Matched J/#psi Eta", kTH1F, {{axisEta}});
        registry.add("hMatchedJpsiPtEta", "Gen Matched J/#psi p_{T} and #eta", kTH2F, {{axisPt}, {axisEta}});
        registry.add("hMatchedJpsiPhi", "Gen Matched J/#psi Phi", kTH1F, {{axisPhi}});
        registry.add("hWeightedJpsiEffPtEta", "Weighted J/#psi PID Efficiency", kTH2F, {{axisPt}, {axisEta}});

        registry.add("hMatchedPromptJpsiPt", "Gen Matched Prompt J/#psi p_{T}", kTH1F, {{axisPt}});
        registry.add("hMatchedPromptJpsiEta", "Gen Matched Prompt J/#psi Eta", kTH1F, {{axisEta}});
        registry.add("hMatchedPromptJpsiPtEta", "Gen Matched Prompt J/#psi p_{T} and #eta", kTH2F, {{axisPt}, {axisEta}});
        registry.add("hMatchedPromptJpsiPhi", "Gen Matched Prompt J/#psi Phi", kTH1F, {{axisPhi}});
        registry.add("hWeightedPromptJpsiEffPtEta", "Weighted Prompt J/#psi PID Efficiency", kTH2F, {{axisPt}, {axisEta}});

        registry.add("hMatchedNonPromptJpsiPt", "Gen Matched Non-Prompt J/#psi p_{T}", kTH1F, {{axisPt}});
        registry.add("hMatchedNonPromptJpsiEta", "Gen Matched Non-Prompt J/#psi Eta", kTH1F, {{axisEta}});
        registry.add("hMatchedNonPromptJpsiPtEta", "Gen Matched Non-Prompt J/#psi p_{T} and #eta", kTH2F, {{axisPt}, {axisEta}});
        registry.add("hMatchedNonPromptJpsiPhi", "Gen Matched Non-Prompt J/#psi Phi", kTH1F, {{axisPhi}});
        registry.add("hWeightedNonPromptJpsiEffPtEta", "Weighted Non-Prompt J/#psi PID Efficiency", kTH2F, {{axisPt}, {axisEta}});

        TFile* fEffMap = TFile::Open(effMapPath.value.c_str());
        if (effMapPath.value.find("toyElectronPidEff") != std::string::npos) { // My artificial map
            hPMElecEffMap = (TH2D*) fEffMap->Get("hElecPidEff");
        }
        else { // Ida's electron maps
            hPMElecEffMap = (TH2D*) fEffMap->Get("effMap");
        }
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
        registry.fill(HIST("hEvtCounter"), 0);
        int i = 0;
        // for (auto& event : events) {
        int pionCounter = 0;
        registry.fill(HIST("hEvtCounter2"), 0);
        if (test.value){
            if (i % 10 == 0) {
                // std::cout << "Event ID: " << event.globalIndex() << "\n";
                // std::cout << "Event i: " << i << "\n";
            }
            i++;
            if (i>Nevts.value) {
                // std::cout << "Reached maximum number of events for testing: " << Nevts.value << ". Stopping.\n";
                return;
            };
        }
        for (auto& track : tracks) {
            registry.fill(HIST("hPtEtaTracks"), track.pt(), track.eta());
            if (track.has_reducedMCTrack()){
                registry.fill(HIST("hPtEtaTracksWithMC"), track.pt(), track.eta());
                auto mcPart = track.reducedMCTrack();
                if (mcPart.pdgCode() == 11) { // electron. for electron or positron, use std::abs.
                    registry.fill(HIST("hPtPMElec"), track.pt());
                    registry.fill(HIST("hEtaPMElec"), track.eta());
                    registry.fill(HIST("hPhiPMElec"), track.phi());
                    registry.fill(HIST("ptResolution"), track.pt(), track.pt() - mcPart.pt());
                    registry.fill(HIST("hPtEtaPMElec"), track.pt(), track.eta());
                    if (mcPart.has_mothers()){
                        auto mother = mcPart.mothers_first_as<aod::ReducedMCTracks>();
                        // auto mother2Id = track.mothersIds()[0];
                        // auto mother = mcStack.rawIteratorAt(motherId);
                        // Other possible ways to get the mother:
                        // 1 (dqEffAssoc): 
                        // currentMCParticle = currentMCParticle.template mothers_first_as<ReducedMCTracks>();
                        // 2 (tableMakerAssoc):
                        // for (auto& m : mctrack.mothersIds()) {
                        // auto aMother = mcTracks.rawIteratorAt(m)
                        // 3 (VarManager):
                        // auto motherId = track.mothersIds()[0];
                        // auto mother = mcStack.rawIteratorAt(motherId);
                        // 4 (MCSignal):
                        // auto mother = currentMCParticle.template mothers_first_as<P>();
                        // auto& mother2Id = mcPart.mothersIds()
                        // auto mother2 = mcTracks.rawIteratorAt(m)

                        if (mother.pdgCode() == 443) { //443 for jpsi, 511 fo B0, 521 fo B+
                            registry.fill(HIST("hPtJpsi"), mother.pt());
                            registry.fill(HIST("hEtaJpsi"), mother.eta());
                            registry.fill(HIST("hPhiJpsi"), mother.phi());
                            registry.fill(HIST("hPtEtaJpsi"), mother.pt(), mother.eta());
                            registry.fill(HIST("hPtRecoElecDaughter"), track.pt());
                            registry.fill(HIST("hPINRecoElecDaughter"), track.tpcInnerParam());
                            registry.fill(HIST("hEtaRecoElecDaughter"), track.eta());
                            if (mother.mcReducedFlags() == 40) {
                                registry.fill(HIST("hPtPromptJpsi"), mother.pt());
                                registry.fill(HIST("hEtaPromptJpsi"), mother.eta());
                                registry.fill(HIST("hPhiPromptJpsi"), mother.phi());
                            }
                            if (mother.mcReducedFlags() == 24) {
                                registry.fill(HIST("hPtNonPromptJpsi"), mother.pt());
                                registry.fill(HIST("hEtaNonPromptJpsi"), mother.eta());
                                registry.fill(HIST("hPhiNonPromptJpsi"), mother.phi());
                            }
                            auto genDaughters = mother.daughters_as<aod::ReducedMCTracks>();
                            registry.fill(HIST("hDaughtersSize"), genDaughters.size());
                            for (const auto& daughter : genDaughters) {
                                registry.fill(HIST("hDaughtersPt"), daughter.pt());
                                registry.fill(HIST("hDaughtersEta"), daughter.eta());
                                registry.fill(HIST("hDaughtersPhi"), daughter.phi());
                                int deltaDaughterTrack = std::abs(daughter.globalIndex() - mcPart.globalIndex());
                                registry.fill(HIST("hDeltaDaughterTrackId"), deltaDaughterTrack);
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
                                                if (genDaughterCand.globalIndex() == daughter.globalIndex()){ //Match!
                                                    registry.fill(HIST("hGenPositPt"), genDaughterCand.pt()); // To compare pt cut with hPositDaughtersPt
                                                    registry.fill(HIST("hGenElecPt"), mcPart.pt()); // To compare pt cut with hElecDaughtersPt
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
                                                    registry.fill(HIST("hMatchedJpsiEta"), mother.eta());
                                                    registry.fill(HIST("hMatchedJpsiPtEta"), mother.pt(), mother.eta());
                                                    registry.fill(HIST("hMatchedJpsiPhi"), mother.phi());
                                                    registry.fill(HIST("hWeightedJpsiEffPtEta"), mother.pt(), mother.eta(), jpsiEff);
                                                    if (mother.mcReducedFlags() == 40) { // Prompt
                                                        registry.fill(HIST("hMatchedPromptJpsiPt"), mother.pt());
                                                        registry.fill(HIST("hMatchedPromptJpsiEta"), mother.eta());
                                                        registry.fill(HIST("hMatchedPromptJpsiPtEta"), mother.pt(), mother.eta());
                                                        registry.fill(HIST("hMatchedPromptJpsiPhi"), mother.phi());
                                                        registry.fill(HIST("hWeightedPromptJpsiEffPtEta"), mother.pt(), mother.eta(), jpsiEff);
                                                    }
                                                    if (mother.mcReducedFlags() == 24) { // Non-prompt
                                                        registry.fill(HIST("hMatchedNonPromptJpsiPt"), mother.pt());
                                                        registry.fill(HIST("hMatchedNonPromptJpsiEta"), mother.eta());
                                                        registry.fill(HIST("hMatchedNonPromptJpsiPtEta"), mother.pt(), mother.eta());
                                                        registry.fill(HIST("hMatchedNonPromptJpsiPhi"), mother.phi());
                                                        registry.fill(HIST("hWeightedNonPromptJpsiEffPtEta"), mother.pt(), mother.eta(), jpsiEff);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if (test.value){
                    switch (mcPart.pdgCode()) {
                        case 211: //pion +
                            pionCounter++;
                            registry.fill(HIST("hGenPionPlusPt"), mcPart.pt());
                            registry.fill(HIST("hGenPionPlusPtEta"), mcPart.pt(), mcPart.eta());
                            if (pionCounter % 10 == 0) {
                                registry.fill(HIST("hGenPionPlusPtEtaLowStat"), mcPart.pt(), mcPart.eta());
                            }
                            if (pionCounter % 100 == 0) {
                                registry.fill(HIST("hGenPionPlusPtEtaLowerStat"), mcPart.pt(), mcPart.eta());
                            }
                            break;
                        case -211: //pion -
                            registry.fill(HIST("hGenPionMinusPt"), mcPart.pt());
                            break;
                        case 2212: //proton
                            registry.fill(HIST("hGenProtonPt"), mcPart.pt());
                            registry.fill(HIST("hGenProtonPlusPtEta"), mcPart.pt(), mcPart.eta());
                            break;
                        case -2212: //antiproton
                            registry.fill(HIST("hGenAntiprotonPt"), mcPart.pt());
                            break;
                        default:
                            break;
                    }
                }
            }
        }
        // }
    }
    PROCESS_SWITCH(AnalysisPIDEfficiency, processJpsiPIDEfficiency, "Calculation of Jpsi PID efficiency in terms of electron PID efficiency", false);
};

// struct AnalysisPIDEfficiencyPositronLoop{
//     HistogramRegistry registry{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
//     Configurable<bool> test{"test", false, "Limits number of events"};
//     Configurable<int> Nevts{"NevtsForTest", 100, "Number of events per DF (for testing)"};
//     Configurable<std::string> effMapPath{"effMapPath", "$HOME/alice/Jpsi-Jets-Analysis/JpsiWorkDir/PIDEfficiency/PIDEfficiencyConverter/toyElectronPidEff.root", "Path to root file with PID-efficiency map"};
//     TH2D* hPMElecEffMap = nullptr;
//     void init(o2::framework::InitContext const&)
//     {
//         registry.add("hEvtCounter", "Event Counter", kTH1F, {{2, -0.5, 1.5}});
//         registry.add("ptResolution", "PM Electrons ptResolution", kTH2F, {{axisPt}, {2000, -1.f, 1.f, "#Delta(p_{T}) (GeV/c)"}});
//         registry.add("hPtPMElec", "PM Electrons p_{T}",kTH1F, {{axisPt}});
//         registry.add("hEtaPMElec", "PM Electrons Eta", kTH1F, {{30, -1.5, +1.5, "#eta"}});
//         registry.add("hPtJpsi", "Gen Jpsi p_{T}", kTH1F, {{axisPt}});
//         registry.add("hEtaJpsi", "Gen Jpsi Eta", kTH1F, {{30, -1.5, +1.5, "#eta"}});
//         registry.add("hPhiJpsi", "Gen Jpsi Phi", kTH1F, {{axisPhi}});
//         registry.add("hDaughtersPt", "Gen Daughters p_{T}",kTH1F, {{axisPt}});
//         registry.add("hDaughtersEta", "Gen Daughters Eta",kTH1F, {{30, -1.5, +1.5, "#eta"}});
//         registry.add("hDaughtersPhi", "Gen Daughters Phi",kTH1F, {{axisPhi}});
//         registry.add("hElecDaughtersPt", "Gen ElecDaughters p_{T}",kTH1F, {{axisPt}});
//         registry.add("hElecDaughtersEta", "Gen ElecDaughters Eta",kTH1F, {{30, -1.5, +1.5, "#eta"}});
//         registry.add("hElecDaughtersPhi", "Gen ElecDaughters Phi",kTH1F, {{axisPhi}});
//         registry.add("hPositDaughtersPt", "Gen PositDaughters p_{T}",kTH1F, {{axisPt}});
//         registry.add("hPositDaughtersEta", "Gen PositDaughters Eta",kTH1F, {{30, -1.5, +1.5, "#eta"}});
//         registry.add("hPositDaughtersPhi", "Gen PositDaughters Phi",kTH1F, {{axisPhi}});

//         registry.add("hGenElecPt", "Gen Elec p_{T}",kTH1F, {{axisPt}});
//         registry.add("hGenPositPt", "Gen Posit p_{T}",kTH1F, {{axisPt}});
//         registry.add("hGenPionPlusPt", "Gen #pi^{+} p_{T}",kTH1F, {{axisPt}});
//         registry.add("hGenPionMinusPt", "Gen #pi^{-} p_{T}",kTH1F, {{axisPt}});
//         registry.add("hGenProtonPt", "Gen Proton p_{T}",kTH1F, {{axisPt}});
//         registry.add("hGenAntiprotonPt", "Gen Antiproton p_{T}",kTH1F, {{axisPt}});

//         registry.add("hDeltaDaughterTrackId", "Delta Track ID", kTH1F, {{5000, -0.5, 5000.5, "(DaughterId - GenElecId)"}});
//         registry.add("hDaughtersSize", "Number of Daughters",kTH1F, {{10, -0.5, 9.5}});
//         registry.add("hPtRecoElecDaughter", "Electrons Daughter pT (All with J/psi mother)", kTH1F, {{axisPt}});
//         registry.add("hPINRecoElecDaughter", "Electrons Daughter pIN (All with J/psi mother)", kTH1F, {{axisP_IN}});
//         registry.add("hEtaRecoElecDaughter", "Electrons Daughter Eta (All with J/psi mother)", kTH1F, {{30, -1.5, +1.5, "#eta"}});
//         registry.add("hPtRecoPositDaughter", "Positrons Daughter pT", kTH1F, {{axisPt}});
//         registry.add("hPINRecoPositDaughter", "Positrons Daughter pIN", kTH1F, {{axisP_IN}});
//         registry.add("hEtaRecoPositDaughter", "Positrons Daughter Eta", kTH1F, {{30, -1.5, +1.5, "#eta"}});
//         registry.add("hPtRecoMatchedElecDaughter", "Electrons Daughter pT (Only With Matched Positron)", kTH1F, {{axisPt}});
//         registry.add("hEtaRecoMatchedElecDaughter", "Electrons Daughter Eta (Only With Matched Positron)", kTH1F, {{30, -1.5, +1.5, "#eta"}});
//         registry.add("hPINRecoMatchedElecDaughter", "Electrons Daughter pIN (Only With Matched Positron)", kTH1F, {{axisP_IN}});
//         registry.add("hPtRecoMatchedPositDaughter", "Electrons Daughter pT (Only With Matched Positron)", kTH1F, {{axisPt}});
//         registry.add("hPINRecoMatchedPositDaughter", "Electrons Daughter pIN (Only With Matched Positron)", kTH1F, {{axisP_IN}});
//         registry.add("hEtaRecoMatchedPositDaughter", "Electrons Daughter Eta (Only With Matched Positron)", kTH1F, {{30, -1.5, +1.5, "#eta"}});
//         registry.add("hElecEffDist", "Electron PID Efficiency", kTH1F, {{axisPIDEff}});
//         registry.add("hPositEffDist", "Positron PID Efficiency", kTH1F, {{axisPIDEff}});
//         registry.add("hJpsiEffDist", "Jpsi PID Efficiency", kTH1F, {{axisPIDEff}});
//         registry.add("hMatchedJpsiPt", "Matched-J/#psi p_{T}", kTH1F, {{axisPt}});
//         registry.add("hMatchedJpsiPtEta", "Matched-J/#psi p_{T} and #eta", kTH2F, {{axisPt}, {axisEta}});
//         registry.add("hWeightedJpsiEffPt", "Weighted J/#psi PID Efficiency", kTH1F, {{axisPt}});
//         registry.add("hWeightedJpsiEffPtEta", "Weighted J/#psi PID Efficiency", kTH2F, {{axisPt}, {axisEta}});
//         // registry.add("hNormJpsiEffpT", "Normalized J/#psi PID Efficiency", kTH1F, {{axisPt}});
//         // registry.add("hNormJpsiEffpTEta", "Normalized J/#psi PID Efficiency", kTH2F, {{axisPt}, {axisEta}});
//         // registry.add("hJpsiEff2D", "J/#psi PID Efficiency", kTH2F, {{axisPt}, {axisEta}});
//         TFile* fEffMap = TFile::Open(effMapPath.value.c_str());
//         hPMElecEffMap = (TH2D*) fEffMap->Get("hElecPidEff");
//         if (!hPMElecEffMap) {
//             LOG(fatal) << "Electron PID map not found";
//         }
//     }

//     void processDummy(MySimpleEvents const&)
//     {
//     }
//     PROCESS_SWITCH(AnalysisPIDEfficiencyPositronLoop, processDummy, "dummy task", true);

//     void processJpsiPIDEfficiencyPositronLoop(MySimpleEvents const& events, MySimpleBarrelTracks const& tracks, ReducedMCEvents const&, ReducedMCTracks const&)
//     {
//         std::cout << "Starting processJpsiPIDEfficiencyPositronLoop\n";
//         int i = 0;
//         for (auto& event : events) {
//             registry.fill(HIST("hEvtCounter"), 0);
//             if (test.value){
//                 if (i % 10 == 0) {
//                     std::cout << "Event ID: " << event.globalIndex() << "\n";
//                     std::cout << "Event i: " << i << "\n";
//                 }
//                 i++;
//                 if (i>Nevts.value) {
//                     std::cout << "Reached maximum number of events for testing: " << Nevts.value << ". Stopping.\n";
//                     return;
//                 };
//             }
//             for (auto& track : tracks) {
//                 if (track.has_reducedMCTrack()){
//                     auto mcPart = track.reducedMCTrack();
//                     if (mcPart.pdgCode() == -11) { // Positron. for electron or positron, use std::abs. \\ TODO: think if this way, some reconstruction bias interfere in the efficiency
//                         registry.fill(HIST("ptResolution"), track.pt(), track.pt() - mcPart.pt());
//                         registry.fill(HIST("hEtaPMElec"), track.eta());
//                         registry.fill(HIST("hPtPMElec"), track.pt());
//                         if (mcPart.has_mothers()){
//                             auto mother = mcPart.mothers_first_as<aod::ReducedMCTracks>();
//                             if (mother.pdgCode() == 443) {
//                                 registry.fill(HIST("hPtJpsi"), mother.pt());
//                                 registry.fill(HIST("hEtaJpsi"), mother.eta());
//                                 registry.fill(HIST("hPhiJpsi"), mother.phi());
//                                 registry.fill(HIST("hPtRecoPositDaughter"), track.pt());
//                                 registry.fill(HIST("hPINRecoPositDaughter"), track.tpcInnerParam());
//                                 registry.fill(HIST("hEtaRecoPositDaughter"), track.eta());
//                                 auto genDaughters = mother.daughters_as<aod::ReducedMCTracks>();
//                                 registry.fill(HIST("hDaughtersSize"), genDaughters.size());
//                                 for (const auto& daughter : genDaughters) {
//                                     registry.fill(HIST("hDaughtersPt"), daughter.pt());
//                                     registry.fill(HIST("hDaughtersEta"), daughter.eta());
//                                     registry.fill(HIST("hDaughtersPhi"), daughter.phi());
//                                     int deltaDaughterTrack = std::abs(daughter.globalIndex() - mcPart.globalIndex());
//                                     registry.fill(HIST("hDeltaDaughterTrackId"), deltaDaughterTrack);
//                                     if (daughter.pdgCode() == -11) {
//                                         registry.fill(HIST("hPositDaughtersPt"),  daughter.pt());
//                                         registry.fill(HIST("hPositDaughtersEta"), daughter.eta());
//                                         registry.fill(HIST("hPositDaughtersPhi"), daughter.phi());
//                                     }
//                                     if (daughter.pdgCode() == 11) { //Electron
//                                         registry.fill(HIST("hElecDaughtersPt"),  daughter.pt());
//                                         registry.fill(HIST("hElecDaughtersEta"), daughter.eta());
//                                         registry.fill(HIST("hElecDaughtersPhi"), daughter.phi());
//                                         for (auto& elecRecoDaughterCand : tracks) { // Looking for reco-level electron daughter. Needed because O2 stores links from reco to gen, not viceversa
//                                             if (elecRecoDaughterCand.has_reducedMCTrack()){
//                                                 auto genDaughterCand = elecRecoDaughterCand.reducedMCTrack();
//                                                 if (genDaughterCand.pdgCode() == 11) { //Electron
//                                                     if (genDaughterCand.globalIndex() == daughter.globalIndex()){
//                                                         registry.fill(HIST("hGenElecPt"), genDaughterCand.pt()); // To compare pt cut with hElecDaughtersPt
//                                                         registry.fill(HIST("hGenPositPt"), mcPart.pt()); // To compare pt cut with hPositDaughtersPt
//                                                         registry.fill(HIST("hPtRecoElecDaughter"), elecRecoDaughterCand.pt());
//                                                         registry.fill(HIST("hPINRecoElecDaughter"), elecRecoDaughterCand.tpcInnerParam());
//                                                         registry.fill(HIST("hEtaRecoElecDaughter"), elecRecoDaughterCand.eta());
//                                                         registry.fill(HIST("hPtRecoMatchedPositDaughter"), track.pt());
//                                                         registry.fill(HIST("hPINRecoMatchedPositDaughter"), track.tpcInnerParam());
//                                                         registry.fill(HIST("hEtaRecoMatchedPositDaughter"), track.eta());
//                                                         int elecMapBin = hPMElecEffMap->FindBin(elecRecoDaughterCand.tpcInnerParam(), elecRecoDaughterCand.eta());
//                                                         double elecEff = hPMElecEffMap->GetBinContent(elecMapBin);
//                                                         int positMapBin = hPMElecEffMap->FindBin(track.tpcInnerParam(), track.eta());
//                                                         double positEff = hPMElecEffMap->GetBinContent(positMapBin);
//                                                         double jpsiEff = positEff * elecEff;
//                                                         registry.fill(HIST("hElecEffDist"), elecEff);
//                                                         registry.fill(HIST("hPositEffDist"), positEff);
//                                                         registry.fill(HIST("hJpsiEffDist"), jpsiEff);
//                                                         registry.fill(HIST("hMatchedJpsiPt"), mother.pt());
//                                                         registry.fill(HIST("hMatchedJpsiPtEta"), mother.pt(), mother.eta());
//                                                         registry.fill(HIST("hWeightedJpsiEffPt"), mother.pt(), jpsiEff);
//                                                         registry.fill(HIST("hWeightedJpsiEffPtEta"), mother.pt(), mother.eta(), jpsiEff);
//                                                     }
//                                                 }
//                                             }
//                                         }
//                                     }
//                                 }
//                             }
//                         }
//                     }
//                     if (test.value){
//                         switch (mcPart.pdgCode()) {
//                             case 211: //pion +
//                                 registry.fill(HIST("hGenPionPlusPt"), mcPart.pt());
//                                 break;
//                             case -211: //pion -
//                                 registry.fill(HIST("hGenPionMinusPt"), mcPart.pt());
//                                 break;
//                             case 2212: //proton
//                                 registry.fill(HIST("hGenProtonPt"), mcPart.pt());
//                                 break;
//                             case -2212: //antiproton
//                                 registry.fill(HIST("hGenAntiprotonPt"), mcPart.pt());
//                                 break;
//                             default:
//                                 break;
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     PROCESS_SWITCH(AnalysisPIDEfficiencyPositronLoop, processJpsiPIDEfficiencyPositronLoop, "Calculation of Jpsi PID efficiency in terms of electron PID efficiency looping over positrons in order to find Jpsi", false);
// };

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    // adaptAnalysisTask<AnalysisPIDEfficiencyPositronLoop>(cfgc),
    adaptAnalysisTask<AnalysisPIDEfficiency>(cfgc)};
}