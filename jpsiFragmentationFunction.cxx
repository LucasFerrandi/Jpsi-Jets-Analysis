#include <vector>

#include "TVector3.h"
// #include <format>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
// #include "PWGJE/Core/JetDQUtilities.h" ?

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"


struct JPsiFragmentationFunctionTask{
  HistogramRegistry registry{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ConfigurableAxis binDielPt{"binDielPt", {1000, 0.f, 30.f}, "Binning of the pT axis of particles"};
  ConfigurableAxis binJetPt{"binJetPt", {1000, 0.f, 150.f}, "Binning of the pT axis of jets"};
  ConfigurableAxis binEta{"binEta", {200, -1.0f, 1.0f}, "Binning of the eta axis"};
  ConfigurableAxis binPhi{"phiAxis", {200, -1.0f, 7.0f}, "Binning of the phi axis"};
  ConfigurableAxis binVtxZ{"binVtxZ", {200, -20.f, 20.f}, ""};
  ConfigurableAxis binDielMass{"binDielMass", {10000, 0.f, 35.f}, "Binning of the dielectron mass axis"};
  ConfigurableAxis binZ{"binZ", {20, 0.f, 1.f}, "Binning of the fragmentation function"};
  ConfigurableAxis binPtRatio{"binPtRatio", {700, 0.f, 1.f}, "Binning of the pT fraction"};

  Configurable<std::vector<float>> ptBounds{"ptBounds", {5.0f, 7.0f, 15.0f, 35.0f}, "Bounds of jet pT ranges"};
  Configurable<std::vector<float>> EBounds{"EBounds", {44.0f, 56.0f, 68.0f, 80.0f, 92.0f, 104.0f, 116.0f, 128.0f, 140.0f}, "Bounds of jet E ranges"};

  // using mcTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
  using mcTracks = soa::Join<aod::TracksIU, aod::McTrackLabels>;

  inline static std::vector<o2::framework::HistPtr> hDielZVsMassHists;
  inline static std::vector<o2::framework::HistPtr> hDielZVsMassHistsE;

  void init(InitContext const&) 
  {
    const AxisSpec dielPtAxis{binDielPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec jetPtAxis{binJetPt, "#it{p}_{T}^{jet} (GeV/#it{c})"};
    const AxisSpec etaAxis{binEta, "#eta"};
    const AxisSpec phiAxis{binPhi, "#phi"};
    const AxisSpec vtxZAxis{binVtxZ, "z vertex (cm)"};
    const AxisSpec dielMassAxis{binDielMass, "Dielectron mass (GeV/#it{c}^{2})"};
    const AxisSpec zAxis{binZ, "z"};
    const AxisSpec pTRatioAxis{binPtRatio, "z"};
    const AxisSpec counterAxis{1, 0.5, 1.5, ""};
    const AxisSpec deltaPtAxis{100, -1.0, +1.0, "#Delta(p_{T})"};
  
    registry.add("h_coll_evtsel", "Collision Event Selection;Collision Event Selection;Entries", {HistType::kTH1F, {{60, 200, 260}}});
    registry.add("h_coll_z", "Collision Z;Collision Z;Entries", {HistType::kTH1F, {{vtxZAxis}}});
    registry.add("h_coll_mult", "Collision Multiplicity;Collision Multiplicity;Entries", {HistType::kTH1F, {{500, -10., 3000.}}});
    registry.add("h_diel_mass", "Candidates Mass Distribution;Dielectron mass;Entries", {HistType::kTH1F, {{dielMassAxis}}});
    registry.add("h_diel_pt", "Dielectron p_{T};Dielectron p_{T};Entries", {HistType::kTH1F, {{dielPtAxis}}});
    registry.add("h_diel_eta", "Dielectron #eta;Dielectron #eta;Entries", {HistType::kTH1F, {{etaAxis}}});
    registry.add("h_diel_y", "Dielectron y;Dielectron y;Entries", {HistType::kTH1F, {{etaAxis}}});
    registry.add("h_diel_phi", "Dielectron #phi;Dielectron #phi;Entries", {HistType::kTH1F, {{phiAxis}}});
    registry.add("h_diel_sign", "Dielectron Sum of Signs;Dielectron sign;Entries", {HistType::kTH1F, {{200, -3, 3}}});
    registry.add("h_Jet_pt", "Jet p_{T};Jet p_{T};Entries", {HistType::kTH1F, {{jetPtAxis}}});
    registry.add("h_Jet_eta", "Jet #eta;Jet #eta;Entries", {HistType::kTH1F, {{etaAxis}}});
    registry.add("h_Jet_y", "Jet y;Jet y;Entries", {HistType::kTH1F, {{etaAxis}}});
    registry.add("h_Jet_phi", "Jet #phi;Jet #phi;Entries", {HistType::kTH1F, {{phiAxis}}});
    registry.add("h_Jet_area", "Jet Area;Jet area;Entries", {HistType::kTH1F, {{200, 0., 1.}}});
    registry.add("h_dielSize_per_jet", "Candidates per Jet;Diels per jets;Entries", {HistType::kTH1F, {{10, 0, 10}}});
    registry.add("h_jets_per_coll", "Jets Per Collision;Jets per collision;Entries", {HistType::kTH1F, {{10, 0, 10}}});
    registry.add("h_diel_jet_projection", "Candidate Momentum Projection to Jet;Dielectron jet projection", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_diel_jet_rejection", "Candidate Momentum Rejection to Jet;Dielectron jet projection", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_diel_pt_projection", "Candidate p_{T} ratio;Dielectron p_{T} projection", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_tracks_per_jets", "TracksId Size per Jet;TracksId size per jet; Entries", {HistType::kTH1F, {{40, 0, 40}}});
    registry.add("h_tracks_per_coll", "Tracks per Collision;Tracks per coll; Entries", {HistType::kTH1F, {{100, 0, 100}}});
    registry.add("h_diels_per_coll", "Candidates per Collision;Diels per coll; Entries", {HistType::kTH1F, {{10, 0, 10}}});
    registry.add("h_area_jet_vs_pT", "Jet Area vs p_{T};p_{T}^{jet} (GeV);A^{jet}", {HistType::kTH2F, {{jetPtAxis}, {50, 0., 1.}}});

    std::string h_diel_z_vs_mass_pTInclusive_name = Form("h_diel_z_vs_mass_pTSemiInclusive_%.0f_to_%.0fMeV", 1000*ptBounds.value[0], 1000*ptBounds.value.back());
    auto hDielZVsMasspTInclusiveHist = registry.add(h_diel_z_vs_mass_pTInclusive_name.c_str(), Form("Mass Distributions in z Bins, %.0f < p_{T}^{jet} < %.0f GeV;Mass (GeV);z", ptBounds.value[0], ptBounds.value.back()), {HistType::kTH2F, {{dielMassAxis}, {zAxis}}});
    for (size_t i = 0; i < ptBounds.value.size() - 1; i++) { //A 2D histogram for each pT range
      std::string h_diel_z_vs_mass_name = Form("h_diel_z_vs_mass_pT_%.0f_to_%.0fMeV", 1000*ptBounds.value[i], 1000*ptBounds.value[i+1]);
      auto hDielZVsMassHist = registry.add(h_diel_z_vs_mass_name.c_str(), Form("Mass Distributions in z bins, %.0f < p_{T}^{jet} < %.0f GeV;Mass (GeV);z", ptBounds.value[i], ptBounds.value[i+1]), {HistType::kTH2F, {{dielMassAxis}, {zAxis}}});
      hDielZVsMassHists.push_back(hDielZVsMassHist);
    }
    registry.add("h_Jet_energy", "Jet Energy;E_{jet};dN/dE_{jet}", {HistType::kTH1F, {{600, 0., 300.}}});
    registry.add("h_diel_jet_distance", "Candidate Angular Distance to Jet Axis;#DeltaR_{cand}^{jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 2.}}});
    registry.add("h_diel_jet_distance_vs_projection", "Candidate p_{T} Ratio vs Angular Distance to Jet Axis;#DeltaR_{cand}^{jet};z", {HistType::kTH2F, {{700, 0., 2.}, {pTRatioAxis}}});
    registry.add("h_diel_otherJet_distance", "Candidate Angular Distance to Same-Event Jets Axis;#DeltaR_{cand}^{jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 5.}}});
    registry.add("h_diel_otherJet_distance_vs_projection", "Candidate p_{T} Ratio vs Angular Distance to Same-Event Jets Axis;#DeltaR_{cand}^{jet};z", {HistType::kTH2F, {{700, 0., 5.}, {pTRatioAxis}}});
    registry.add("h_jet_otherJet_distance", "Candidate-Jet Angular Distance to Same-Event Jets Axis;#DeltaR_{cand}^{jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 5.}}});
    registry.add("h_jet_otherJet_distance_vs_projection", "Candidate p_{T} Ratio vs Candidate-Jet Angular Distance to Same-Event Jets Axis;#DeltaR_{cand}^{jet};z", {HistType::kTH2F, {{700, 0., 5.}, {pTRatioAxis}}});
    registry.add("h_diel_otherJet_Azimutdistance", "Candidate Azimuthal Distance to Same-Event Jets Axis;#Delta#phi_{cand}^{jet};dN/d(#Delta#phi)", {HistType::kTH1F, {{1000, -1., 4.}}});
    registry.add("h_diel_otherJet_Azimutdistance_vs_projection", "Candidate p_{T} Ratio vs Azimuthal Distance to Same-Event Jets Axis;#Delta#phi_{cand}^{jet};z", {HistType::kTH2F, {{700, -1., 4.}, {pTRatioAxis}}});
    registry.add("h_jet_otherJet_Azimutdistance", "Candidate-Jet Azimuthal Distance to Same-Event Jets Axis;#Delta#phi_{cand}^{jet};dN/d(#Delta#phi)", {HistType::kTH1F, {{1000, -1., 4.}}});
    registry.add("h_jet_otherJet_Azimutdistance_vs_projection", "Candidate p_{T} Ratio vs Candidate-Jet Azimuthal Distance to Same-Event Jets Axis;#Delta#phi_{cand}^{jet};z", {HistType::kTH2F, {{700, -1., 4.}, {pTRatioAxis}}});
    
    std::string h_diel_z_vs_mass_EInclusive_name = Form("h_diel_z_vs_mass_EnergySemiInclusive_%.0f_to_%.0fMeV", 1000*EBounds.value[0], 1000*EBounds.value.back());
    auto hDielZVsMassEInclusiveHist = registry.add(h_diel_z_vs_mass_EInclusive_name.c_str(), Form("Mass Distributions in z Bins, %.0f < E_{jet} < %.0f GeV;Mass (GeV);z", EBounds.value[0], EBounds.value.back()), {HistType::kTH2F, {{dielMassAxis}, {zAxis}}});
    for (size_t i = 0; i < EBounds.value.size() - 1; i++) { //A 2D histogram for each E bin, aiming analysis of Xi(E,z)
      std::string h_diel_z_vs_mass_E_name = Form("h_diel_z_vs_mass_Energy_%.0f_to_%.0fMeV", 1000*EBounds.value[i], 1000*EBounds.value[i+1]);
      auto hDielZVsMassEHist = registry.add(h_diel_z_vs_mass_E_name.c_str(), Form("Mass Distributions in z bins, %.0f < E_{jet} < %.0f GeV;Mass (GeV);z", EBounds.value[i], EBounds.value[i+1]), {HistType::kTH2F, {{dielMassAxis}, {zAxis}}});
      hDielZVsMassHistsE.push_back(hDielZVsMassEHist);
    }
    hDielZVsMassHists.push_back(hDielZVsMasspTInclusiveHist);
    hDielZVsMassHistsE.push_back(hDielZVsMassEInclusiveHist);

    // MC Histos
    registry.add("h_mcpeventCounter", "MCP Events Counter", kTH1F, {counterAxis});
    registry.add("h_mcdeventCounter", "MCD Events Counter", kTH1F, {counterAxis});
    registry.add("h_etaHistogram_all", "Eta, All Tracks", kTH1F, {etaAxis});
    registry.add("h_ptHistogram_all", "pT, All Tracks", kTH1F, {dielPtAxis});
    registry.add("h_ptResolution_all", "pT Resolution, All Tracks", kTH2F, {dielPtAxis, deltaPtAxis});
    registry.add("h_numberOfRecoCollisions", "numberOfRecoCollisions", kTH1F, {{10,-0.5f, 9.5f}});
    registry.add("h_multiplicityCorrelation", "multiplicityCorrelations", kTH2F, {{100, -0.5f, 99.5f}, {100,-0.5f, 99.5f}});
    // registry.add("h_pdg_code_negvalues", "PDG code;PDG code;Entries", {HistType::kTH1F, {{10000, -10000, 0}}});
    registry.add("h_pdg_code_smallvalues", "PDG code;PDG code;Entries", {HistType::kTH1F, {{20001, -10000.5, 10000.5}}});
    // registry.add("h_pdg_code_bigvalues", "PDG code;PDG code;Entries", {HistType::kTH1F, {{12000000, 46000000, 58000000}}});
    // registry.add("h_pdg_code_verybigvalues", "PDG code;PDG code;Entries", {HistType::kTH1F, {{8000000, 1000000000, 1008000000}}});
    // registry.add("h_pdg_code_all", "PDG code;PDG code;Entries", {HistType::kTH1F, {{10000, 0, 1008000000}}});
    registry.add("h_pdg_code_mothers", "PDG code, All Mothers;PDG code;Entries", {HistType::kTH1F, {{20001, -10000.5, 10000.5}}});
    // registry.add("h_pdg_code_daughters", "PDG code;PDG code;Entries", {HistType::kTH1F, {{20001, -10000, 10000}}});

    // MC Match and Resolution Histos
    registry.add("h_mcp_coll_z", "Collision Z;Collision Z (cm);Entries", {HistType::kTH1F, {{vtxZAxis}}});
    registry.add("h_mcd_coll_z", "Collision Z;Collision Z (cm);Entries", {HistType::kTH1F, {{vtxZAxis}}});
    registry.add("h_match_coll_z", "MCP vs MCD Collision Z;MCP Collision Z (cm);MCD Collision Z (cm)", {HistType::kTH2F, {{vtxZAxis}, {vtxZAxis}}});
    registry.add("h_resolution_coll_z", "Collision Z Resolution;MCP posZ (cm);#Delta posZ (cm)", {HistType::kTH2F, {{vtxZAxis}, {vtxZAxis}}});
    registry.add("h_mcp_jpsipt", "J/Psi p_{T};J/Psi p_{T} (GeV);Entries", {HistType::kTH1F, {{dielPtAxis}}});
    registry.add("h_mcd_jpsipt", "J/Psi p_{T};J/Psi p_{T} (GeV);Entries", {HistType::kTH1F, {{dielPtAxis}}});
    registry.add("h_match_jpsipt", "MCP vs MCD J/Psi p_{T};MCP J/Psi p_{T} (GeV);MCD J/Psi p_{T} (GeV)", {HistType::kTH2F, {{dielPtAxis}, {dielPtAxis}}});
    registry.add("h_resolution_jpsipt", "J/Psi p_{T} Resolution;MCP p_{T} (GeV);#Delta p_{T} (GeV)", {HistType::kTH2F, {{dielPtAxis}, {dielPtAxis}}});
    registry.add("h_mothers_size", "Number of Mothers;Number of Mothers;Entries", {HistType::kTH1F, {{6, -0.5, 5.5}}});
    // TODO: Add eta and phi histos here
  
  } 
  int jet_globalid = 0;

  void processDummy(JetCollision const&)
  {
  }
  PROCESS_SWITCH(JPsiFragmentationFunctionTask, processDummy, "dummy task", true);

  void processJPsiJetsData(JetCollision const& collision, soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents> const& jets, aod::CandidatesDielectronData const& cands, JetTracks const& tracks)
  {
    registry.fill(HIST("h_coll_z"), collision.posZ());
    registry.fill(HIST("h_coll_mult"), collision.multiplicity());
    registry.fill(HIST("h_coll_evtsel"), collision.eventSel());
    int dielsInColl = 0;
    for (auto& cand : cands) {
      dielsInColl++;
    }
    registry.fill(HIST("h_diels_per_coll"), dielsInColl);
    int trackInColl = 0;
    for (auto& track : tracks) {
      trackInColl++;
    }
    registry.fill(HIST("h_tracks_per_coll"), trackInColl);
    int jetInColl = 0;
    for (auto& jet : jets) {
      bool jetHasJPsi = false;
      registry.fill(HIST("h_Jet_pt"), jet.pt());
      registry.fill(HIST("h_Jet_eta"), jet.eta());
      registry.fill(HIST("h_Jet_y"), jet.y());
      registry.fill(HIST("h_Jet_phi"), jet.phi());
      registry.fill(HIST("h_Jet_area"), jet.area());
      registry.fill(HIST("h_Jet_energy"), jet.energy());
      if (jet.pt() > 10){
        registry.fill(HIST("h_area_jet_vs_pT"), jet.pt(), jet.area());
      }
      TVector3 jetVector(jet.px(), jet.py(), jet.pz());
      for (auto& jpsiCandidate : jet.candidates_as<aod::CandidatesDielectronData>()) {
        jetHasJPsi = true;
        registry.fill(HIST("h_diel_mass"), jpsiCandidate.mass());
        registry.fill(HIST("h_diel_pt"), jpsiCandidate.pt());
        registry.fill(HIST("h_diel_eta"), jpsiCandidate.eta());
        registry.fill(HIST("h_diel_y"), jpsiCandidate.y());
        registry.fill(HIST("h_diel_phi"), jpsiCandidate.phi());
        registry.fill(HIST("h_diel_sign"), jpsiCandidate.sign());
        TVector3 jpsiVector(jpsiCandidate.px(), jpsiCandidate.py(), jpsiCandidate.pz());
        double z_parallel = (jetVector * jpsiVector) / (jetVector * jetVector);
        double z_perpendicular = (jetVector.Cross(jpsiVector)).Mag() / jetVector.Mag();
        double pt_ratio = jpsiCandidate.pt() / jet.pt();
        registry.fill(HIST("h_diel_jet_projection"), z_parallel);
        registry.fill(HIST("h_diel_jet_rejection"), z_perpendicular);
        registry.fill(HIST("h_diel_pt_projection"), pt_ratio);
        double candDistToAxis = jetutilities::deltaR(jpsiCandidate, jet);
        registry.fill(HIST("h_diel_jet_distance"), candDistToAxis);
        registry.fill(HIST("h_diel_jet_distance_vs_projection"), candDistToAxis, z_parallel);
        int jetCounter = 0;
        for (auto& otherJet : jets) { //Distances JPsi-jet
          double candJetDistToJets = jetutilities::deltaR(jet, otherJet);
          if (candJetDistToJets < 0.001){
            continue;
          }
          registry.fill(HIST("h_jet_otherJet_distance"), candJetDistToJets);
          registry.fill(HIST("h_jet_otherJet_distance_vs_projection"), candJetDistToJets, z_parallel);
          double candDistToJets = jetutilities::deltaR(jpsiCandidate, otherJet); //otherJet
          registry.fill(HIST("h_diel_otherJet_distance"), candDistToJets);
          registry.fill(HIST("h_diel_otherJet_distance_vs_projection"), candDistToJets, z_parallel);
          double candAzimutDistToJets = RecoDecay::constrainAngle(std::abs(jpsiCandidate.phi() - otherJet.phi()), 0, 2);
          registry.fill(HIST("h_diel_otherJet_Azimutdistance"), candAzimutDistToJets);
          registry.fill(HIST("h_diel_otherJet_Azimutdistance_vs_projection"), candAzimutDistToJets, z_parallel);
          double candJetAzimutDistToJets = RecoDecay::constrainAngle(std::abs(jet.phi() - otherJet.phi()), 0, 2);
          registry.fill(HIST("h_jet_otherJet_Azimutdistance"), candJetAzimutDistToJets);
          registry.fill(HIST("h_jet_otherJet_Azimutdistance_vs_projection"), candJetAzimutDistToJets, z_parallel);
        }
      }
      registry.fill(HIST("h_dielSize_per_jet"), jet.candidatesIds().size()); //Always = 1 if doCandidateJetFinding = true
      registry.fill(HIST("h_tracks_per_jets"), jet.tracksIds().size());
      jet_globalid++;
      jetInColl++;
      for (size_t i = 0; i < hDielZVsMassHists.size(); i++) { //Checks in which pT range the jet is
        if (ptBounds.value[i] < jet.pt() && jet.pt() < ptBounds.value[i+1]) {
          for (auto& jpsiCandidate : jet.candidates_as<aod::CandidatesDielectronData>()) {
            TVector3 jpsiVector(jpsiCandidate.px(), jpsiCandidate.py(), jpsiCandidate.pz());
            double pt_ratio = jpsiCandidate.pt() / jet.pt();
            // std::string z_mass_histname ="h_diel_z_vs_mass_" + std::format("{:.2f}", ptBounds.value[i]) + "_to_" + std::format("{:.2f}", ptBounds.value[i+1]) + "GeV";
            // registry.fill(HIST(h_diel_z_vs_mass_names[i].c_str()), pt_ratio, jpsiCandidate.mass());
            // registry.fill(HIST(hDielZVsMassHists[i]), pt_ratio, jpsiCandidate.mass());
            // hDielZVsMassHists[i]->Fill(pt_ratio, jpsiCandidate.mass());
            // std::get<std::shared_ptr<TH2F>>(hDielZVsMassHists[i])->Fill(pt_ratio, jpsiCandidate.mass());
            std::get<std::shared_ptr<TH2>>(hDielZVsMassHists[i])->Fill(jpsiCandidate.mass(), pt_ratio);
            std::get<std::shared_ptr<TH2>>(hDielZVsMassHists.back())->Fill(jpsiCandidate.mass(), pt_ratio);
            // registry.fill(HIST(h_diel_z_vs_mass_pTInclusive_name.c_str()), jpsiCandidate.mass(), pt_ratio);
            // registry.get<TH2F>(hDielZVsMassHists[i])->Fill(pt_ratio, jpsiCandidate.mass());
          }
        }
      }
      for (size_t i = 0; i < hDielZVsMassHistsE.size(); i++) { //Checks in which E range the jet is
        if (EBounds.value[i] < jet.energy() && jet.energy() < EBounds.value[i+1]) {
          for (auto& jpsiCandidate : jet.candidates_as<aod::CandidatesDielectronData>()) {
            TVector3 jpsiVector(jpsiCandidate.px(), jpsiCandidate.py(), jpsiCandidate.pz());
            double pt_ratio = jpsiCandidate.pt() / jet.pt();
            std::get<std::shared_ptr<TH2>>(hDielZVsMassHistsE[i])->Fill(jpsiCandidate.mass(), pt_ratio);
            std::get<std::shared_ptr<TH2>>(hDielZVsMassHistsE.back())->Fill(jpsiCandidate.mass(), pt_ratio);
            // registry.fill(HIST(h_diel_z_vs_mass_EInclusive_name.c_str()), jpsiCandidate.mass(), pt_ratio);
          }
        }
      }
    }
  registry.fill(HIST("h_jets_per_coll"), jetInColl);
  // registry.fill(HIST("h_coll_id"), coll_id);
  // coll_id++;
  }
  PROCESS_SWITCH(JPsiFragmentationFunctionTask, processJPsiJetsData, "JPsi Fragmentation Function Data Task", false);

  int mcPartCounter = 0;
  int trackCounter = 0;
  void processJPsiMCValidation(aod::Collision const& collision, mcTracks const& tracks, aod::McParticles const&) // mcParticle is accessed by the tracks
  {
    registry.fill(HIST("h_mcdeventCounter"), 1.0);
    for (const auto& track : tracks) {
      trackCounter++;
      // if( track.tpcNClsCrossedRows() < 70 ) continue; //badly tracked
      registry.fill(HIST("h_etaHistogram_all"), track.eta());
      registry.fill(HIST("h_ptHistogram_all"), track.pt());
      if(track.has_mcParticle()){ 
        mcPartCounter++;
        auto mcParticle = track.mcParticle();
        registry.fill(HIST("h_ptResolution_all"), track.pt(), track.pt() - mcParticle.pt());
        // registry.fill(HIST("h_pdg_code_all"), mcParticle.pdgCode());
        // if(mcParticle.pdgCode() < 0) {
          // registry.fill(HIST("h_pdg_code_negvalues"), mcParticle.pdgCode());
          // std::cout << "Negative PDG code: " << mcParticle.pdgCode() << std::endl;
        // } else if (mcParticle.pdgCode() < 10000) {
        if (abs(mcParticle.pdgCode()) < 10000) { 
          registry.fill(HIST("h_pdg_code_smallvalues"), mcParticle.pdgCode());
          // std::cout << "Small PDG code: " << mcParticle.pdgCode() << std::endl;
        // } else if (mcParticle.pdgCode() < 58000000) {
        //   registry.fill(HIST("h_pdg_code_bigvalues"), mcParticle.pdgCode());
          // std::cout << "Big PDG code: " << mcParticle.pdgCode() << std::endl;
        // } else if (mcParticle.pdgCode() < 1008000000) {
          // registry.fill(HIST("h_pdg_code_verybigvalues"), mcParticle.pdgCode());
          // std::cout << "Very big PDG code: " << mcParticle.pdgCode() << std::endl;
        // } else if (mcParticle.pdgCode() > 1008000000) {
          // std::cout << "PDG code out of range: " << mcParticle.pdgCode() << std::endl;
        }
        registry.fill(HIST("h_mothers_size"), mcParticle.mothersIds().size());
        if(mcParticle.has_mothers()){
          // Check first mother
          auto const& mother = mcParticle.mothers_first_as<aod::McParticles>();
          // LOGF(info, "First mother: %d has pdg code %d", mother.globalIndex(), mother.pdgCode());
          if(mother.pdgCode() > -10000 && mother.pdgCode() < 10000) {
            registry.fill(HIST("h_pdg_code_mothers"), mother.pdgCode());
          }
          // if(mcParticle.mothersIds().size() > 1){
          //   std::cout << "Mother 1: " << mcParticle.mothersIds()[0].pdgCode() <<  std::endl;
          //   std::cout << "Mother 2: " << mcParticle.mothersIds()[1].pdgCode() <<  std::endl;
          }
          if(mother.pdgCode() == 443){
            registry.fill(HIST("h_mcp_jpsipt"), mother.pt()); // This way duplicate JPsis needs to be excluded!
            // registry.fill(HIST("h_mcp_jpsipt"), mcParticle.pt());
            // registry.fill(HIST("h_mcp_coll_z"), collision.posZ());
            // // std::cout << "JPsi found with PDG code: " << mother.pdgCode() << " and globalIndex: " << mother.globalIndex() << std::endl;
            // for (size_t i = 0; i < hDielZVsMassHists.size(); i++) { //Checks in which pT range the JPsi is
            //   if (ptBounds.value[i] < mcParticle.pt() && mcParticle.pt() < ptBounds.value[i+1]) {
            //     std::get<std::shared_ptr<TH2>>(hDielZVsMassHists[i])->Fill(mcParticle.mass(), 1.0); //z=1.0 for JPsi
            //     std::get<std::shared_ptr<TH2>>(hDielZVsMassHists.back())->Fill(mcParticle.mass(), 1.0);
            //   }
            // }
            // for (size_t i = 0; i < hDielZVsMassHistsE.size(); i++) { //Checks in which E range the JPsi is
            //   if (EBounds.value[i] < mcParticle.energy() && mcParticle.energy() < EBounds.value[i+1]) {
            //     std::get<std::shared_ptr<TH2>>(hDielZVsMassHistsE[i])->Fill(mcParticle.mass(), 1.0); //z=1.0 for JPsi
            //     std::get<std::shared_ptr<TH2>>(hDielZVsMassHistsE.back())->Fill(mcParticle.mass(), 1.0);
            //   }
            // }
          }
        }
          // Access mothers (apparently max 2)
          // int motherCounter = 0;
          // for (auto& m : mcParticle.mothers_as<aod::McParticles>()) {
          //   std::cout << "Mother " << motherCounter << " PDG code and globalIndex: " << m.pdgCode() << ", "<< m.globalIndex() << std::endl;
          //   motherCounter++;
          // }
        // if (mcParticle.has_daughters()) {
        //   auto const& daughter = mcParticle.daughters_first_as<aod::McParticles>();
        //   LOGF(info, "First daughter: %d has pdg code %d", daughter.globalIndex(), daughter.pdgCode());
        //   if(daughter.pdgCode() > -10000 && daughter.pdgCode() < 10000) {
        //     registry.fill(HIST("h_pdg_code_daughters"), daughter.pdgCode());
        //   }
        //   for (auto& d : mcParticle.daughters_as<aod::McParticles>()) {
        //     LOGF(info, "D2 %d %d", mcParticle.globalIndex(), d.globalIndex());
        //   }
        // }
      }
    }
    // std::cout << "Number of MC particles: " << mcPartCounter << std::endl;
    // std::cout << "Number of tracks: " << trackCounter << std::endl;
  } 
  PROCESS_SWITCH(JPsiFragmentationFunctionTask, processJPsiMCValidation, "JPsi Fragmentation Function MC Validation", false);


  //   void processJPsiJetsMC(JetCollision const& collision, soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents> const& jets, aod::CandidatesDielectronMCD const& cands, JetTracks const& tracks) {
  //     registry.fill(HIST("h_coll_z"), collision.posZ());
  //     registry.fill(HIST("h_coll_mult"), collision.multiplicity());
  //     registry.fill(HIST("h_coll_evtsel"), collision.eventSel());
  //     int dielsInColl = 0;
  //     for (auto& cand : cands) {
  //       dielsInColl++;
  //     }
  //     registry.fill(HIST("h_diels_per_coll"), dielsInColl);
  //     int trackInColl = 0;
  //     for (auto& track : tracks) {
  //       trackInColl++;
  //     }
  //     registry.fill(HIST("h_tracks_per_coll"), trackInColl);
  //     int jetInColl = 0;
  //     for (auto& jet : jets) {
  //       bool jetHasJPsi = false;
  //       registry.fill(HIST("h_Jet_pt"), jet.pt());
  //       registry.fill(HIST("h_Jet_eta"), jet.eta());
  //       registry.fill(HIST("h_Jet_y"), jet.y());
  //       registry.fill(HIST("h_Jet_phi"), jet.phi());
  //       registry.fill(HIST("h_Jet_area"), jet.area());
  //       registry.fill(HIST("h_Jet_energy"), jet.energy());
  //       if (jet.pt() > 10){
  //         registry.fill(HIST("h_area_jet_vs_pT"), jet.pt(), jet.area());
  //       }
  //       TVector3 jetVector(jet.px(), jet.py(), jet.pz());
  //       for (auto& jpsiCandidate : jet.candidates_as<aod::CandidatesDielectronMCD>()) {
  //         jetHasJPsi = true;
  //         registry.fill(HIST("h_diel_mass"), jpsiCandidate.mass());
  //         registry.fill(HIST("h_diel_pt"), jpsiCandidate.pt());
  //         registry.fill(HIST("h_diel_eta"), jpsiCandidate.eta());
  //         registry.fill(HIST("h_diel_y"), jpsiCandidate.y());
  //         registry.fill(HIST("h_diel_phi"), jpsiCandidate.phi());
  //         registry.fill(HIST("h_diel_sign"), jpsiCandidate.sign());
  //         TVector3 jpsiVector(jpsiCandidate.px(), j
  // psiCandidate.py(), jpsiCandidate.pz());
  //         double z_parallel = (jetVector * jpsiVector) / (jetVector * jetVector);
  //         double z_perpendicular = (jetVector.Cross(jpsiVector)).Mag() / jetVector.Mag();
  //         double pt_ratio = jpsiCandidate.pt() / jet.pt();
  //         registry.fill(HIST("h_diel_jet_projection"), z_parallel);
  //         registry.fill(HIST("h_diel_jet_rejection"), z_perpendicular);
  //         registry.fill(HIST("h_diel_pt_projection"), pt_ratio);
  //         double candDistToAxis = jetutilities::deltaR(jpsiCandidate, jet);
  //         registry.fill(HIST("h_diel_jet_distance"), candDistToAxis);
  //         registry.fill(HIST("h_diel_jet_distance_vs_projection"), candDistToAxis, z_parallel);
  //         int jetCounter = 0;
  //         for (auto& otherJet : jets) { //Distances JPsi-jet
  //           double candJetDistToJets = jetutilities::deltaR(jet, otherJet);
  //           if (candJetDistToJets < 0.001){
  //             continue;
  //           }
  //           registry.fill(HIST("h_jet_otherJet_distance"), candJetDistToJets);
  //           registry.fill(HIST("h_jet_otherJet_distance_vs_projection"), candJetDistToJets, z_parallel);
  //           double candDistToJets = jetutilities::deltaR(jpsiCandidate, otherJet); //otherJet
  //           registry.fill(HIST("h_diel_otherJet_distance"), candDistToJets);
  //           registry.fill(HIST("h_diel_otherJet_distance_vs_projection"), candDistToJets, z_parallel);
  //           double candAzimutDistToJets = RecoDecay::constrainAngle(std::abs(jpsiCandidate.phi() - otherJet.phi()), 0, 2);
  //           registry.fill(HIST("h_diel_otherJet_Azimutdistance"), candAzimutDistToJets);
  //           registry.fill(HIST("h_diel_otherJet_Azimutdistance_vs_projection"), candAzimutDistToJets, z_parallel);
  //           double candJetAzimutDistToJets = RecoDecay::constrainAngle(std::abs(jet.phi() - otherJet.phi()), 0, 2);
  //           registry.fill(HIST("h_jet_otherJet_Azimutdistance"), candJetAzimutDistToJets);
  //           registry.fill(HIST("h_jet_otherJet_Azimutdistance_vs
  // _projection"), candJetAzimutDistToJets, z_parallel);
  //         }
  //       }
  //       registry.fill(HIST("h_dielSize_per_jet"), jet.candidatesIds().size()); //Always = 1 if doCandidateJetFinding = true
  //       registry.fill(HIST("h_tracks_per_jets"), jet.tracksIds().size());
  //       jet_globalid++;
  //       jetInColl++;
  //       for (size_t i = 0; i < hDielZVsMassHists.size(); i++) { //Checks in which pT range the jet is
  //         if (ptBounds.value[i] < jet.pt() && jet.pt() < ptBounds.value[i+1]) {
  //           for (auto& jpsiCandidate : jet.candidates_as<aod::CandidatesDielectronMCD>()) {
  //             TVector3 jpsiVector(jpsiCandidate.px(), jpsiCandidate.py(), jpsiCandidate.pz());
  //             double pt_ratio = jpsiCandidate.pt() / jet.pt();
  //             std::get<std::shared_ptr<TH2>>(hDielZVsMassHists[i])->Fill(jpsiCandidate.mass(), pt_ratio);
  //             std::get<std::shared_ptr<TH2>>(hDielZVsMassHists.back())->Fill(jpsiCandidate.mass(), pt_ratio);
  //           }
  //         }
  //       }
  //       for (size_t i = 0; i < hDielZVsMassHistsE.size(); i++) { //Checks in which E range the jet is
  //         if (EBounds.value[i] < jet.energy() && jet.energy() < EBounds.value[i+1]) {
  //           for (auto& jpsiCandidate : jet.candidates_as<aod::CandidatesDielectronMCD>()) {
  //             TVector3 jpsiVector(jpsiCandidate.px(), jpsiCandidate.py(), jpsiCandidate.pz());
  //             double pt_ratio = jpsiCandidate.pt() / jet.pt();
  //             std::get<std::shared_ptr<TH2>>(hDielZVsMassHistsE[i])->Fill(jpsiCandidate.mass(), pt_ratio);
  //             std::get<std::shared_ptr<TH2>>(hDielZVsMassHistsE.back())->Fill(jpsiCandidate.mass(), pt_ratio);
  //           }
  //         }
  //       }
  //     }
  //     registry.fill(HIST("h_jets_per_coll"), jetInColl);
  //   }

  //   PROCESS_SWITCH(JPsiFragmentationFunctionTask, processJPsiJetsMC, "JPsi Fragmentation Function MC Task", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JPsiFragmentationFunctionTask>(cfgc)};
}
