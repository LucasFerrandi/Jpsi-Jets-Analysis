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
  Configurable<std::vector<float>> ptBounds{"ptBounds", {5.0f, 7.0f, 15.0f, 35.0f}, "Bounds of jet pT ranges"};
  Configurable<std::vector<float>> EBounds{"EBounds", {44.0f, 56.0f, 68.0f, 80.0f, 92.0f, 104.0f, 116.0f, 128.0f, 140.0f}, "Bounds of jet E ranges"};
  inline static std::vector<o2::framework::HistPtr> hDielZVsMassHists;
  inline static std::vector<o2::framework::HistPtr> hDielZVsMassHistsE;
  void init(InitContext const&) {
    registry.add("h_coll_evtsel", "Collision Event Selection;Collision Event Selection;Entries", {HistType::kTH1F, {{60, 200, 260}}});
    registry.add("h_coll_z", "Collision Z;Collision Z;Entries", {HistType::kTH1F, {{100, -20., 20.}}});
    registry.add("h_coll_mult", "Collision Multiplicity;Collision Multiplicity;Entries", {HistType::kTH1F, {{500, -10., 3000.}}});
    // registry.add("h_coll_ambig", ";Collision is_ambiguous", {HistType::kTH1F, {{4, -2., 2.}}});
    registry.add("h_diel_mass", "Candidates Mass Distribution;Dielectron mass;Entries", {HistType::kTH1F, {{10000, 0., 18.}}});
    registry.add("h_diel_pt", "Dielectron p_{T};Dielectron p_{T};Entries", {HistType::kTH1F, {{1000, 0., 30.}}});
    registry.add("h_diel_eta", "Dielectron #eta;Dielectron #eta;Entries", {HistType::kTH1F, {{200, -2., 2.}}});
    registry.add("h_diel_y", "Dielectron y;Dielectron y;Entries", {HistType::kTH1F, {{200, -2., 2.}}});
    registry.add("h_diel_phi", "Dielectron #phi;Dielectron #phi;Entries", {HistType::kTH1F, {{200, -1., 7.}}});
    registry.add("h_diel_sign", "Dielectron Sum of Signs;Dielectron sign;Entries", {HistType::kTH1F, {{200, -3, 3}}});
    registry.add("h_Jet_pt", "Jet p_{T};Jet p_{T};Entries", {HistType::kTH1F, {{1000, 0., 150.}}});
    registry.add("h_Jet_eta", "Jet #eta;Jet #eta;Entries", {HistType::kTH1F, {{200, -1., 1.}}});
    registry.add("h_Jet_y", "Jet y;Jet y;Entries", {HistType::kTH1F, {{200, -1., 1.}}});
    registry.add("h_Jet_phi", "Jet #phi;Jet #phi;Entries", {HistType::kTH1F, {{200, -1., 7.}}});
    registry.add("h_Jet_area", "Jet Area;Jet area;Entries", {HistType::kTH1F, {{200, 0., 1.}}});
    registry.add("h_dielSize_per_jet", "Candidates per Jet;Diels per jets;Entries", {HistType::kTH1F, {{10, 0, 10}}});
    registry.add("h_jets_per_coll", "Jets Per Collision;Jets per collision;Entries", {HistType::kTH1F, {{10, 0, 10}}});
    registry.add("h_diel_jet_projection", "Candidate Momentum Projection to Jet;Dielectron jet projection", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_diel_jet_rejection", "Candidate Momentum Rejection to Jet;Dielectron jet projection", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_diel_pt_projection", "Candidate p_{T} ratio;Dielectron p_{T} projection", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_tracks_per_jets", "TracksId Size per Jet;TracksId size per jet; Entries", {HistType::kTH1F, {{40, 0, 40}}});
    registry.add("h_tracks_per_coll", "Tracks per Collision;Tracks per coll; Entries", {HistType::kTH1F, {{100, 0, 100}}});
    registry.add("h_diels_per_coll", "Candidates per Collision;Diels per coll; Entries", {HistType::kTH1F, {{10, 0, 10}}});
    registry.add("h_area_jet_vs_pT", "Jet Area vs p_{T};p_{T}^{jet} (GeV);A^{jet}", {HistType::kTH2F, {{500, 0., 100.}, {50, 0., 1.}}});
    std::string h_diel_z_vs_mass_pTInclusive_name = Form("h_diel_z_vs_mass_pTSemiInclusive_%.0f_to_%.0fMeV", 1000*ptBounds.value[0], 1000*ptBounds.value.back());
    auto hDielZVsMasspTInclusiveHist = registry.add(h_diel_z_vs_mass_pTInclusive_name.c_str(), Form("Mass Distributions in z Bins, %.0f < p_{T}^{jet} < %.0f GeV;Mass (GeV);z", ptBounds.value[0], ptBounds.value.back()), {HistType::kTH2F, {{10000, 0., 35.}, {20, 0, 1}}});
    for (size_t i = 0; i < ptBounds.value.size() - 1; i++) { //A 2D histogram for each pT range
      // registry.add(Form("h_diel_z_vs_mass_%s_to_%sGeV", ptBounds->at(i).data(), ptBounds->at(i+1).data()), ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});

      std::string h_diel_z_vs_mass_name = Form("h_diel_z_vs_mass_pT_%.0f_to_%.0fMeV", 1000*ptBounds.value[i], 1000*ptBounds.value[i+1]);
      // std::string h_diel_z_vs_mass_name = std::format("h_diel_z_vs_mass_{:.2f}_to_{:.2f}GeV", ptBounds.value[i], ptBounds.value[i + 1]);
      // registry.add(h_diel_z_vs_mass_name.c_str(), ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});
      // h_diel_z_vs_mass_names.push_back(h_diel_z_vs_mass_name);
      // auto h_diel_z_vs_mass = registry.add(h_diel_z_vs_mass_name.c_str(), ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});
      auto hDielZVsMassHist = registry.add(h_diel_z_vs_mass_name.c_str(), Form("Mass Distributions in z bins, %.0f < p_{T}^{jet} < %.0f GeV;Mass (GeV);z", ptBounds.value[i], ptBounds.value[i+1]), {HistType::kTH2F, {{2000, 0., 35.}, {20, 0, 1}}});
      hDielZVsMassHists.push_back(hDielZVsMassHist);
      // hDielZVsMassHists.push_back(h_diel_z_vs_mass);
      // registry.add(h_diel_z_vs_mass_name, ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});
      // h_diel_z_vs_mass_names.push_back(h_diel_z_vs_mass_name);

      // std::string z_mass_histname ="h_diel_z_vs_mass_" + std::format("{:.2f}", ptBounds->at(i)) + "_to_" + std::format("{:.2f}", ptBounds->at(i + 1)) + "GeV";
      // registry.add(z_mass_histname.c_str(), ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});

      // hP[prefixInd] = registry.add<TH1>(Form("%s/hP", histPrefixes[prefixInd].data()), "#it{p};#it{p} (GeV/#it{c})", HistType::kTH1F, {{500, 0., 6.}});
    }
    registry.add("h_Jet_energy", "Jet Energy;E_{jet};dN/dE_{jet}", {HistType::kTH1F, {{600, 0., 300.}}});
    registry.add("h_diel_jet_distance", "Candidate Angular Distance to Jet Axis;#DeltaR_{cand}^{jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 2.}}});
    registry.add("h_diel_jet_distance_vs_projection", "Candidate p_{T} Ratio vs Angular Distance to Jet Axis;#DeltaR_{cand}^{jet};z", {HistType::kTH2F, {{700, 0., 2.}, {700, 0., 1.}}});
    registry.add("h_diel_otherJet_distance", "Candidate Angular Distance to Same-Event Jets Axis;#DeltaR_{cand}^{jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 5.}}});
    registry.add("h_diel_otherJet_distance_vs_projection", "Candidate p_{T} Ratio vs Angular Distance to Same-Event Jets Axis;#DeltaR_{cand}^{jet};z", {HistType::kTH2F, {{700, 0., 5.}, {700, 0., 1.}}});
    registry.add("h_jet_otherJet_distance", "Candidate-Jet Angular Distance to Same-Event Jets Axis;#DeltaR_{cand}^{jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 5.}}});
    registry.add("h_jet_otherJet_distance_vs_projection", "Candidate p_{T} Ratio vs Candidate-Jet Angular Distance to Same-Event Jets Axis;#DeltaR_{cand}^{jet};z", {HistType::kTH2F, {{700, 0., 5.}, {700, 0., 1.}}});
    registry.add("h_diel_otherJet_Azimutdistance", "Candidate Azimuthal Distance to Same-Event Jets Axis;#Delta#phi_{cand}^{jet};dN/d(#Delta#phi)", {HistType::kTH1F, {{1000, -1., 4.}}});
    registry.add("h_diel_otherJet_Azimutdistance_vs_projection", "Candidate p_{T} Ratio vs Azimuthal Distance to Same-Event Jets Axis;#Delta#phi_{cand}^{jet};z", {HistType::kTH2F, {{700, -1., 4.}, {700, 0., 1.}}});
    registry.add("h_jet_otherJet_Azimutdistance", "Candidate-Jet Azimuthal Distance to Same-Event Jets Axis;#Delta#phi_{cand}^{jet};dN/d(#Delta#phi)", {HistType::kTH1F, {{1000, -1., 4.}}});
    registry.add("h_jet_otherJet_Azimutdistance_vs_projection", "Candidate p_{T} Ratio vs Candidate-Jet Azimuthal Distance to Same-Event Jets Axis;#Delta#phi_{cand}^{jet};z", {HistType::kTH2F, {{700, -1., 4.}, {700, 0., 1.}}});
    std::string h_diel_z_vs_mass_EInclusive_name = Form("h_diel_z_vs_mass_EnergySemiInclusive_%.0f_to_%.0fMeV", 1000*EBounds.value[0], 1000*EBounds.value.back());
    auto hDielZVsMassEInclusiveHist = registry.add(h_diel_z_vs_mass_EInclusive_name.c_str(), Form("Mass Distributions in z Bins, %.0f < E_{jet} < %.0f GeV;Mass (GeV);z", EBounds.value[0], EBounds.value.back()), {HistType::kTH2F, {{10000, 0., 35.}, {20, 0, 1}}});
    for (size_t i = 0; i < EBounds.value.size() - 1; i++) { //A 2D histogram for each E bin, aiming analysis of Xi(E,z)
      std::string h_diel_z_vs_mass_E_name = Form("h_diel_z_vs_mass_Energy_%.0f_to_%.0fMeV", 1000*EBounds.value[i], 1000*EBounds.value[i+1]);
      auto hDielZVsMassEHist = registry.add(h_diel_z_vs_mass_E_name.c_str(), Form("Mass Distributions in z bins, %.0f < E_{jet} < %.0f GeV;Mass (GeV);z", EBounds.value[i], EBounds.value[i+1]), {HistType::kTH2F, {{2000, 0., 35.}, {20, 0, 1}}});
      hDielZVsMassHistsE.push_back(hDielZVsMassEHist);
    }
    hDielZVsMassHists.push_back(hDielZVsMasspTInclusiveHist);
    hDielZVsMassHistsE.push_back(hDielZVsMassEInclusiveHist);
    // std::string z_mass_histname1 ="h_diel_z_vs_mass_" + std::format("{:.2f}", ptBounds->at(0)) + "_to_" + std::format("{:.2f}", ptBounds->at(0 + 1)) + "GeV";
    // registry.add("h_diel_z_vs_mass_5_to_7GeV", ";z; Mass", {HistType::kTH2F, {{10error: no match for 'operator=' (operand types are '__gnu_cxx::__alloc_traits<std::allocator<std::shared_ptr<TH2F> >, std::shared_ptr<TH2F> >::value_type' {aka 'std::shared_ptr<TH2F>'} and 'std::shared_ptr<TH2>')
    // registry.add("h_diel_z_vs_mass_7_to_15GeV", ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});
    // std::string z_mass_histname3 ="h_diel_z_vs_mass_" + std::format("{:.2f}", ptBounds->at(2)) + "_to_" + std::format("{:.2f}", ptBounds->at(2 + 1)) + "GeV";
    // registry.add("h_diel_z_vs_mass_15_to_35GeV", ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});
  }
  int jet_globalid = 0;
  void process(JetCollision const& collision, soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents> const& jets, aod::CandidatesDielectronData const& cands, JetTracks const& tracks) {
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JPsiFragmentationFunctionTask>(cfgc)};
}
