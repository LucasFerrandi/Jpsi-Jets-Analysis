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
  // inline static std::vector<std::shared_ptr<TH2F>> hDielZVsMassHists;
  inline static std::vector<o2::framework::HistPtr> hDielZVsMassHists;


  //TEST
  // static std::vector<std::shared_ptr<TH2F>> hDielZVsMassHists;
  // static std::vector<std::string> h_diel_z_vs_mass_names; 
  // static std::vector<decltype(registry.add(""))> hDielZVsMassHists;
  // static std::vector<decltype(std::declval<HistogramRegistry>().add(""))> hDielZVsMassHists;
  // std::array<std::shared_ptr<TH2>, nCharges> hTPCSigvsP;

  void init(InitContext const&) {
    // registry.add("h_coll_id", ";Coll Count", {HistType::kTH1F, {{10000, 0, 10000}}});
    registry.add("h_coll_evtsel", ";Collision Event Selection", {HistType::kTH1F, {{60, 200, 260}}});
    registry.add("h_coll_z", ";Collision Z", {HistType::kTH1F, {{100, -20., 20.}}});
    registry.add("h_coll_mult", ";Collision Multiplicity", {HistType::kTH1F, {{500, -10., 3000.}}});
    // registry.add("h_coll_ambig", ";Collision is_ambiguous", {HistType::kTH1F, {{4, -2., 2.}}});
    registry.add("h_diel_mass", ";Dielectron mass", {HistType::kTH1F, {{10000, 0., 18.}}});
    registry.add("h_diel_pt", ";Dielectron p_T", {HistType::kTH1F, {{1000, 0., 30.}}});
    registry.add("h_diel_eta", ";Dielectron #eta", {HistType::kTH1F, {{200, -2., 2.}}});
    registry.add("h_diel_y", ";Dielectron y", {HistType::kTH1F, {{200, -2., 2.}}});
    registry.add("h_diel_phi", ";Dielectron #phi", {HistType::kTH1F, {{200, -1., 7.}}});
    registry.add("h_diel_sign", ";Dielectron sign", {HistType::kTH1F, {{200, -3, 3}}});
    registry.add("h_Jet_pt", ";Jet p_T", {HistType::kTH1F, {{1000, 0., 150.}}});
    registry.add("h_Jet_eta", ";Jet #eta", {HistType::kTH1F, {{200, -1., 1.}}});
    registry.add("h_Jet_y", ";Jet y", {HistType::kTH1F, {{200, -1., 1.}}});
    registry.add("h_Jet_phi", ";Jet #phi", {HistType::kTH1F, {{200, -1., 7.}}});
    registry.add("h_Jet_area", ";Jet area", {HistType::kTH1F, {{200, 0., 1.}}});
    registry.add("h_dielSize_per_jet", ";Diels per jets", {HistType::kTH1F, {{10, 0, 10}}});
    registry.add("h_jets_per_coll", ";Jets per collision", {HistType::kTH1F, {{10, 0, 10}}});
    registry.add("h_diel_jet_projection", ";Dielectron jet projection", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_diel_jet_rejection", ";Dielectron jet projection", {HistType::kTH1F, {{100, 0., 1.}}});
    registry.add("h_diel_pt_projection", ";Dielectron p_T projection", {HistType::kTH1F, {{100, 0., 1.5}}});
    // registry.add("h_jets_clusters", ";Clusters per jet", {HistType::kTH1F, {{200, 0, 200}}});
    registry.add("h_tracks_per_jets", ";TracksId size per jets; Entries", {HistType::kTH1F, {{40, 0, 40}}});
    registry.add("h_tracks_per_coll", ";Tracks per coll; Entries", {HistType::kTH1F, {{100, 0, 100}}});
    registry.add("h_diels_per_coll", ";Diels per coll; Entries", {HistType::kTH1F, {{10, 0, 10}}});
    registry.add("h_area_jet_vs_pT", ";#p_{T}^{jet} (GeV);#A^{jet}", {HistType::kTH2F, {{500, 0., 100.}, {50, 0., 1.}}});
    registry.add("h_diel_z_vs_mass_pTInclusive", ";Mass (GeV);z", {HistType::kTH2F, {{10000, 0., 35.}, {20, 0, 1}}});
    for (size_t i = 0; i < ptBounds.value.size() - 1; i++) { //A 2D histogram for each pT range
      // registry.add(Form("h_diel_z_vs_mass_%s_to_%sGeV", ptBounds->at(i).data(), ptBounds->at(i+1).data()), ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});


      std::string h_diel_z_vs_mass_name = Form("h_diel_z_vs_mass_%.0f_to_%.0fMeV", 1000*ptBounds.value[i], 1000*ptBounds.value[i+1]);
      // std::string h_diel_z_vs_mass_name = std::format("h_diel_z_vs_mass_{:.2f}_to_{:.2f}GeV", ptBounds.value[i], ptBounds.value[i + 1]);
      // registry.add(h_diel_z_vs_mass_name.c_str(), ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});
      // h_diel_z_vs_mass_names.push_back(h_diel_z_vs_mass_name);
      // auto h_diel_z_vs_mass = registry.add(h_diel_z_vs_mass_name.c_str(), ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});
      auto hDielZVsMassHist = registry.add(h_diel_z_vs_mass_name.c_str(), ";Mass (GeV);z", {HistType::kTH2F, {{2000, 0., 35.}, {20, 0, 1}}});
      hDielZVsMassHists.push_back(hDielZVsMassHist);
      // hDielZVsMassHists.push_back(h_diel_z_vs_mass);
      // registry.add(h_diel_z_vs_mass_name, ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});
      // h_diel_z_vs_mass_names.push_back(h_diel_z_vs_mass_name);

      // std::string z_mass_histname ="h_diel_z_vs_mass_" + std::format("{:.2f}", ptBounds->at(i)) + "_to_" + std::format("{:.2f}", ptBounds->at(i + 1)) + "GeV";
      // registry.add(z_mass_histname.c_str(), ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});

      // hP[prefixInd] = registry.add<TH1>(Form("%s/hP", histPrefixes[prefixInd].data()), "#it{p};#it{p} (GeV/#it{c})", HistType::kTH1F, {{500, 0., 6.}});
    }
    registry.add("h_Jet_energy", ";#E_{jet}", {HistType::kTH1F, {{600, 0., 300.}}});
    registry.add("h_diel_jet_distance", ";#DeltaR_{Diel,jet};dN/d(#DeltaR)", {HistType::kTH1F, {{1000, 0., 2.}}});
    registry.add("h_diel_jet_distance_vs_projection", ";#DeltaR_{Diel,jet};z^{Diel,jet}_{||}", {HistType::kTH2F, {{700, 0., 2.}, {700, 0., 1.}}});
    // registry.add("h_xi", ";#E_{jet}; #Xi(E_{jet}, 0.425) in 8 GeV bins", {HistType::kTH1F, {{10, 0, 10}}});
      
    // std::string z_mass_histname1 ="h_diel_z_vs_mass_" + std::format("{:.2f}", ptBounds->at(0)) + "_to_" + std::format("{:.2f}", ptBounds->at(0 + 1)) + "GeV";
    // registry.add("h_diel_z_vs_mass_5_to_7GeV", ";z; Mass", {HistType::kTH2F, {{10error: no match for 'operator=' (operand types are '__gnu_cxx::__alloc_traits<std::allocator<std::shared_ptr<TH2F> >, std::shared_ptr<TH2F> >::value_type' {aka 'std::shared_ptr<TH2F>'} and 'std::shared_ptr<TH2>')
    // registry.add("h_diel_z_vs_mass_7_to_15GeV", ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});
    // std::string z_mass_histname3 ="h_diel_z_vs_mass_" + std::format("{:.2f}", ptBounds->at(2)) + "_to_" + std::format("{:.2f}", ptBounds->at(2 + 1)) + "GeV";
    // registry.add("h_diel_z_vs_mass_15_to_35GeV", ";z; Mass", {HistType::kTH2F, {{10, 0, 1}, {10000, 0., 35.}}});
  }

  // int coll_id = 0;
  int jet_globalid = 0;
  void process(JetCollision const& collision, soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents> const& jets, aod::CandidatesDielectronData const& cands, JetTracks const& tracks) {
    registry.fill(HIST("h_coll_z"), collision.posZ());
    registry.fill(HIST("h_coll_mult"), collision.multiplicity());
    registry.fill(HIST("h_coll_evtsel"), collision.eventSel());
    // registry.fill(HIST("h_coll_ambig"), collision.isAmbiguous()); //got "no member named 'isAmbiguous'"

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
        if (pt_ratio == 1){
          pt_ratio = 1.4;
        }
        registry.fill(HIST("h_diel_pt_projection"), pt_ratio);
        double axisDistance = jetutilities::deltaR(jet, jpsiCandidate);
        registry.fill(HIST("h_diel_jet_distance_vs_projection"), axisDistance, z_parallel);
        registry.fill(HIST("h_diel_jet_distance"), axisDistance);
      }
      registry.fill(HIST("h_dielSize_per_jet"), jet.candidatesIds().size());
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
            // registry.get<TH2F>(hDielZVsMassHists[i])->Fill(pt_ratio, jpsiCandidate.mass());
            

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
