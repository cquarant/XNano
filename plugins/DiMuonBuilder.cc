#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"

#include <vector>
#include <string>
#include "TLorentzVector.h"

#include "helper.h"
#include "KinVtxFitter.h"

constexpr bool debug = false;

class DiMuonBuilder : public edm::global::EDProducer<> {

public:

  typedef std::vector<pat::Muon> MuonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit DiMuonBuilder(const edm::ParameterSet &cfg):
    l1_selection_{cfg.getParameter<std::string>("lep1Selection")},
    l2_selection_{cfg.getParameter<std::string>("lep2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    src_{consumes<MuonCollection>( cfg.getParameter<edm::InputTag>("src") )},
    ttracks_src_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracksSrc") )},
    beamSpotSrc_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )}{
      produces<pat::CompositeCandidateCollection>("SelectedDiMuons");
    }

  ~DiMuonBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::Muon> l1_selection_;    
  const StringCutObjectSelector<pat::Muon> l2_selection_;    
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_;
  const edm::EDGetTokenT<MuonCollection> src_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
};

void DiMuonBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {
  
  // input
  edm::Handle<MuonCollection> muons;
  evt.getByToken(src_, muons);
  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_src_, ttracks);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  evt.getByToken(beamSpotSrc_, beamSpotHandle);
  const reco::BeamSpot& beamSpot = *beamSpotHandle;

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());

  if(debug){
    std::cout << "Passing l1 selection (pt > 2.5 GeV): ";
    for(size_t l0_idx = 0; l0_idx < muons->size(); ++l0_idx) {
      edm::Ptr<pat::Muon> l0_ptr(muons, l0_idx);
      std::cout << l1_selection_(*l0_ptr) << "  ";
    }
    std::cout << std::endl;
    std::cout << "Passing l2 selection (pt > 2.5 GeV): ";
    for(size_t l0_idx = 0; l0_idx < muons->size(); ++l0_idx) {
      edm::Ptr<pat::Muon> l0_ptr(muons, l0_idx);
      std::cout << l2_selection_(*l0_ptr) << "  ";
    }
    std::cout << std::endl;
  }

  for(size_t l1_idx = 0; l1_idx < muons->size(); ++l1_idx) {
    edm::Ptr<pat::Muon> l1_ptr(muons, l1_idx);
    if(!l1_selection_(*l1_ptr)) continue; 
    
    for(size_t l2_idx = l1_idx + 1; l2_idx < muons->size(); ++l2_idx) {
      edm::Ptr<pat::Muon> l2_ptr(muons, l2_idx);
      if(!l2_selection_(*l2_ptr)) continue;


      // Muons must be different
      if (l1_idx==l2_idx) continue;

      pat::CompositeCandidate muon_pair;
      muon_pair.setP4(l1_ptr->p4() + l2_ptr->p4());
      muon_pair.setCharge(l1_ptr->charge() + l2_ptr->charge());
      muon_pair.addUserFloat("lep_deltaR", reco::deltaR(*l1_ptr, *l2_ptr));
      
      // Put the lepton passing the corresponding selection
      muon_pair.addUserInt("l1_idx", l1_idx );
      muon_pair.addUserInt("l2_idx", l2_idx );


      // Use UserCands as they should not use memory but keep the Ptr itself
      muon_pair.addUserCand("l1", l1_ptr );
      muon_pair.addUserCand("l2", l2_ptr );

      if(debug){
	std::cout << std::endl;
	std::cout << "Muon pair: " << l1_idx << "   " << l2_idx << std::endl;
	std::cout << "Prefit vertex selection [var (cut): observed val] " << std::endl;
	std::cout << "|Vz(l1) - Vz(l2)| (<1): "   << abs(l1_ptr->vz() - l2_ptr->vz()) << std::endl;
	std::cout << "pair mass (in [1,5]): "     << muon_pair.mass() << std::endl;
	std::cout << "charge (==0): "             << muon_pair.charge() << std::endl;
	std::cout << "lep deltaR (>0.02): "       << muon_pair.userFloat("lep_deltaR") << std::endl;
	std::cout << "PASSED: " << pre_vtx_selection_(muon_pair) << std::endl;
      }

      if( !pre_vtx_selection_(muon_pair) ) continue; // before making the SV, cut on the info we have
      
      KinVtxFitter fitter(
			  {ttracks->at(l1_idx), ttracks->at(l2_idx)},
			  {l1_ptr->mass(), l2_ptr->mass()},
			  {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
			  );
      std::cout << "Fitter success: " << fitter.success() << std::endl;
      if ( !fitter.success() ) continue;

      // save quantities after fit needed for selection and to be saved in the final ntuples
      KinematicState fitted_cand = fitter.fitted_candidate();
      RefCountedKinematicVertex fitted_vtx = fitter.fitted_refvtx();

      muon_pair.addUserFloat("sv_prob", fitter.prob());
      muon_pair.addUserFloat("fitted_mass", fitter.success() ? fitted_cand.mass() : -1);

      // cut on the SV info
      if( !post_vtx_selection_(muon_pair) ) continue;


      // DCA between the two muons (used at HLT)
      float DCA = 10.;
      TrajectoryStateClosestToPoint mu1TS = (ttracks->at(l1_idx)).impactPointTSCP();
      TrajectoryStateClosestToPoint mu2TS = (ttracks->at(l2_idx)).impactPointTSCP();
      if (mu1TS.isValid() && mu2TS.isValid()) {
	ClosestApproachInRPhi cApp;
	cApp.calculate(mu1TS.theState(), mu2TS.theState());
	if (cApp.status()) DCA = cApp.distance();
      }
      muon_pair.addUserFloat("DCA", DCA); 

      // Lxy (used at HLT)  
      math::XYZVector pperp(l1_ptr->px() + l2_ptr->px(), l1_ptr->py() + l2_ptr->py(), 0.);
      GlobalError fitted_vtx_err( fitted_vtx->error().cxx(), fitted_vtx->error().cyx(), fitted_vtx->error().cyy(), fitted_vtx->error().czx(), fitted_vtx->error().czy(), fitted_vtx->error().czz() );
      GlobalPoint dispFromBS( -1*( (beamSpot.x0() - fitted_vtx->position().x()) + (fitted_vtx->position().z() - beamSpot.z0()) * beamSpot.dxdz()), -1*((beamSpot.y0() - fitted_vtx->position().y()) + (fitted_vtx->position().z() - beamSpot.z0()) * beamSpot.dydz()), 0);
      float lxy = dispFromBS.perp();
      float lxyerr = sqrt(fitted_vtx_err.rerr(dispFromBS));
      float lxySign = lxy/lxyerr;
      muon_pair.addUserFloat("LxySign", lxySign); 

      // CosAlpha (used at HLT)  
      reco::Vertex::Point vperp(dispFromBS.x(),dispFromBS.y(),0.);
      float cosAlpha = vperp.Dot(pperp)/(vperp.R()*pperp.R());
      muon_pair.addUserFloat("cosAlpha", cosAlpha); 


      // save further quantities, to be saved in the final ntuples: JPsi infos after fit
      TVector3 B_J(fitted_cand.globalMomentum().x(),
		   fitted_cand.globalMomentum().y(),
		   fitted_cand.globalMomentum().z());
      muon_pair.addUserFloat("fitted_pt",  B_J.Pt());
      muon_pair.addUserFloat("fitted_eta", B_J.Eta());
      muon_pair.addUserFloat("fitted_phi", B_J.Phi());

      // save further quantities, to be saved in the final ntuples: JPsi vertex after fit
      muon_pair.addUserFloat("fitted_vtxX",  fitted_vtx->position().x());
      muon_pair.addUserFloat("fitted_vtxY",  fitted_vtx->position().y());
      muon_pair.addUserFloat("fitted_vtxZ",  fitted_vtx->position().z());
      muon_pair.addUserFloat("fitted_vtxEx", fitted_vtx->error().cxx());
      muon_pair.addUserFloat("fitted_vtxEy", fitted_vtx->error().cyy());
      muon_pair.addUserFloat("fitted_vtxEz", fitted_vtx->error().czz());

      // save further quantities, to be saved in the final ntuples: muons before fit
      // Muons post fit are saved only after the very final B fit
      muon_pair.addUserFloat("mu1_pt",  l1_ptr->pt());
      muon_pair.addUserFloat("mu1_eta", l1_ptr->eta());
      muon_pair.addUserFloat("mu1_phi", l1_ptr->phi());
      muon_pair.addUserFloat("mu1_dr",  l1_ptr->userFloat("dr"));
      muon_pair.addUserFloat("mu2_pt",  l2_ptr->pt());
      muon_pair.addUserFloat("mu2_eta", l2_ptr->eta());
      muon_pair.addUserFloat("mu2_phi", l2_ptr->phi());
      muon_pair.addUserFloat("mu2_dr",  l2_ptr->userFloat("dr"));   

      // save further quantities, to be saved in the final ntuples: fired paths

      // 2017 and 2018 interesting triggers
      // muon_pair.addUserInt("mu1_fired_Dimuon25_Jpsi",      l1_ptr->userInt("HLT_Dimuon25_Jpsi"));
      // muon_pair.addUserInt("mu1_fired_Dimuon18_PsiPrime",  l1_ptr->userInt("HLT_Dimuon18_PsiPrime"));
      // muon_pair.addUserInt("mu1_fired_DoubleMu4_JpsiTrk_Displaced",     l1_ptr->userInt("HLT_DoubleMu4_JpsiTrk_Displaced"));
      // muon_pair.addUserInt("mu1_fired_DoubleMu4_PsiPrimeTrk_Displaced", l1_ptr->userInt("HLT_DoubleMu4_PsiPrimeTrk_Displaced"));
      // muon_pair.addUserInt("mu1_fired_DoubleMu4_JpsiTrkTrk_Displaced",  l1_ptr->userInt("HLT_DoubleMu4_JpsiTrkTrk_Displaced"));  

      // All 2018 Charmonium triggers mu1 fired
      // muon_pair.addUserInt("mu1_fired_Dimuon0_Jpsi3p5_Muon2",                  l1_ptr->userInt("HLT_Dimuon0_Jpsi3p5_Muon2"));
      // muon_pair.addUserInt("mu1_fired_Dimuon0_Jpsi_L1_4R_0er1p5R",             l1_ptr->userInt("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R"));
      // muon_pair.addUserInt("mu1_fired_Dimuon0_Jpsi_L1_NoOS",                   l1_ptr->userInt("HLT_Dimuon0_Jpsi_L1_NoOS"));
      // muon_pair.addUserInt("mu1_fired_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", l1_ptr->userInt("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R"));
      // muon_pair.addUserInt("mu1_fired_Dimuon0_Jpsi_NoVertexing_NoOS",          l1_ptr->userInt("HLT_Dimuon0_Jpsi_NoVertexing_NoOS"));
      // muon_pair.addUserInt("mu1_fired_Dimuon0_Jpsi_NoVertexing",               l1_ptr->userInt("HLT_Dimuon0_Jpsi_NoVertexing"));
      // muon_pair.addUserInt("mu1_fired_Dimuon0_Jpsi",                           l1_ptr->userInt("HLT_Dimuon0_Jpsi"));
      // muon_pair.addUserInt("mu1_fired_Dimuon0_LowMass_L1_0er1p5R", l1_ptr->userInt("HLT_Dimuon0_LowMass_L1_0er1p5R"));
      // muon_pair.addUserInt("mu1_fired_Dimuon0_LowMass_L1_0er1p5",  l1_ptr->userInt("HLT_Dimuon0_LowMass_L1_0er1p5"));
      // muon_pair.addUserInt("mu1_fired_Dimuon0_LowMass_L1_4R",      l1_ptr->userInt("HLT_Dimuon0_LowMass_L1_4R"));
      // muon_pair.addUserInt("mu1_fired_Dimuon0_LowMass_L1_4",       l1_ptr->userInt("HLT_Dimuon0_LowMass_L1_4"));
      // muon_pair.addUserInt("mu1_fired_Dimuon0_LowMass",            l1_ptr->userInt("HLT_Dimuon0_LowMass"));
      // muon_pair.addUserInt("mu1_fired_Dimuon10_PsiPrime_Barrel_Seagulls", l1_ptr->userInt("HLT_Dimuon10_PsiPrime_Barrel_Seagulls"));
      // muon_pair.addUserInt("mu1_fired_Dimuon18_PsiPrime_noCorrL1",        l1_ptr->userInt("HLT_Dimuon18_PsiPrime_noCorrL1"));
      // muon_pair.addUserInt("mu1_fired_Dimuon18_PsiPrime",                 l1_ptr->userInt("HLT_Dimuon18_PsiPrime"));
      // muon_pair.addUserInt("mu1_fired_Dimuon20_Jpsi_Barrel_Seagulls", l1_ptr->userInt("HLT_Dimuon20_Jpsi_Barrel_Seagulls"));
      // muon_pair.addUserInt("mu1_fired_Dimuon25_Jpsi_noCorrL1",        l1_ptr->userInt("HLT_Dimuon25_Jpsi_noCorrL1"));
      // muon_pair.addUserInt("mu1_fired_Dimuon25_Jpsi",                 l1_ptr->userInt("HLT_Dimuon25_Jpsi"));
      // muon_pair.addUserInt("mu1_fired_DoubleMu2_Jpsi_DoubleTkMu0_Phi",    l1_ptr->userInt("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi"));
      // muon_pair.addUserInt("mu1_fired_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05", l1_ptr->userInt("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05"));
      // muon_pair.addUserInt("mu1_fired_DoubleMu4_3_Bs",                 l1_ptr->userInt("HLT_DoubleMu4_3_Bs"));
      // muon_pair.addUserInt("mu1_fired_DoubleMu4_3_Jpsi",               l1_ptr->userInt("HLT_DoubleMu4_3_Jpsi"));
      // muon_pair.addUserInt("mu1_fired_DoubleMu4_JpsiTrkTrk_Displaced", l1_ptr->userInt("HLT_DoubleMu4_JpsiTrkTrk_Displaced"));
      // muon_pair.addUserInt("mu1_fired_DoubleMu4_JpsiTrk_Displaced",    l1_ptr->userInt("HLT_DoubleMu4_JpsiTrk_Displaced"));
      // muon_pair.addUserInt("mu1_fired_DoubleMu4_Jpsi_NoVertexing",     l1_ptr->userInt("HLT_DoubleMu4_Jpsi_NoVertexing"));
      // muon_pair.addUserInt("mu1_fired_DoubleMu4_PsiPrimeTrk_Displaced",l1_ptr->userInt("HLT_DoubleMu4_PsiPrimeTrk_Displaced"));
      // muon_pair.addUserInt("mu1_fired_Mu30_TkMu0_Psi",      l1_ptr->userInt("HLT_Mu30_TkMu0_Psi"));
      // muon_pair.addUserInt("mu1_fired_Mu7p5_L2Mu2_Jpsi",    l1_ptr->userInt("HLT_Mu7p5_L2Mu2_Jpsi"));
      // muon_pair.addUserInt("mu1_fired_Mu7p5_Track2_Jpsi",   l1_ptr->userInt("HLT_Mu7p5_Track2_Jpsi"));
      // muon_pair.addUserInt("mu1_fired_Mu7p5_Track3p5_Jpsi", l1_ptr->userInt("HLT_Mu7p5_Track3p5_Jpsi"));
      // muon_pair.addUserInt("mu1_fired_Mu7p5_Track7_Jpsi",   l1_ptr->userInt("HLT_Mu7p5_Track7_Jpsi"));

      // BParking triggers mu1 fired
      muon_pair.addUserInt("mu1_fired_Mu7_IP4",      l1_ptr->userInt("HLT_Mu7_IP4"));
      muon_pair.addUserInt("mu1_fired_Mu8_IP6",      l1_ptr->userInt("HLT_Mu8_IP6"));
      muon_pair.addUserInt("mu1_fired_Mu8_IP5",      l1_ptr->userInt("HLT_Mu8_IP5"));
      muon_pair.addUserInt("mu1_fired_Mu8_IP3",      l1_ptr->userInt("HLT_Mu8_IP3"));
      muon_pair.addUserInt("mu1_fired_Mu8p5_IP3p5",  l1_ptr->userInt("HLT_Mu8p5_IP3p5"));
      muon_pair.addUserInt("mu1_fired_Mu9_IP6",      l1_ptr->userInt("HLT_Mu9_IP6"));
      muon_pair.addUserInt("mu1_fired_Mu9_IP5",      l1_ptr->userInt("HLT_Mu9_IP5"));
      muon_pair.addUserInt("mu1_fired_Mu9_IP4",      l1_ptr->userInt("HLT_Mu9_IP4"));
      muon_pair.addUserInt("mu1_fired_Mu10p5_IP3p5", l1_ptr->userInt("HLT_Mu10p5_IP3p5"));
      muon_pair.addUserInt("mu1_fired_Mu12_IP6",     l1_ptr->userInt("HLT_Mu12_IP6"));

      // 2017 and 2018 interesting Jpsi triggers mu2 fired
      // muon_pair.addUserInt("mu2_fired_Dimuon25_Jpsi",      l2_ptr->userInt("HLT_Dimuon25_Jpsi"));
      // muon_pair.addUserInt("mu2_fired_Dimuon18_PsiPrime",  l2_ptr->userInt("HLT_Dimuon18_PsiPrime"));
      // muon_pair.addUserInt("mu2_fired_DoubleMu4_JpsiTrk_Displaced",     l2_ptr->userInt("HLT_DoubleMu4_JpsiTrk_Displaced"));
      // muon_pair.addUserInt("mu2_fired_DoubleMu4_PsiPrimeTrk_Displaced", l2_ptr->userInt("HLT_DoubleMu4_PsiPrimeTrk_Displaced"));
      // muon_pair.addUserInt("mu2_fired_DoubleMu4_JpsiTrkTrk_Displaced",  l2_ptr->userInt("HLT_DoubleMu4_JpsiTrkTrk_Displaced"));   

      // All 2018 Charmonium triggers mu2 fired
      // muon_pair.addUserInt("mu2_fired_Dimuon0_Jpsi3p5_Muon2",                  l2_ptr->userInt("HLT_Dimuon0_Jpsi3p5_Muon2"));
      // muon_pair.addUserInt("mu2_fired_Dimuon0_Jpsi_L1_4R_0er1p5R",             l2_ptr->userInt("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R"));
      // muon_pair.addUserInt("mu2_fired_Dimuon0_Jpsi_L1_NoOS",                   l2_ptr->userInt("HLT_Dimuon0_Jpsi_L1_NoOS"));
      // muon_pair.addUserInt("mu2_fired_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", l2_ptr->userInt("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R"));
      // muon_pair.addUserInt("mu2_fired_Dimuon0_Jpsi_NoVertexing_NoOS",          l2_ptr->userInt("HLT_Dimuon0_Jpsi_NoVertexing_NoOS"));
      // muon_pair.addUserInt("mu2_fired_Dimuon0_Jpsi_NoVertexing",               l2_ptr->userInt("HLT_Dimuon0_Jpsi_NoVertexing"));
      // muon_pair.addUserInt("mu2_fired_Dimuon0_Jpsi",                           l2_ptr->userInt("HLT_Dimuon0_Jpsi"));
      // muon_pair.addUserInt("mu2_fired_Dimuon0_LowMass_L1_0er1p5R", l2_ptr->userInt("HLT_Dimuon0_LowMass_L1_0er1p5R"));
      // muon_pair.addUserInt("mu2_fired_Dimuon0_LowMass_L1_0er1p5",  l2_ptr->userInt("HLT_Dimuon0_LowMass_L1_0er1p5"));
      // muon_pair.addUserInt("mu2_fired_Dimuon0_LowMass_L1_4R",      l2_ptr->userInt("HLT_Dimuon0_LowMass_L1_4R"));
      // muon_pair.addUserInt("mu2_fired_Dimuon0_LowMass_L1_4",       l2_ptr->userInt("HLT_Dimuon0_LowMass_L1_4"));
      // muon_pair.addUserInt("mu2_fired_Dimuon0_LowMass",            l2_ptr->userInt("HLT_Dimuon0_LowMass"));
      // muon_pair.addUserInt("mu2_fired_Dimuon10_PsiPrime_Barrel_Seagulls", l2_ptr->userInt("HLT_Dimuon10_PsiPrime_Barrel_Seagulls"));
      // muon_pair.addUserInt("mu2_fired_Dimuon18_PsiPrime_noCorrL1",        l2_ptr->userInt("HLT_Dimuon18_PsiPrime_noCorrL1"));
      // muon_pair.addUserInt("mu2_fired_Dimuon18_PsiPrime",                 l2_ptr->userInt("HLT_Dimuon18_PsiPrime"));
      // muon_pair.addUserInt("mu2_fired_Dimuon20_Jpsi_Barrel_Seagulls", l2_ptr->userInt("HLT_Dimuon20_Jpsi_Barrel_Seagulls"));
      // muon_pair.addUserInt("mu2_fired_Dimuon25_Jpsi_noCorrL1",        l2_ptr->userInt("HLT_Dimuon25_Jpsi_noCorrL1"));
      // muon_pair.addUserInt("mu2_fired_Dimuon25_Jpsi",                 l2_ptr->userInt("HLT_Dimuon25_Jpsi"));
      // muon_pair.addUserInt("mu2_fired_DoubleMu2_Jpsi_DoubleTkMu0_Phi",    l2_ptr->userInt("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi"));
      // muon_pair.addUserInt("mu2_fired_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05", l2_ptr->userInt("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05"));
      // muon_pair.addUserInt("mu2_fired_DoubleMu4_3_Bs",                 l2_ptr->userInt("HLT_DoubleMu4_3_Bs"));
      // muon_pair.addUserInt("mu2_fired_DoubleMu4_3_Jpsi",               l2_ptr->userInt("HLT_DoubleMu4_3_Jpsi"));
      // muon_pair.addUserInt("mu2_fired_DoubleMu4_JpsiTrkTrk_Displaced", l2_ptr->userInt("HLT_DoubleMu4_JpsiTrkTrk_Displaced"));
      // muon_pair.addUserInt("mu2_fired_DoubleMu4_JpsiTrk_Displaced",    l2_ptr->userInt("HLT_DoubleMu4_JpsiTrk_Displaced"));
      // muon_pair.addUserInt("mu2_fired_DoubleMu4_Jpsi_NoVertexing",     l2_ptr->userInt("HLT_DoubleMu4_Jpsi_NoVertexing"));
      // muon_pair.addUserInt("mu2_fired_DoubleMu4_PsiPrimeTrk_Displaced",l2_ptr->userInt("HLT_DoubleMu4_PsiPrimeTrk_Displaced"));
      // muon_pair.addUserInt("mu2_fired_Mu30_TkMu0_Psi",      l2_ptr->userInt("HLT_Mu30_TkMu0_Psi"));
      // muon_pair.addUserInt("mu2_fired_Mu7p5_L2Mu2_Jpsi",    l2_ptr->userInt("HLT_Mu7p5_L2Mu2_Jpsi"));
      // muon_pair.addUserInt("mu2_fired_Mu7p5_Track2_Jpsi",   l2_ptr->userInt("HLT_Mu7p5_Track2_Jpsi"));
      // muon_pair.addUserInt("mu2_fired_Mu7p5_Track3p5_Jpsi", l2_ptr->userInt("HLT_Mu7p5_Track3p5_Jpsi"));
      // muon_pair.addUserInt("mu2_fired_Mu7p5_Track7_Jpsi",   l2_ptr->userInt("HLT_Mu7p5_Track7_Jpsi"));

      // BParking triggers mu2 fired
      muon_pair.addUserInt("mu2_fired_Mu7_IP4",      l2_ptr->userInt("HLT_Mu7_IP4"));
      muon_pair.addUserInt("mu2_fired_Mu8_IP6",      l2_ptr->userInt("HLT_Mu8_IP6"));
      muon_pair.addUserInt("mu2_fired_Mu8_IP5",      l2_ptr->userInt("HLT_Mu8_IP5"));
      muon_pair.addUserInt("mu2_fired_Mu8_IP3",      l2_ptr->userInt("HLT_Mu8_IP3"));
      muon_pair.addUserInt("mu2_fired_Mu8p5_IP3p5",  l2_ptr->userInt("HLT_Mu8p5_IP3p5"));
      muon_pair.addUserInt("mu2_fired_Mu9_IP6",      l2_ptr->userInt("HLT_Mu9_IP6"));
      muon_pair.addUserInt("mu2_fired_Mu9_IP5",      l2_ptr->userInt("HLT_Mu9_IP5"));
      muon_pair.addUserInt("mu2_fired_Mu9_IP4",      l2_ptr->userInt("HLT_Mu9_IP4"));
      muon_pair.addUserInt("mu2_fired_Mu10p5_IP3p5", l2_ptr->userInt("HLT_Mu10p5_IP3p5"));
      muon_pair.addUserInt("mu2_fired_Mu12_IP6",     l2_ptr->userInt("HLT_Mu12_IP6"));

      // push in the event
      ret_value->push_back(muon_pair);
    }
  }

  evt.put(std::move(ret_value),  "SelectedDiMuons");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DiMuonBuilder);
