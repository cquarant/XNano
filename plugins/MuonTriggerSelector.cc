// class to produce 2 pat::MuonCollections

// one matched to the selected triggers
// another fitered wrt the selected triggers

// chiara: 
// for the moment track purity is not required. To be activated when switching to 10.6.x

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <TLorentzVector.h>
#include "helper.h"

using namespace std;

constexpr bool debug = false;

class MuonTriggerSelector : public edm::EDProducer {
  
public:
    
  explicit MuonTriggerSelector(const edm::ParameterSet &iConfig);
    
  ~MuonTriggerSelector() override {};
  
private:
  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
  edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
  
  
  // for trigger match
  std::vector<std::string> HLTPaths_;
  const double drForTriggerMatch_;
  
  // Offline muons selection  
  const double ptMin_;           
  const double absEtaMax_;       
};


MuonTriggerSelector::MuonTriggerSelector(const edm::ParameterSet &iConfig):
  muonSrc_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  //
  beamSpotSrc_(consumes<reco::BeamSpot>( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),
  HLTPaths_(iConfig.getParameter<std::vector<std::string>>("HLTPaths")),
  drForTriggerMatch_(iConfig.getParameter<double>("drForTriggerMatch")),    
  //
  ptMin_(iConfig.getParameter<double>("ptMin")),
  absEtaMax_(iConfig.getParameter<double>("absEtaMax"))
{
  // produce the SelectedMuons collection (all muons passing the preselection)
  produces<pat::MuonCollection>("SelectedMuons");
  produces<TransientTrackCollection>("SelectedTransientMuons");  
  // produce the triggering muons collection
  produces<pat::MuonCollection>("trgMuons");
}

void MuonTriggerSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
  const reco::BeamSpot& beamSpot = *beamSpotHandle;

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonSrc_, muons);

  // Outputs
  std::unique_ptr<pat::MuonCollection>      muons_out      ( new pat::MuonCollection );
  std::unique_ptr<TransientTrackCollection> trans_muons_out( new TransientTrackCollection );
  std::unique_ptr<pat::MuonCollection>      trgmuons_out   ( new pat::MuonCollection );
  
  std::vector<int> muonIsTrigger(muons->size(), 0);
  std::vector<std::vector<int>> fires;
  std::vector<std::vector<float>> matcher; 
  std::vector<std::vector<float>> DR;
  

  // ----------------------------------------------------------
  // debug
  /*
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl;
  std::cout << "\n == TRIGGER PATHS= " << std::endl;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    if (triggerBits->accept(i)) 
      std::cout << "Trigger " << names.triggerName(i) 
		<< ": Pass = " << (triggerBits->accept(i)) 
		<< ", Was Run = " << (triggerBits->wasrun(i))
		<< std::endl;
  }

  std::cout << std::endl;
  std::cout << "\n TRIGGER OBJECTS " << std::endl;
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    obj.unpackPathNames(names);
    std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
    // Print trigger object collection and type
    std::cout << "\t   Collection: " << obj.collection() << std::endl;
    std::cout << "\t   Type IDs:   ";
    std::cout << std::endl;
    std::vector<string> pathNamesAll  = obj.pathNames(false);
    std::vector<string> pathNamesLast = obj.pathNames(true);
    // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
    // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
    // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
    std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
    for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
      bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
      bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );      //false = last; true = L3 accepted ==> isL3 always true
      bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
      bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
      std::cout << "   " << pathNamesAll[h];
      if (isBoth) std::cout << "(L,3)";
      if (isL3 && !isBoth) std::cout << "(*,3)";
      if (isLF && !isBoth) std::cout << "(L,*)";
      if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
    }
    std::cout << std::endl;
  }
  // debug end
  */
  // ----------------------------------------------------------


  // Loop over reconstructed muons
  for(const pat::Muon &muon : *muons){
    
    if(debug) std::cout << "Slimmed muons: muon Pt = " << muon.pt() 
			<< " Eta = " << muon.eta() << " Phi = " << muon.phi()  <<endl;
    
    // These vectors have one entry per HLT path
    std::vector<int> frs(HLTPaths_.size(),0);              
    std::vector<float> temp_matched_to(HLTPaths_.size(),1000.);
    std::vector<float> temp_DR(HLTPaths_.size(),1000.);

    // Loop over trigger paths
    int ipath=-1;
    for (const std::string path: HLTPaths_){
      
      if(debug) std::cout << "ipath = " << ipath << ", path = " << path << std::endl;
      ipath++;
      
      // Here we loop over trigger objects
      float minDr = 1000.;
      float minPt = 1000.;
      
      const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
	obj.unpackPathNames(triggerNames);
	if(debug) std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;

	// Is this HLT object matching the muon? 
	bool hltMatchOffline = false;
	TVector3 p3triggerObj;
	p3triggerObj.SetPtEtaPhi(obj.pt(), obj.eta(), obj.phi());  
	TVector3 p3Muon;                                           
	p3Muon.SetPtEtaPhi(muon.pt(), muon.eta(), muon.phi());    
	float dr = p3triggerObj.DeltaR(p3Muon);
	if (dr<drForTriggerMatch_) hltMatchOffline = true;

	// Is this HLT object matching this path?
	bool hltMatchThisPath = false;
	if (hltMatchOffline==true) { 
	  char cstr[ (path+"*").size() + 1]; 
	  strcpy( cstr, (path+"*").c_str() ); 
	  bool isBoth = obj.hasPathName( cstr, true, true );
	  if (isBoth) {
	    if(debug)  std::cout << "Matching path " << path << std::endl;;
	    hltMatchThisPath = true;
	  }
	}

	// Look for minDR between reco muon and matched HLT object for this path
	if (hltMatchThisPath==true && hltMatchOffline==true) {
	  frs[ipath]=1;
	  if (dr<minDr) {
	    minDr = dr;
	    minPt = obj.pt();
	  }
	}
	
      } // Loop over trigger objects
    
      // Here we store the minimum between reco muon and all its matched HLT objects for this HLT path
      temp_DR[ipath]=minDr;
      temp_matched_to[ipath]=minPt;
      
    } // Loop over paths
    
    // One vector per muon. Each vector : one element per path (corresponding to the closest HLT object
    fires.push_back(frs);                 // This is used in order to see if a reco muon fired a Trigger (1) or not (0).
    matcher.push_back(temp_matched_to);   // PT of the reco muon matched to HLT object
    DR.push_back(temp_DR);
    
  } // Loop over reco muons


  // Now check for different reco muons that are matched to the same HLTObject.
  for(unsigned int path=0; path<HLTPaths_.size(); path++){
    for(unsigned int iMuo=0; iMuo<muons->size(); iMuo++){
      for(unsigned int im=(iMuo+1); im<muons->size(); im++){
	if(matcher[iMuo][path]!=1000. && matcher[iMuo][path]==matcher[im][path]){
	  if(DR[iMuo][path]<DR[im][path]){     // Keep the one that has the minimum DR with the HLT object
	    fires[im][path]=0;
	    matcher[im][path]=1000.;
	    DR[im][path]=1000.;                       
	  }
	  else{
	    fires[iMuo][path]=0;
	    matcher[iMuo][path]=1000.;
	    DR[iMuo][path]=1000.;                       
	  }
	}              
      }
      if(matcher[iMuo][path]!=1000.) muonIsTrigger[iMuo]=1;
    }
  }
  // 

  if(debug) std::cout << "number of Muons=" <<muons->size() << endl;

  // Create a collection with all trg muons
  for(const pat::Muon & muon : *muons){
    unsigned int iMuo(&muon -&(muons->at(0)));
    if(muonIsTrigger[iMuo]==1){
      pat::Muon recoTriggerMuonCand(muon);
      trgmuons_out->emplace_back(recoTriggerMuonCand);
    }
  }

  // Save all reco muons passing the selection
  for(const pat::Muon & muon : *muons){
    unsigned int iMuo(&muon - &(muons->at(0)) );
    
    if( muon.pt()<ptMin_ ) continue;
    if( fabs(muon.eta())>absEtaMax_ ) continue;
    // if( !((muon.track())->highPurity()) ) continue;      // chiara: to-be-done
    // if (!muon.passed(13)) continue;   // 13 = "SoftCutBasedId"; from DataFormats/MuonReco/interface/Muon.h. Boh, does not work, to be fixed

    const reco::TransientTrack muonTT((*(muon.bestTrack())),&(*bFieldHandle));  
    if(!muonTT.isValid()) continue; 
    
    muons_out->emplace_back(muon);

    int isPFcand = (int) muon.isPFMuon();
    int isGlobal = (int) muon.isGlobalMuon();
    int isTracker = (int) muon.isTrackerMuon();
    int isLoose = (int)muon.isLooseMuon();
    muons_out->back().addUserInt("isPFcand", isPFcand);    
    muons_out->back().addUserInt("isGlobal", isGlobal);    
    muons_out->back().addUserInt("isTracker", isTracker);    
    muons_out->back().addUserInt("looseId", isLoose);
    muons_out->back().addUserInt("charge", muon.charge());

    // dr cut (same quantity as in HLTMuonDimuonL3Filter, to emulate HLT)
    float mudr = fabs( (- (muon.vx()-beamSpot.x0()) * muon.py() + (muon.vy()-beamSpot.y0()) * muon.px() ) / muon.pt() );
    muons_out->back().addUserFloat("dr", mudr);  

    for(unsigned int i=0; i<HLTPaths_.size(); i++){
      muons_out->back().addUserInt(HLTPaths_[i],fires[iMuo][i]);
    }
    trans_muons_out->emplace_back(muonTT);
  }

  iEvent.put(std::move(muons_out),       "SelectedMuons");
  iEvent.put(std::move(trans_muons_out), "SelectedTransientMuons");
  iEvent.put(std::move(trgmuons_out),    "trgMuons");
}



DEFINE_FWK_MODULE(MuonTriggerSelector);
