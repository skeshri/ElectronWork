// -*- C++ -*-
//
// Package:    ElectronWork/ElectronNtupler
// Class:      PhotonNtuplerAOD
// 
/**\class PhotonNtuplerAOD PhotonNtuplerAOD.cc ElectronWork/PhotonNtuplerAOD/plugins/PhotonNtuplerAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ilya Kravchenko
//         Created:  Thu, 10 Jul 2014 09:54:13 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"

// #include "RecoEgamma/PhotonIdentification/interface/GEDPhoIDTools.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"



//
// class declaration
//

class PhotonNtuplerAOD : public edm::EDAnalyzer {
   public:
      explicit PhotonNtuplerAOD(const edm::ParameterSet&);
      ~PhotonNtuplerAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<edm::View<reco::Photon> > photonCollectionToken_;
      edm::EDGetTokenT<double> rhoToken_;

  TTree *photonTree_;

  Int_t nPV_;        // number of reconsrtucted primary vertices
  Float_t rho_;      // the rho variable


  // all photon variables
  Int_t nPhotons_;

  std::vector<Float_t> pt_;
  std::vector<Float_t> eta_;
  std::vector<Float_t> phi_;

  std::vector<Float_t> full5x5_sigmaIetaIeta_;
  std::vector<Float_t> hOverE_;
  std::vector<Float_t> hasPixelSeed_;

  std::vector<Float_t> isoChargedHadrons_;
  std::vector<Float_t> isoNeutralHadrons_;
  std::vector<Float_t> isoPhotons_;



};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PhotonNtuplerAOD::PhotonNtuplerAOD(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  photonCollectionToken_(consumes<edm::View<reco::Photon> >(iConfig.getParameter<edm::InputTag>("photons"))),
  rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
{

  edm::Service<TFileService> fs;
  photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");
  
  photonTree_->Branch("nPV"        ,  &nPV_     , "nPV/I");
  photonTree_->Branch("rho"        ,  &rho_ , "rho/F");

  photonTree_->Branch("nPho",  &nPhotons_ , "nPho/I");
  photonTree_->Branch("pt"  ,  &pt_    );
  photonTree_->Branch("eta" ,  &eta_ );
  photonTree_->Branch("phi" ,  &phi_ );

  photonTree_->Branch("hOverE",  &hOverE_);
  photonTree_->Branch("hasPixelSeed"           ,  &hasPixelSeed_);
  photonTree_->Branch("full5x5_sigmaIetaIeta"  , &full5x5_sigmaIetaIeta_);
  photonTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_);
  photonTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_);
  photonTree_->Branch("isoPhotons"             , &isoPhotons_);
 
}


PhotonNtuplerAOD::~PhotonNtuplerAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhotonNtuplerAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  
  // // An object needed for isolation calculations
  // GEDPhoIDTools *GEDIdTool = new GEDPhoIDTools(iEvent);

  // Get photon collection
  edm::Handle<edm::View<reco::Photon> > collection;
  iEvent.getByToken(photonCollectionToken_, collection);

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  //const reco::Vertex &pv = vertices->front();
  nPV_    = vertices->size();

  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  for (VertexCollection::const_iterator vtx = vertices->begin(); 
       vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    if (  /*!vtx->isFake() &&*/ 
        !(vtx->chi2()==0 && vtx->ndof()==0) 
        &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
        && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }

  if ( firstGoodVertex==vertices->end() )
    return; // skip event if there are no good PVs

  // Get rho
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  // Clear vectors
  nPhotons_ = 0;
  pt_.clear();
  eta_.clear();
  phi_.clear();

  // Loop over photons
  
  const auto& pho_refs = collection->refVector();
  
  for( const auto& pho : pho_refs ) {
    
    // Kinematics
    if( pho->pt() < 10 ) 
      continue;
    
    nPhotons_++;
    pt_  .push_back( pho->pt() );
    eta_ .push_back( pho->superCluster()->eta() );
    phi_ .push_back( pho->superCluster()->phi() );

    full5x5_sigmaIetaIeta_ .push_back( pho->full5x5_sigmaIetaIeta() );
    hOverE_                .push_back( pho->hadTowOverEm() );
    hasPixelSeed_          .push_back( pho->hasPixelSeed() );

    // reco::PhotonRef recophoRef(collection, pho);
    // GEDIdTool->setPhotonP4(recophoRef, firstGoodVertex);
    // isoChargedHadrons_ .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::h));
    // isoPhotons_        .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::gamma));
    // isoNeutralHadrons_ .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::h0));

    isoChargedHadrons_ .push_back(-999);
    isoPhotons_        .push_back(-999);
    isoNeutralHadrons_ .push_back(-999);


   }
   
  // Save the info
  photonTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
PhotonNtuplerAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonNtuplerAOD::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
PhotonNtuplerAOD::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
PhotonNtuplerAOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
PhotonNtuplerAOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
PhotonNtuplerAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonNtuplerAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonNtuplerAOD);
