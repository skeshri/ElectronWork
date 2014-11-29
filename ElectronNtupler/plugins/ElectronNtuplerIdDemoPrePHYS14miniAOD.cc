// -*- C++ -*-
//
// Package:    ElectronWork/ElectronNtupler
// Class:      ElectronNtuplerIdDemoPrePHYS14miniAOD
// 
/**\class ElectronNtuplerIdDemoPrePHYS14miniAOD ElectronNtuplerIdDemoPrePHYS14miniAOD.cc ElectronWork/ElectronNtuplerIdDemoPrePHYS14miniAOD/plugins/ElectronNtuplerIdDemoPrePHYS14miniAOD.cc

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

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"



//
// class declaration
//

class ElectronNtuplerIdDemoPrePHYS14miniAOD : public edm::EDAnalyzer {
   public:
      explicit ElectronNtuplerIdDemoPrePHYS14miniAOD(const edm::ParameterSet&);
      ~ElectronNtuplerIdDemoPrePHYS14miniAOD();

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
      edm::EDGetTokenT<edm::View<reco::GsfElectron> > electronCollectionToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;

  TTree *electronTree_;

  // all electron variables
  Int_t nElectrons_;

  std::vector<Float_t> pt_;
  std::vector<Float_t> etaSC_;
  std::vector<Float_t> phiSC_;

  std::vector<Int_t>   passVetoId_;     
  std::vector<Int_t>   passTightId_;     

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
ElectronNtuplerIdDemoPrePHYS14miniAOD::ElectronNtuplerIdDemoPrePHYS14miniAOD(const edm::ParameterSet& iConfig):
  electronCollectionToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
  electronVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"))),
  electronTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap")))
{

  edm::Service<TFileService> fs;
  electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");
  
  electronTree_->Branch("nEle"    ,  &nElectrons_ , "nEle/I");

  electronTree_->Branch("pt"    ,  &pt_    );
  electronTree_->Branch("etaSC" ,  &etaSC_ );
  electronTree_->Branch("phiSC" ,  &phiSC_ );

  electronTree_->Branch("passVetoId", &passVetoId_);
  electronTree_->Branch("passTightId", &passTightId_);
 
}


ElectronNtuplerIdDemoPrePHYS14miniAOD::~ElectronNtuplerIdDemoPrePHYS14miniAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronNtuplerIdDemoPrePHYS14miniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  // Get the electron collection
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByToken(electronCollectionToken_, electrons);

  // Get the electron ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);

  // Loop over electrons
  nElectrons_ = 0;
  pt_.clear();
  etaSC_.clear();
  phiSC_.clear();

  passVetoId_.clear();     
  passTightId_.clear();     
  
  // Loop over electrons
  // for (const pat::Electron &el : *electrons) {    
  const auto& ele_refs = electrons->refVector();
  for( const auto& el : ele_refs ) {

    // Kinematics
    if( el->pt() < 10 ) 
      continue;

    nElectrons_++;
    pt_.push_back( el->pt() );
    etaSC_.push_back( el->superCluster()->eta() );
    phiSC_.push_back( el->superCluster()->phi() );
     
    // Look up the ID decision for this electron in 
    // the ValueMap object and store it
    bool isPassVeto  = (*veto_id_decisions)[el];
    bool isPassTight = (*tight_id_decisions)[el];
    passVetoId_.push_back( isPassVeto );
    passTightId_.push_back( isPassTight );

  }
   
  // Save the info for this event
  electronTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
ElectronNtuplerIdDemoPrePHYS14miniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronNtuplerIdDemoPrePHYS14miniAOD::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
ElectronNtuplerIdDemoPrePHYS14miniAOD::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ElectronNtuplerIdDemoPrePHYS14miniAOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ElectronNtuplerIdDemoPrePHYS14miniAOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ElectronNtuplerIdDemoPrePHYS14miniAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronNtuplerIdDemoPrePHYS14miniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronNtuplerIdDemoPrePHYS14miniAOD);
