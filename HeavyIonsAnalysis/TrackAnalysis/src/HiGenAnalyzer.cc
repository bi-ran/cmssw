// -*- C++ -*-
//
// Package:    HiGenAnalyzer
// Class:      HiGenAnalyzer
//
/**\class HiGenAnalyzer HiGenAnalyzer.cc

   Description: Analyzer that studies (HI) gen event info

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Yetkin Yilmaz, Frank Ma
//         Created:  Tue Dec 18 09:44:41 EST 2007
// $Id: HiGenAnalyzer.cc,v 1.9 2012/09/28 20:10:52 yilmaz Exp $
//
//


// system include files
#include <string>
#include <fstream>

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

// root include file
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include <vector>

//
// class decleration
//

struct GenEvent {
  int event;

  int mult;

  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<int> pdg;
  std::vector<int> chg;
  std::vector<int> sube;
};

class HiGenAnalyzer : public edm::EDAnalyzer {
public:
  explicit HiGenAnalyzer(const edm::ParameterSet&);
  ~HiGenAnalyzer();

private:
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  TTree* t_;
  GenEvent p_;

  Double_t etaMax_;
  Double_t ptMin_;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticleSrc_;

  edm::Service<TFileService> f;
};

//
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

HiGenAnalyzer::HiGenAnalyzer(const edm::ParameterSet& iConfig)
{
  etaMax_ = iConfig.getUntrackedParameter<Double_t>("etaMax");
  ptMin_ = iConfig.getUntrackedParameter<Double_t>("ptMin");

  genParticleSrc_ = consumes<reco::GenParticleCollection>(
     iConfig.getParameter<edm::InputTag>("genParticleSrc"));
}

HiGenAnalyzer::~HiGenAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called to for each event  ------------
void
HiGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  p_.pt.clear();
  p_.eta.clear();
  p_.phi.clear();
  p_.pdg.clear();
  p_.chg.clear();
  p_.sube.clear();

  p_.event = iEvent.id().event();
  p_.mult = 0;

  edm::Handle<reco::GenParticleCollection> particles;
  iEvent.getByToken(genParticleSrc_, particles);

  for (auto const& p : *particles) {
     if (p.pt() < ptMin_) { continue; }
     if (std::abs(p.eta()) > etaMax_) { continue; }
     if (p.status() != 1) { continue; }

     p_.pt.push_back(p.pt());
     p_.eta.push_back(p.eta());
     p_.phi.push_back(p.phi());
     p_.pdg.push_back(p.pdgId());
     p_.chg.push_back(p.charge());
     p_.sube.push_back(p.collisionId());

     ++(p_.mult);
  }

  t_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
HiGenAnalyzer::beginRun(const edm::Run&, const edm::EventSetup& iSetup)
{
}

void
HiGenAnalyzer::beginJob()
{
  t_ = f->make<TTree>("hi", "");

  t_->Branch("event", &p_.event, "event/I");

  t_->Branch("mult", &p_.mult, "mult/I");

  t_->Branch("pt", &p_.pt);
  t_->Branch("eta", &p_.eta);
  t_->Branch("phi", &p_.phi);
  t_->Branch("pdg", &p_.pdg);
  t_->Branch("chg", &p_.chg);
  t_->Branch("sube", &p_.sube);
}

// ------------ method called once each job just after ending the event loop  ------------
void
HiGenAnalyzer::endJob() { }

//define this as a plug-in
DEFINE_FWK_MODULE(HiGenAnalyzer);
