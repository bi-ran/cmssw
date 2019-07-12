// -*- C++ -*-
//
// Package:    TrackAnalyzer
// Class:      TrackAnalyzer
//
/**\class TrackAnalyzer TrackAnalyzer.cc MitHig/TrackAnalyzer/src/TrackAnalyzer.cc

Description: <one line class summary>

Implementation:
Prepare the Track Tree for analysis
*/
//
// Original Author:  Yilmaz Yetkin, Yen-Jie Lee
// Updated: Frank Ma, Matt Nguyen
//         Created:  Tue Sep 30 15:14:28 CEST 2008
// $Id: TrackAnalyzer.cc,v 1.55 2013/06/11 20:58:09 yjlee Exp $

// system include files
#include <vector>
#include <string>
#include <map>

// CMSSW user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// Root include files
#include "TTree.h"

//
// class decleration
//

#define MAXTRACKS 60000

struct TrackEvent {

  // event information
  int nRun;
  int nEv;
  int nLumi;

  int nTrk;

  // Vertex information
  float xVtx;
  float yVtx;
  float zVtx;
  float xVtxErr;
  float yVtxErr;
  float zVtxErr;

  // -- rec tracks --
  int trkCharge[MAXTRACKS];
  float trkPt[MAXTRACKS];
  float trkEta[MAXTRACKS];
  float trkPhi[MAXTRACKS];
  float trkPtError[MAXTRACKS];
  float trkChi2[MAXTRACKS];
  unsigned char trkNHit[MAXTRACKS];
  unsigned char trkNlayer[MAXTRACKS];
  unsigned char trkNdof[MAXTRACKS];
  float trkDz1[MAXTRACKS];
  float trkDzError1[MAXTRACKS];
  float trkDxy1[MAXTRACKS];
  float trkDxyError1[MAXTRACKS];
  unsigned char trkAlgo[MAXTRACKS];
  unsigned char trkOriginalAlgo[MAXTRACKS];
  float trkMVA[MAXTRACKS];

  int pfType[MAXTRACKS];
  float pfCandPt[MAXTRACKS];
  float pfEcal[MAXTRACKS];
  float pfHcal[MAXTRACKS];
};

class TrackAnalyzer : public edm::EDAnalyzer {

  public:
    explicit TrackAnalyzer(const edm::ParameterSet&);
    ~TrackAnalyzer();

  private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    void fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);

    // ----------member data ---------------------------
    double trackPtMin_;

    edm::EDGetTokenT<std::vector<reco::Track>> trackSrc_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
    edm::EDGetTokenT<std::vector<float>> mvaSrc_;
    edm::EDGetTokenT<reco::PFCandidateCollection> pfCandSrc_;

    edm::Service<TFileService> fs;

    TTree* trackTree_;

    TrackEvent p_;
};

//--------------------------------------------------------------------------------------------------
TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig)
{
  trackPtMin_ = iConfig.getParameter<double>("trackPtMin");

  vertexSrc_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"));
  trackSrc_ = consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("trackSrc"));
  mvaSrc_ = consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("mvaSrc"));
  pfCandSrc_ = consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandSrc"));
}

//--------------------------------------------------------------------------------------------------
TrackAnalyzer::~TrackAnalyzer()
{
}

//--------------------------------------------------------------------------------------------------
  void
TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Get tracker geometry
  //  cout <<"StartFill"<<endl;

  p_.nEv = (int)iEvent.id().event();
  p_.nRun = (int)iEvent.id().run();
  p_.nLumi = (int)iEvent.luminosityBlock();

  p_.nTrk = 0;

  fill(iEvent, iSetup);

  trackTree_->Fill();
}

//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::fill(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexSrc_, vertices);

  edm::Handle<std::vector<reco::Track>> tracks;
  iEvent.getByToken(trackSrc_, tracks);

  edm::Handle<std::vector<float>> mvavals;
  iEvent.getByToken(mvaSrc_, mvavals);

  edm::Handle<reco::PFCandidateCollection> pfcands;
  iEvent.getByToken(pfCandSrc_, pfcands);

  reco::Vertex vertex(math::XYZPoint(0, 0, -999), math::Error<3>::type());

  if (!vertices->empty()) {
    vertex = (*vertices)[0];

    p_.xVtx = vertex.position().x();
    p_.yVtx = vertex.position().y();
    p_.zVtx = vertex.position().z();

    p_.xVtxErr = vertex.xError();
    p_.yVtxErr = vertex.yError();
    p_.zVtxErr = vertex.zError();
  }

  auto track_pfcand_map = std::unordered_map<int, int>();
  for (std::size_t i = 0; i < pfcands->size(); ++i) {
    auto const& cand = (*pfcands)[i];

    int type = cand.particleId();
    // only charged hadrons and leptons can be asscociated with a track
    if (!(type == reco::PFCandidate::h ||
          type == reco::PFCandidate::e ||
          type == reco::PFCandidate::mu)
       ) { continue; }

    auto key = cand.trackRef().key();
    track_pfcand_map[key] = i;
  }

  for (uint64_t it = 0; it < tracks->size(); ++it) {
    reco::Track const& trk = (*tracks)[it];

    if (trk.quality(reco::TrackBase::highPurity) != 1) { continue; }
    if (trk.pt() < trackPtMin_) { continue; }

    reco::TrackRef trackref = reco::TrackRef(tracks, it);

    p_.trkCharge[p_.nTrk] = trk.charge();
    p_.trkPt[p_.nTrk] = trk.pt();
    p_.trkEta[p_.nTrk] = trk.eta();
    p_.trkPhi[p_.nTrk] = trk.phi();
    p_.trkPtError[p_.nTrk] = trk.ptError();
    p_.trkChi2[p_.nTrk] = trk.chi2();
    p_.trkNHit[p_.nTrk] = trk.numberOfValidHits();
    p_.trkNdof[p_.nTrk] = trk.ndof();

    p_.trkDz1[p_.nTrk] = trk.dz(vertex.position());
    p_.trkDzError1[p_.nTrk] = std::sqrt(trk.dzError() * trk.dzError()
      + p_.zVtxErr * p_.zVtxErr);
    p_.trkDxy1[p_.nTrk] = trk.dxy(vertex.position());
    p_.trkDxyError1[p_.nTrk] = std::sqrt(trk.dxyError() * trk.dxyError()
      + p_.xVtxErr * p_.yVtxErr);

    p_.trkAlgo[p_.nTrk] = trk.algo();
    p_.trkOriginalAlgo[p_.nTrk] = trk.originalAlgo();

    p_.trkMVA[p_.nTrk] = (trk.algo() == 11) ? 1. : (*mvavals)[it];

    auto index = trackref.key();
    if (track_pfcand_map.find(index) != std::end(track_pfcand_map)) {
      auto const& cand = (*pfcands)[track_pfcand_map[index]];

      p_.pfType[p_.nTrk] = cand.particleId();
      p_.pfCandPt[p_.nTrk] = cand.pt();
      p_.pfEcal[p_.nTrk] = cand.ecalEnergy();
      p_.pfHcal[p_.nTrk] = cand.hcalEnergy();
    } else {
      p_.pfType[p_.nTrk] = -1;
      p_.pfCandPt[p_.nTrk] = -999;
      p_.pfEcal[p_.nTrk] = -999;
      p_.pfHcal[p_.nTrk] = -999;
    }

    p_.nTrk++;
  }
}

// ------------ method called once each job just before starting event loop  ------------
  void
TrackAnalyzer::beginJob()
{
  trackTree_ = fs->make<TTree>("trackTree", "v1");

  // event
  trackTree_->Branch("nEv", &p_.nEv, "nEv/I");
  trackTree_->Branch("nLumi", &p_.nLumi, "nLumi/I");
  trackTree_->Branch("nRun", &p_.nRun, "nRun/I");

  trackTree_->Branch("nTrk", &p_.nTrk, "nTrk/I");

  // vertex
  trackTree_->Branch("xVtx", &p_.xVtx, "xVtx/F");
  trackTree_->Branch("yVtx", &p_.yVtx, "yVtx/F");
  trackTree_->Branch("zVtx", &p_.zVtx, "zVtx/F");
  trackTree_->Branch("xVtxErr", &p_.xVtxErr, "xVtxErr/F");
  trackTree_->Branch("yVtxErr", &p_.yVtxErr, "yVtxErr/F");
  trackTree_->Branch("zVtxErr", &p_.zVtxErr, "zVtxErr/F");

  // Tracks
  trackTree_->Branch("trkCharge", &p_.trkCharge, "trkCharge[nTrk]/I");
  trackTree_->Branch("trkPt", &p_.trkPt, "trkPt[nTrk]/F");
  trackTree_->Branch("trkEta", &p_.trkEta, "trkEta[nTrk]/F");
  trackTree_->Branch("trkPhi", &p_.trkPhi, "trkPhi[nTrk]/F");
  trackTree_->Branch("trkPtError", &p_.trkPtError, "trkPtError[nTrk]/F");
  trackTree_->Branch("trkChi2", &p_.trkChi2, "trkChi2[nTrk]/F");
  trackTree_->Branch("trkNHit", &p_.trkNHit, "trkNHit[nTrk]/b");
  trackTree_->Branch("trkNlayer", &p_.trkNlayer, "trkNlayer[nTrk]/b");
  trackTree_->Branch("trkNdof", &p_.trkNdof, "trkNdof[nTrk]/b");
  trackTree_->Branch("trkAlgo", &p_.trkAlgo, "trkAlgo[nTrk]/b");
  trackTree_->Branch("trkOriginalAlgo", &p_.trkOriginalAlgo, "trkOriginalAlgo[nTrk]/b");
  trackTree_->Branch("trkDxy1", &p_.trkDxy1, "trkDxy1[nTrk]/F");
  trackTree_->Branch("trkDxyError1", &p_.trkDxyError1, "trkDxyError1[nTrk]/F");
  trackTree_->Branch("trkDz1", &p_.trkDz1, "trkDz1[nTrk]/F");
  trackTree_->Branch("trkDzError1", &p_.trkDzError1, "trkDzError1[nTrk]/F");

  trackTree_->Branch("trkMVA", &p_.trkMVA, "trkMVA[nTrk]/F");

  trackTree_->Branch("pfType", &p_.pfType, "pfType[nTrk]/I");
  trackTree_->Branch("pfCandPt", &p_.pfCandPt, "pfCandPt[nTrk]/F");
  trackTree_->Branch("pfEcal", &p_.pfEcal, "pfEcal[nTrk]/F");
  trackTree_->Branch("pfHcal", &p_.pfHcal, "pfHcal[nTrk]/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void
TrackAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);
