/*
  Based on the jet response analyzer
  Modified by Matt Nguyen, November 2010

*/

#include "HeavyIonsAnalysis/JetAnalysis/interface/HiInclusiveJetAnalyzer.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace edm;
using namespace reco;

HiInclusiveJetAnalyzer::HiInclusiveJetAnalyzer(const edm::ParameterSet& iConfig)
{
  jetTagLabel_ = iConfig.getParameter<InputTag>("jetTag");
  jetTag_ = consumes<reco::JetView>(jetTagLabel_);
  jetTagPat_ = consumes<pat::JetCollection>(jetTagLabel_);
  
  pfCandidateLabel_ = consumes<reco::PFCandidateCollection>(
    iConfig.getUntrackedParameter<edm::InputTag>("pfCandidateLabel"));

  jetPtMin_ = iConfig.getParameter<double>("jetPtMin");
  jetAbsEtaMax_ = iConfig.getUntrackedParameter<double>("jetAbsEtaMax", 5.0);

  isMC_ = iConfig.getUntrackedParameter<bool>("isMC", false);

  if (isMC_) {
    genjetTag_ = consumes<edm::View<reco::GenJet>>(
      iConfig.getParameter<InputTag>("genjetTag"));
    genPtMin_ = iConfig.getUntrackedParameter<double>("genPtMin",10);

    eventGenInfoTag_ = consumes<GenEventInfoProduct>(
      iConfig.getParameter<InputTag>("eventInfoTag"));
  }
}

HiInclusiveJetAnalyzer::~HiInclusiveJetAnalyzer() { }

void
HiInclusiveJetAnalyzer::beginRun(const edm::Run& run,
                                 const edm::EventSetup & es) {}

void
HiInclusiveJetAnalyzer::beginJob() {
  t = f_->make<TTree>("t", "jets");

  t->Branch("run", &jets_.run, "run/I");
  t->Branch("evt", &jets_.evt, "evt/I");
  t->Branch("lumi", &jets_.lumi, "lumi/I");

  t->Branch("nref", &jets_.nref, "nref/I");
  t->Branch("rawpt", jets_.rawpt, "rawpt[nref]/F");
  t->Branch("jtpt", jets_.jtpt, "jtpt[nref]/F");
  t->Branch("jteta", jets_.jteta, "jteta[nref]/F");
  t->Branch("jtphi", jets_.jtphi, "jtphi[nref]/F");
  t->Branch("jtpu", jets_.jtpu, "jtpu[nref]/F");
  t->Branch("jtm", jets_.jtm, "jtm[nref]/F");
  t->Branch("jtarea", jets_.jtarea, "jtarea[nref]/F");

  t->Branch("WTAeta", jets_.WTAeta, "WTAeta[nref]/F");
  t->Branch("WTAphi", jets_.WTAphi, "WTAphi[nref]/F");

  t->Branch("jtPfCHF", jets_.jtPfCHF, "jtPfCHF[nref]/F");
  t->Branch("jtPfNHF", jets_.jtPfNHF, "jtPfNHF[nref]/F");
  t->Branch("jtPfCEF", jets_.jtPfCEF, "jtPfCEF[nref]/F");
  t->Branch("jtPfNEF", jets_.jtPfNEF, "jtPfNEF[nref]/F");
  t->Branch("jtPfMUF", jets_.jtPfMUF, "jtPfMUF[nref]/F");

  t->Branch("jtPfCHM", jets_.jtPfCHM, "jtPfCHM[nref]/I");
  t->Branch("jtPfNHM", jets_.jtPfNHM, "jtPfNHM[nref]/I");
  t->Branch("jtPfCEM", jets_.jtPfCEM, "jtPfCEM[nref]/I");
  t->Branch("jtPfNEM", jets_.jtPfNEM, "jtPfNEM[nref]/I");
  t->Branch("jtPfMUM", jets_.jtPfMUM, "jtPfMUM[nref]/I");

  if (isMC_) {
    t->Branch("pthat", &jets_.pthat, "pthat/F");

    t->Branch("refpt", jets_.refpt, "refpt[nref]/F");
    t->Branch("refeta", jets_.refeta, "refeta[nref]/F");
    t->Branch("refphi", jets_.refphi, "refphi[nref]/F");
    t->Branch("refm", jets_.refm, "refm[nref]/F");
    t->Branch("refarea", jets_.refarea, "refarea[nref]/F");
    t->Branch("subid", jets_.subid, "subid[nref]/I");

    t->Branch("ngen", &jets_.ngen, "ngen/I");

    t->Branch("genpt", jets_.genpt, "genpt[ngen]/F");
    t->Branch("geneta", jets_.geneta, "geneta[ngen]/F");
    t->Branch("genphi", jets_.genphi, "genphi[ngen]/F");
    t->Branch("genm", jets_.genm, "genm[ngen]/F");
    t->Branch("gensubid", jets_.gensubid, "gensubid[ngen]/I");

    t->Branch("WTAgeneta", jets_.WTAgeneta, "WTAgeneta[ngen]/F");
    t->Branch("WTAgenphi", jets_.WTAgenphi, "WTAgenphi[ngen]/F");
  }
}

void
HiInclusiveJetAnalyzer::analyze(const Event& iEvent,
                                const EventSetup& iSetup) {
  int event = iEvent.id().event();
  int run = iEvent.id().run();
  int lumi = iEvent.id().luminosityBlock();

  jets_.run = run;
  jets_.evt = event;
  jets_.lumi = lumi;

  edm::Handle<pat::JetCollection> patjets;
  iEvent.getByToken(jetTagPat_, patjets);

  edm::Handle<reco::JetView> jets;
  iEvent.getByToken(jetTag_, jets);

  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByToken(pfCandidateLabel_, pfCandidates);

  jets_.nref = 0;

  for (unsigned int j = 0; j < patjets->size(); ++j) {
    auto const& jet = (*patjets)[j];

    if (jet.pt() < jetPtMin_) continue;
    if (std::abs(jet.eta()) > jetAbsEtaMax_) continue;

    jets_.rawpt[jets_.nref]=(*patjets)[j].correctedJet("Uncorrected").pt();
    jets_.jtpt[jets_.nref] = jet.pt();
    jets_.jteta[jets_.nref] = jet.eta();
    jets_.jtphi[jets_.nref] = jet.phi();
    jets_.jtpu[jets_.nref] = jet.pileup();
    jets_.jtm[jets_.nref] = jet.mass();
    jets_.jtarea[jets_.nref] = jet.jetArea();

    if (jet.isPFJet()) {
      std::vector<fastjet::PseudoJet> candidates;
      for (auto const& it : jet.getJetConstituents()) {
	candidates.push_back(fastjet::PseudoJet(
	  (*it).px(), (*it).py(), (*it).pz(), (*it).energy()));
      }
      auto cs = new fastjet::ClusterSequence(candidates, WTAjtDef);
      std::vector<fastjet::PseudoJet> wtajt = fastjet::sorted_by_pt(cs->inclusive_jets(0));

      jets_.WTAeta[jets_.nref] = (!wtajt.empty()) ? wtajt[0].eta() : -999;
      jets_.WTAphi[jets_.nref] = (!wtajt.empty()) ? wtajt[0].phi_std() : -999;
      delete cs;
    }

    if ((*patjets)[j].isPFJet()) {
      jets_.jtPfCHF[jets_.nref] = (*patjets)[j].chargedHadronEnergyFraction();
      jets_.jtPfNHF[jets_.nref] = (*patjets)[j].neutralHadronEnergyFraction();
      jets_.jtPfCEF[jets_.nref] = (*patjets)[j].chargedEmEnergyFraction();
      jets_.jtPfNEF[jets_.nref] = (*patjets)[j].neutralEmEnergyFraction();
      jets_.jtPfMUF[jets_.nref] = (*patjets)[j].muonEnergyFraction();

      jets_.jtPfCHM[jets_.nref] = (*patjets)[j].chargedHadronMultiplicity();
      jets_.jtPfNHM[jets_.nref] = (*patjets)[j].neutralHadronMultiplicity();
      jets_.jtPfCEM[jets_.nref] = (*patjets)[j].electronMultiplicity();
      jets_.jtPfNEM[jets_.nref] = (*patjets)[j].photonMultiplicity();
      jets_.jtPfMUM[jets_.nref] = (*patjets)[j].muonMultiplicity();
    } else {
      jets_.jtPfCHF[jets_.nref] = -999;
      jets_.jtPfNHF[jets_.nref] = -999;
      jets_.jtPfCEF[jets_.nref] = -999;
      jets_.jtPfNEF[jets_.nref] = -999;
      jets_.jtPfMUF[jets_.nref] = -999;

      jets_.jtPfCHM[jets_.nref] = -999;
      jets_.jtPfNHM[jets_.nref] = -999;
      jets_.jtPfCEM[jets_.nref] = -999;
      jets_.jtPfNEM[jets_.nref] = -999;
      jets_.jtPfMUM[jets_.nref] = -999;
    }

    if (isMC_) {
      const reco::GenJet* genjet = (*patjets)[j].genJet();

      if (genjet) {
        jets_.refpt[jets_.nref] = genjet->pt();
        jets_.refeta[jets_.nref] = genjet->eta();
        jets_.refphi[jets_.nref] = genjet->phi();
        jets_.refm[jets_.nref] = genjet->mass();
        jets_.refarea[jets_.nref] = genjet->jetArea();

        const GenParticle* gencon = genjet->getGenConstituent(0);
        jets_.subid[jets_.nref] = gencon->collisionId();
      } else {
        jets_.refpt[jets_.nref] = -999;
        jets_.refeta[jets_.nref] = -999;
        jets_.refphi[jets_.nref] = -999;
        jets_.refm[jets_.nref] = -999;
        jets_.refarea[jets_.nref] = -999;
        jets_.subid[jets_.nref] = -999;
      }
    }

    jets_.nref++;
  }

  if (isMC_) {
    edm::Handle<GenEventInfoProduct> hEventInfo;
    iEvent.getByToken(eventGenInfoTag_, hEventInfo);

    jets_.pthat = hEventInfo->qScale();

    edm::Handle<edm::View<reco::GenJet>> genjets;
    iEvent.getByToken(genjetTag_, genjets);
    
    jets_.ngen = 0;

    for (auto const& genjet : *genjets) {
      std::vector<fastjet::PseudoJet> candidates;
      for (auto const& it : genjet.getJetConstituents()) {
	candidates.push_back(fastjet::PseudoJet(
	  (*it).px(), (*it).py(), (*it).pz(), (*it).energy()));
      }
      auto cs = new fastjet::ClusterSequence(candidates, WTAjtDef);
      std::vector<fastjet::PseudoJet> wtajt = fastjet::sorted_by_pt(cs->inclusive_jets(0));

      jets_.WTAgeneta[jets_.ngen] = (!wtajt.empty()) ? wtajt[0].eta() : -999;
      jets_.WTAgenphi[jets_.ngen] = (!wtajt.empty()) ? wtajt[0].phi_std() : -999;
      delete cs;

      // threshold to reduce size of output in minbias PbPb
      if (genjet.pt() > genPtMin_) {
        jets_.genpt [jets_.ngen] = genjet.pt();
        jets_.geneta[jets_.ngen] = genjet.eta();
        jets_.genphi[jets_.ngen] = genjet.phi();
        jets_.genm  [jets_.ngen] = genjet.mass();

        const GenParticle* gencon = genjet.getGenConstituent(0);
        jets_.gensubid[jets_.ngen] = gencon->collisionId();

        jets_.ngen++;
      }
    }
  }

  t->Fill();

  memset(&jets_, 0, sizeof(jets_));
}

DEFINE_FWK_MODULE(HiInclusiveJetAnalyzer);
