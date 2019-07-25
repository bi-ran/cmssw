#ifndef MNguyen_HiInclusiveJetAnalyzer_inclusiveJetAnalyzer_
#define MNguyen_HiInclusiveJetAnalyzer_inclusiveJetAnalyzer_

// system include files
#include <string>

// ROOT headers
#include "TTree.h"

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// external headers
#include <fastjet/JetDefinition.hh>

//

/**\class HiInclusiveJetAnalyzer

   \author Matt Nguyen
   \date   November 2010
*/

class HiInclusiveJetAnalyzer : public edm::EDAnalyzer {
public:

  explicit HiInclusiveJetAnalyzer(const edm::ParameterSet&);
  ~HiInclusiveJetAnalyzer();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginRun(const edm::Run& r, const edm::EventSetup& c);
  virtual void beginJob();

private:

  // for reWTA reclustering-----------------------
  fastjet::JetDefinition const WTAjtDef = fastjet::JetDefinition(
     fastjet::JetAlgorithm::antikt_algorithm, 2, fastjet::WTA_pt_scheme);
  //--------------------------------------------

  edm::InputTag jetTagLabel_;
  edm::EDGetTokenT<reco::JetView> jetTag_;
  edm::EDGetTokenT<pat::JetCollection> jetTagPat_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidateLabel_;
  edm::EDGetTokenT<edm::View<reco::GenJet>> genjetTag_;
  edm::EDGetTokenT<GenEventInfoProduct> eventGenInfoTag_;
  
  bool isMC_;

  double genPtMin_;
  double jetPtMin_;
  double jetAbsEtaMax_;

  TTree *t;
  edm::Service<TFileService> f_;

  static constexpr int MAXJETS = 1000;

  struct JRA {
    int run;
    int evt;
    int lumi;

    int nref;

    float jtpt[MAXJETS];
    float rawpt[MAXJETS];
    float jteta[MAXJETS];
    float jtphi[MAXJETS];

    float WTAeta[MAXJETS];
    float WTAphi[MAXJETS];

    float jtpu[MAXJETS];
    float jtm[MAXJETS];
    float jtarea[MAXJETS];

    float jtPfCHF[MAXJETS];
    float jtPfNHF[MAXJETS];
    float jtPfCEF[MAXJETS];
    float jtPfNEF[MAXJETS];
    float jtPfMUF[MAXJETS];

    int jtPfCHM[MAXJETS];
    int jtPfNHM[MAXJETS];
    int jtPfCEM[MAXJETS];
    int jtPfNEM[MAXJETS];
    int jtPfMUM[MAXJETS];

    float refpt[MAXJETS];
    float refeta[MAXJETS];
    float refphi[MAXJETS];
    float refm[MAXJETS];
    float refarea[MAXJETS];
    float refflavour[MAXJETS];
    int subid[MAXJETS];

    float pthat;

    int ngen;

    float genpt[MAXJETS];
    float geneta[MAXJETS];
    float genphi[MAXJETS];
    float genm[MAXJETS];
    int gensubid[MAXJETS];

    float WTAgeneta[MAXJETS];
    float WTAgenphi[MAXJETS];
  };

  JRA jets_;
};

#endif
