import FWCore.ParameterSet.Config as cms

inclusiveJetAnalyzer = cms.EDAnalyzer(
    "HiInclusiveJetAnalyzer",
    jetTag = cms.InputTag("ak4PFJets"),
    jetPtMin = cms.double(5.0),
    genjetTag = cms.InputTag("ak4HiGenJets"),
    isMC = cms.untracked.bool(False), 
    # dummy parameters below
    matchTag = cms.untracked.InputTag("akPu4PFpatJets"),
    doSubEvent = cms.untracked.bool(False),
    fillGenJets = cms.untracked.bool(False),
    rParam = cms.double(0.4),
    )
