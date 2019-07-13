import FWCore.ParameterSet.Config as cms

inclusiveJetAnalyzer = cms.EDAnalyzer(
    "HiInclusiveJetAnalyzer",
    jetTag = cms.InputTag("ak4PFJets"),
    jetPtMin = cms.double(5.0),
    genjetTag = cms.InputTag("ak4GenJets"),
    isMC = cms.untracked.bool(False), 
    # dummy parameters below
    matchTag = cms.untracked.InputTag("ak4PFpatJets"),
    doGenTaus = cms.untracked.bool(False),
    doSubEvent = cms.untracked.bool(False),
    fillGenJets = cms.untracked.bool(False),
    rParam = cms.double(0.4),
    )
