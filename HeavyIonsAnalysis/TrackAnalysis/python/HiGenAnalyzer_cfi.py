import FWCore.ParameterSet.Config as cms

HiGenParticleAna = cms.EDAnalyzer(
    'HiGenAnalyzer',
    ptMin = cms.untracked.double(5),
    etaMax = cms.untracked.double(2.5),
    genParticleSrc = cms.InputTag("genParticles"),
    genHiSrc = cms.InputTag("heavyIon"),
    )
