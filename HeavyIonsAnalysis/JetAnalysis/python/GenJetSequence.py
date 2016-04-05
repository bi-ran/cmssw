import FWCore.ParameterSet.Config as cms

from RecoHI.HiJetAlgos.HiGenJets_cff import *
from RecoJets.Configuration.GenJetParticles_cff import *
from HeavyIonsAnalysis.JetAnalysis.akSoftDrop4GenJets_cfi import akSoftDrop4GenJets

akSoftDrop1GenJets = akSoftDrop4GenJets.clone(rParam = 0.1)
akSoftDrop2GenJets = akSoftDrop4GenJets.clone(rParam = 0.2)
akSoftDrop3GenJets = akSoftDrop4GenJets.clone(rParam = 0.3)
akSoftDrop5GenJets = akSoftDrop4GenJets.clone(rParam = 0.5)
akSoftDrop6GenJets = akSoftDrop4GenJets.clone(rParam = 0.6)

akHiGenJets = cms.Sequence(
    genParticlesForJets +
    ak1HiGenJets +
    ak2HiGenJets +
    ak3HiGenJets +
    ak4HiGenJets +
    ak5HiGenJets +
    ak6HiGenJets +
    akSoftDrop4GenJets +
    akSoftDrop5GenJets
)
