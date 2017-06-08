import FWCore.ParameterSet.Config as cms

from RecoHI.HiJetAlgos.HiGenJets_cff import *
from RecoJets.Configuration.GenJetParticles_cff import *
from RecoHI.HiJetAlgos.HiSignalParticleProducer_cfi import hiSignalGenParticles
from HeavyIonsAnalysis.JetAnalysis.akSoftDrop4GenJets_cfi import akSoftDrop4GenJets

genParticlesForJetsSignal = genParticlesForJets.clone(src = cms.InputTag("hiSignalGenParticles"))
akSoftDrop4GenJets.src            = cms.InputTag("genParticlesForJetsSignal")

#add positive beta grooming
akSoftDropZ05B154GenJets = akSoftDrop4GenJets.clone(zcut=cms.double(0.5), beta=cms.double(1.5))

akSoftDrop5GenJets = akSoftDrop4GenJets.clone(rParam = cms.double(0.5) , R0 = cms.double(0.5))
akSoftDropZ05B155GenJets = akSoftDropZ05B154GenJets.clone(rParam = cms.double(0.5) , R0 = cms.double(0.5))

#Njettiness gen jets
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
ak2GenNjettiness = Njettiness.clone(
                    src = cms.InputTag("ak2HiSignalGenJets"),
                    R0  = cms.double( 0.2)
)

ak3GenNjettiness = Njettiness.clone(
                    src = cms.InputTag("ak3HiSignalGenJets"),
                    R0  = cms.double( 0.3)
)
ak4GenNjettiness = Njettiness.clone(
                    src = cms.InputTag("ak4HiSignalGenJets"),
                    R0  = cms.double( 0.4)
)
ak5GenNjettiness = Njettiness.clone(
                    src = cms.InputTag("ak5HiSignalGenJets"),
                    R0  = cms.double( 0.5)
)

akHiGenJets = cms.Sequence(
    genParticlesForJets +
    hiSignalGenParticles +
    genParticlesForJetsSignal +
    ak1HiGenJets +
    ak2HiGenJets +
    ak3HiGenJets +
    ak4HiGenJets +
    ak5HiGenJets +
    ak6HiGenJets +
    akSoftDrop4GenJets +
    akSoftDrop5GenJets +
    akSoftDropZ05B154GenJets +
    akSoftDropZ05B155GenJets
)
