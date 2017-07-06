import FWCore.ParameterSet.Config as cms
from RecoHI.HiJetAlgos.HiRecoJets_cff import akPu2CaloJets, akPu3CaloJets, akPu4CaloJets, ak2CaloJets, ak3CaloJets, ak4CaloJets
from RecoHI.HiJetAlgos.HiRecoPFJets_cff import PFTowers, ak5PFJets, akPu2PFJets, akPu3PFJets, akPu4PFJets
from RecoJets.JetProducers.akCs4PFJets_cfi import akCs4PFJets

#------------------------------------------------------------------------------
#    Create jet producers and change default settings of jet producers
#------------------------------------------------------------------------------
ak5PFJets.doAreaFastjet = cms.bool(True)

ak2PFJets = ak5PFJets.clone(rParam       = cms.double(0.2))
ak3PFJets = ak5PFJets.clone(rParam       = cms.double(0.3))
ak4PFJets = ak5PFJets.clone(rParam       = cms.double(0.4))

akPu2PFJets.jetPtMin = 1
akPu2CaloJets.jetPtMin = 1
akPu3PFJets.jetPtMin = 1
akPu3CaloJets.jetPtMin = 1
akPu4PFJets.jetPtMin = 1
akPu4CaloJets.jetPtMin = 1

ak2PFJets.jetPtMin = 1
ak2PFJets.src = cms.InputTag("particleFlow")
ak2CaloJets.jetPtMin = 1
ak3PFJets.jetPtMin = 1
ak3PFJets.src = cms.InputTag("particleFlow")
ak3CaloJets.jetPtMin = 1
ak4PFJets.jetPtMin = 1
ak4PFJets.src = cms.InputTag("particleFlow")
ak4CaloJets.jetPtMin = 1

akCs4PFJets.jetPtMin = 1
akCs4PFJets.src = cms.InputTag("particleFlow")

PFTowers.src = cms.InputTag("particleFlow")

#------------------------------------------------------------------------------
#    SoftDrop jets - no background subtraction
#------------------------------------------------------------------------------
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
akSoftDrop4PFJets = cms.EDProducer(
    "SoftDropJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0   = cms.double(0.4),
    useOnlyCharged = cms.bool(False),
    useExplicitGhosts = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
akSoftDrop4PFJets.src = cms.InputTag("particleFlow")

akSoftDropZ05B154PFJets = akSoftDrop4PFJets.clone(zcut=cms.double(0.5), beta=cms.double(1.5))

#------------------------------------------------------------------------------
#    KT jets and rho estimators
#------------------------------------------------------------------------------
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
kt4PFJets.src = cms.InputTag('particleFlow')
kt4PFJets.doAreaFastjet = True
kt4PFJets.jetPtMin      = cms.double(0.0)
kt4PFJets.GhostArea     = cms.double(0.005)
from RecoHI.HiJetAlgos.hiFJRhoProducer import hiFJRhoProducer
hiFJRhoProducer.jetSource = cms.InputTag("kt4PFJets")
from RecoHI.HiJetAlgos.hiFJGridEmptyAreaCalculator_cff import hiFJGridEmptyAreaCalculator
hiFJGridEmptyAreaCalculator.jetSource = cms.InputTag("kt4PFJets")
hiFJGridEmptyAreaCalculator.CentralityBinSrc = cms.InputTag("centralityBin","HFtowersPlusTrunc")
hiFJGridEmptyAreaCalculator.pfCandSource = cms.InputTag('particleFlow')

#------------------------------------------------------------------------------
#    CS jets
#------------------------------------------------------------------------------
akCs4PFJets.rho      = cms.InputTag('hiFJGridEmptyAreaCalculator','mapToRhoCorr1Bin')
akCs4PFJets.rhom      = cms.InputTag('hiFJGridEmptyAreaCalculator','mapToRhoMCorr1Bin')
akCs2PFJets = akCs4PFJets.clone(rParam       = cms.double(0.2))
akCs3PFJets = akCs4PFJets.clone(rParam       = cms.double(0.3))

#------------------------------------------------------------------------------
#    SoftDrop jets - after constituent subtraction
#------------------------------------------------------------------------------
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
akCsSoftDrop4PFJets = cms.EDProducer(
    "SoftDropJetProducer",
    PFJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0   = cms.double(0.4),
    useOnlyCharged = cms.bool(False),
    useExplicitGhosts = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
akCsSoftDrop4PFJets.src = cms.InputTag("akCs4PFJets","pfParticlesCs")

akCsSoftDropZ05B154PFJets = akCsSoftDrop4PFJets.clone(zcut=cms.double(0.5), beta=cms.double(1.5))


hiReRecoPFJets = cms.Sequence(
    akPu2PFJets
    +
    akPu3PFJets
    +
    akPu4PFJets
    +
    ak2PFJets
    +
    ak3PFJets
    +
    ak4PFJets
    +
    akCs2PFJets
    +
    akCs3PFJets
    +
    akCs4PFJets
    +
    akSoftDrop4PFJets
    +
    akSoftDropZ05B154PFJets
    +
    akCsSoftDrop4PFJets
    +
    akCsSoftDropZ05B154PFJets
)

hiReRecoCaloJets = cms.Sequence(
    akPu3CaloJets
    +
    akPu4CaloJets
    +
    ak3CaloJets
    +
    ak4CaloJets
)
