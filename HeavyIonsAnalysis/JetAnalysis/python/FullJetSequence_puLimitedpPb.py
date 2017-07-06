import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.jets.HiReRecoJets_pPb_cff import PFTowers, kt4PFJets, hiFJRhoProducer, hiFJGridEmptyAreaCalculator, akPu2CaloJets, akPu2PFJets, akCs2PFJets, akPu3CaloJets, akPu3PFJets, akCs3PFJets, akPu4CaloJets, akPu4PFJets, akCs4PFJets, akCsSoftDrop4PFJets, akCsSoftDropZ05B154PFJets, akSoftDrop4PFJets, akSoftDropZ05B154PFJets, ak2PFJets, ak3PFJets, ak4PFJets, ak2CaloJets, ak3CaloJets, ak4CaloJets

#jet analyzers
from HeavyIonsAnalysis.JetAnalysis.jets.ak3PFJetSequence_pPb_jec_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4PFJetSequence_pPb_jec_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu3PFJetSequence_pPb_jec_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu4PFJetSequence_pPb_jec_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4CaloJetSequence_pPb_jec_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu4CaloJetSequence_pPb_jec_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs3PFJetSequence_pPb_jec_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs4PFJetSequence_pPb_jec_cff import *

#softdrop analyzers
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDrop4PFJetSequence_pPb_jec_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akSoftDropZ05B154PFJetSequence_pPb_jec_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCsSoftDrop4PFJetSequence_pPb_jec_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCsSoftDropZ05B154PFJetSequence_pPb_jec_cff import *

highPurityTracks = cms.EDFilter("TrackSelector",
                                src = cms.InputTag("generalTracks"),
                                cut = cms.string('quality("highPurity")')
)

#the following lines are in the wrong python config
#they should be in a python config handling reco, not analyzers. To be fixed
akPu2PFJets.minimumTowersFraction = cms.double(0.5)
akPu3PFJets.minimumTowersFraction = cms.double(0.5)
akPu4PFJets.minimumTowersFraction = cms.double(0.5)
akPu4CaloJets.minimumTowersFraction = cms.double(0.)

#left commented by Marta since I don't know why we are turning this on here
#ak4PFJetAnalyzer.doSubEvent = True
#ak4CaloJetAnalyzer.doSubEvent = True

jetSequences = cms.Sequence(
    PFTowers +
    kt4PFJets +
    hiFJRhoProducer +
    hiFJGridEmptyAreaCalculator +
    ak3PFJets +
    ak4PFJets +
    akPu3PFJets +
    akPu4PFJets +
    ak4CaloJets +
    akPu4CaloJets +
    akCs3PFJets +
    akCs4PFJets +
    akSoftDrop4PFJets +
    akSoftDropZ05B154PFJets +
    akCsSoftDrop4PFJets +
    akCsSoftDropZ05B154PFJets +
    highPurityTracks +
    ak3PFJetSequence +
    ak4PFJetSequence +
    akPu3PFJetSequence +
    akPu4PFJetSequence +
    ak4CaloJetSequence +
    akPu4CaloJetSequence +
    akCs3PFJetSequence +
    akCs4PFJetSequence +
    akSoftDrop4PFJetSequence +
    akSoftDropZ05B154PFJetSequence +
    akCsSoftDrop4PFJetSequence +
    akCsSoftDropZ05B154PFJetSequence
)
