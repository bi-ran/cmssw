import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.jets.HiReRecoJets_HI_cff import PFTowers, kt4PFJets, hiFJRhoProducer, hiFJGridEmptyAreaCalculator, akPu2CaloJets, akPu2PFJets, akCs2PFJets, akPu3CaloJets, akPu3PFJets, akCs3PFJets, akPu4CaloJets, akPu4PFJets, akCs4PFJets, akCsSoftDrop4PFJets, akCsSoftDropZ05B154PFJets

#jet analyzers
from HeavyIonsAnalysis.JetAnalysis.jets.akPu2CaloJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu2PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs2PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu3CaloJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu3PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs3PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu4CaloJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akPu4PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCs4PFJetSequence_PbPb_mc_cff import *

#SoftDrop analyzers
from HeavyIonsAnalysis.JetAnalysis.jets.akCsSoftDrop4PFJetSequence_PbPb_mc_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.akCsSoftDropZ05B154PFJetSequence_PbPb_mc_cff import *

highPurityTracks = cms.EDFilter("TrackSelector",
                                src = cms.InputTag("hiGeneralTracks"),
                                cut = cms.string('quality("highPurity")'))

from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import *
offlinePrimaryVertices.TrackLabel = 'highPurityTracks'

#the following lines are in the wrong python config
#they should be in a python config handling reco, not analyzers. To be fixed
akPu3PFJets.minimumTowersFraction = cms.double(0.5)
akPu4PFJets.minimumTowersFraction = cms.double(0.5)

akPu3CaloJets.minimumTowersFraction = cms.double(0.)
akPu4CaloJets.minimumTowersFraction = cms.double(0.)


jetSequences = cms.Sequence(
    PFTowers + 
    kt4PFJets +
    hiFJRhoProducer +
    hiFJGridEmptyAreaCalculator + 

    akPu2CaloJets +
    akPu2PFJets +
    akCs2PFJets +

    akPu3CaloJets +
    akPu3PFJets +
    akCs3PFJets +

    akPu4CaloJets +
    akPu4PFJets +
    akCs4PFJets +

    akCsSoftDrop4PFJets +
    akCsSoftDropZ05B154PFJets +
    
    highPurityTracks +
    offlinePrimaryVertices +

    akPu2CaloJetSequence +
    akPu2PFJetSequence +
    akCs2PFJetSequence +

    akPu3CaloJetSequence +
    akPu3PFJetSequence +
    akCs3PFJetSequence +

    akPu4CaloJetSequence +
    akPu4PFJetSequence +
    akCs4PFJetSequence +

    akCsSoftDrop4PFJetSequence +
    akCsSoftDropZ05B154PFJetSequence
)
