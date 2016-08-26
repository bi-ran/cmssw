import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.jets.HiReRecoJets_pp_cff import *

from HeavyIonsAnalysis.JetAnalysis.jets.ak4PFJetSequence_pPb_jec_cff import *
from HeavyIonsAnalysis.JetAnalysis.jets.ak4CaloJetSequence_pPb_jec_cff import *

ak4PFJetAnalyzer.doSubEvent = True
ak4CaloJetAnalyzer.doSubEvent = True

from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
ak5GenJets = ak5GenJets
ak4GenJets = ak5GenJets.clone(rParam = 0.4)
from RecoJets.Configuration.GenJetParticles_cff import *

akGenJets = cms.Sequence(
    genParticlesForJets +
    ak4GenJets
)

from HeavyIonsAnalysis.JetAnalysis.makePartons_cff import *
highPurityTracks = cms.EDFilter("TrackSelector",
                                src = cms.InputTag("generalTracks"),
                                cut = cms.string('quality("highPurity")')
)

jetSequences = cms.Sequence(
    akGenJets +
#    ppReRecoPFJets +
#    ppReRecoCaloJets +
    makePartons +
    highPurityTracks +
    ak4CaloJets +
    ak4PFJetSequence +
    ak4CaloJetSequence 
    )
