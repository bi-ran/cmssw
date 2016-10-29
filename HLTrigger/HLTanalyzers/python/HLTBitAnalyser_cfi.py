import FWCore.ParameterSet.Config as cms

hltbitanalysis = cms.EDAnalyzer("HLTBitAnalyzer",
    ### L1 Legacy and Stage 1 objects
    l1GctHFBitCounts                = cms.InputTag("hltGctDigis"),
    l1GctHFRingSums                 = cms.InputTag("hltGctDigis"),
    l1GtReadoutRecord               = cms.InputTag("hltGtDigis::HLT"),
    l1extramc                       = cms.string('hltL1extraParticles'),
    l1extramu                       = cms.string('hltL1extraParticles'),

    ### L1 Stage 2 objects
    gObjectMapRecord                = cms.InputTag("hltGtStage2ObjectMap"),
    gmtStage2Digis                  = cms.string("hltGmtStage2Digis"),
    caloStage2Digis                 = cms.string("hltCaloStage2Digis"),

    ### HLT
    hltresults                      = cms.InputTag("TriggerResults::HLT"),
    HLTProcessName                  = cms.string("HLT"),
                                
    ### Run parameters
    RunParameters = cms.PSet(
    HistogramFile = cms.untracked.string('hltbitanalysis.root'),
    UseL1Stage2   = cms.bool(True)
    )
                                
)
