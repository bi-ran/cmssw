import FWCore.ParameterSet.Config as cms

ppTrack = cms.EDAnalyzer(
    'TrackAnalyzer',
    trackPtMin = cms.double(0.01),
    vertexSrc = cms.InputTag('offlinePrimaryVertices'),
    trackSrc = cms.InputTag('generalTracks'),
    mvaSrc = cms.InputTag("generalTracks", "MVAValues"),
    pfCandSrc = cms.InputTag('particleFlow'),
)
