
import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akPuSoftDropZ05B153Calomatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akPuSoftDropZ05B153CaloJets"),
    matched = cms.InputTag("ak3HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.3
)

akPuSoftDropZ05B153CalomatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B153HiSignalGenJets"),
    matched = cms.InputTag("ak3HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.3
)

# akPuSoftDropZ05B153Caloparton = patJetPartonMatch.clone(
#     src = cms.InputTag("akPuSoftDropZ05B153CaloJets")
# )

akPuSoftDropZ05B153Calocorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
    # primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akPuSoftDropZ05B153CaloJets"),
    payload = "AKPu3Calo_offline"
)

akPuSoftDropZ05B153CaloJetID = cms.EDProducer(
    'JetIDProducer',
    JetIDParams,
    src = cms.InputTag('akPuSoftDropZ05B153CaloJets'))

# akPuSoftDropZ05B153Caloclean = heavyIonCleanedGenJets.clone(
#     src = cms.InputTag('ak3HiSignalGenJets'))

# akPuSoftDropZ05B153CalobTagger = bTaggers(
#     "akPuSoftDropZ05B153Calo", 0.3)

# create objects locally since they dont load properly otherwise
# akPuSoftDropZ05B153Calomatch = akPuSoftDropZ05B153CalobTagger.match
akPuSoftDropZ05B153Caloparton = patJetPartonMatch.clone(
    src = cms.InputTag("akPuSoftDropZ05B153CaloJets"),
    matched = cms.InputTag("genParticles"))
# akPuSoftDropZ05B153CaloPatJetFlavourAssociationLegacy = akPuSoftDropZ05B153CalobTagger.PatJetFlavourAssociationLegacy
# akPuSoftDropZ05B153CaloPatJetPartons = akPuSoftDropZ05B153CalobTagger.PatJetPartons
# akPuSoftDropZ05B153CaloJetTracksAssociatorAtVertex = akPuSoftDropZ05B153CalobTagger.JetTracksAssociatorAtVertex
# akPuSoftDropZ05B153CaloJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
# akPuSoftDropZ05B153CaloSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B153CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akPuSoftDropZ05B153CaloSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B153CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akPuSoftDropZ05B153CaloCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B153CalobTagger.CombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B153CaloCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B153CalobTagger.CombinedSecondaryVertexV2BJetTags
# akPuSoftDropZ05B153CaloJetBProbabilityBJetTags = akPuSoftDropZ05B153CalobTagger.JetBProbabilityBJetTags
# akPuSoftDropZ05B153CaloSoftPFMuonByPtBJetTags = akPuSoftDropZ05B153CalobTagger.SoftPFMuonByPtBJetTags
# akPuSoftDropZ05B153CaloSoftPFMuonByIP3dBJetTags = akPuSoftDropZ05B153CalobTagger.SoftPFMuonByIP3dBJetTags
# akPuSoftDropZ05B153CaloTrackCountingHighEffBJetTags = akPuSoftDropZ05B153CalobTagger.TrackCountingHighEffBJetTags
# akPuSoftDropZ05B153CaloTrackCountingHighPurBJetTags = akPuSoftDropZ05B153CalobTagger.TrackCountingHighPurBJetTags
# akPuSoftDropZ05B153CaloPatJetPartonAssociationLegacy = akPuSoftDropZ05B153CalobTagger.PatJetPartonAssociationLegacy

# akPuSoftDropZ05B153CaloImpactParameterTagInfos = akPuSoftDropZ05B153CalobTagger.ImpactParameterTagInfos
# akPuSoftDropZ05B153CaloImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akPuSoftDropZ05B153CaloJetProbabilityBJetTags = akPuSoftDropZ05B153CalobTagger.JetProbabilityBJetTags
# 
# akPuSoftDropZ05B153CaloSecondaryVertexTagInfos = akPuSoftDropZ05B153CalobTagger.SecondaryVertexTagInfos
# akPuSoftDropZ05B153CaloSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B153CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akPuSoftDropZ05B153CaloSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B153CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akPuSoftDropZ05B153CaloCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B153CalobTagger.CombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B153CaloCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B153CalobTagger.CombinedSecondaryVertexV2BJetTags
# 
# akPuSoftDropZ05B153CaloSecondaryVertexNegativeTagInfos = akPuSoftDropZ05B153CalobTagger.SecondaryVertexNegativeTagInfos
# akPuSoftDropZ05B153CaloNegativeSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B153CalobTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
# akPuSoftDropZ05B153CaloNegativeSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B153CalobTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
# akPuSoftDropZ05B153CaloNegativeCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B153CalobTagger.NegativeCombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B153CaloPositiveCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B153CalobTagger.PositiveCombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B153CaloNegativeCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B153CalobTagger.NegativeCombinedSecondaryVertexV2BJetTags
# akPuSoftDropZ05B153CaloPositiveCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B153CalobTagger.PositiveCombinedSecondaryVertexV2BJetTags
# 
# akPuSoftDropZ05B153CaloSoftPFMuonsTagInfos = akPuSoftDropZ05B153CalobTagger.SoftPFMuonsTagInfos
# akPuSoftDropZ05B153CaloSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akPuSoftDropZ05B153CaloSoftPFMuonBJetTags = akPuSoftDropZ05B153CalobTagger.SoftPFMuonBJetTags
# akPuSoftDropZ05B153CaloSoftPFMuonByIP3dBJetTags = akPuSoftDropZ05B153CalobTagger.SoftPFMuonByIP3dBJetTags
# akPuSoftDropZ05B153CaloSoftPFMuonByPtBJetTags = akPuSoftDropZ05B153CalobTagger.SoftPFMuonByPtBJetTags
# akPuSoftDropZ05B153CaloNegativeSoftPFMuonByPtBJetTags = akPuSoftDropZ05B153CalobTagger.NegativeSoftPFMuonByPtBJetTags
# akPuSoftDropZ05B153CaloPositiveSoftPFMuonByPtBJetTags = akPuSoftDropZ05B153CalobTagger.PositiveSoftPFMuonByPtBJetTags
# akPuSoftDropZ05B153CaloPatJetFlavourIdLegacy = cms.Sequence(akPuSoftDropZ05B153CaloPatJetPartonAssociationLegacy*akPuSoftDropZ05B153CaloPatJetFlavourAssociationLegacy)
# Not working with our PU sub, but keep it here for reference
# akPuSoftDropZ05B153CaloPatJetFlavourAssociation = akPuSoftDropZ05B153CalobTagger.PatJetFlavourAssociation
# akPuSoftDropZ05B153CaloPatJetFlavourId = cms.Sequence(
#     akPuSoftDropZ05B153CaloPatJetPartons
#     *
#     akPuSoftDropZ05B153CaloPatJetFlavourAssociation)

# akPuSoftDropZ05B153CaloJetBtaggingIP  = cms.Sequence(
#     akPuSoftDropZ05B153CaloImpactParameterTagInfos *
#     akPuSoftDropZ05B153CaloTrackCountingHighEffBJetTags +
#     akPuSoftDropZ05B153CaloTrackCountingHighPurBJetTags +
#     akPuSoftDropZ05B153CaloJetProbabilityBJetTags +
#     akPuSoftDropZ05B153CaloJetBProbabilityBJetTags 
# )

# akPuSoftDropZ05B153CaloJetBtaggingSV = cms.Sequence(
#     akPuSoftDropZ05B153CaloImpactParameterTagInfos
#     *
#     akPuSoftDropZ05B153CaloSecondaryVertexTagInfos
#     *
#     akPuSoftDropZ05B153CaloSimpleSecondaryVertexHighEffBJetTags
#     +
#     akPuSoftDropZ05B153CaloSimpleSecondaryVertexHighPurBJetTags
#     +
#     akPuSoftDropZ05B153CaloCombinedSecondaryVertexBJetTags
#     +
#     akPuSoftDropZ05B153CaloCombinedSecondaryVertexV2BJetTags
# )
# 
# akPuSoftDropZ05B153CaloJetBtaggingNegSV = cms.Sequence(
#     akPuSoftDropZ05B153CaloImpactParameterTagInfos
#     *
#     akPuSoftDropZ05B153CaloSecondaryVertexNegativeTagInfos
#     *
#     akPuSoftDropZ05B153CaloNegativeSimpleSecondaryVertexHighEffBJetTags
#     +
#     akPuSoftDropZ05B153CaloNegativeSimpleSecondaryVertexHighPurBJetTags
#     +
#     akPuSoftDropZ05B153CaloNegativeCombinedSecondaryVertexBJetTags
#     +
#     akPuSoftDropZ05B153CaloPositiveCombinedSecondaryVertexBJetTags
#     +
#     akPuSoftDropZ05B153CaloNegativeCombinedSecondaryVertexV2BJetTags
#     +
#     akPuSoftDropZ05B153CaloPositiveCombinedSecondaryVertexV2BJetTags
# )
# 
# akPuSoftDropZ05B153CaloJetBtaggingMu = cms.Sequence(
#     akPuSoftDropZ05B153CaloSoftPFMuonsTagInfos
#     *
#     akPuSoftDropZ05B153CaloSoftPFMuonBJetTags
#     +
#     akPuSoftDropZ05B153CaloSoftPFMuonByIP3dBJetTags
#     +
#     akPuSoftDropZ05B153CaloSoftPFMuonByPtBJetTags
#     +
#     akPuSoftDropZ05B153CaloNegativeSoftPFMuonByPtBJetTags
#     +
#     akPuSoftDropZ05B153CaloPositiveSoftPFMuonByPtBJetTags
# )
# 
# akPuSoftDropZ05B153CaloJetBtagging = cms.Sequence(
#     akPuSoftDropZ05B153CaloJetBtaggingIP
#     *
#     akPuSoftDropZ05B153CaloJetBtaggingSV
#     *
#     akPuSoftDropZ05B153CaloJetBtaggingNegSV
#     # *
#     # akPuSoftDropZ05B153CaloJetBtaggingMu
# )

akPuSoftDropZ05B153CalopatJetsWithBtagging = patJets.clone(
    jetSource              = cms.InputTag("akPuSoftDropZ05B153CaloJets"),
    genJetMatch            = cms.InputTag("akPuSoftDropZ05B153Calomatch"),
    genPartonMatch         = cms.InputTag("akPuSoftDropZ05B153Caloparton"),
    jetCorrFactorsSource   = cms.VInputTag(cms.InputTag("akPuSoftDropZ05B153Calocorr")),
#     JetPartonMapSource     = cms.InputTag("akPuSoftDropZ05B153CaloPatJetFlavourAssociationLegacy"),
#     JetFlavourInfoSource   = cms.InputTag("akPuSoftDropZ05B153CaloPatJetFlavourAssociation"),
#     trackAssociationSource = cms.InputTag("akPuSoftDropZ05B153CaloJetTracksAssociatorAtVertex"),
#     useLegacyJetMCFlavour  = True,
#     discriminatorSources   = cms.VInputTag(
#         cms.InputTag("akPuSoftDropZ05B153CaloSimpleSecondaryVertexHighEffBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B153CaloSimpleSecondaryVertexHighPurBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B153CaloCombinedSecondaryVertexBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B153CaloCombinedSecondaryVertexV2BJetTags"),
#         cms.InputTag("akPuSoftDropZ05B153CaloJetBProbabilityBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B153CaloJetProbabilityBJetTags"),
#         # cms.InputTag("akPuSoftDropZ05B153CaloSoftPFMuonByPtBJetTags"),
#         # cms.InputTag("akPuSoftDropZ05B153CaloSoftPFMuonByIP3dBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B153CaloTrackCountingHighEffBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B153CaloTrackCountingHighPurBJetTags"),
#     ),
    jetIDMap = cms.InputTag("akPuSoftDropZ05B153CaloJetID"),
    addBTagInfo         = False,
    addTagInfos         = False,
    addDiscriminators   = False,
    addAssociatedTracks = False,
    addJetCharge        = False,
    addJetID            = False,
    getJetMCFlavour     = False,
    addGenPartonMatch   = True,
    addGenJetMatch      = True,
    embedGenJetMatch    = True,
    embedGenPartonMatch = True,
    # embedCaloTowers   = False,
    # embedPFCandidates = False
)

akPuSoftDropZ05B153CaloNjettiness = Njettiness.clone(
    src = cms.InputTag("akPuSoftDropZ05B153CaloJets"),
    R0  = cms.double( 0.3)
)

# akPuSoftDropZ05B153CalopatJetsWithBtagging.userData.userFloats.src += [
#     'akPuSoftDropZ05B153CaloNjettiness:tau1',
#     'akPuSoftDropZ05B153CaloNjettiness:tau2',
#     'akPuSoftDropZ05B153CaloNjettiness:tau3'
# ]

akPuSoftDropZ05B153CaloJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("akPuSoftDropZ05B153CalopatJetsWithBtagging"),
    genjetTag = 'ak3HiSignalGenJets',
    rParam = 0.3,
    matchJets = cms.untracked.bool(False),
    matchTag = 'patJetsWithBtagging',
    pfCandidateLabel = cms.untracked.InputTag('particleFlow'),
    trackTag = cms.InputTag("generalTracks"),
    fillGenJets = True,
    isMC = True,
    doSubEvent = True,
    useHepMC = cms.untracked.bool(False),
    genParticles = cms.untracked.InputTag("genParticles"),
    eventInfoTag = cms.InputTag("generator"),
    # doLifeTimeTagging = cms.untracked.bool(True),
    # doLifeTimeTaggingExtras = cms.untracked.bool(False),
    # bTagJetName = cms.untracked.string("akPuSoftDropZ05B153Calo"),
    jetName = cms.untracked.string("akPuSoftDropZ05B153Calo"),
    genPtMin = cms.untracked.double(5),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
    doTower = cms.untracked.bool(True),
    doSubJets = cms.untracked.bool(True),
    doGenSubJets = cms.untracked.bool(True),     
    subjetGenTag = cms.untracked.InputTag("akSoftDropZ05B153GenJets"),
    doGenTaus = cms.untracked.bool(False),
    genTau1 = cms.InputTag("akSoftDropZ05B153GenNjettiness","tau1"),
    genTau2 = cms.InputTag("akSoftDropZ05B153GenNjettiness","tau2"),
    genTau3 = cms.InputTag("akSoftDropZ05B153GenNjettiness","tau3"),
    doGenSym = cms.untracked.bool(True),
    genSym = cms.InputTag("akSoftDropZ05B153GenJets","sym"),
    genDroppedBranches = cms.InputTag("akSoftDropZ05B153GenJets","droppedBranches")
)

akPuSoftDropZ05B153CaloJetSequence_mc = cms.Sequence(
    # akPuSoftDropZ05B153Caloclean
    # *
    akPuSoftDropZ05B153Calomatch
    *
    # akPuSoftDropZ05B153CalomatchGroomed
    # *
    akPuSoftDropZ05B153Caloparton
    *
    akPuSoftDropZ05B153Calocorr
    *
    # akPuSoftDropZ05B153CaloJetID
    # *
    # akPuSoftDropZ05B153CaloPatJetFlavourIdLegacy
    # *
    # Use legacy algo till PU implemented
    # akPuSoftDropZ05B153CaloPatJetFlavourId
    # *
    # akPuSoftDropZ05B153CaloJetTracksAssociatorAtVertex
    # *
    # akPuSoftDropZ05B153CaloJetBtagging
    # *
    # No constituents for calo jets in pp. Must be removed for pp calo jets but
    # I'm not sure how to do this transparently (Marta)
    akPuSoftDropZ05B153CaloNjettiness
    *
    akPuSoftDropZ05B153CalopatJetsWithBtagging
    *
    akPuSoftDropZ05B153CaloJetAnalyzer
)

akPuSoftDropZ05B153CaloJetSequence_data = cms.Sequence(
    akPuSoftDropZ05B153Calocorr
    *
    # akPuSoftDropZ05B153CaloJetID
    # *
    # akPuSoftDropZ05B153CaloJetTracksAssociatorAtVertex
    # *
    # akPuSoftDropZ05B153CaloJetBtagging
    # *
    akPuSoftDropZ05B153CaloNjettiness 
    *
    akPuSoftDropZ05B153CalopatJetsWithBtagging
    *
    akPuSoftDropZ05B153CaloJetAnalyzer
)

akPuSoftDropZ05B153CaloJetSequence_jec = cms.Sequence(
    akPuSoftDropZ05B153CaloJetSequence_mc)
akPuSoftDropZ05B153CaloJetSequence_mb = cms.Sequence(
    akPuSoftDropZ05B153CaloJetSequence_mc)

akPuSoftDropZ05B153CaloJetSequence = cms.Sequence(
    akPuSoftDropZ05B153CaloJetSequence_mc)
akPuSoftDropZ05B153CalopatJetsWithBtagging.userData.userFloats.src += ['akPuSoftDropZ05B153CaloJets:sym']
akPuSoftDropZ05B153CalopatJetsWithBtagging.userData.userInts.src += ['akPuSoftDropZ05B153CaloJets:droppedBranches']
