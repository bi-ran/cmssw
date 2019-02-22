
import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akSoftDropZ05B154Calomatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B154CaloJets"),
    matched = cms.InputTag("ak4HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
)

akSoftDropZ05B154CalomatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B154HiSignalGenJets"),
    matched = cms.InputTag("ak4HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
)

# akSoftDropZ05B154Caloparton = patJetPartonMatch.clone(
#     src = cms.InputTag("akSoftDropZ05B154CaloJets")
# )

akSoftDropZ05B154Calocorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
    # primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akSoftDropZ05B154CaloJets"),
    payload = "AK4Calo_offline"
)

akSoftDropZ05B154CaloJetID = cms.EDProducer(
    'JetIDProducer',
    JetIDParams,
    src = cms.InputTag('akSoftDropZ05B154CaloJets'))

# akSoftDropZ05B154Caloclean = heavyIonCleanedGenJets.clone(
#     src = cms.InputTag('ak4HiSignalGenJets'))

# akSoftDropZ05B154CalobTagger = bTaggers(
#     "akSoftDropZ05B154Calo", 0.4)

# create objects locally since they dont load properly otherwise
# akSoftDropZ05B154Calomatch = akSoftDropZ05B154CalobTagger.match
akSoftDropZ05B154Caloparton = patJetPartonMatch.clone(
    src = cms.InputTag("akSoftDropZ05B154CaloJets"),
    matched = cms.InputTag("hiSignalGenParticles"))
# akSoftDropZ05B154CaloPatJetFlavourAssociationLegacy = akSoftDropZ05B154CalobTagger.PatJetFlavourAssociationLegacy
# akSoftDropZ05B154CaloPatJetPartons = akSoftDropZ05B154CalobTagger.PatJetPartons
# akSoftDropZ05B154CaloJetTracksAssociatorAtVertex = akSoftDropZ05B154CalobTagger.JetTracksAssociatorAtVertex
# akSoftDropZ05B154CaloJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
# akSoftDropZ05B154CaloSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B154CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akSoftDropZ05B154CaloSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B154CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akSoftDropZ05B154CaloCombinedSecondaryVertexBJetTags = akSoftDropZ05B154CalobTagger.CombinedSecondaryVertexBJetTags
# akSoftDropZ05B154CaloCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B154CalobTagger.CombinedSecondaryVertexV2BJetTags
# akSoftDropZ05B154CaloJetBProbabilityBJetTags = akSoftDropZ05B154CalobTagger.JetBProbabilityBJetTags
# akSoftDropZ05B154CaloSoftPFMuonByPtBJetTags = akSoftDropZ05B154CalobTagger.SoftPFMuonByPtBJetTags
# akSoftDropZ05B154CaloSoftPFMuonByIP3dBJetTags = akSoftDropZ05B154CalobTagger.SoftPFMuonByIP3dBJetTags
# akSoftDropZ05B154CaloTrackCountingHighEffBJetTags = akSoftDropZ05B154CalobTagger.TrackCountingHighEffBJetTags
# akSoftDropZ05B154CaloTrackCountingHighPurBJetTags = akSoftDropZ05B154CalobTagger.TrackCountingHighPurBJetTags
# akSoftDropZ05B154CaloPatJetPartonAssociationLegacy = akSoftDropZ05B154CalobTagger.PatJetPartonAssociationLegacy

# akSoftDropZ05B154CaloImpactParameterTagInfos = akSoftDropZ05B154CalobTagger.ImpactParameterTagInfos
# akSoftDropZ05B154CaloImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akSoftDropZ05B154CaloJetProbabilityBJetTags = akSoftDropZ05B154CalobTagger.JetProbabilityBJetTags
# 
# akSoftDropZ05B154CaloSecondaryVertexTagInfos = akSoftDropZ05B154CalobTagger.SecondaryVertexTagInfos
# akSoftDropZ05B154CaloSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B154CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akSoftDropZ05B154CaloSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B154CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akSoftDropZ05B154CaloCombinedSecondaryVertexBJetTags = akSoftDropZ05B154CalobTagger.CombinedSecondaryVertexBJetTags
# akSoftDropZ05B154CaloCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B154CalobTagger.CombinedSecondaryVertexV2BJetTags
# 
# akSoftDropZ05B154CaloSecondaryVertexNegativeTagInfos = akSoftDropZ05B154CalobTagger.SecondaryVertexNegativeTagInfos
# akSoftDropZ05B154CaloNegativeSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B154CalobTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
# akSoftDropZ05B154CaloNegativeSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B154CalobTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
# akSoftDropZ05B154CaloNegativeCombinedSecondaryVertexBJetTags = akSoftDropZ05B154CalobTagger.NegativeCombinedSecondaryVertexBJetTags
# akSoftDropZ05B154CaloPositiveCombinedSecondaryVertexBJetTags = akSoftDropZ05B154CalobTagger.PositiveCombinedSecondaryVertexBJetTags
# akSoftDropZ05B154CaloNegativeCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B154CalobTagger.NegativeCombinedSecondaryVertexV2BJetTags
# akSoftDropZ05B154CaloPositiveCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B154CalobTagger.PositiveCombinedSecondaryVertexV2BJetTags
# 
# akSoftDropZ05B154CaloSoftPFMuonsTagInfos = akSoftDropZ05B154CalobTagger.SoftPFMuonsTagInfos
# akSoftDropZ05B154CaloSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akSoftDropZ05B154CaloSoftPFMuonBJetTags = akSoftDropZ05B154CalobTagger.SoftPFMuonBJetTags
# akSoftDropZ05B154CaloSoftPFMuonByIP3dBJetTags = akSoftDropZ05B154CalobTagger.SoftPFMuonByIP3dBJetTags
# akSoftDropZ05B154CaloSoftPFMuonByPtBJetTags = akSoftDropZ05B154CalobTagger.SoftPFMuonByPtBJetTags
# akSoftDropZ05B154CaloNegativeSoftPFMuonByPtBJetTags = akSoftDropZ05B154CalobTagger.NegativeSoftPFMuonByPtBJetTags
# akSoftDropZ05B154CaloPositiveSoftPFMuonByPtBJetTags = akSoftDropZ05B154CalobTagger.PositiveSoftPFMuonByPtBJetTags
# akSoftDropZ05B154CaloPatJetFlavourIdLegacy = cms.Sequence(akSoftDropZ05B154CaloPatJetPartonAssociationLegacy*akSoftDropZ05B154CaloPatJetFlavourAssociationLegacy)
# Not working with our PU sub, but keep it here for reference
# akSoftDropZ05B154CaloPatJetFlavourAssociation = akSoftDropZ05B154CalobTagger.PatJetFlavourAssociation
# akSoftDropZ05B154CaloPatJetFlavourId = cms.Sequence(
#     akSoftDropZ05B154CaloPatJetPartons
#     *
#     akSoftDropZ05B154CaloPatJetFlavourAssociation)

# akSoftDropZ05B154CaloJetBtaggingIP  = cms.Sequence(
#     akSoftDropZ05B154CaloImpactParameterTagInfos *
#     akSoftDropZ05B154CaloTrackCountingHighEffBJetTags +
#     akSoftDropZ05B154CaloTrackCountingHighPurBJetTags +
#     akSoftDropZ05B154CaloJetProbabilityBJetTags +
#     akSoftDropZ05B154CaloJetBProbabilityBJetTags 
# )

# akSoftDropZ05B154CaloJetBtaggingSV = cms.Sequence(
#     akSoftDropZ05B154CaloImpactParameterTagInfos
#     *
#     akSoftDropZ05B154CaloSecondaryVertexTagInfos
#     *
#     akSoftDropZ05B154CaloSimpleSecondaryVertexHighEffBJetTags
#     +
#     akSoftDropZ05B154CaloSimpleSecondaryVertexHighPurBJetTags
#     +
#     akSoftDropZ05B154CaloCombinedSecondaryVertexBJetTags
#     +
#     akSoftDropZ05B154CaloCombinedSecondaryVertexV2BJetTags
# )
# 
# akSoftDropZ05B154CaloJetBtaggingNegSV = cms.Sequence(
#     akSoftDropZ05B154CaloImpactParameterTagInfos
#     *
#     akSoftDropZ05B154CaloSecondaryVertexNegativeTagInfos
#     *
#     akSoftDropZ05B154CaloNegativeSimpleSecondaryVertexHighEffBJetTags
#     +
#     akSoftDropZ05B154CaloNegativeSimpleSecondaryVertexHighPurBJetTags
#     +
#     akSoftDropZ05B154CaloNegativeCombinedSecondaryVertexBJetTags
#     +
#     akSoftDropZ05B154CaloPositiveCombinedSecondaryVertexBJetTags
#     +
#     akSoftDropZ05B154CaloNegativeCombinedSecondaryVertexV2BJetTags
#     +
#     akSoftDropZ05B154CaloPositiveCombinedSecondaryVertexV2BJetTags
# )
# 
# akSoftDropZ05B154CaloJetBtaggingMu = cms.Sequence(
#     akSoftDropZ05B154CaloSoftPFMuonsTagInfos
#     *
#     akSoftDropZ05B154CaloSoftPFMuonBJetTags
#     +
#     akSoftDropZ05B154CaloSoftPFMuonByIP3dBJetTags
#     +
#     akSoftDropZ05B154CaloSoftPFMuonByPtBJetTags
#     +
#     akSoftDropZ05B154CaloNegativeSoftPFMuonByPtBJetTags
#     +
#     akSoftDropZ05B154CaloPositiveSoftPFMuonByPtBJetTags
# )
# 
# akSoftDropZ05B154CaloJetBtagging = cms.Sequence(
#     akSoftDropZ05B154CaloJetBtaggingIP
#     *
#     akSoftDropZ05B154CaloJetBtaggingSV
#     *
#     akSoftDropZ05B154CaloJetBtaggingNegSV
#     # *
#     # akSoftDropZ05B154CaloJetBtaggingMu
# )

akSoftDropZ05B154CalopatJetsWithBtagging = patJets.clone(
    jetSource              = cms.InputTag("akSoftDropZ05B154CaloJets"),
    genJetMatch            = cms.InputTag("akSoftDropZ05B154Calomatch"),
    genPartonMatch         = cms.InputTag("akSoftDropZ05B154Caloparton"),
    jetCorrFactorsSource   = cms.VInputTag(cms.InputTag("akSoftDropZ05B154Calocorr")),
#     JetPartonMapSource     = cms.InputTag("akSoftDropZ05B154CaloPatJetFlavourAssociationLegacy"),
#     JetFlavourInfoSource   = cms.InputTag("akSoftDropZ05B154CaloPatJetFlavourAssociation"),
#     trackAssociationSource = cms.InputTag("akSoftDropZ05B154CaloJetTracksAssociatorAtVertex"),
#     useLegacyJetMCFlavour  = True,
#     discriminatorSources   = cms.VInputTag(
#         cms.InputTag("akSoftDropZ05B154CaloSimpleSecondaryVertexHighEffBJetTags"),
#         cms.InputTag("akSoftDropZ05B154CaloSimpleSecondaryVertexHighPurBJetTags"),
#         cms.InputTag("akSoftDropZ05B154CaloCombinedSecondaryVertexBJetTags"),
#         cms.InputTag("akSoftDropZ05B154CaloCombinedSecondaryVertexV2BJetTags"),
#         cms.InputTag("akSoftDropZ05B154CaloJetBProbabilityBJetTags"),
#         cms.InputTag("akSoftDropZ05B154CaloJetProbabilityBJetTags"),
#         # cms.InputTag("akSoftDropZ05B154CaloSoftPFMuonByPtBJetTags"),
#         # cms.InputTag("akSoftDropZ05B154CaloSoftPFMuonByIP3dBJetTags"),
#         cms.InputTag("akSoftDropZ05B154CaloTrackCountingHighEffBJetTags"),
#         cms.InputTag("akSoftDropZ05B154CaloTrackCountingHighPurBJetTags"),
#     ),
    jetIDMap = cms.InputTag("akSoftDropZ05B154CaloJetID"),
    addBTagInfo         = False,
    addTagInfos         = False,
    addDiscriminators   = False,
    addAssociatedTracks = False,
    addJetCharge        = False,
    addJetID            = False,
    getJetMCFlavour     = False,
    addGenPartonMatch   = False,
    addGenJetMatch      = False,
    embedGenJetMatch    = False,
    embedGenPartonMatch = False,
    # embedCaloTowers   = False,
    # embedPFCandidates = False
)

akSoftDropZ05B154CaloNjettiness = Njettiness.clone(
    src = cms.InputTag("akSoftDropZ05B154CaloJets"),
    R0  = cms.double( 0.4)
)

# akSoftDropZ05B154CalopatJetsWithBtagging.userData.userFloats.src += [
#     'akSoftDropZ05B154CaloNjettiness:tau1',
#     'akSoftDropZ05B154CaloNjettiness:tau2',
#     'akSoftDropZ05B154CaloNjettiness:tau3'
# ]

akSoftDropZ05B154CaloJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("akSoftDropZ05B154CalopatJetsWithBtagging"),
    genjetTag = 'ak4HiSignalGenJets',
    rParam = 0.4,
    matchJets = cms.untracked.bool(False),
    matchTag = 'patJetsWithBtagging',
    pfCandidateLabel = cms.untracked.InputTag('particleFlowTmp'),
    trackTag = cms.InputTag("hiGeneralTracks"),
    fillGenJets = False,
    isMC = False,
    doSubEvent = False,
    useHepMC = cms.untracked.bool(False),
    genParticles = cms.untracked.InputTag("genParticles"),
    eventInfoTag = cms.InputTag("generator"),
    # doLifeTimeTagging = cms.untracked.bool(True),
    # doLifeTimeTaggingExtras = cms.untracked.bool(False),
    # bTagJetName = cms.untracked.string("akSoftDropZ05B154Calo"),
    jetName = cms.untracked.string("akSoftDropZ05B154Calo"),
    genPtMin = cms.untracked.double(5),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
    doTower = cms.untracked.bool(True),
    doSubJets = cms.untracked.bool(True),
    doGenSubJets = cms.untracked.bool(False),     
    subjetGenTag = cms.untracked.InputTag("akSoftDropZ05B154GenJets"),
    doGenTaus = cms.untracked.bool(False),
    genTau1 = cms.InputTag("akSoftDropZ05B154GenNjettiness","tau1"),
    genTau2 = cms.InputTag("akSoftDropZ05B154GenNjettiness","tau2"),
    genTau3 = cms.InputTag("akSoftDropZ05B154GenNjettiness","tau3"),
    doGenSym = cms.untracked.bool(False),
    genSym = cms.InputTag("akSoftDropZ05B154GenJets","sym"),
    genDroppedBranches = cms.InputTag("akSoftDropZ05B154GenJets","droppedBranches")
)

akSoftDropZ05B154CaloJetSequence_mc = cms.Sequence(
    # akSoftDropZ05B154Caloclean
    # *
    akSoftDropZ05B154Calomatch
    *
    # akSoftDropZ05B154CalomatchGroomed
    # *
    akSoftDropZ05B154Caloparton
    *
    akSoftDropZ05B154Calocorr
    *
    # akSoftDropZ05B154CaloJetID
    # *
    # akSoftDropZ05B154CaloPatJetFlavourIdLegacy
    # *
    # Use legacy algo till PU implemented
    # akSoftDropZ05B154CaloPatJetFlavourId
    # *
    # akSoftDropZ05B154CaloJetTracksAssociatorAtVertex
    # *
    # akSoftDropZ05B154CaloJetBtagging
    # *
    # No constituents for calo jets in pp. Must be removed for pp calo jets but
    # I'm not sure how to do this transparently (Marta)
    akSoftDropZ05B154CaloNjettiness
    *
    akSoftDropZ05B154CalopatJetsWithBtagging
    *
    akSoftDropZ05B154CaloJetAnalyzer
)

akSoftDropZ05B154CaloJetSequence_data = cms.Sequence(
    akSoftDropZ05B154Calocorr
    *
    # akSoftDropZ05B154CaloJetID
    # *
    # akSoftDropZ05B154CaloJetTracksAssociatorAtVertex
    # *
    # akSoftDropZ05B154CaloJetBtagging
    # *
    akSoftDropZ05B154CaloNjettiness 
    *
    akSoftDropZ05B154CalopatJetsWithBtagging
    *
    akSoftDropZ05B154CaloJetAnalyzer
)

akSoftDropZ05B154CaloJetSequence_jec = cms.Sequence(
    akSoftDropZ05B154CaloJetSequence_mc)
akSoftDropZ05B154CaloJetSequence_mb = cms.Sequence(
    akSoftDropZ05B154CaloJetSequence_mc)

akSoftDropZ05B154CaloJetSequence = cms.Sequence(
    akSoftDropZ05B154CaloJetSequence_data)
akSoftDropZ05B154CalopatJetsWithBtagging.userData.userFloats.src += ['akSoftDropZ05B154CaloJets:sym']
akSoftDropZ05B154CalopatJetsWithBtagging.userData.userInts.src += ['akSoftDropZ05B154CaloJets:droppedBranches']
