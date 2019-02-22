
import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akPuSoftDropZ05B154Calomatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akPuSoftDropZ05B154CaloJets"),
    matched = cms.InputTag("ak4HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
)

akPuSoftDropZ05B154CalomatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B154HiSignalGenJets"),
    matched = cms.InputTag("ak4HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
)

# akPuSoftDropZ05B154Caloparton = patJetPartonMatch.clone(
#     src = cms.InputTag("akPuSoftDropZ05B154CaloJets")
# )

akPuSoftDropZ05B154Calocorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
    # primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akPuSoftDropZ05B154CaloJets"),
    payload = "AKPu4Calo_offline"
)

akPuSoftDropZ05B154CaloJetID = cms.EDProducer(
    'JetIDProducer',
    JetIDParams,
    src = cms.InputTag('akPuSoftDropZ05B154CaloJets'))

# akPuSoftDropZ05B154Caloclean = heavyIonCleanedGenJets.clone(
#     src = cms.InputTag('ak4HiSignalGenJets'))

# akPuSoftDropZ05B154CalobTagger = bTaggers(
#     "akPuSoftDropZ05B154Calo", 0.4)

# create objects locally since they dont load properly otherwise
# akPuSoftDropZ05B154Calomatch = akPuSoftDropZ05B154CalobTagger.match
akPuSoftDropZ05B154Caloparton = patJetPartonMatch.clone(
    src = cms.InputTag("akPuSoftDropZ05B154CaloJets"),
    matched = cms.InputTag("hiSignalGenParticles"))
# akPuSoftDropZ05B154CaloPatJetFlavourAssociationLegacy = akPuSoftDropZ05B154CalobTagger.PatJetFlavourAssociationLegacy
# akPuSoftDropZ05B154CaloPatJetPartons = akPuSoftDropZ05B154CalobTagger.PatJetPartons
# akPuSoftDropZ05B154CaloJetTracksAssociatorAtVertex = akPuSoftDropZ05B154CalobTagger.JetTracksAssociatorAtVertex
# akPuSoftDropZ05B154CaloJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
# akPuSoftDropZ05B154CaloSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B154CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akPuSoftDropZ05B154CaloSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B154CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akPuSoftDropZ05B154CaloCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B154CalobTagger.CombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B154CaloCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B154CalobTagger.CombinedSecondaryVertexV2BJetTags
# akPuSoftDropZ05B154CaloJetBProbabilityBJetTags = akPuSoftDropZ05B154CalobTagger.JetBProbabilityBJetTags
# akPuSoftDropZ05B154CaloSoftPFMuonByPtBJetTags = akPuSoftDropZ05B154CalobTagger.SoftPFMuonByPtBJetTags
# akPuSoftDropZ05B154CaloSoftPFMuonByIP3dBJetTags = akPuSoftDropZ05B154CalobTagger.SoftPFMuonByIP3dBJetTags
# akPuSoftDropZ05B154CaloTrackCountingHighEffBJetTags = akPuSoftDropZ05B154CalobTagger.TrackCountingHighEffBJetTags
# akPuSoftDropZ05B154CaloTrackCountingHighPurBJetTags = akPuSoftDropZ05B154CalobTagger.TrackCountingHighPurBJetTags
# akPuSoftDropZ05B154CaloPatJetPartonAssociationLegacy = akPuSoftDropZ05B154CalobTagger.PatJetPartonAssociationLegacy

# akPuSoftDropZ05B154CaloImpactParameterTagInfos = akPuSoftDropZ05B154CalobTagger.ImpactParameterTagInfos
# akPuSoftDropZ05B154CaloImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akPuSoftDropZ05B154CaloJetProbabilityBJetTags = akPuSoftDropZ05B154CalobTagger.JetProbabilityBJetTags
# 
# akPuSoftDropZ05B154CaloSecondaryVertexTagInfos = akPuSoftDropZ05B154CalobTagger.SecondaryVertexTagInfos
# akPuSoftDropZ05B154CaloSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B154CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akPuSoftDropZ05B154CaloSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B154CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akPuSoftDropZ05B154CaloCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B154CalobTagger.CombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B154CaloCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B154CalobTagger.CombinedSecondaryVertexV2BJetTags
# 
# akPuSoftDropZ05B154CaloSecondaryVertexNegativeTagInfos = akPuSoftDropZ05B154CalobTagger.SecondaryVertexNegativeTagInfos
# akPuSoftDropZ05B154CaloNegativeSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B154CalobTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
# akPuSoftDropZ05B154CaloNegativeSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B154CalobTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
# akPuSoftDropZ05B154CaloNegativeCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B154CalobTagger.NegativeCombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B154CaloPositiveCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B154CalobTagger.PositiveCombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B154CaloNegativeCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B154CalobTagger.NegativeCombinedSecondaryVertexV2BJetTags
# akPuSoftDropZ05B154CaloPositiveCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B154CalobTagger.PositiveCombinedSecondaryVertexV2BJetTags
# 
# akPuSoftDropZ05B154CaloSoftPFMuonsTagInfos = akPuSoftDropZ05B154CalobTagger.SoftPFMuonsTagInfos
# akPuSoftDropZ05B154CaloSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akPuSoftDropZ05B154CaloSoftPFMuonBJetTags = akPuSoftDropZ05B154CalobTagger.SoftPFMuonBJetTags
# akPuSoftDropZ05B154CaloSoftPFMuonByIP3dBJetTags = akPuSoftDropZ05B154CalobTagger.SoftPFMuonByIP3dBJetTags
# akPuSoftDropZ05B154CaloSoftPFMuonByPtBJetTags = akPuSoftDropZ05B154CalobTagger.SoftPFMuonByPtBJetTags
# akPuSoftDropZ05B154CaloNegativeSoftPFMuonByPtBJetTags = akPuSoftDropZ05B154CalobTagger.NegativeSoftPFMuonByPtBJetTags
# akPuSoftDropZ05B154CaloPositiveSoftPFMuonByPtBJetTags = akPuSoftDropZ05B154CalobTagger.PositiveSoftPFMuonByPtBJetTags
# akPuSoftDropZ05B154CaloPatJetFlavourIdLegacy = cms.Sequence(akPuSoftDropZ05B154CaloPatJetPartonAssociationLegacy*akPuSoftDropZ05B154CaloPatJetFlavourAssociationLegacy)
# Not working with our PU sub, but keep it here for reference
# akPuSoftDropZ05B154CaloPatJetFlavourAssociation = akPuSoftDropZ05B154CalobTagger.PatJetFlavourAssociation
# akPuSoftDropZ05B154CaloPatJetFlavourId = cms.Sequence(
#     akPuSoftDropZ05B154CaloPatJetPartons
#     *
#     akPuSoftDropZ05B154CaloPatJetFlavourAssociation)

# akPuSoftDropZ05B154CaloJetBtaggingIP  = cms.Sequence(
#     akPuSoftDropZ05B154CaloImpactParameterTagInfos *
#     akPuSoftDropZ05B154CaloTrackCountingHighEffBJetTags +
#     akPuSoftDropZ05B154CaloTrackCountingHighPurBJetTags +
#     akPuSoftDropZ05B154CaloJetProbabilityBJetTags +
#     akPuSoftDropZ05B154CaloJetBProbabilityBJetTags 
# )

# akPuSoftDropZ05B154CaloJetBtaggingSV = cms.Sequence(
#     akPuSoftDropZ05B154CaloImpactParameterTagInfos
#     *
#     akPuSoftDropZ05B154CaloSecondaryVertexTagInfos
#     *
#     akPuSoftDropZ05B154CaloSimpleSecondaryVertexHighEffBJetTags
#     +
#     akPuSoftDropZ05B154CaloSimpleSecondaryVertexHighPurBJetTags
#     +
#     akPuSoftDropZ05B154CaloCombinedSecondaryVertexBJetTags
#     +
#     akPuSoftDropZ05B154CaloCombinedSecondaryVertexV2BJetTags
# )
# 
# akPuSoftDropZ05B154CaloJetBtaggingNegSV = cms.Sequence(
#     akPuSoftDropZ05B154CaloImpactParameterTagInfos
#     *
#     akPuSoftDropZ05B154CaloSecondaryVertexNegativeTagInfos
#     *
#     akPuSoftDropZ05B154CaloNegativeSimpleSecondaryVertexHighEffBJetTags
#     +
#     akPuSoftDropZ05B154CaloNegativeSimpleSecondaryVertexHighPurBJetTags
#     +
#     akPuSoftDropZ05B154CaloNegativeCombinedSecondaryVertexBJetTags
#     +
#     akPuSoftDropZ05B154CaloPositiveCombinedSecondaryVertexBJetTags
#     +
#     akPuSoftDropZ05B154CaloNegativeCombinedSecondaryVertexV2BJetTags
#     +
#     akPuSoftDropZ05B154CaloPositiveCombinedSecondaryVertexV2BJetTags
# )
# 
# akPuSoftDropZ05B154CaloJetBtaggingMu = cms.Sequence(
#     akPuSoftDropZ05B154CaloSoftPFMuonsTagInfos
#     *
#     akPuSoftDropZ05B154CaloSoftPFMuonBJetTags
#     +
#     akPuSoftDropZ05B154CaloSoftPFMuonByIP3dBJetTags
#     +
#     akPuSoftDropZ05B154CaloSoftPFMuonByPtBJetTags
#     +
#     akPuSoftDropZ05B154CaloNegativeSoftPFMuonByPtBJetTags
#     +
#     akPuSoftDropZ05B154CaloPositiveSoftPFMuonByPtBJetTags
# )
# 
# akPuSoftDropZ05B154CaloJetBtagging = cms.Sequence(
#     akPuSoftDropZ05B154CaloJetBtaggingIP
#     *
#     akPuSoftDropZ05B154CaloJetBtaggingSV
#     *
#     akPuSoftDropZ05B154CaloJetBtaggingNegSV
#     # *
#     # akPuSoftDropZ05B154CaloJetBtaggingMu
# )

akPuSoftDropZ05B154CalopatJetsWithBtagging = patJets.clone(
    jetSource              = cms.InputTag("akPuSoftDropZ05B154CaloJets"),
    genJetMatch            = cms.InputTag("akPuSoftDropZ05B154Calomatch"),
    genPartonMatch         = cms.InputTag("akPuSoftDropZ05B154Caloparton"),
    jetCorrFactorsSource   = cms.VInputTag(cms.InputTag("akPuSoftDropZ05B154Calocorr")),
#     JetPartonMapSource     = cms.InputTag("akPuSoftDropZ05B154CaloPatJetFlavourAssociationLegacy"),
#     JetFlavourInfoSource   = cms.InputTag("akPuSoftDropZ05B154CaloPatJetFlavourAssociation"),
#     trackAssociationSource = cms.InputTag("akPuSoftDropZ05B154CaloJetTracksAssociatorAtVertex"),
#     useLegacyJetMCFlavour  = True,
#     discriminatorSources   = cms.VInputTag(
#         cms.InputTag("akPuSoftDropZ05B154CaloSimpleSecondaryVertexHighEffBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B154CaloSimpleSecondaryVertexHighPurBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B154CaloCombinedSecondaryVertexBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B154CaloCombinedSecondaryVertexV2BJetTags"),
#         cms.InputTag("akPuSoftDropZ05B154CaloJetBProbabilityBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B154CaloJetProbabilityBJetTags"),
#         # cms.InputTag("akPuSoftDropZ05B154CaloSoftPFMuonByPtBJetTags"),
#         # cms.InputTag("akPuSoftDropZ05B154CaloSoftPFMuonByIP3dBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B154CaloTrackCountingHighEffBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B154CaloTrackCountingHighPurBJetTags"),
#     ),
    jetIDMap = cms.InputTag("akPuSoftDropZ05B154CaloJetID"),
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

akPuSoftDropZ05B154CaloNjettiness = Njettiness.clone(
    src = cms.InputTag("akPuSoftDropZ05B154CaloJets"),
    R0  = cms.double( 0.4)
)

# akPuSoftDropZ05B154CalopatJetsWithBtagging.userData.userFloats.src += [
#     'akPuSoftDropZ05B154CaloNjettiness:tau1',
#     'akPuSoftDropZ05B154CaloNjettiness:tau2',
#     'akPuSoftDropZ05B154CaloNjettiness:tau3'
# ]

akPuSoftDropZ05B154CaloJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("akPuSoftDropZ05B154CalopatJetsWithBtagging"),
    genjetTag = 'ak4HiSignalGenJets',
    rParam = 0.4,
    matchJets = cms.untracked.bool(False),
    matchTag = 'patJetsWithBtagging',
    pfCandidateLabel = cms.untracked.InputTag('particleFlowTmp'),
    trackTag = cms.InputTag("hiGeneralTracks"),
    fillGenJets = True,
    isMC = True,
    doSubEvent = True,
    useHepMC = cms.untracked.bool(False),
    genParticles = cms.untracked.InputTag("genParticles"),
    eventInfoTag = cms.InputTag("generator"),
    # doLifeTimeTagging = cms.untracked.bool(True),
    # doLifeTimeTaggingExtras = cms.untracked.bool(False),
    # bTagJetName = cms.untracked.string("akPuSoftDropZ05B154Calo"),
    jetName = cms.untracked.string("akPuSoftDropZ05B154Calo"),
    genPtMin = cms.untracked.double(5),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
    doTower = cms.untracked.bool(True),
    doSubJets = cms.untracked.bool(True),
    doGenSubJets = cms.untracked.bool(True),     
    subjetGenTag = cms.untracked.InputTag("akSoftDropZ05B154GenJets"),
    doGenTaus = cms.untracked.bool(False),
    genTau1 = cms.InputTag("akSoftDropZ05B154GenNjettiness","tau1"),
    genTau2 = cms.InputTag("akSoftDropZ05B154GenNjettiness","tau2"),
    genTau3 = cms.InputTag("akSoftDropZ05B154GenNjettiness","tau3"),
    doGenSym = cms.untracked.bool(True),
    genSym = cms.InputTag("akSoftDropZ05B154GenJets","sym"),
    genDroppedBranches = cms.InputTag("akSoftDropZ05B154GenJets","droppedBranches")
)

akPuSoftDropZ05B154CaloJetSequence_mc = cms.Sequence(
    # akPuSoftDropZ05B154Caloclean
    # *
    akPuSoftDropZ05B154Calomatch
    *
    # akPuSoftDropZ05B154CalomatchGroomed
    # *
    akPuSoftDropZ05B154Caloparton
    *
    akPuSoftDropZ05B154Calocorr
    *
    # akPuSoftDropZ05B154CaloJetID
    # *
    # akPuSoftDropZ05B154CaloPatJetFlavourIdLegacy
    # *
    # Use legacy algo till PU implemented
    # akPuSoftDropZ05B154CaloPatJetFlavourId
    # *
    # akPuSoftDropZ05B154CaloJetTracksAssociatorAtVertex
    # *
    # akPuSoftDropZ05B154CaloJetBtagging
    # *
    # No constituents for calo jets in pp. Must be removed for pp calo jets but
    # I'm not sure how to do this transparently (Marta)
    akPuSoftDropZ05B154CaloNjettiness
    *
    akPuSoftDropZ05B154CalopatJetsWithBtagging
    *
    akPuSoftDropZ05B154CaloJetAnalyzer
)

akPuSoftDropZ05B154CaloJetSequence_data = cms.Sequence(
    akPuSoftDropZ05B154Calocorr
    *
    # akPuSoftDropZ05B154CaloJetID
    # *
    # akPuSoftDropZ05B154CaloJetTracksAssociatorAtVertex
    # *
    # akPuSoftDropZ05B154CaloJetBtagging
    # *
    akPuSoftDropZ05B154CaloNjettiness 
    *
    akPuSoftDropZ05B154CalopatJetsWithBtagging
    *
    akPuSoftDropZ05B154CaloJetAnalyzer
)

akPuSoftDropZ05B154CaloJetSequence_jec = cms.Sequence(
    akPuSoftDropZ05B154CaloJetSequence_mc)
akPuSoftDropZ05B154CaloJetSequence_mb = cms.Sequence(
    akPuSoftDropZ05B154CaloJetSequence_mc)

akPuSoftDropZ05B154CaloJetSequence = cms.Sequence(
    akPuSoftDropZ05B154CaloJetSequence_jec)
akPuSoftDropZ05B154CaloJetAnalyzer.genPtMin = cms.untracked.double(1)
akPuSoftDropZ05B154CaloJetAnalyzer.jetPtMin = cms.double(1)
akPuSoftDropZ05B154CalopatJetsWithBtagging.userData.userFloats.src += ['akPuSoftDropZ05B154CaloJets:sym']
akPuSoftDropZ05B154CalopatJetsWithBtagging.userData.userInts.src += ['akPuSoftDropZ05B154CaloJets:droppedBranches']
