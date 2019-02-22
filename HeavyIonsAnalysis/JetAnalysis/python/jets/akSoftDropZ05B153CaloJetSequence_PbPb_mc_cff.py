
import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akSoftDropZ05B153Calomatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B153CaloJets"),
    matched = cms.InputTag("ak3HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.3
)

akSoftDropZ05B153CalomatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B153HiSignalGenJets"),
    matched = cms.InputTag("ak3HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.3
)

# akSoftDropZ05B153Caloparton = patJetPartonMatch.clone(
#     src = cms.InputTag("akSoftDropZ05B153CaloJets")
# )

akSoftDropZ05B153Calocorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
    # primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akSoftDropZ05B153CaloJets"),
    payload = "AK3Calo_offline"
)

akSoftDropZ05B153CaloJetID = cms.EDProducer(
    'JetIDProducer',
    JetIDParams,
    src = cms.InputTag('akSoftDropZ05B153CaloJets'))

# akSoftDropZ05B153Caloclean = heavyIonCleanedGenJets.clone(
#     src = cms.InputTag('ak3HiSignalGenJets'))

# akSoftDropZ05B153CalobTagger = bTaggers(
#     "akSoftDropZ05B153Calo", 0.3)

# create objects locally since they dont load properly otherwise
# akSoftDropZ05B153Calomatch = akSoftDropZ05B153CalobTagger.match
akSoftDropZ05B153Caloparton = patJetPartonMatch.clone(
    src = cms.InputTag("akSoftDropZ05B153CaloJets"),
    matched = cms.InputTag("hiSignalGenParticles"))
# akSoftDropZ05B153CaloPatJetFlavourAssociationLegacy = akSoftDropZ05B153CalobTagger.PatJetFlavourAssociationLegacy
# akSoftDropZ05B153CaloPatJetPartons = akSoftDropZ05B153CalobTagger.PatJetPartons
# akSoftDropZ05B153CaloJetTracksAssociatorAtVertex = akSoftDropZ05B153CalobTagger.JetTracksAssociatorAtVertex
# akSoftDropZ05B153CaloJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
# akSoftDropZ05B153CaloSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B153CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akSoftDropZ05B153CaloSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B153CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akSoftDropZ05B153CaloCombinedSecondaryVertexBJetTags = akSoftDropZ05B153CalobTagger.CombinedSecondaryVertexBJetTags
# akSoftDropZ05B153CaloCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B153CalobTagger.CombinedSecondaryVertexV2BJetTags
# akSoftDropZ05B153CaloJetBProbabilityBJetTags = akSoftDropZ05B153CalobTagger.JetBProbabilityBJetTags
# akSoftDropZ05B153CaloSoftPFMuonByPtBJetTags = akSoftDropZ05B153CalobTagger.SoftPFMuonByPtBJetTags
# akSoftDropZ05B153CaloSoftPFMuonByIP3dBJetTags = akSoftDropZ05B153CalobTagger.SoftPFMuonByIP3dBJetTags
# akSoftDropZ05B153CaloTrackCountingHighEffBJetTags = akSoftDropZ05B153CalobTagger.TrackCountingHighEffBJetTags
# akSoftDropZ05B153CaloTrackCountingHighPurBJetTags = akSoftDropZ05B153CalobTagger.TrackCountingHighPurBJetTags
# akSoftDropZ05B153CaloPatJetPartonAssociationLegacy = akSoftDropZ05B153CalobTagger.PatJetPartonAssociationLegacy

# akSoftDropZ05B153CaloImpactParameterTagInfos = akSoftDropZ05B153CalobTagger.ImpactParameterTagInfos
# akSoftDropZ05B153CaloImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akSoftDropZ05B153CaloJetProbabilityBJetTags = akSoftDropZ05B153CalobTagger.JetProbabilityBJetTags
# 
# akSoftDropZ05B153CaloSecondaryVertexTagInfos = akSoftDropZ05B153CalobTagger.SecondaryVertexTagInfos
# akSoftDropZ05B153CaloSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B153CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akSoftDropZ05B153CaloSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B153CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akSoftDropZ05B153CaloCombinedSecondaryVertexBJetTags = akSoftDropZ05B153CalobTagger.CombinedSecondaryVertexBJetTags
# akSoftDropZ05B153CaloCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B153CalobTagger.CombinedSecondaryVertexV2BJetTags
# 
# akSoftDropZ05B153CaloSecondaryVertexNegativeTagInfos = akSoftDropZ05B153CalobTagger.SecondaryVertexNegativeTagInfos
# akSoftDropZ05B153CaloNegativeSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B153CalobTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
# akSoftDropZ05B153CaloNegativeSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B153CalobTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
# akSoftDropZ05B153CaloNegativeCombinedSecondaryVertexBJetTags = akSoftDropZ05B153CalobTagger.NegativeCombinedSecondaryVertexBJetTags
# akSoftDropZ05B153CaloPositiveCombinedSecondaryVertexBJetTags = akSoftDropZ05B153CalobTagger.PositiveCombinedSecondaryVertexBJetTags
# akSoftDropZ05B153CaloNegativeCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B153CalobTagger.NegativeCombinedSecondaryVertexV2BJetTags
# akSoftDropZ05B153CaloPositiveCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B153CalobTagger.PositiveCombinedSecondaryVertexV2BJetTags
# 
# akSoftDropZ05B153CaloSoftPFMuonsTagInfos = akSoftDropZ05B153CalobTagger.SoftPFMuonsTagInfos
# akSoftDropZ05B153CaloSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akSoftDropZ05B153CaloSoftPFMuonBJetTags = akSoftDropZ05B153CalobTagger.SoftPFMuonBJetTags
# akSoftDropZ05B153CaloSoftPFMuonByIP3dBJetTags = akSoftDropZ05B153CalobTagger.SoftPFMuonByIP3dBJetTags
# akSoftDropZ05B153CaloSoftPFMuonByPtBJetTags = akSoftDropZ05B153CalobTagger.SoftPFMuonByPtBJetTags
# akSoftDropZ05B153CaloNegativeSoftPFMuonByPtBJetTags = akSoftDropZ05B153CalobTagger.NegativeSoftPFMuonByPtBJetTags
# akSoftDropZ05B153CaloPositiveSoftPFMuonByPtBJetTags = akSoftDropZ05B153CalobTagger.PositiveSoftPFMuonByPtBJetTags
# akSoftDropZ05B153CaloPatJetFlavourIdLegacy = cms.Sequence(akSoftDropZ05B153CaloPatJetPartonAssociationLegacy*akSoftDropZ05B153CaloPatJetFlavourAssociationLegacy)
# Not working with our PU sub, but keep it here for reference
# akSoftDropZ05B153CaloPatJetFlavourAssociation = akSoftDropZ05B153CalobTagger.PatJetFlavourAssociation
# akSoftDropZ05B153CaloPatJetFlavourId = cms.Sequence(
#     akSoftDropZ05B153CaloPatJetPartons
#     *
#     akSoftDropZ05B153CaloPatJetFlavourAssociation)

# akSoftDropZ05B153CaloJetBtaggingIP  = cms.Sequence(
#     akSoftDropZ05B153CaloImpactParameterTagInfos *
#     akSoftDropZ05B153CaloTrackCountingHighEffBJetTags +
#     akSoftDropZ05B153CaloTrackCountingHighPurBJetTags +
#     akSoftDropZ05B153CaloJetProbabilityBJetTags +
#     akSoftDropZ05B153CaloJetBProbabilityBJetTags 
# )

# akSoftDropZ05B153CaloJetBtaggingSV = cms.Sequence(
#     akSoftDropZ05B153CaloImpactParameterTagInfos
#     *
#     akSoftDropZ05B153CaloSecondaryVertexTagInfos
#     *
#     akSoftDropZ05B153CaloSimpleSecondaryVertexHighEffBJetTags
#     +
#     akSoftDropZ05B153CaloSimpleSecondaryVertexHighPurBJetTags
#     +
#     akSoftDropZ05B153CaloCombinedSecondaryVertexBJetTags
#     +
#     akSoftDropZ05B153CaloCombinedSecondaryVertexV2BJetTags
# )
# 
# akSoftDropZ05B153CaloJetBtaggingNegSV = cms.Sequence(
#     akSoftDropZ05B153CaloImpactParameterTagInfos
#     *
#     akSoftDropZ05B153CaloSecondaryVertexNegativeTagInfos
#     *
#     akSoftDropZ05B153CaloNegativeSimpleSecondaryVertexHighEffBJetTags
#     +
#     akSoftDropZ05B153CaloNegativeSimpleSecondaryVertexHighPurBJetTags
#     +
#     akSoftDropZ05B153CaloNegativeCombinedSecondaryVertexBJetTags
#     +
#     akSoftDropZ05B153CaloPositiveCombinedSecondaryVertexBJetTags
#     +
#     akSoftDropZ05B153CaloNegativeCombinedSecondaryVertexV2BJetTags
#     +
#     akSoftDropZ05B153CaloPositiveCombinedSecondaryVertexV2BJetTags
# )
# 
# akSoftDropZ05B153CaloJetBtaggingMu = cms.Sequence(
#     akSoftDropZ05B153CaloSoftPFMuonsTagInfos
#     *
#     akSoftDropZ05B153CaloSoftPFMuonBJetTags
#     +
#     akSoftDropZ05B153CaloSoftPFMuonByIP3dBJetTags
#     +
#     akSoftDropZ05B153CaloSoftPFMuonByPtBJetTags
#     +
#     akSoftDropZ05B153CaloNegativeSoftPFMuonByPtBJetTags
#     +
#     akSoftDropZ05B153CaloPositiveSoftPFMuonByPtBJetTags
# )
# 
# akSoftDropZ05B153CaloJetBtagging = cms.Sequence(
#     akSoftDropZ05B153CaloJetBtaggingIP
#     *
#     akSoftDropZ05B153CaloJetBtaggingSV
#     *
#     akSoftDropZ05B153CaloJetBtaggingNegSV
#     # *
#     # akSoftDropZ05B153CaloJetBtaggingMu
# )

akSoftDropZ05B153CalopatJetsWithBtagging = patJets.clone(
    jetSource              = cms.InputTag("akSoftDropZ05B153CaloJets"),
    genJetMatch            = cms.InputTag("akSoftDropZ05B153Calomatch"),
    genPartonMatch         = cms.InputTag("akSoftDropZ05B153Caloparton"),
    jetCorrFactorsSource   = cms.VInputTag(cms.InputTag("akSoftDropZ05B153Calocorr")),
#     JetPartonMapSource     = cms.InputTag("akSoftDropZ05B153CaloPatJetFlavourAssociationLegacy"),
#     JetFlavourInfoSource   = cms.InputTag("akSoftDropZ05B153CaloPatJetFlavourAssociation"),
#     trackAssociationSource = cms.InputTag("akSoftDropZ05B153CaloJetTracksAssociatorAtVertex"),
#     useLegacyJetMCFlavour  = True,
#     discriminatorSources   = cms.VInputTag(
#         cms.InputTag("akSoftDropZ05B153CaloSimpleSecondaryVertexHighEffBJetTags"),
#         cms.InputTag("akSoftDropZ05B153CaloSimpleSecondaryVertexHighPurBJetTags"),
#         cms.InputTag("akSoftDropZ05B153CaloCombinedSecondaryVertexBJetTags"),
#         cms.InputTag("akSoftDropZ05B153CaloCombinedSecondaryVertexV2BJetTags"),
#         cms.InputTag("akSoftDropZ05B153CaloJetBProbabilityBJetTags"),
#         cms.InputTag("akSoftDropZ05B153CaloJetProbabilityBJetTags"),
#         # cms.InputTag("akSoftDropZ05B153CaloSoftPFMuonByPtBJetTags"),
#         # cms.InputTag("akSoftDropZ05B153CaloSoftPFMuonByIP3dBJetTags"),
#         cms.InputTag("akSoftDropZ05B153CaloTrackCountingHighEffBJetTags"),
#         cms.InputTag("akSoftDropZ05B153CaloTrackCountingHighPurBJetTags"),
#     ),
    jetIDMap = cms.InputTag("akSoftDropZ05B153CaloJetID"),
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

akSoftDropZ05B153CaloNjettiness = Njettiness.clone(
    src = cms.InputTag("akSoftDropZ05B153CaloJets"),
    R0  = cms.double( 0.3)
)

# akSoftDropZ05B153CalopatJetsWithBtagging.userData.userFloats.src += [
#     'akSoftDropZ05B153CaloNjettiness:tau1',
#     'akSoftDropZ05B153CaloNjettiness:tau2',
#     'akSoftDropZ05B153CaloNjettiness:tau3'
# ]

akSoftDropZ05B153CaloJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("akSoftDropZ05B153CalopatJetsWithBtagging"),
    genjetTag = 'ak3HiSignalGenJets',
    rParam = 0.3,
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
    # bTagJetName = cms.untracked.string("akSoftDropZ05B153Calo"),
    jetName = cms.untracked.string("akSoftDropZ05B153Calo"),
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

akSoftDropZ05B153CaloJetSequence_mc = cms.Sequence(
    # akSoftDropZ05B153Caloclean
    # *
    akSoftDropZ05B153Calomatch
    *
    # akSoftDropZ05B153CalomatchGroomed
    # *
    akSoftDropZ05B153Caloparton
    *
    akSoftDropZ05B153Calocorr
    *
    # akSoftDropZ05B153CaloJetID
    # *
    # akSoftDropZ05B153CaloPatJetFlavourIdLegacy
    # *
    # Use legacy algo till PU implemented
    # akSoftDropZ05B153CaloPatJetFlavourId
    # *
    # akSoftDropZ05B153CaloJetTracksAssociatorAtVertex
    # *
    # akSoftDropZ05B153CaloJetBtagging
    # *
    # No constituents for calo jets in pp. Must be removed for pp calo jets but
    # I'm not sure how to do this transparently (Marta)
    akSoftDropZ05B153CaloNjettiness
    *
    akSoftDropZ05B153CalopatJetsWithBtagging
    *
    akSoftDropZ05B153CaloJetAnalyzer
)

akSoftDropZ05B153CaloJetSequence_data = cms.Sequence(
    akSoftDropZ05B153Calocorr
    *
    # akSoftDropZ05B153CaloJetID
    # *
    # akSoftDropZ05B153CaloJetTracksAssociatorAtVertex
    # *
    # akSoftDropZ05B153CaloJetBtagging
    # *
    akSoftDropZ05B153CaloNjettiness 
    *
    akSoftDropZ05B153CalopatJetsWithBtagging
    *
    akSoftDropZ05B153CaloJetAnalyzer
)

akSoftDropZ05B153CaloJetSequence_jec = cms.Sequence(
    akSoftDropZ05B153CaloJetSequence_mc)
akSoftDropZ05B153CaloJetSequence_mb = cms.Sequence(
    akSoftDropZ05B153CaloJetSequence_mc)

akSoftDropZ05B153CaloJetSequence = cms.Sequence(
    akSoftDropZ05B153CaloJetSequence_mc)
akSoftDropZ05B153CalopatJetsWithBtagging.userData.userFloats.src += ['akSoftDropZ05B153CaloJets:sym']
akSoftDropZ05B153CalopatJetsWithBtagging.userData.userInts.src += ['akSoftDropZ05B153CaloJets:droppedBranches']
