
import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akSoftDropZ05B152Calomatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B152CaloJets"),
    matched = cms.InputTag("ak2HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.2
)

akSoftDropZ05B152CalomatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B152HiSignalGenJets"),
    matched = cms.InputTag("ak2HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.2
)

# akSoftDropZ05B152Caloparton = patJetPartonMatch.clone(
#     src = cms.InputTag("akSoftDropZ05B152CaloJets")
# )

akSoftDropZ05B152Calocorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
    # primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akSoftDropZ05B152CaloJets"),
    payload = "AK2Calo_offline"
)

akSoftDropZ05B152CaloJetID = cms.EDProducer(
    'JetIDProducer',
    JetIDParams,
    src = cms.InputTag('akSoftDropZ05B152CaloJets'))

# akSoftDropZ05B152Caloclean = heavyIonCleanedGenJets.clone(
#     src = cms.InputTag('ak2HiSignalGenJets'))

# akSoftDropZ05B152CalobTagger = bTaggers(
#     "akSoftDropZ05B152Calo", 0.2)

# create objects locally since they dont load properly otherwise
# akSoftDropZ05B152Calomatch = akSoftDropZ05B152CalobTagger.match
akSoftDropZ05B152Caloparton = patJetPartonMatch.clone(
    src = cms.InputTag("akSoftDropZ05B152CaloJets"),
    matched = cms.InputTag("genParticles"))
# akSoftDropZ05B152CaloPatJetFlavourAssociationLegacy = akSoftDropZ05B152CalobTagger.PatJetFlavourAssociationLegacy
# akSoftDropZ05B152CaloPatJetPartons = akSoftDropZ05B152CalobTagger.PatJetPartons
# akSoftDropZ05B152CaloJetTracksAssociatorAtVertex = akSoftDropZ05B152CalobTagger.JetTracksAssociatorAtVertex
# akSoftDropZ05B152CaloJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
# akSoftDropZ05B152CaloSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B152CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akSoftDropZ05B152CaloSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B152CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akSoftDropZ05B152CaloCombinedSecondaryVertexBJetTags = akSoftDropZ05B152CalobTagger.CombinedSecondaryVertexBJetTags
# akSoftDropZ05B152CaloCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B152CalobTagger.CombinedSecondaryVertexV2BJetTags
# akSoftDropZ05B152CaloJetBProbabilityBJetTags = akSoftDropZ05B152CalobTagger.JetBProbabilityBJetTags
# akSoftDropZ05B152CaloSoftPFMuonByPtBJetTags = akSoftDropZ05B152CalobTagger.SoftPFMuonByPtBJetTags
# akSoftDropZ05B152CaloSoftPFMuonByIP3dBJetTags = akSoftDropZ05B152CalobTagger.SoftPFMuonByIP3dBJetTags
# akSoftDropZ05B152CaloTrackCountingHighEffBJetTags = akSoftDropZ05B152CalobTagger.TrackCountingHighEffBJetTags
# akSoftDropZ05B152CaloTrackCountingHighPurBJetTags = akSoftDropZ05B152CalobTagger.TrackCountingHighPurBJetTags
# akSoftDropZ05B152CaloPatJetPartonAssociationLegacy = akSoftDropZ05B152CalobTagger.PatJetPartonAssociationLegacy

# akSoftDropZ05B152CaloImpactParameterTagInfos = akSoftDropZ05B152CalobTagger.ImpactParameterTagInfos
# akSoftDropZ05B152CaloImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akSoftDropZ05B152CaloJetProbabilityBJetTags = akSoftDropZ05B152CalobTagger.JetProbabilityBJetTags
# 
# akSoftDropZ05B152CaloSecondaryVertexTagInfos = akSoftDropZ05B152CalobTagger.SecondaryVertexTagInfos
# akSoftDropZ05B152CaloSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B152CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akSoftDropZ05B152CaloSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B152CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akSoftDropZ05B152CaloCombinedSecondaryVertexBJetTags = akSoftDropZ05B152CalobTagger.CombinedSecondaryVertexBJetTags
# akSoftDropZ05B152CaloCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B152CalobTagger.CombinedSecondaryVertexV2BJetTags
# 
# akSoftDropZ05B152CaloSecondaryVertexNegativeTagInfos = akSoftDropZ05B152CalobTagger.SecondaryVertexNegativeTagInfos
# akSoftDropZ05B152CaloNegativeSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B152CalobTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
# akSoftDropZ05B152CaloNegativeSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B152CalobTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
# akSoftDropZ05B152CaloNegativeCombinedSecondaryVertexBJetTags = akSoftDropZ05B152CalobTagger.NegativeCombinedSecondaryVertexBJetTags
# akSoftDropZ05B152CaloPositiveCombinedSecondaryVertexBJetTags = akSoftDropZ05B152CalobTagger.PositiveCombinedSecondaryVertexBJetTags
# akSoftDropZ05B152CaloNegativeCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B152CalobTagger.NegativeCombinedSecondaryVertexV2BJetTags
# akSoftDropZ05B152CaloPositiveCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B152CalobTagger.PositiveCombinedSecondaryVertexV2BJetTags
# 
# akSoftDropZ05B152CaloSoftPFMuonsTagInfos = akSoftDropZ05B152CalobTagger.SoftPFMuonsTagInfos
# akSoftDropZ05B152CaloSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akSoftDropZ05B152CaloSoftPFMuonBJetTags = akSoftDropZ05B152CalobTagger.SoftPFMuonBJetTags
# akSoftDropZ05B152CaloSoftPFMuonByIP3dBJetTags = akSoftDropZ05B152CalobTagger.SoftPFMuonByIP3dBJetTags
# akSoftDropZ05B152CaloSoftPFMuonByPtBJetTags = akSoftDropZ05B152CalobTagger.SoftPFMuonByPtBJetTags
# akSoftDropZ05B152CaloNegativeSoftPFMuonByPtBJetTags = akSoftDropZ05B152CalobTagger.NegativeSoftPFMuonByPtBJetTags
# akSoftDropZ05B152CaloPositiveSoftPFMuonByPtBJetTags = akSoftDropZ05B152CalobTagger.PositiveSoftPFMuonByPtBJetTags
# akSoftDropZ05B152CaloPatJetFlavourIdLegacy = cms.Sequence(akSoftDropZ05B152CaloPatJetPartonAssociationLegacy*akSoftDropZ05B152CaloPatJetFlavourAssociationLegacy)
# Not working with our PU sub, but keep it here for reference
# akSoftDropZ05B152CaloPatJetFlavourAssociation = akSoftDropZ05B152CalobTagger.PatJetFlavourAssociation
# akSoftDropZ05B152CaloPatJetFlavourId = cms.Sequence(
#     akSoftDropZ05B152CaloPatJetPartons
#     *
#     akSoftDropZ05B152CaloPatJetFlavourAssociation)

# akSoftDropZ05B152CaloJetBtaggingIP  = cms.Sequence(
#     akSoftDropZ05B152CaloImpactParameterTagInfos *
#     akSoftDropZ05B152CaloTrackCountingHighEffBJetTags +
#     akSoftDropZ05B152CaloTrackCountingHighPurBJetTags +
#     akSoftDropZ05B152CaloJetProbabilityBJetTags +
#     akSoftDropZ05B152CaloJetBProbabilityBJetTags 
# )

# akSoftDropZ05B152CaloJetBtaggingSV = cms.Sequence(
#     akSoftDropZ05B152CaloImpactParameterTagInfos
#     *
#     akSoftDropZ05B152CaloSecondaryVertexTagInfos
#     *
#     akSoftDropZ05B152CaloSimpleSecondaryVertexHighEffBJetTags
#     +
#     akSoftDropZ05B152CaloSimpleSecondaryVertexHighPurBJetTags
#     +
#     akSoftDropZ05B152CaloCombinedSecondaryVertexBJetTags
#     +
#     akSoftDropZ05B152CaloCombinedSecondaryVertexV2BJetTags
# )
# 
# akSoftDropZ05B152CaloJetBtaggingNegSV = cms.Sequence(
#     akSoftDropZ05B152CaloImpactParameterTagInfos
#     *
#     akSoftDropZ05B152CaloSecondaryVertexNegativeTagInfos
#     *
#     akSoftDropZ05B152CaloNegativeSimpleSecondaryVertexHighEffBJetTags
#     +
#     akSoftDropZ05B152CaloNegativeSimpleSecondaryVertexHighPurBJetTags
#     +
#     akSoftDropZ05B152CaloNegativeCombinedSecondaryVertexBJetTags
#     +
#     akSoftDropZ05B152CaloPositiveCombinedSecondaryVertexBJetTags
#     +
#     akSoftDropZ05B152CaloNegativeCombinedSecondaryVertexV2BJetTags
#     +
#     akSoftDropZ05B152CaloPositiveCombinedSecondaryVertexV2BJetTags
# )
# 
# akSoftDropZ05B152CaloJetBtaggingMu = cms.Sequence(
#     akSoftDropZ05B152CaloSoftPFMuonsTagInfos
#     *
#     akSoftDropZ05B152CaloSoftPFMuonBJetTags
#     +
#     akSoftDropZ05B152CaloSoftPFMuonByIP3dBJetTags
#     +
#     akSoftDropZ05B152CaloSoftPFMuonByPtBJetTags
#     +
#     akSoftDropZ05B152CaloNegativeSoftPFMuonByPtBJetTags
#     +
#     akSoftDropZ05B152CaloPositiveSoftPFMuonByPtBJetTags
# )
# 
# akSoftDropZ05B152CaloJetBtagging = cms.Sequence(
#     akSoftDropZ05B152CaloJetBtaggingIP
#     *
#     akSoftDropZ05B152CaloJetBtaggingSV
#     *
#     akSoftDropZ05B152CaloJetBtaggingNegSV
#     # *
#     # akSoftDropZ05B152CaloJetBtaggingMu
# )

akSoftDropZ05B152CalopatJetsWithBtagging = patJets.clone(
    jetSource              = cms.InputTag("akSoftDropZ05B152CaloJets"),
    genJetMatch            = cms.InputTag("akSoftDropZ05B152Calomatch"),
    genPartonMatch         = cms.InputTag("akSoftDropZ05B152Caloparton"),
    jetCorrFactorsSource   = cms.VInputTag(cms.InputTag("akSoftDropZ05B152Calocorr")),
#     JetPartonMapSource     = cms.InputTag("akSoftDropZ05B152CaloPatJetFlavourAssociationLegacy"),
#     JetFlavourInfoSource   = cms.InputTag("akSoftDropZ05B152CaloPatJetFlavourAssociation"),
#     trackAssociationSource = cms.InputTag("akSoftDropZ05B152CaloJetTracksAssociatorAtVertex"),
#     useLegacyJetMCFlavour  = True,
#     discriminatorSources   = cms.VInputTag(
#         cms.InputTag("akSoftDropZ05B152CaloSimpleSecondaryVertexHighEffBJetTags"),
#         cms.InputTag("akSoftDropZ05B152CaloSimpleSecondaryVertexHighPurBJetTags"),
#         cms.InputTag("akSoftDropZ05B152CaloCombinedSecondaryVertexBJetTags"),
#         cms.InputTag("akSoftDropZ05B152CaloCombinedSecondaryVertexV2BJetTags"),
#         cms.InputTag("akSoftDropZ05B152CaloJetBProbabilityBJetTags"),
#         cms.InputTag("akSoftDropZ05B152CaloJetProbabilityBJetTags"),
#         # cms.InputTag("akSoftDropZ05B152CaloSoftPFMuonByPtBJetTags"),
#         # cms.InputTag("akSoftDropZ05B152CaloSoftPFMuonByIP3dBJetTags"),
#         cms.InputTag("akSoftDropZ05B152CaloTrackCountingHighEffBJetTags"),
#         cms.InputTag("akSoftDropZ05B152CaloTrackCountingHighPurBJetTags"),
#     ),
    jetIDMap = cms.InputTag("akSoftDropZ05B152CaloJetID"),
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

akSoftDropZ05B152CaloNjettiness = Njettiness.clone(
    src = cms.InputTag("akSoftDropZ05B152CaloJets"),
    R0  = cms.double( 0.2)
)

# akSoftDropZ05B152CalopatJetsWithBtagging.userData.userFloats.src += [
#     'akSoftDropZ05B152CaloNjettiness:tau1',
#     'akSoftDropZ05B152CaloNjettiness:tau2',
#     'akSoftDropZ05B152CaloNjettiness:tau3'
# ]

akSoftDropZ05B152CaloJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("akSoftDropZ05B152CalopatJetsWithBtagging"),
    genjetTag = 'ak2HiSignalGenJets',
    rParam = 0.2,
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
    # bTagJetName = cms.untracked.string("akSoftDropZ05B152Calo"),
    jetName = cms.untracked.string("akSoftDropZ05B152Calo"),
    genPtMin = cms.untracked.double(5),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
    doTower = cms.untracked.bool(True),
    doSubJets = cms.untracked.bool(True),
    doGenSubJets = cms.untracked.bool(True),     
    subjetGenTag = cms.untracked.InputTag("akSoftDropZ05B152GenJets"),
    doGenTaus = cms.untracked.bool(False),
    genTau1 = cms.InputTag("akSoftDropZ05B152GenNjettiness","tau1"),
    genTau2 = cms.InputTag("akSoftDropZ05B152GenNjettiness","tau2"),
    genTau3 = cms.InputTag("akSoftDropZ05B152GenNjettiness","tau3"),
    doGenSym = cms.untracked.bool(True),
    genSym = cms.InputTag("akSoftDropZ05B152GenJets","sym"),
    genDroppedBranches = cms.InputTag("akSoftDropZ05B152GenJets","droppedBranches")
)

akSoftDropZ05B152CaloJetSequence_mc = cms.Sequence(
    # akSoftDropZ05B152Caloclean
    # *
    akSoftDropZ05B152Calomatch
    *
    # akSoftDropZ05B152CalomatchGroomed
    # *
    akSoftDropZ05B152Caloparton
    *
    akSoftDropZ05B152Calocorr
    *
    # akSoftDropZ05B152CaloJetID
    # *
    # akSoftDropZ05B152CaloPatJetFlavourIdLegacy
    # *
    # Use legacy algo till PU implemented
    # akSoftDropZ05B152CaloPatJetFlavourId
    # *
    # akSoftDropZ05B152CaloJetTracksAssociatorAtVertex
    # *
    # akSoftDropZ05B152CaloJetBtagging
    # *
    # No constituents for calo jets in pp. Must be removed for pp calo jets but
    # I'm not sure how to do this transparently (Marta)
    akSoftDropZ05B152CaloNjettiness
    *
    akSoftDropZ05B152CalopatJetsWithBtagging
    *
    akSoftDropZ05B152CaloJetAnalyzer
)

akSoftDropZ05B152CaloJetSequence_data = cms.Sequence(
    akSoftDropZ05B152Calocorr
    *
    # akSoftDropZ05B152CaloJetID
    # *
    # akSoftDropZ05B152CaloJetTracksAssociatorAtVertex
    # *
    # akSoftDropZ05B152CaloJetBtagging
    # *
    akSoftDropZ05B152CaloNjettiness 
    *
    akSoftDropZ05B152CalopatJetsWithBtagging
    *
    akSoftDropZ05B152CaloJetAnalyzer
)

akSoftDropZ05B152CaloJetSequence_jec = cms.Sequence(
    akSoftDropZ05B152CaloJetSequence_mc)
akSoftDropZ05B152CaloJetSequence_mb = cms.Sequence(
    akSoftDropZ05B152CaloJetSequence_mc)

akSoftDropZ05B152CaloJetSequence = cms.Sequence(
    akSoftDropZ05B152CaloJetSequence_mc)
akSoftDropZ05B152CalopatJetsWithBtagging.userData.userFloats.src += ['akSoftDropZ05B152CaloJets:sym']
akSoftDropZ05B152CalopatJetsWithBtagging.userData.userInts.src += ['akSoftDropZ05B152CaloJets:droppedBranches']
