
import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akPuSoftDropZ05B152Calomatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akPuSoftDropZ05B152CaloJets"),
    matched = cms.InputTag("ak2HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.2
)

akPuSoftDropZ05B152CalomatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B152HiSignalGenJets"),
    matched = cms.InputTag("ak2HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.2
)

# akPuSoftDropZ05B152Caloparton = patJetPartonMatch.clone(
#     src = cms.InputTag("akPuSoftDropZ05B152CaloJets")
# )

akPuSoftDropZ05B152Calocorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
    # primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akPuSoftDropZ05B152CaloJets"),
    payload = "AKPu2Calo_offline"
)

akPuSoftDropZ05B152CaloJetID = cms.EDProducer(
    'JetIDProducer',
    JetIDParams,
    src = cms.InputTag('akPuSoftDropZ05B152CaloJets'))

# akPuSoftDropZ05B152Caloclean = heavyIonCleanedGenJets.clone(
#     src = cms.InputTag('ak2HiSignalGenJets'))

# akPuSoftDropZ05B152CalobTagger = bTaggers(
#     "akPuSoftDropZ05B152Calo", 0.2)

# create objects locally since they dont load properly otherwise
# akPuSoftDropZ05B152Calomatch = akPuSoftDropZ05B152CalobTagger.match
akPuSoftDropZ05B152Caloparton = patJetPartonMatch.clone(
    src = cms.InputTag("akPuSoftDropZ05B152CaloJets"),
    matched = cms.InputTag("genParticles"))
# akPuSoftDropZ05B152CaloPatJetFlavourAssociationLegacy = akPuSoftDropZ05B152CalobTagger.PatJetFlavourAssociationLegacy
# akPuSoftDropZ05B152CaloPatJetPartons = akPuSoftDropZ05B152CalobTagger.PatJetPartons
# akPuSoftDropZ05B152CaloJetTracksAssociatorAtVertex = akPuSoftDropZ05B152CalobTagger.JetTracksAssociatorAtVertex
# akPuSoftDropZ05B152CaloJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
# akPuSoftDropZ05B152CaloSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B152CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akPuSoftDropZ05B152CaloSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B152CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akPuSoftDropZ05B152CaloCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B152CalobTagger.CombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B152CaloCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B152CalobTagger.CombinedSecondaryVertexV2BJetTags
# akPuSoftDropZ05B152CaloJetBProbabilityBJetTags = akPuSoftDropZ05B152CalobTagger.JetBProbabilityBJetTags
# akPuSoftDropZ05B152CaloSoftPFMuonByPtBJetTags = akPuSoftDropZ05B152CalobTagger.SoftPFMuonByPtBJetTags
# akPuSoftDropZ05B152CaloSoftPFMuonByIP3dBJetTags = akPuSoftDropZ05B152CalobTagger.SoftPFMuonByIP3dBJetTags
# akPuSoftDropZ05B152CaloTrackCountingHighEffBJetTags = akPuSoftDropZ05B152CalobTagger.TrackCountingHighEffBJetTags
# akPuSoftDropZ05B152CaloTrackCountingHighPurBJetTags = akPuSoftDropZ05B152CalobTagger.TrackCountingHighPurBJetTags
# akPuSoftDropZ05B152CaloPatJetPartonAssociationLegacy = akPuSoftDropZ05B152CalobTagger.PatJetPartonAssociationLegacy

# akPuSoftDropZ05B152CaloImpactParameterTagInfos = akPuSoftDropZ05B152CalobTagger.ImpactParameterTagInfos
# akPuSoftDropZ05B152CaloImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akPuSoftDropZ05B152CaloJetProbabilityBJetTags = akPuSoftDropZ05B152CalobTagger.JetProbabilityBJetTags
# 
# akPuSoftDropZ05B152CaloSecondaryVertexTagInfos = akPuSoftDropZ05B152CalobTagger.SecondaryVertexTagInfos
# akPuSoftDropZ05B152CaloSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B152CalobTagger.SimpleSecondaryVertexHighEffBJetTags
# akPuSoftDropZ05B152CaloSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B152CalobTagger.SimpleSecondaryVertexHighPurBJetTags
# akPuSoftDropZ05B152CaloCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B152CalobTagger.CombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B152CaloCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B152CalobTagger.CombinedSecondaryVertexV2BJetTags
# 
# akPuSoftDropZ05B152CaloSecondaryVertexNegativeTagInfos = akPuSoftDropZ05B152CalobTagger.SecondaryVertexNegativeTagInfos
# akPuSoftDropZ05B152CaloNegativeSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B152CalobTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
# akPuSoftDropZ05B152CaloNegativeSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B152CalobTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
# akPuSoftDropZ05B152CaloNegativeCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B152CalobTagger.NegativeCombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B152CaloPositiveCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B152CalobTagger.PositiveCombinedSecondaryVertexBJetTags
# akPuSoftDropZ05B152CaloNegativeCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B152CalobTagger.NegativeCombinedSecondaryVertexV2BJetTags
# akPuSoftDropZ05B152CaloPositiveCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B152CalobTagger.PositiveCombinedSecondaryVertexV2BJetTags
# 
# akPuSoftDropZ05B152CaloSoftPFMuonsTagInfos = akPuSoftDropZ05B152CalobTagger.SoftPFMuonsTagInfos
# akPuSoftDropZ05B152CaloSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
# akPuSoftDropZ05B152CaloSoftPFMuonBJetTags = akPuSoftDropZ05B152CalobTagger.SoftPFMuonBJetTags
# akPuSoftDropZ05B152CaloSoftPFMuonByIP3dBJetTags = akPuSoftDropZ05B152CalobTagger.SoftPFMuonByIP3dBJetTags
# akPuSoftDropZ05B152CaloSoftPFMuonByPtBJetTags = akPuSoftDropZ05B152CalobTagger.SoftPFMuonByPtBJetTags
# akPuSoftDropZ05B152CaloNegativeSoftPFMuonByPtBJetTags = akPuSoftDropZ05B152CalobTagger.NegativeSoftPFMuonByPtBJetTags
# akPuSoftDropZ05B152CaloPositiveSoftPFMuonByPtBJetTags = akPuSoftDropZ05B152CalobTagger.PositiveSoftPFMuonByPtBJetTags
# akPuSoftDropZ05B152CaloPatJetFlavourIdLegacy = cms.Sequence(akPuSoftDropZ05B152CaloPatJetPartonAssociationLegacy*akPuSoftDropZ05B152CaloPatJetFlavourAssociationLegacy)
# Not working with our PU sub, but keep it here for reference
# akPuSoftDropZ05B152CaloPatJetFlavourAssociation = akPuSoftDropZ05B152CalobTagger.PatJetFlavourAssociation
# akPuSoftDropZ05B152CaloPatJetFlavourId = cms.Sequence(
#     akPuSoftDropZ05B152CaloPatJetPartons
#     *
#     akPuSoftDropZ05B152CaloPatJetFlavourAssociation)

# akPuSoftDropZ05B152CaloJetBtaggingIP  = cms.Sequence(
#     akPuSoftDropZ05B152CaloImpactParameterTagInfos *
#     akPuSoftDropZ05B152CaloTrackCountingHighEffBJetTags +
#     akPuSoftDropZ05B152CaloTrackCountingHighPurBJetTags +
#     akPuSoftDropZ05B152CaloJetProbabilityBJetTags +
#     akPuSoftDropZ05B152CaloJetBProbabilityBJetTags 
# )

# akPuSoftDropZ05B152CaloJetBtaggingSV = cms.Sequence(
#     akPuSoftDropZ05B152CaloImpactParameterTagInfos
#     *
#     akPuSoftDropZ05B152CaloSecondaryVertexTagInfos
#     *
#     akPuSoftDropZ05B152CaloSimpleSecondaryVertexHighEffBJetTags
#     +
#     akPuSoftDropZ05B152CaloSimpleSecondaryVertexHighPurBJetTags
#     +
#     akPuSoftDropZ05B152CaloCombinedSecondaryVertexBJetTags
#     +
#     akPuSoftDropZ05B152CaloCombinedSecondaryVertexV2BJetTags
# )
# 
# akPuSoftDropZ05B152CaloJetBtaggingNegSV = cms.Sequence(
#     akPuSoftDropZ05B152CaloImpactParameterTagInfos
#     *
#     akPuSoftDropZ05B152CaloSecondaryVertexNegativeTagInfos
#     *
#     akPuSoftDropZ05B152CaloNegativeSimpleSecondaryVertexHighEffBJetTags
#     +
#     akPuSoftDropZ05B152CaloNegativeSimpleSecondaryVertexHighPurBJetTags
#     +
#     akPuSoftDropZ05B152CaloNegativeCombinedSecondaryVertexBJetTags
#     +
#     akPuSoftDropZ05B152CaloPositiveCombinedSecondaryVertexBJetTags
#     +
#     akPuSoftDropZ05B152CaloNegativeCombinedSecondaryVertexV2BJetTags
#     +
#     akPuSoftDropZ05B152CaloPositiveCombinedSecondaryVertexV2BJetTags
# )
# 
# akPuSoftDropZ05B152CaloJetBtaggingMu = cms.Sequence(
#     akPuSoftDropZ05B152CaloSoftPFMuonsTagInfos
#     *
#     akPuSoftDropZ05B152CaloSoftPFMuonBJetTags
#     +
#     akPuSoftDropZ05B152CaloSoftPFMuonByIP3dBJetTags
#     +
#     akPuSoftDropZ05B152CaloSoftPFMuonByPtBJetTags
#     +
#     akPuSoftDropZ05B152CaloNegativeSoftPFMuonByPtBJetTags
#     +
#     akPuSoftDropZ05B152CaloPositiveSoftPFMuonByPtBJetTags
# )
# 
# akPuSoftDropZ05B152CaloJetBtagging = cms.Sequence(
#     akPuSoftDropZ05B152CaloJetBtaggingIP
#     *
#     akPuSoftDropZ05B152CaloJetBtaggingSV
#     *
#     akPuSoftDropZ05B152CaloJetBtaggingNegSV
#     # *
#     # akPuSoftDropZ05B152CaloJetBtaggingMu
# )

akPuSoftDropZ05B152CalopatJetsWithBtagging = patJets.clone(
    jetSource              = cms.InputTag("akPuSoftDropZ05B152CaloJets"),
    genJetMatch            = cms.InputTag("akPuSoftDropZ05B152Calomatch"),
    genPartonMatch         = cms.InputTag("akPuSoftDropZ05B152Caloparton"),
    jetCorrFactorsSource   = cms.VInputTag(cms.InputTag("akPuSoftDropZ05B152Calocorr")),
#     JetPartonMapSource     = cms.InputTag("akPuSoftDropZ05B152CaloPatJetFlavourAssociationLegacy"),
#     JetFlavourInfoSource   = cms.InputTag("akPuSoftDropZ05B152CaloPatJetFlavourAssociation"),
#     trackAssociationSource = cms.InputTag("akPuSoftDropZ05B152CaloJetTracksAssociatorAtVertex"),
#     useLegacyJetMCFlavour  = True,
#     discriminatorSources   = cms.VInputTag(
#         cms.InputTag("akPuSoftDropZ05B152CaloSimpleSecondaryVertexHighEffBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B152CaloSimpleSecondaryVertexHighPurBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B152CaloCombinedSecondaryVertexBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B152CaloCombinedSecondaryVertexV2BJetTags"),
#         cms.InputTag("akPuSoftDropZ05B152CaloJetBProbabilityBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B152CaloJetProbabilityBJetTags"),
#         # cms.InputTag("akPuSoftDropZ05B152CaloSoftPFMuonByPtBJetTags"),
#         # cms.InputTag("akPuSoftDropZ05B152CaloSoftPFMuonByIP3dBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B152CaloTrackCountingHighEffBJetTags"),
#         cms.InputTag("akPuSoftDropZ05B152CaloTrackCountingHighPurBJetTags"),
#     ),
    jetIDMap = cms.InputTag("akPuSoftDropZ05B152CaloJetID"),
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

akPuSoftDropZ05B152CaloNjettiness = Njettiness.clone(
    src = cms.InputTag("akPuSoftDropZ05B152CaloJets"),
    R0  = cms.double( 0.2)
)

# akPuSoftDropZ05B152CalopatJetsWithBtagging.userData.userFloats.src += [
#     'akPuSoftDropZ05B152CaloNjettiness:tau1',
#     'akPuSoftDropZ05B152CaloNjettiness:tau2',
#     'akPuSoftDropZ05B152CaloNjettiness:tau3'
# ]

akPuSoftDropZ05B152CaloJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("akPuSoftDropZ05B152CalopatJetsWithBtagging"),
    genjetTag = 'ak2HiSignalGenJets',
    rParam = 0.2,
    matchJets = cms.untracked.bool(False),
    matchTag = 'patJetsWithBtagging',
    pfCandidateLabel = cms.untracked.InputTag('particleFlow'),
    trackTag = cms.InputTag("generalTracks"),
    fillGenJets = False,
    isMC = False,
    doSubEvent = False,
    useHepMC = cms.untracked.bool(False),
    genParticles = cms.untracked.InputTag("genParticles"),
    eventInfoTag = cms.InputTag("generator"),
    # doLifeTimeTagging = cms.untracked.bool(True),
    # doLifeTimeTaggingExtras = cms.untracked.bool(False),
    # bTagJetName = cms.untracked.string("akPuSoftDropZ05B152Calo"),
    jetName = cms.untracked.string("akPuSoftDropZ05B152Calo"),
    genPtMin = cms.untracked.double(5),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
    doTower = cms.untracked.bool(True),
    doSubJets = cms.untracked.bool(True),
    doGenSubJets = cms.untracked.bool(False),     
    subjetGenTag = cms.untracked.InputTag("akSoftDropZ05B152GenJets"),
    doGenTaus = cms.untracked.bool(False),
    genTau1 = cms.InputTag("akSoftDropZ05B152GenNjettiness","tau1"),
    genTau2 = cms.InputTag("akSoftDropZ05B152GenNjettiness","tau2"),
    genTau3 = cms.InputTag("akSoftDropZ05B152GenNjettiness","tau3"),
    doGenSym = cms.untracked.bool(False),
    genSym = cms.InputTag("akSoftDropZ05B152GenJets","sym"),
    genDroppedBranches = cms.InputTag("akSoftDropZ05B152GenJets","droppedBranches")
)

akPuSoftDropZ05B152CaloJetSequence_mc = cms.Sequence(
    # akPuSoftDropZ05B152Caloclean
    # *
    akPuSoftDropZ05B152Calomatch
    *
    # akPuSoftDropZ05B152CalomatchGroomed
    # *
    akPuSoftDropZ05B152Caloparton
    *
    akPuSoftDropZ05B152Calocorr
    *
    # akPuSoftDropZ05B152CaloJetID
    # *
    # akPuSoftDropZ05B152CaloPatJetFlavourIdLegacy
    # *
    # Use legacy algo till PU implemented
    # akPuSoftDropZ05B152CaloPatJetFlavourId
    # *
    # akPuSoftDropZ05B152CaloJetTracksAssociatorAtVertex
    # *
    # akPuSoftDropZ05B152CaloJetBtagging
    # *
    # No constituents for calo jets in pp. Must be removed for pp calo jets but
    # I'm not sure how to do this transparently (Marta)
    akPuSoftDropZ05B152CaloNjettiness
    *
    akPuSoftDropZ05B152CalopatJetsWithBtagging
    *
    akPuSoftDropZ05B152CaloJetAnalyzer
)

akPuSoftDropZ05B152CaloJetSequence_data = cms.Sequence(
    akPuSoftDropZ05B152Calocorr
    *
    # akPuSoftDropZ05B152CaloJetID
    # *
    # akPuSoftDropZ05B152CaloJetTracksAssociatorAtVertex
    # *
    # akPuSoftDropZ05B152CaloJetBtagging
    # *
    akPuSoftDropZ05B152CaloNjettiness 
    *
    akPuSoftDropZ05B152CalopatJetsWithBtagging
    *
    akPuSoftDropZ05B152CaloJetAnalyzer
)

akPuSoftDropZ05B152CaloJetSequence_jec = cms.Sequence(
    akPuSoftDropZ05B152CaloJetSequence_mc)
akPuSoftDropZ05B152CaloJetSequence_mb = cms.Sequence(
    akPuSoftDropZ05B152CaloJetSequence_mc)

akPuSoftDropZ05B152CaloJetSequence = cms.Sequence(
    akPuSoftDropZ05B152CaloJetSequence_data)
akPuSoftDropZ05B152CalopatJetsWithBtagging.userData.userFloats.src += ['akPuSoftDropZ05B152CaloJets:sym']
akPuSoftDropZ05B152CalopatJetsWithBtagging.userData.userInts.src += ['akPuSoftDropZ05B152CaloJets:droppedBranches']
