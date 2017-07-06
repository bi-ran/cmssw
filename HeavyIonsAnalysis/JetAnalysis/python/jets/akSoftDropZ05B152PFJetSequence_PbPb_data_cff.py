

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akSoftDropZ05B152PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B152PFJets"),
    matched = cms.InputTag("ak2HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.2
    )

akSoftDropZ05B152PFmatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B152HiSignalGenJets"),
    matched = cms.InputTag("ak2HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.2
    )

akSoftDropZ05B152PFparton = patJetPartonMatch.clone(src = cms.InputTag("akSoftDropZ05B152PFJets")
                                                        )

akSoftDropZ05B152PFcorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akSoftDropZ05B152PFJets"),
    payload = "AK2PF_offline"
    )

akSoftDropZ05B152PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akSoftDropZ05B152CaloJets'))

#akSoftDropZ05B152PFclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak2HiSignalGenJets'))

akSoftDropZ05B152PFbTagger = bTaggers("akSoftDropZ05B152PF",0.2)

#create objects locally since they dont load properly otherwise
#akSoftDropZ05B152PFmatch = akSoftDropZ05B152PFbTagger.match
akSoftDropZ05B152PFparton = patJetPartonMatch.clone(src = cms.InputTag("akSoftDropZ05B152PFJets"), matched = cms.InputTag("hiSignalGenParticles"))
akSoftDropZ05B152PFPatJetFlavourAssociationLegacy = akSoftDropZ05B152PFbTagger.PatJetFlavourAssociationLegacy
akSoftDropZ05B152PFPatJetPartons = akSoftDropZ05B152PFbTagger.PatJetPartons
akSoftDropZ05B152PFJetTracksAssociatorAtVertex = akSoftDropZ05B152PFbTagger.JetTracksAssociatorAtVertex
akSoftDropZ05B152PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akSoftDropZ05B152PFSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B152PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akSoftDropZ05B152PFSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B152PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akSoftDropZ05B152PFCombinedSecondaryVertexBJetTags = akSoftDropZ05B152PFbTagger.CombinedSecondaryVertexBJetTags
akSoftDropZ05B152PFCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B152PFbTagger.CombinedSecondaryVertexV2BJetTags
akSoftDropZ05B152PFJetBProbabilityBJetTags = akSoftDropZ05B152PFbTagger.JetBProbabilityBJetTags
akSoftDropZ05B152PFSoftPFMuonByPtBJetTags = akSoftDropZ05B152PFbTagger.SoftPFMuonByPtBJetTags
akSoftDropZ05B152PFSoftPFMuonByIP3dBJetTags = akSoftDropZ05B152PFbTagger.SoftPFMuonByIP3dBJetTags
akSoftDropZ05B152PFTrackCountingHighEffBJetTags = akSoftDropZ05B152PFbTagger.TrackCountingHighEffBJetTags
akSoftDropZ05B152PFTrackCountingHighPurBJetTags = akSoftDropZ05B152PFbTagger.TrackCountingHighPurBJetTags
akSoftDropZ05B152PFPatJetPartonAssociationLegacy = akSoftDropZ05B152PFbTagger.PatJetPartonAssociationLegacy

akSoftDropZ05B152PFImpactParameterTagInfos = akSoftDropZ05B152PFbTagger.ImpactParameterTagInfos
akSoftDropZ05B152PFImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akSoftDropZ05B152PFJetProbabilityBJetTags = akSoftDropZ05B152PFbTagger.JetProbabilityBJetTags

akSoftDropZ05B152PFSecondaryVertexTagInfos = akSoftDropZ05B152PFbTagger.SecondaryVertexTagInfos
akSoftDropZ05B152PFSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B152PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akSoftDropZ05B152PFSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B152PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akSoftDropZ05B152PFCombinedSecondaryVertexBJetTags = akSoftDropZ05B152PFbTagger.CombinedSecondaryVertexBJetTags
akSoftDropZ05B152PFCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B152PFbTagger.CombinedSecondaryVertexV2BJetTags

akSoftDropZ05B152PFSecondaryVertexNegativeTagInfos = akSoftDropZ05B152PFbTagger.SecondaryVertexNegativeTagInfos
akSoftDropZ05B152PFNegativeSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B152PFbTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akSoftDropZ05B152PFNegativeSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B152PFbTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akSoftDropZ05B152PFNegativeCombinedSecondaryVertexBJetTags = akSoftDropZ05B152PFbTagger.NegativeCombinedSecondaryVertexBJetTags
akSoftDropZ05B152PFPositiveCombinedSecondaryVertexBJetTags = akSoftDropZ05B152PFbTagger.PositiveCombinedSecondaryVertexBJetTags
akSoftDropZ05B152PFNegativeCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B152PFbTagger.NegativeCombinedSecondaryVertexV2BJetTags
akSoftDropZ05B152PFPositiveCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B152PFbTagger.PositiveCombinedSecondaryVertexV2BJetTags

akSoftDropZ05B152PFSoftPFMuonsTagInfos = akSoftDropZ05B152PFbTagger.SoftPFMuonsTagInfos
akSoftDropZ05B152PFSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akSoftDropZ05B152PFSoftPFMuonBJetTags = akSoftDropZ05B152PFbTagger.SoftPFMuonBJetTags
akSoftDropZ05B152PFSoftPFMuonByIP3dBJetTags = akSoftDropZ05B152PFbTagger.SoftPFMuonByIP3dBJetTags
akSoftDropZ05B152PFSoftPFMuonByPtBJetTags = akSoftDropZ05B152PFbTagger.SoftPFMuonByPtBJetTags
akSoftDropZ05B152PFNegativeSoftPFMuonByPtBJetTags = akSoftDropZ05B152PFbTagger.NegativeSoftPFMuonByPtBJetTags
akSoftDropZ05B152PFPositiveSoftPFMuonByPtBJetTags = akSoftDropZ05B152PFbTagger.PositiveSoftPFMuonByPtBJetTags
akSoftDropZ05B152PFPatJetFlavourIdLegacy = cms.Sequence(akSoftDropZ05B152PFPatJetPartonAssociationLegacy*akSoftDropZ05B152PFPatJetFlavourAssociationLegacy)
#Not working with our PU sub, but keep it here for reference
#akSoftDropZ05B152PFPatJetFlavourAssociation = akSoftDropZ05B152PFbTagger.PatJetFlavourAssociation
#akSoftDropZ05B152PFPatJetFlavourId = cms.Sequence(akSoftDropZ05B152PFPatJetPartons*akSoftDropZ05B152PFPatJetFlavourAssociation)

akSoftDropZ05B152PFJetBtaggingIP       = cms.Sequence(akSoftDropZ05B152PFImpactParameterTagInfos *
            (akSoftDropZ05B152PFTrackCountingHighEffBJetTags +
             akSoftDropZ05B152PFTrackCountingHighPurBJetTags +
             akSoftDropZ05B152PFJetProbabilityBJetTags +
             akSoftDropZ05B152PFJetBProbabilityBJetTags 
            )
            )

akSoftDropZ05B152PFJetBtaggingSV = cms.Sequence(akSoftDropZ05B152PFImpactParameterTagInfos
            *
            akSoftDropZ05B152PFSecondaryVertexTagInfos
            * (akSoftDropZ05B152PFSimpleSecondaryVertexHighEffBJetTags+
                akSoftDropZ05B152PFSimpleSecondaryVertexHighPurBJetTags+
                akSoftDropZ05B152PFCombinedSecondaryVertexBJetTags+
                akSoftDropZ05B152PFCombinedSecondaryVertexV2BJetTags
              )
            )

akSoftDropZ05B152PFJetBtaggingNegSV = cms.Sequence(akSoftDropZ05B152PFImpactParameterTagInfos
            *
            akSoftDropZ05B152PFSecondaryVertexNegativeTagInfos
            * (akSoftDropZ05B152PFNegativeSimpleSecondaryVertexHighEffBJetTags+
                akSoftDropZ05B152PFNegativeSimpleSecondaryVertexHighPurBJetTags+
                akSoftDropZ05B152PFNegativeCombinedSecondaryVertexBJetTags+
                akSoftDropZ05B152PFPositiveCombinedSecondaryVertexBJetTags+
                akSoftDropZ05B152PFNegativeCombinedSecondaryVertexV2BJetTags+
                akSoftDropZ05B152PFPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akSoftDropZ05B152PFJetBtaggingMu = cms.Sequence(akSoftDropZ05B152PFSoftPFMuonsTagInfos * (akSoftDropZ05B152PFSoftPFMuonBJetTags
                +
                akSoftDropZ05B152PFSoftPFMuonByIP3dBJetTags
                +
                akSoftDropZ05B152PFSoftPFMuonByPtBJetTags
                +
                akSoftDropZ05B152PFNegativeSoftPFMuonByPtBJetTags
                +
                akSoftDropZ05B152PFPositiveSoftPFMuonByPtBJetTags
              )
            )

akSoftDropZ05B152PFJetBtagging = cms.Sequence(akSoftDropZ05B152PFJetBtaggingIP
            *akSoftDropZ05B152PFJetBtaggingSV
            *akSoftDropZ05B152PFJetBtaggingNegSV
#            *akSoftDropZ05B152PFJetBtaggingMu
            )

akSoftDropZ05B152PFpatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akSoftDropZ05B152PFJets"),
        genJetMatch          = cms.InputTag("akSoftDropZ05B152PFmatch"),
        genPartonMatch       = cms.InputTag("akSoftDropZ05B152PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akSoftDropZ05B152PFcorr")),
        JetPartonMapSource   = cms.InputTag("akSoftDropZ05B152PFPatJetFlavourAssociationLegacy"),
	JetFlavourInfoSource   = cms.InputTag("akSoftDropZ05B152PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akSoftDropZ05B152PFJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = True,
        discriminatorSources = cms.VInputTag(cms.InputTag("akSoftDropZ05B152PFSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akSoftDropZ05B152PFSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akSoftDropZ05B152PFCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akSoftDropZ05B152PFCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akSoftDropZ05B152PFJetBProbabilityBJetTags"),
            cms.InputTag("akSoftDropZ05B152PFJetProbabilityBJetTags"),
            #cms.InputTag("akSoftDropZ05B152PFSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akSoftDropZ05B152PFSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akSoftDropZ05B152PFTrackCountingHighEffBJetTags"),
            cms.InputTag("akSoftDropZ05B152PFTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akSoftDropZ05B152PFJetID"),
        addBTagInfo = True,
        addTagInfos = True,
        addDiscriminators = True,
        addAssociatedTracks = True,
        addJetCharge = False,
        addJetID = False,
        getJetMCFlavour = False,
        addGenPartonMatch = False,
        addGenJetMatch = False,
        embedGenJetMatch = False,
        embedGenPartonMatch = False,
        # embedCaloTowers = False,
        # embedPFCandidates = True
        )

akSoftDropZ05B152PFNjettiness = Njettiness.clone(
		    src = cms.InputTag("akSoftDropZ05B152PFJets"),
           	    R0  = cms.double( 0.2)
)
akSoftDropZ05B152PFpatJetsWithBtagging.userData.userFloats.src += ['akSoftDropZ05B152PFNjettiness:tau1','akSoftDropZ05B152PFNjettiness:tau2','akSoftDropZ05B152PFNjettiness:tau3']

akSoftDropZ05B152PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akSoftDropZ05B152PFpatJetsWithBtagging"),
                                                             genjetTag = 'ak2HiSignalGenJets',
                                                             rParam = 0.2,
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
                                                             doLifeTimeTagging = cms.untracked.bool(True),
                                                             doLifeTimeTaggingExtras = cms.untracked.bool(False),
                                                             bTagJetName = cms.untracked.string("akSoftDropZ05B152PF"),
                                                             jetName = cms.untracked.string("akSoftDropZ05B152PF"),
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

akSoftDropZ05B152PFJetSequence_mc = cms.Sequence(
                                                  #akSoftDropZ05B152PFclean
                                                  #*
                                                  akSoftDropZ05B152PFmatch
                                                  #*
                                                  #akSoftDropZ05B152PFmatchGroomed
                                                  *
                                                  akSoftDropZ05B152PFparton
                                                  *
                                                  akSoftDropZ05B152PFcorr
                                                  *
                                                  #akSoftDropZ05B152PFJetID
                                                  #*
                                                  akSoftDropZ05B152PFPatJetFlavourIdLegacy
                                                  #*
			                          #akSoftDropZ05B152PFPatJetFlavourId  # Use legacy algo till PU implemented
                                                  *
                                                  akSoftDropZ05B152PFJetTracksAssociatorAtVertex
                                                  *
                                                  akSoftDropZ05B152PFJetBtagging
                                                  *
                                                  akSoftDropZ05B152PFNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akSoftDropZ05B152PFpatJetsWithBtagging
                                                  *
                                                  akSoftDropZ05B152PFJetAnalyzer
                                                  )

akSoftDropZ05B152PFJetSequence_data = cms.Sequence(akSoftDropZ05B152PFcorr
                                                    *
                                                    #akSoftDropZ05B152PFJetID
                                                    #*
                                                    akSoftDropZ05B152PFJetTracksAssociatorAtVertex
                                                    *
                                                    akSoftDropZ05B152PFJetBtagging
                                                    *
                                                    akSoftDropZ05B152PFNjettiness 
                                                    *
                                                    akSoftDropZ05B152PFpatJetsWithBtagging
                                                    *
                                                    akSoftDropZ05B152PFJetAnalyzer
                                                    )

akSoftDropZ05B152PFJetSequence_jec = cms.Sequence(akSoftDropZ05B152PFJetSequence_mc)
akSoftDropZ05B152PFJetSequence_mb = cms.Sequence(akSoftDropZ05B152PFJetSequence_mc)

akSoftDropZ05B152PFJetSequence = cms.Sequence(akSoftDropZ05B152PFJetSequence_data)
akSoftDropZ05B152PFpatJetsWithBtagging.userData.userFloats.src += ['akSoftDropZ05B152PFJets:sym']
akSoftDropZ05B152PFpatJetsWithBtagging.userData.userInts.src += ['akSoftDropZ05B152PFJets:droppedBranches']
