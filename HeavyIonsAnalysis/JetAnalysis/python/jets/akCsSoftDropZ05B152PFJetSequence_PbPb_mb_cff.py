

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akCsSoftDropZ05B152PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akCsSoftDropZ05B152PFJets"),
    matched = cms.InputTag("ak2HiCleanedGenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.2
    )

akCsSoftDropZ05B152PFmatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B152HiGenJets"),
    matched = cms.InputTag("ak2HiCleanedGenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.2
    )

akCsSoftDropZ05B152PFparton = patJetPartonMatch.clone(src = cms.InputTag("akCsSoftDropZ05B152PFJets")
                                                        )

akCsSoftDropZ05B152PFcorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akCsSoftDropZ05B152PFJets"),
    payload = "AK2PF_offline"
    )

akCsSoftDropZ05B152PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akCsSoftDropZ05B152CaloJets'))

#akCsSoftDropZ05B152PFclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak2HiCleanedGenJets'))

akCsSoftDropZ05B152PFbTagger = bTaggers("akCsSoftDropZ05B152PF",0.2)

#create objects locally since they dont load properly otherwise
#akCsSoftDropZ05B152PFmatch = akCsSoftDropZ05B152PFbTagger.match
akCsSoftDropZ05B152PFparton = patJetPartonMatch.clone(src = cms.InputTag("akCsSoftDropZ05B152PFJets"), matched = cms.InputTag("selectedPartons"))
akCsSoftDropZ05B152PFPatJetFlavourAssociationLegacy = akCsSoftDropZ05B152PFbTagger.PatJetFlavourAssociationLegacy
akCsSoftDropZ05B152PFPatJetPartons = akCsSoftDropZ05B152PFbTagger.PatJetPartons
akCsSoftDropZ05B152PFJetTracksAssociatorAtVertex = akCsSoftDropZ05B152PFbTagger.JetTracksAssociatorAtVertex
akCsSoftDropZ05B152PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akCsSoftDropZ05B152PFSimpleSecondaryVertexHighEffBJetTags = akCsSoftDropZ05B152PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akCsSoftDropZ05B152PFSimpleSecondaryVertexHighPurBJetTags = akCsSoftDropZ05B152PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akCsSoftDropZ05B152PFCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B152PFbTagger.CombinedSecondaryVertexBJetTags
akCsSoftDropZ05B152PFCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B152PFbTagger.CombinedSecondaryVertexV2BJetTags
akCsSoftDropZ05B152PFJetBProbabilityBJetTags = akCsSoftDropZ05B152PFbTagger.JetBProbabilityBJetTags
akCsSoftDropZ05B152PFSoftPFMuonByPtBJetTags = akCsSoftDropZ05B152PFbTagger.SoftPFMuonByPtBJetTags
akCsSoftDropZ05B152PFSoftPFMuonByIP3dBJetTags = akCsSoftDropZ05B152PFbTagger.SoftPFMuonByIP3dBJetTags
akCsSoftDropZ05B152PFTrackCountingHighEffBJetTags = akCsSoftDropZ05B152PFbTagger.TrackCountingHighEffBJetTags
akCsSoftDropZ05B152PFTrackCountingHighPurBJetTags = akCsSoftDropZ05B152PFbTagger.TrackCountingHighPurBJetTags
akCsSoftDropZ05B152PFPatJetPartonAssociationLegacy = akCsSoftDropZ05B152PFbTagger.PatJetPartonAssociationLegacy

akCsSoftDropZ05B152PFImpactParameterTagInfos = akCsSoftDropZ05B152PFbTagger.ImpactParameterTagInfos
akCsSoftDropZ05B152PFImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akCsSoftDropZ05B152PFJetProbabilityBJetTags = akCsSoftDropZ05B152PFbTagger.JetProbabilityBJetTags

akCsSoftDropZ05B152PFSecondaryVertexTagInfos = akCsSoftDropZ05B152PFbTagger.SecondaryVertexTagInfos
akCsSoftDropZ05B152PFSimpleSecondaryVertexHighEffBJetTags = akCsSoftDropZ05B152PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akCsSoftDropZ05B152PFSimpleSecondaryVertexHighPurBJetTags = akCsSoftDropZ05B152PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akCsSoftDropZ05B152PFCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B152PFbTagger.CombinedSecondaryVertexBJetTags
akCsSoftDropZ05B152PFCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B152PFbTagger.CombinedSecondaryVertexV2BJetTags

akCsSoftDropZ05B152PFSecondaryVertexNegativeTagInfos = akCsSoftDropZ05B152PFbTagger.SecondaryVertexNegativeTagInfos
akCsSoftDropZ05B152PFNegativeSimpleSecondaryVertexHighEffBJetTags = akCsSoftDropZ05B152PFbTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akCsSoftDropZ05B152PFNegativeSimpleSecondaryVertexHighPurBJetTags = akCsSoftDropZ05B152PFbTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akCsSoftDropZ05B152PFNegativeCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B152PFbTagger.NegativeCombinedSecondaryVertexBJetTags
akCsSoftDropZ05B152PFPositiveCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B152PFbTagger.PositiveCombinedSecondaryVertexBJetTags
akCsSoftDropZ05B152PFNegativeCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B152PFbTagger.NegativeCombinedSecondaryVertexV2BJetTags
akCsSoftDropZ05B152PFPositiveCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B152PFbTagger.PositiveCombinedSecondaryVertexV2BJetTags

akCsSoftDropZ05B152PFSoftPFMuonsTagInfos = akCsSoftDropZ05B152PFbTagger.SoftPFMuonsTagInfos
akCsSoftDropZ05B152PFSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akCsSoftDropZ05B152PFSoftPFMuonBJetTags = akCsSoftDropZ05B152PFbTagger.SoftPFMuonBJetTags
akCsSoftDropZ05B152PFSoftPFMuonByIP3dBJetTags = akCsSoftDropZ05B152PFbTagger.SoftPFMuonByIP3dBJetTags
akCsSoftDropZ05B152PFSoftPFMuonByPtBJetTags = akCsSoftDropZ05B152PFbTagger.SoftPFMuonByPtBJetTags
akCsSoftDropZ05B152PFNegativeSoftPFMuonByPtBJetTags = akCsSoftDropZ05B152PFbTagger.NegativeSoftPFMuonByPtBJetTags
akCsSoftDropZ05B152PFPositiveSoftPFMuonByPtBJetTags = akCsSoftDropZ05B152PFbTagger.PositiveSoftPFMuonByPtBJetTags
akCsSoftDropZ05B152PFPatJetFlavourIdLegacy = cms.Sequence(akCsSoftDropZ05B152PFPatJetPartonAssociationLegacy*akCsSoftDropZ05B152PFPatJetFlavourAssociationLegacy)
#Not working with our PU sub, but keep it here for reference
#akCsSoftDropZ05B152PFPatJetFlavourAssociation = akCsSoftDropZ05B152PFbTagger.PatJetFlavourAssociation
#akCsSoftDropZ05B152PFPatJetFlavourId = cms.Sequence(akCsSoftDropZ05B152PFPatJetPartons*akCsSoftDropZ05B152PFPatJetFlavourAssociation)

akCsSoftDropZ05B152PFJetBtaggingIP       = cms.Sequence(akCsSoftDropZ05B152PFImpactParameterTagInfos *
            (akCsSoftDropZ05B152PFTrackCountingHighEffBJetTags +
             akCsSoftDropZ05B152PFTrackCountingHighPurBJetTags +
             akCsSoftDropZ05B152PFJetProbabilityBJetTags +
             akCsSoftDropZ05B152PFJetBProbabilityBJetTags 
            )
            )

akCsSoftDropZ05B152PFJetBtaggingSV = cms.Sequence(akCsSoftDropZ05B152PFImpactParameterTagInfos
            *
            akCsSoftDropZ05B152PFSecondaryVertexTagInfos
            * (akCsSoftDropZ05B152PFSimpleSecondaryVertexHighEffBJetTags+
                akCsSoftDropZ05B152PFSimpleSecondaryVertexHighPurBJetTags+
                akCsSoftDropZ05B152PFCombinedSecondaryVertexBJetTags+
                akCsSoftDropZ05B152PFCombinedSecondaryVertexV2BJetTags
              )
            )

akCsSoftDropZ05B152PFJetBtaggingNegSV = cms.Sequence(akCsSoftDropZ05B152PFImpactParameterTagInfos
            *
            akCsSoftDropZ05B152PFSecondaryVertexNegativeTagInfos
            * (akCsSoftDropZ05B152PFNegativeSimpleSecondaryVertexHighEffBJetTags+
                akCsSoftDropZ05B152PFNegativeSimpleSecondaryVertexHighPurBJetTags+
                akCsSoftDropZ05B152PFNegativeCombinedSecondaryVertexBJetTags+
                akCsSoftDropZ05B152PFPositiveCombinedSecondaryVertexBJetTags+
                akCsSoftDropZ05B152PFNegativeCombinedSecondaryVertexV2BJetTags+
                akCsSoftDropZ05B152PFPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akCsSoftDropZ05B152PFJetBtaggingMu = cms.Sequence(akCsSoftDropZ05B152PFSoftPFMuonsTagInfos * (akCsSoftDropZ05B152PFSoftPFMuonBJetTags
                +
                akCsSoftDropZ05B152PFSoftPFMuonByIP3dBJetTags
                +
                akCsSoftDropZ05B152PFSoftPFMuonByPtBJetTags
                +
                akCsSoftDropZ05B152PFNegativeSoftPFMuonByPtBJetTags
                +
                akCsSoftDropZ05B152PFPositiveSoftPFMuonByPtBJetTags
              )
            )

akCsSoftDropZ05B152PFJetBtagging = cms.Sequence(akCsSoftDropZ05B152PFJetBtaggingIP
            *akCsSoftDropZ05B152PFJetBtaggingSV
            *akCsSoftDropZ05B152PFJetBtaggingNegSV
#            *akCsSoftDropZ05B152PFJetBtaggingMu
            )

akCsSoftDropZ05B152PFpatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akCsSoftDropZ05B152PFJets"),
        genJetMatch          = cms.InputTag("akCsSoftDropZ05B152PFmatch"),
        genPartonMatch       = cms.InputTag("akCsSoftDropZ05B152PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akCsSoftDropZ05B152PFcorr")),
        JetPartonMapSource   = cms.InputTag("akCsSoftDropZ05B152PFPatJetFlavourAssociationLegacy"),
	JetFlavourInfoSource   = cms.InputTag("akCsSoftDropZ05B152PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akCsSoftDropZ05B152PFJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = True,
        discriminatorSources = cms.VInputTag(cms.InputTag("akCsSoftDropZ05B152PFSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akCsSoftDropZ05B152PFSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akCsSoftDropZ05B152PFCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akCsSoftDropZ05B152PFCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akCsSoftDropZ05B152PFJetBProbabilityBJetTags"),
            cms.InputTag("akCsSoftDropZ05B152PFJetProbabilityBJetTags"),
            #cms.InputTag("akCsSoftDropZ05B152PFSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akCsSoftDropZ05B152PFSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akCsSoftDropZ05B152PFTrackCountingHighEffBJetTags"),
            cms.InputTag("akCsSoftDropZ05B152PFTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akCsSoftDropZ05B152PFJetID"),
        addBTagInfo = True,
        addTagInfos = True,
        addDiscriminators = True,
        addAssociatedTracks = True,
        addJetCharge = False,
        addJetID = False,
        getJetMCFlavour = True,
        addGenPartonMatch = True,
        addGenJetMatch = True,
        embedGenJetMatch = True,
        embedGenPartonMatch = True,
        # embedCaloTowers = False,
        # embedPFCandidates = True
        )

akCsSoftDropZ05B152PFNjettiness = Njettiness.clone(
		    src = cms.InputTag("akCsSoftDropZ05B152PFJets"),
           	    R0  = cms.double( 0.2)
)
akCsSoftDropZ05B152PFpatJetsWithBtagging.userData.userFloats.src += ['akCsSoftDropZ05B152PFNjettiness:tau1','akCsSoftDropZ05B152PFNjettiness:tau2','akCsSoftDropZ05B152PFNjettiness:tau3']

akCsSoftDropZ05B152PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akCsSoftDropZ05B152PFpatJetsWithBtagging"),
                                                             genjetTag = 'ak2HiGenJets',
                                                             rParam = 0.2,
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
                                                             doLifeTimeTagging = cms.untracked.bool(True),
                                                             doLifeTimeTaggingExtras = cms.untracked.bool(False),
                                                             bTagJetName = cms.untracked.string("akCsSoftDropZ05B152PF"),
                                                             jetName = cms.untracked.string("akCsSoftDropZ05B152PF"),
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

akCsSoftDropZ05B152PFJetSequence_mc = cms.Sequence(
                                                  #akCsSoftDropZ05B152PFclean
                                                  #*
                                                  akCsSoftDropZ05B152PFmatch
                                                  #*
                                                  #akCsSoftDropZ05B152PFmatchGroomed
                                                  *
                                                  akCsSoftDropZ05B152PFparton
                                                  *
                                                  akCsSoftDropZ05B152PFcorr
                                                  *
                                                  #akCsSoftDropZ05B152PFJetID
                                                  #*
                                                  akCsSoftDropZ05B152PFPatJetFlavourIdLegacy
                                                  #*
			                          #akCsSoftDropZ05B152PFPatJetFlavourId  # Use legacy algo till PU implemented
                                                  *
                                                  akCsSoftDropZ05B152PFJetTracksAssociatorAtVertex
                                                  *
                                                  akCsSoftDropZ05B152PFJetBtagging
                                                  *
                                                  akCsSoftDropZ05B152PFNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akCsSoftDropZ05B152PFpatJetsWithBtagging
                                                  *
                                                  akCsSoftDropZ05B152PFJetAnalyzer
                                                  )

akCsSoftDropZ05B152PFJetSequence_data = cms.Sequence(akCsSoftDropZ05B152PFcorr
                                                    *
                                                    #akCsSoftDropZ05B152PFJetID
                                                    #*
                                                    akCsSoftDropZ05B152PFJetTracksAssociatorAtVertex
                                                    *
                                                    akCsSoftDropZ05B152PFJetBtagging
                                                    *
                                                    akCsSoftDropZ05B152PFNjettiness 
                                                    *
                                                    akCsSoftDropZ05B152PFpatJetsWithBtagging
                                                    *
                                                    akCsSoftDropZ05B152PFJetAnalyzer
                                                    )

akCsSoftDropZ05B152PFJetSequence_jec = cms.Sequence(akCsSoftDropZ05B152PFJetSequence_mc)
akCsSoftDropZ05B152PFJetSequence_mb = cms.Sequence(akCsSoftDropZ05B152PFJetSequence_mc)

akCsSoftDropZ05B152PFJetSequence = cms.Sequence(akCsSoftDropZ05B152PFJetSequence_mb)
akCsSoftDropZ05B152PFpatJetsWithBtagging.userData.userFloats.src += ['akCsSoftDropZ05B152PFJets:sym']
akCsSoftDropZ05B152PFpatJetsWithBtagging.userData.userInts.src += ['akCsSoftDropZ05B152PFJets:droppedBranches']
