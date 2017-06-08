

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akCsSoftDropZ05B154PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akCsSoftDropZ05B154PFJets"),
    matched = cms.InputTag("ak4HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.4
    )

akCsSoftDropZ05B154PFmatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B154HiSignalGenJets"),
    matched = cms.InputTag("ak4HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.4
    )

akCsSoftDropZ05B154PFparton = patJetPartonMatch.clone(src = cms.InputTag("akCsSoftDropZ05B154PFJets")
                                                        )

akCsSoftDropZ05B154PFcorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akCsSoftDropZ05B154PFJets"),
    payload = "AK4PF_offline"
    )

akCsSoftDropZ05B154PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akCsSoftDropZ05B154CaloJets'))

#akCsSoftDropZ05B154PFclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak4HiSignalGenJets'))

akCsSoftDropZ05B154PFbTagger = bTaggers("akCsSoftDropZ05B154PF",0.4)

#create objects locally since they dont load properly otherwise
#akCsSoftDropZ05B154PFmatch = akCsSoftDropZ05B154PFbTagger.match
akCsSoftDropZ05B154PFparton = patJetPartonMatch.clone(src = cms.InputTag("akCsSoftDropZ05B154PFJets"), matched = cms.InputTag("hiSignalGenParticles"))
akCsSoftDropZ05B154PFPatJetFlavourAssociationLegacy = akCsSoftDropZ05B154PFbTagger.PatJetFlavourAssociationLegacy
akCsSoftDropZ05B154PFPatJetPartons = akCsSoftDropZ05B154PFbTagger.PatJetPartons
akCsSoftDropZ05B154PFJetTracksAssociatorAtVertex = akCsSoftDropZ05B154PFbTagger.JetTracksAssociatorAtVertex
akCsSoftDropZ05B154PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akCsSoftDropZ05B154PFSimpleSecondaryVertexHighEffBJetTags = akCsSoftDropZ05B154PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akCsSoftDropZ05B154PFSimpleSecondaryVertexHighPurBJetTags = akCsSoftDropZ05B154PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akCsSoftDropZ05B154PFCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B154PFbTagger.CombinedSecondaryVertexBJetTags
akCsSoftDropZ05B154PFCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B154PFbTagger.CombinedSecondaryVertexV2BJetTags
akCsSoftDropZ05B154PFJetBProbabilityBJetTags = akCsSoftDropZ05B154PFbTagger.JetBProbabilityBJetTags
akCsSoftDropZ05B154PFSoftPFMuonByPtBJetTags = akCsSoftDropZ05B154PFbTagger.SoftPFMuonByPtBJetTags
akCsSoftDropZ05B154PFSoftPFMuonByIP3dBJetTags = akCsSoftDropZ05B154PFbTagger.SoftPFMuonByIP3dBJetTags
akCsSoftDropZ05B154PFTrackCountingHighEffBJetTags = akCsSoftDropZ05B154PFbTagger.TrackCountingHighEffBJetTags
akCsSoftDropZ05B154PFTrackCountingHighPurBJetTags = akCsSoftDropZ05B154PFbTagger.TrackCountingHighPurBJetTags
akCsSoftDropZ05B154PFPatJetPartonAssociationLegacy = akCsSoftDropZ05B154PFbTagger.PatJetPartonAssociationLegacy

akCsSoftDropZ05B154PFImpactParameterTagInfos = akCsSoftDropZ05B154PFbTagger.ImpactParameterTagInfos
akCsSoftDropZ05B154PFImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akCsSoftDropZ05B154PFJetProbabilityBJetTags = akCsSoftDropZ05B154PFbTagger.JetProbabilityBJetTags

akCsSoftDropZ05B154PFSecondaryVertexTagInfos = akCsSoftDropZ05B154PFbTagger.SecondaryVertexTagInfos
akCsSoftDropZ05B154PFSimpleSecondaryVertexHighEffBJetTags = akCsSoftDropZ05B154PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akCsSoftDropZ05B154PFSimpleSecondaryVertexHighPurBJetTags = akCsSoftDropZ05B154PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akCsSoftDropZ05B154PFCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B154PFbTagger.CombinedSecondaryVertexBJetTags
akCsSoftDropZ05B154PFCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B154PFbTagger.CombinedSecondaryVertexV2BJetTags

akCsSoftDropZ05B154PFSecondaryVertexNegativeTagInfos = akCsSoftDropZ05B154PFbTagger.SecondaryVertexNegativeTagInfos
akCsSoftDropZ05B154PFNegativeSimpleSecondaryVertexHighEffBJetTags = akCsSoftDropZ05B154PFbTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akCsSoftDropZ05B154PFNegativeSimpleSecondaryVertexHighPurBJetTags = akCsSoftDropZ05B154PFbTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akCsSoftDropZ05B154PFNegativeCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B154PFbTagger.NegativeCombinedSecondaryVertexBJetTags
akCsSoftDropZ05B154PFPositiveCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B154PFbTagger.PositiveCombinedSecondaryVertexBJetTags
akCsSoftDropZ05B154PFNegativeCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B154PFbTagger.NegativeCombinedSecondaryVertexV2BJetTags
akCsSoftDropZ05B154PFPositiveCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B154PFbTagger.PositiveCombinedSecondaryVertexV2BJetTags

akCsSoftDropZ05B154PFSoftPFMuonsTagInfos = akCsSoftDropZ05B154PFbTagger.SoftPFMuonsTagInfos
akCsSoftDropZ05B154PFSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akCsSoftDropZ05B154PFSoftPFMuonBJetTags = akCsSoftDropZ05B154PFbTagger.SoftPFMuonBJetTags
akCsSoftDropZ05B154PFSoftPFMuonByIP3dBJetTags = akCsSoftDropZ05B154PFbTagger.SoftPFMuonByIP3dBJetTags
akCsSoftDropZ05B154PFSoftPFMuonByPtBJetTags = akCsSoftDropZ05B154PFbTagger.SoftPFMuonByPtBJetTags
akCsSoftDropZ05B154PFNegativeSoftPFMuonByPtBJetTags = akCsSoftDropZ05B154PFbTagger.NegativeSoftPFMuonByPtBJetTags
akCsSoftDropZ05B154PFPositiveSoftPFMuonByPtBJetTags = akCsSoftDropZ05B154PFbTagger.PositiveSoftPFMuonByPtBJetTags
akCsSoftDropZ05B154PFPatJetFlavourIdLegacy = cms.Sequence(akCsSoftDropZ05B154PFPatJetPartonAssociationLegacy*akCsSoftDropZ05B154PFPatJetFlavourAssociationLegacy)
#Not working with our PU sub, but keep it here for reference
#akCsSoftDropZ05B154PFPatJetFlavourAssociation = akCsSoftDropZ05B154PFbTagger.PatJetFlavourAssociation
#akCsSoftDropZ05B154PFPatJetFlavourId = cms.Sequence(akCsSoftDropZ05B154PFPatJetPartons*akCsSoftDropZ05B154PFPatJetFlavourAssociation)

akCsSoftDropZ05B154PFJetBtaggingIP       = cms.Sequence(akCsSoftDropZ05B154PFImpactParameterTagInfos *
            (akCsSoftDropZ05B154PFTrackCountingHighEffBJetTags +
             akCsSoftDropZ05B154PFTrackCountingHighPurBJetTags +
             akCsSoftDropZ05B154PFJetProbabilityBJetTags +
             akCsSoftDropZ05B154PFJetBProbabilityBJetTags 
            )
            )

akCsSoftDropZ05B154PFJetBtaggingSV = cms.Sequence(akCsSoftDropZ05B154PFImpactParameterTagInfos
            *
            akCsSoftDropZ05B154PFSecondaryVertexTagInfos
            * (akCsSoftDropZ05B154PFSimpleSecondaryVertexHighEffBJetTags+
                akCsSoftDropZ05B154PFSimpleSecondaryVertexHighPurBJetTags+
                akCsSoftDropZ05B154PFCombinedSecondaryVertexBJetTags+
                akCsSoftDropZ05B154PFCombinedSecondaryVertexV2BJetTags
              )
            )

akCsSoftDropZ05B154PFJetBtaggingNegSV = cms.Sequence(akCsSoftDropZ05B154PFImpactParameterTagInfos
            *
            akCsSoftDropZ05B154PFSecondaryVertexNegativeTagInfos
            * (akCsSoftDropZ05B154PFNegativeSimpleSecondaryVertexHighEffBJetTags+
                akCsSoftDropZ05B154PFNegativeSimpleSecondaryVertexHighPurBJetTags+
                akCsSoftDropZ05B154PFNegativeCombinedSecondaryVertexBJetTags+
                akCsSoftDropZ05B154PFPositiveCombinedSecondaryVertexBJetTags+
                akCsSoftDropZ05B154PFNegativeCombinedSecondaryVertexV2BJetTags+
                akCsSoftDropZ05B154PFPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akCsSoftDropZ05B154PFJetBtaggingMu = cms.Sequence(akCsSoftDropZ05B154PFSoftPFMuonsTagInfos * (akCsSoftDropZ05B154PFSoftPFMuonBJetTags
                +
                akCsSoftDropZ05B154PFSoftPFMuonByIP3dBJetTags
                +
                akCsSoftDropZ05B154PFSoftPFMuonByPtBJetTags
                +
                akCsSoftDropZ05B154PFNegativeSoftPFMuonByPtBJetTags
                +
                akCsSoftDropZ05B154PFPositiveSoftPFMuonByPtBJetTags
              )
            )

akCsSoftDropZ05B154PFJetBtagging = cms.Sequence(akCsSoftDropZ05B154PFJetBtaggingIP
            *akCsSoftDropZ05B154PFJetBtaggingSV
            *akCsSoftDropZ05B154PFJetBtaggingNegSV
#            *akCsSoftDropZ05B154PFJetBtaggingMu
            )

akCsSoftDropZ05B154PFpatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akCsSoftDropZ05B154PFJets"),
        genJetMatch          = cms.InputTag("akCsSoftDropZ05B154PFmatch"),
        genPartonMatch       = cms.InputTag("akCsSoftDropZ05B154PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akCsSoftDropZ05B154PFcorr")),
        JetPartonMapSource   = cms.InputTag("akCsSoftDropZ05B154PFPatJetFlavourAssociationLegacy"),
	JetFlavourInfoSource   = cms.InputTag("akCsSoftDropZ05B154PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akCsSoftDropZ05B154PFJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = True,
        discriminatorSources = cms.VInputTag(cms.InputTag("akCsSoftDropZ05B154PFSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akCsSoftDropZ05B154PFSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akCsSoftDropZ05B154PFCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akCsSoftDropZ05B154PFCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akCsSoftDropZ05B154PFJetBProbabilityBJetTags"),
            cms.InputTag("akCsSoftDropZ05B154PFJetProbabilityBJetTags"),
            #cms.InputTag("akCsSoftDropZ05B154PFSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akCsSoftDropZ05B154PFSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akCsSoftDropZ05B154PFTrackCountingHighEffBJetTags"),
            cms.InputTag("akCsSoftDropZ05B154PFTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akCsSoftDropZ05B154PFJetID"),
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

akCsSoftDropZ05B154PFNjettiness = Njettiness.clone(
		    src = cms.InputTag("akCsSoftDropZ05B154PFJets"),
           	    R0  = cms.double( 0.4)
)
akCsSoftDropZ05B154PFpatJetsWithBtagging.userData.userFloats.src += ['akCsSoftDropZ05B154PFNjettiness:tau1','akCsSoftDropZ05B154PFNjettiness:tau2','akCsSoftDropZ05B154PFNjettiness:tau3']

akCsSoftDropZ05B154PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akCsSoftDropZ05B154PFpatJetsWithBtagging"),
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
                                                             doLifeTimeTagging = cms.untracked.bool(True),
                                                             doLifeTimeTaggingExtras = cms.untracked.bool(False),
                                                             bTagJetName = cms.untracked.string("akCsSoftDropZ05B154PF"),
                                                             jetName = cms.untracked.string("akCsSoftDropZ05B154PF"),
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

akCsSoftDropZ05B154PFJetSequence_mc = cms.Sequence(
                                                  #akCsSoftDropZ05B154PFclean
                                                  #*
                                                  akCsSoftDropZ05B154PFmatch
                                                  #*
                                                  #akCsSoftDropZ05B154PFmatchGroomed
                                                  *
                                                  akCsSoftDropZ05B154PFparton
                                                  *
                                                  akCsSoftDropZ05B154PFcorr
                                                  *
                                                  #akCsSoftDropZ05B154PFJetID
                                                  #*
                                                  akCsSoftDropZ05B154PFPatJetFlavourIdLegacy
                                                  #*
			                          #akCsSoftDropZ05B154PFPatJetFlavourId  # Use legacy algo till PU implemented
                                                  *
                                                  akCsSoftDropZ05B154PFJetTracksAssociatorAtVertex
                                                  *
                                                  akCsSoftDropZ05B154PFJetBtagging
                                                  *
                                                  akCsSoftDropZ05B154PFNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akCsSoftDropZ05B154PFpatJetsWithBtagging
                                                  *
                                                  akCsSoftDropZ05B154PFJetAnalyzer
                                                  )

akCsSoftDropZ05B154PFJetSequence_data = cms.Sequence(akCsSoftDropZ05B154PFcorr
                                                    *
                                                    #akCsSoftDropZ05B154PFJetID
                                                    #*
                                                    akCsSoftDropZ05B154PFJetTracksAssociatorAtVertex
                                                    *
                                                    akCsSoftDropZ05B154PFJetBtagging
                                                    *
                                                    akCsSoftDropZ05B154PFNjettiness 
                                                    *
                                                    akCsSoftDropZ05B154PFpatJetsWithBtagging
                                                    *
                                                    akCsSoftDropZ05B154PFJetAnalyzer
                                                    )

akCsSoftDropZ05B154PFJetSequence_jec = cms.Sequence(akCsSoftDropZ05B154PFJetSequence_mc)
akCsSoftDropZ05B154PFJetSequence_mb = cms.Sequence(akCsSoftDropZ05B154PFJetSequence_mc)

akCsSoftDropZ05B154PFJetSequence = cms.Sequence(akCsSoftDropZ05B154PFJetSequence_mc)
akCsSoftDropZ05B154PFpatJetsWithBtagging.userData.userFloats.src += ['akCsSoftDropZ05B154PFJets:sym']
akCsSoftDropZ05B154PFpatJetsWithBtagging.userData.userInts.src += ['akCsSoftDropZ05B154PFJets:droppedBranches']
