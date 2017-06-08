

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akPuSoftDropZ05B154PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akPuSoftDropZ05B154PFJets"),
    matched = cms.InputTag("ak4HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
    )

akPuSoftDropZ05B154PFmatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B154HiSignalGenJets"),
    matched = cms.InputTag("ak4HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.4
    )

akPuSoftDropZ05B154PFparton = patJetPartonMatch.clone(src = cms.InputTag("akPuSoftDropZ05B154PFJets")
                                                        )

akPuSoftDropZ05B154PFcorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akPuSoftDropZ05B154PFJets"),
    payload = "AKPu4PF_offline"
    )

akPuSoftDropZ05B154PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akPuSoftDropZ05B154CaloJets'))

#akPuSoftDropZ05B154PFclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak4HiSignalGenJets'))

akPuSoftDropZ05B154PFbTagger = bTaggers("akPuSoftDropZ05B154PF",0.4)

#create objects locally since they dont load properly otherwise
#akPuSoftDropZ05B154PFmatch = akPuSoftDropZ05B154PFbTagger.match
akPuSoftDropZ05B154PFparton = patJetPartonMatch.clone(src = cms.InputTag("akPuSoftDropZ05B154PFJets"), matched = cms.InputTag("hiSignalGenParticles"))
akPuSoftDropZ05B154PFPatJetFlavourAssociationLegacy = akPuSoftDropZ05B154PFbTagger.PatJetFlavourAssociationLegacy
akPuSoftDropZ05B154PFPatJetPartons = akPuSoftDropZ05B154PFbTagger.PatJetPartons
akPuSoftDropZ05B154PFJetTracksAssociatorAtVertex = akPuSoftDropZ05B154PFbTagger.JetTracksAssociatorAtVertex
akPuSoftDropZ05B154PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akPuSoftDropZ05B154PFSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B154PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akPuSoftDropZ05B154PFSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B154PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akPuSoftDropZ05B154PFCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B154PFbTagger.CombinedSecondaryVertexBJetTags
akPuSoftDropZ05B154PFCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B154PFbTagger.CombinedSecondaryVertexV2BJetTags
akPuSoftDropZ05B154PFJetBProbabilityBJetTags = akPuSoftDropZ05B154PFbTagger.JetBProbabilityBJetTags
akPuSoftDropZ05B154PFSoftPFMuonByPtBJetTags = akPuSoftDropZ05B154PFbTagger.SoftPFMuonByPtBJetTags
akPuSoftDropZ05B154PFSoftPFMuonByIP3dBJetTags = akPuSoftDropZ05B154PFbTagger.SoftPFMuonByIP3dBJetTags
akPuSoftDropZ05B154PFTrackCountingHighEffBJetTags = akPuSoftDropZ05B154PFbTagger.TrackCountingHighEffBJetTags
akPuSoftDropZ05B154PFTrackCountingHighPurBJetTags = akPuSoftDropZ05B154PFbTagger.TrackCountingHighPurBJetTags
akPuSoftDropZ05B154PFPatJetPartonAssociationLegacy = akPuSoftDropZ05B154PFbTagger.PatJetPartonAssociationLegacy

akPuSoftDropZ05B154PFImpactParameterTagInfos = akPuSoftDropZ05B154PFbTagger.ImpactParameterTagInfos
akPuSoftDropZ05B154PFImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akPuSoftDropZ05B154PFJetProbabilityBJetTags = akPuSoftDropZ05B154PFbTagger.JetProbabilityBJetTags

akPuSoftDropZ05B154PFSecondaryVertexTagInfos = akPuSoftDropZ05B154PFbTagger.SecondaryVertexTagInfos
akPuSoftDropZ05B154PFSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B154PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akPuSoftDropZ05B154PFSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B154PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akPuSoftDropZ05B154PFCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B154PFbTagger.CombinedSecondaryVertexBJetTags
akPuSoftDropZ05B154PFCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B154PFbTagger.CombinedSecondaryVertexV2BJetTags

akPuSoftDropZ05B154PFSecondaryVertexNegativeTagInfos = akPuSoftDropZ05B154PFbTagger.SecondaryVertexNegativeTagInfos
akPuSoftDropZ05B154PFNegativeSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B154PFbTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akPuSoftDropZ05B154PFNegativeSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B154PFbTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akPuSoftDropZ05B154PFNegativeCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B154PFbTagger.NegativeCombinedSecondaryVertexBJetTags
akPuSoftDropZ05B154PFPositiveCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B154PFbTagger.PositiveCombinedSecondaryVertexBJetTags
akPuSoftDropZ05B154PFNegativeCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B154PFbTagger.NegativeCombinedSecondaryVertexV2BJetTags
akPuSoftDropZ05B154PFPositiveCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B154PFbTagger.PositiveCombinedSecondaryVertexV2BJetTags

akPuSoftDropZ05B154PFSoftPFMuonsTagInfos = akPuSoftDropZ05B154PFbTagger.SoftPFMuonsTagInfos
akPuSoftDropZ05B154PFSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akPuSoftDropZ05B154PFSoftPFMuonBJetTags = akPuSoftDropZ05B154PFbTagger.SoftPFMuonBJetTags
akPuSoftDropZ05B154PFSoftPFMuonByIP3dBJetTags = akPuSoftDropZ05B154PFbTagger.SoftPFMuonByIP3dBJetTags
akPuSoftDropZ05B154PFSoftPFMuonByPtBJetTags = akPuSoftDropZ05B154PFbTagger.SoftPFMuonByPtBJetTags
akPuSoftDropZ05B154PFNegativeSoftPFMuonByPtBJetTags = akPuSoftDropZ05B154PFbTagger.NegativeSoftPFMuonByPtBJetTags
akPuSoftDropZ05B154PFPositiveSoftPFMuonByPtBJetTags = akPuSoftDropZ05B154PFbTagger.PositiveSoftPFMuonByPtBJetTags
akPuSoftDropZ05B154PFPatJetFlavourIdLegacy = cms.Sequence(akPuSoftDropZ05B154PFPatJetPartonAssociationLegacy*akPuSoftDropZ05B154PFPatJetFlavourAssociationLegacy)
#Not working with our PU sub, but keep it here for reference
#akPuSoftDropZ05B154PFPatJetFlavourAssociation = akPuSoftDropZ05B154PFbTagger.PatJetFlavourAssociation
#akPuSoftDropZ05B154PFPatJetFlavourId = cms.Sequence(akPuSoftDropZ05B154PFPatJetPartons*akPuSoftDropZ05B154PFPatJetFlavourAssociation)

akPuSoftDropZ05B154PFJetBtaggingIP       = cms.Sequence(akPuSoftDropZ05B154PFImpactParameterTagInfos *
            (akPuSoftDropZ05B154PFTrackCountingHighEffBJetTags +
             akPuSoftDropZ05B154PFTrackCountingHighPurBJetTags +
             akPuSoftDropZ05B154PFJetProbabilityBJetTags +
             akPuSoftDropZ05B154PFJetBProbabilityBJetTags 
            )
            )

akPuSoftDropZ05B154PFJetBtaggingSV = cms.Sequence(akPuSoftDropZ05B154PFImpactParameterTagInfos
            *
            akPuSoftDropZ05B154PFSecondaryVertexTagInfos
            * (akPuSoftDropZ05B154PFSimpleSecondaryVertexHighEffBJetTags+
                akPuSoftDropZ05B154PFSimpleSecondaryVertexHighPurBJetTags+
                akPuSoftDropZ05B154PFCombinedSecondaryVertexBJetTags+
                akPuSoftDropZ05B154PFCombinedSecondaryVertexV2BJetTags
              )
            )

akPuSoftDropZ05B154PFJetBtaggingNegSV = cms.Sequence(akPuSoftDropZ05B154PFImpactParameterTagInfos
            *
            akPuSoftDropZ05B154PFSecondaryVertexNegativeTagInfos
            * (akPuSoftDropZ05B154PFNegativeSimpleSecondaryVertexHighEffBJetTags+
                akPuSoftDropZ05B154PFNegativeSimpleSecondaryVertexHighPurBJetTags+
                akPuSoftDropZ05B154PFNegativeCombinedSecondaryVertexBJetTags+
                akPuSoftDropZ05B154PFPositiveCombinedSecondaryVertexBJetTags+
                akPuSoftDropZ05B154PFNegativeCombinedSecondaryVertexV2BJetTags+
                akPuSoftDropZ05B154PFPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akPuSoftDropZ05B154PFJetBtaggingMu = cms.Sequence(akPuSoftDropZ05B154PFSoftPFMuonsTagInfos * (akPuSoftDropZ05B154PFSoftPFMuonBJetTags
                +
                akPuSoftDropZ05B154PFSoftPFMuonByIP3dBJetTags
                +
                akPuSoftDropZ05B154PFSoftPFMuonByPtBJetTags
                +
                akPuSoftDropZ05B154PFNegativeSoftPFMuonByPtBJetTags
                +
                akPuSoftDropZ05B154PFPositiveSoftPFMuonByPtBJetTags
              )
            )

akPuSoftDropZ05B154PFJetBtagging = cms.Sequence(akPuSoftDropZ05B154PFJetBtaggingIP
            *akPuSoftDropZ05B154PFJetBtaggingSV
            *akPuSoftDropZ05B154PFJetBtaggingNegSV
#            *akPuSoftDropZ05B154PFJetBtaggingMu
            )

akPuSoftDropZ05B154PFpatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akPuSoftDropZ05B154PFJets"),
        genJetMatch          = cms.InputTag("akPuSoftDropZ05B154PFmatch"),
        genPartonMatch       = cms.InputTag("akPuSoftDropZ05B154PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPuSoftDropZ05B154PFcorr")),
        JetPartonMapSource   = cms.InputTag("akPuSoftDropZ05B154PFPatJetFlavourAssociationLegacy"),
	JetFlavourInfoSource   = cms.InputTag("akPuSoftDropZ05B154PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akPuSoftDropZ05B154PFJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = True,
        discriminatorSources = cms.VInputTag(cms.InputTag("akPuSoftDropZ05B154PFSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akPuSoftDropZ05B154PFSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akPuSoftDropZ05B154PFCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akPuSoftDropZ05B154PFCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akPuSoftDropZ05B154PFJetBProbabilityBJetTags"),
            cms.InputTag("akPuSoftDropZ05B154PFJetProbabilityBJetTags"),
            #cms.InputTag("akPuSoftDropZ05B154PFSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akPuSoftDropZ05B154PFSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akPuSoftDropZ05B154PFTrackCountingHighEffBJetTags"),
            cms.InputTag("akPuSoftDropZ05B154PFTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akPuSoftDropZ05B154PFJetID"),
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

akPuSoftDropZ05B154PFNjettiness = Njettiness.clone(
		    src = cms.InputTag("akPuSoftDropZ05B154PFJets"),
           	    R0  = cms.double( 0.4)
)
akPuSoftDropZ05B154PFpatJetsWithBtagging.userData.userFloats.src += ['akPuSoftDropZ05B154PFNjettiness:tau1','akPuSoftDropZ05B154PFNjettiness:tau2','akPuSoftDropZ05B154PFNjettiness:tau3']

akPuSoftDropZ05B154PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akPuSoftDropZ05B154PFpatJetsWithBtagging"),
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
                                                             bTagJetName = cms.untracked.string("akPuSoftDropZ05B154PF"),
                                                             jetName = cms.untracked.string("akPuSoftDropZ05B154PF"),
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

akPuSoftDropZ05B154PFJetSequence_mc = cms.Sequence(
                                                  #akPuSoftDropZ05B154PFclean
                                                  #*
                                                  akPuSoftDropZ05B154PFmatch
                                                  #*
                                                  #akPuSoftDropZ05B154PFmatchGroomed
                                                  *
                                                  akPuSoftDropZ05B154PFparton
                                                  *
                                                  akPuSoftDropZ05B154PFcorr
                                                  *
                                                  #akPuSoftDropZ05B154PFJetID
                                                  #*
                                                  akPuSoftDropZ05B154PFPatJetFlavourIdLegacy
                                                  #*
			                          #akPuSoftDropZ05B154PFPatJetFlavourId  # Use legacy algo till PU implemented
                                                  *
                                                  akPuSoftDropZ05B154PFJetTracksAssociatorAtVertex
                                                  *
                                                  akPuSoftDropZ05B154PFJetBtagging
                                                  *
                                                  akPuSoftDropZ05B154PFNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akPuSoftDropZ05B154PFpatJetsWithBtagging
                                                  *
                                                  akPuSoftDropZ05B154PFJetAnalyzer
                                                  )

akPuSoftDropZ05B154PFJetSequence_data = cms.Sequence(akPuSoftDropZ05B154PFcorr
                                                    *
                                                    #akPuSoftDropZ05B154PFJetID
                                                    #*
                                                    akPuSoftDropZ05B154PFJetTracksAssociatorAtVertex
                                                    *
                                                    akPuSoftDropZ05B154PFJetBtagging
                                                    *
                                                    akPuSoftDropZ05B154PFNjettiness 
                                                    *
                                                    akPuSoftDropZ05B154PFpatJetsWithBtagging
                                                    *
                                                    akPuSoftDropZ05B154PFJetAnalyzer
                                                    )

akPuSoftDropZ05B154PFJetSequence_jec = cms.Sequence(akPuSoftDropZ05B154PFJetSequence_mc)
akPuSoftDropZ05B154PFJetSequence_mb = cms.Sequence(akPuSoftDropZ05B154PFJetSequence_mc)

akPuSoftDropZ05B154PFJetSequence = cms.Sequence(akPuSoftDropZ05B154PFJetSequence_jec)
akPuSoftDropZ05B154PFJetAnalyzer.genPtMin = cms.untracked.double(1)
akPuSoftDropZ05B154PFJetAnalyzer.jetPtMin = cms.double(1)
akPuSoftDropZ05B154PFpatJetsWithBtagging.userData.userFloats.src += ['akPuSoftDropZ05B154PFJets:sym']
akPuSoftDropZ05B154PFpatJetsWithBtagging.userData.userInts.src += ['akPuSoftDropZ05B154PFJets:droppedBranches']
