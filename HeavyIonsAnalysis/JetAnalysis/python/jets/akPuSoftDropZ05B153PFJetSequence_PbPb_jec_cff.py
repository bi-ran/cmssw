

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akPuSoftDropZ05B153PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akPuSoftDropZ05B153PFJets"),
    matched = cms.InputTag("ak3HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.3
    )

akPuSoftDropZ05B153PFmatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B153HiSignalGenJets"),
    matched = cms.InputTag("ak3HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.3
    )

akPuSoftDropZ05B153PFparton = patJetPartonMatch.clone(src = cms.InputTag("akPuSoftDropZ05B153PFJets")
                                                        )

akPuSoftDropZ05B153PFcorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akPuSoftDropZ05B153PFJets"),
    payload = "AKPu3PF_offline"
    )

akPuSoftDropZ05B153PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akPuSoftDropZ05B153CaloJets'))

#akPuSoftDropZ05B153PFclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak3HiSignalGenJets'))

akPuSoftDropZ05B153PFbTagger = bTaggers("akPuSoftDropZ05B153PF",0.3)

#create objects locally since they dont load properly otherwise
#akPuSoftDropZ05B153PFmatch = akPuSoftDropZ05B153PFbTagger.match
akPuSoftDropZ05B153PFparton = patJetPartonMatch.clone(src = cms.InputTag("akPuSoftDropZ05B153PFJets"), matched = cms.InputTag("hiSignalGenParticles"))
akPuSoftDropZ05B153PFPatJetFlavourAssociationLegacy = akPuSoftDropZ05B153PFbTagger.PatJetFlavourAssociationLegacy
akPuSoftDropZ05B153PFPatJetPartons = akPuSoftDropZ05B153PFbTagger.PatJetPartons
akPuSoftDropZ05B153PFJetTracksAssociatorAtVertex = akPuSoftDropZ05B153PFbTagger.JetTracksAssociatorAtVertex
akPuSoftDropZ05B153PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akPuSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akPuSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akPuSoftDropZ05B153PFCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B153PFbTagger.CombinedSecondaryVertexBJetTags
akPuSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B153PFbTagger.CombinedSecondaryVertexV2BJetTags
akPuSoftDropZ05B153PFJetBProbabilityBJetTags = akPuSoftDropZ05B153PFbTagger.JetBProbabilityBJetTags
akPuSoftDropZ05B153PFSoftPFMuonByPtBJetTags = akPuSoftDropZ05B153PFbTagger.SoftPFMuonByPtBJetTags
akPuSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags = akPuSoftDropZ05B153PFbTagger.SoftPFMuonByIP3dBJetTags
akPuSoftDropZ05B153PFTrackCountingHighEffBJetTags = akPuSoftDropZ05B153PFbTagger.TrackCountingHighEffBJetTags
akPuSoftDropZ05B153PFTrackCountingHighPurBJetTags = akPuSoftDropZ05B153PFbTagger.TrackCountingHighPurBJetTags
akPuSoftDropZ05B153PFPatJetPartonAssociationLegacy = akPuSoftDropZ05B153PFbTagger.PatJetPartonAssociationLegacy

akPuSoftDropZ05B153PFImpactParameterTagInfos = akPuSoftDropZ05B153PFbTagger.ImpactParameterTagInfos
akPuSoftDropZ05B153PFImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akPuSoftDropZ05B153PFJetProbabilityBJetTags = akPuSoftDropZ05B153PFbTagger.JetProbabilityBJetTags

akPuSoftDropZ05B153PFSecondaryVertexTagInfos = akPuSoftDropZ05B153PFbTagger.SecondaryVertexTagInfos
akPuSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akPuSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akPuSoftDropZ05B153PFCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B153PFbTagger.CombinedSecondaryVertexBJetTags
akPuSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B153PFbTagger.CombinedSecondaryVertexV2BJetTags

akPuSoftDropZ05B153PFSecondaryVertexNegativeTagInfos = akPuSoftDropZ05B153PFbTagger.SecondaryVertexNegativeTagInfos
akPuSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighEffBJetTags = akPuSoftDropZ05B153PFbTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akPuSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighPurBJetTags = akPuSoftDropZ05B153PFbTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akPuSoftDropZ05B153PFNegativeCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B153PFbTagger.NegativeCombinedSecondaryVertexBJetTags
akPuSoftDropZ05B153PFPositiveCombinedSecondaryVertexBJetTags = akPuSoftDropZ05B153PFbTagger.PositiveCombinedSecondaryVertexBJetTags
akPuSoftDropZ05B153PFNegativeCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B153PFbTagger.NegativeCombinedSecondaryVertexV2BJetTags
akPuSoftDropZ05B153PFPositiveCombinedSecondaryVertexV2BJetTags = akPuSoftDropZ05B153PFbTagger.PositiveCombinedSecondaryVertexV2BJetTags

akPuSoftDropZ05B153PFSoftPFMuonsTagInfos = akPuSoftDropZ05B153PFbTagger.SoftPFMuonsTagInfos
akPuSoftDropZ05B153PFSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akPuSoftDropZ05B153PFSoftPFMuonBJetTags = akPuSoftDropZ05B153PFbTagger.SoftPFMuonBJetTags
akPuSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags = akPuSoftDropZ05B153PFbTagger.SoftPFMuonByIP3dBJetTags
akPuSoftDropZ05B153PFSoftPFMuonByPtBJetTags = akPuSoftDropZ05B153PFbTagger.SoftPFMuonByPtBJetTags
akPuSoftDropZ05B153PFNegativeSoftPFMuonByPtBJetTags = akPuSoftDropZ05B153PFbTagger.NegativeSoftPFMuonByPtBJetTags
akPuSoftDropZ05B153PFPositiveSoftPFMuonByPtBJetTags = akPuSoftDropZ05B153PFbTagger.PositiveSoftPFMuonByPtBJetTags
akPuSoftDropZ05B153PFPatJetFlavourIdLegacy = cms.Sequence(akPuSoftDropZ05B153PFPatJetPartonAssociationLegacy*akPuSoftDropZ05B153PFPatJetFlavourAssociationLegacy)
#Not working with our PU sub, but keep it here for reference
#akPuSoftDropZ05B153PFPatJetFlavourAssociation = akPuSoftDropZ05B153PFbTagger.PatJetFlavourAssociation
#akPuSoftDropZ05B153PFPatJetFlavourId = cms.Sequence(akPuSoftDropZ05B153PFPatJetPartons*akPuSoftDropZ05B153PFPatJetFlavourAssociation)

akPuSoftDropZ05B153PFJetBtaggingIP       = cms.Sequence(akPuSoftDropZ05B153PFImpactParameterTagInfos *
            (akPuSoftDropZ05B153PFTrackCountingHighEffBJetTags +
             akPuSoftDropZ05B153PFTrackCountingHighPurBJetTags +
             akPuSoftDropZ05B153PFJetProbabilityBJetTags +
             akPuSoftDropZ05B153PFJetBProbabilityBJetTags 
            )
            )

akPuSoftDropZ05B153PFJetBtaggingSV = cms.Sequence(akPuSoftDropZ05B153PFImpactParameterTagInfos
            *
            akPuSoftDropZ05B153PFSecondaryVertexTagInfos
            * (akPuSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags+
                akPuSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags+
                akPuSoftDropZ05B153PFCombinedSecondaryVertexBJetTags+
                akPuSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags
              )
            )

akPuSoftDropZ05B153PFJetBtaggingNegSV = cms.Sequence(akPuSoftDropZ05B153PFImpactParameterTagInfos
            *
            akPuSoftDropZ05B153PFSecondaryVertexNegativeTagInfos
            * (akPuSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighEffBJetTags+
                akPuSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighPurBJetTags+
                akPuSoftDropZ05B153PFNegativeCombinedSecondaryVertexBJetTags+
                akPuSoftDropZ05B153PFPositiveCombinedSecondaryVertexBJetTags+
                akPuSoftDropZ05B153PFNegativeCombinedSecondaryVertexV2BJetTags+
                akPuSoftDropZ05B153PFPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akPuSoftDropZ05B153PFJetBtaggingMu = cms.Sequence(akPuSoftDropZ05B153PFSoftPFMuonsTagInfos * (akPuSoftDropZ05B153PFSoftPFMuonBJetTags
                +
                akPuSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags
                +
                akPuSoftDropZ05B153PFSoftPFMuonByPtBJetTags
                +
                akPuSoftDropZ05B153PFNegativeSoftPFMuonByPtBJetTags
                +
                akPuSoftDropZ05B153PFPositiveSoftPFMuonByPtBJetTags
              )
            )

akPuSoftDropZ05B153PFJetBtagging = cms.Sequence(akPuSoftDropZ05B153PFJetBtaggingIP
            *akPuSoftDropZ05B153PFJetBtaggingSV
            *akPuSoftDropZ05B153PFJetBtaggingNegSV
#            *akPuSoftDropZ05B153PFJetBtaggingMu
            )

akPuSoftDropZ05B153PFpatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akPuSoftDropZ05B153PFJets"),
        genJetMatch          = cms.InputTag("akPuSoftDropZ05B153PFmatch"),
        genPartonMatch       = cms.InputTag("akPuSoftDropZ05B153PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akPuSoftDropZ05B153PFcorr")),
        JetPartonMapSource   = cms.InputTag("akPuSoftDropZ05B153PFPatJetFlavourAssociationLegacy"),
	JetFlavourInfoSource   = cms.InputTag("akPuSoftDropZ05B153PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akPuSoftDropZ05B153PFJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = True,
        discriminatorSources = cms.VInputTag(cms.InputTag("akPuSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akPuSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akPuSoftDropZ05B153PFCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akPuSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akPuSoftDropZ05B153PFJetBProbabilityBJetTags"),
            cms.InputTag("akPuSoftDropZ05B153PFJetProbabilityBJetTags"),
            #cms.InputTag("akPuSoftDropZ05B153PFSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akPuSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akPuSoftDropZ05B153PFTrackCountingHighEffBJetTags"),
            cms.InputTag("akPuSoftDropZ05B153PFTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akPuSoftDropZ05B153PFJetID"),
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

akPuSoftDropZ05B153PFNjettiness = Njettiness.clone(
		    src = cms.InputTag("akPuSoftDropZ05B153PFJets"),
           	    R0  = cms.double( 0.3)
)
akPuSoftDropZ05B153PFpatJetsWithBtagging.userData.userFloats.src += ['akPuSoftDropZ05B153PFNjettiness:tau1','akPuSoftDropZ05B153PFNjettiness:tau2','akPuSoftDropZ05B153PFNjettiness:tau3']

akPuSoftDropZ05B153PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akPuSoftDropZ05B153PFpatJetsWithBtagging"),
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
                                                             doLifeTimeTagging = cms.untracked.bool(True),
                                                             doLifeTimeTaggingExtras = cms.untracked.bool(False),
                                                             bTagJetName = cms.untracked.string("akPuSoftDropZ05B153PF"),
                                                             jetName = cms.untracked.string("akPuSoftDropZ05B153PF"),
                                                             genPtMin = cms.untracked.double(5),
                                                             hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
							     doTower = cms.untracked.bool(True),
							     doSubJets = cms.untracked.bool(True),
                                                             doGenSubJets = cms.untracked.bool(False),     
                                                             subjetGenTag = cms.untracked.InputTag("akSoftDropZ05B153GenJets"),
                                                             doGenTaus = cms.untracked.bool(False),
                                                             genTau1 = cms.InputTag("akSoftDropZ05B153GenNjettiness","tau1"),
                                                             genTau2 = cms.InputTag("akSoftDropZ05B153GenNjettiness","tau2"),
                                                             genTau3 = cms.InputTag("akSoftDropZ05B153GenNjettiness","tau3"),
                                                             doGenSym = cms.untracked.bool(False),
                                                             genSym = cms.InputTag("akSoftDropZ05B153GenJets","sym"),
                                                             genDroppedBranches = cms.InputTag("akSoftDropZ05B153GenJets","droppedBranches")
                                                             )

akPuSoftDropZ05B153PFJetSequence_mc = cms.Sequence(
                                                  #akPuSoftDropZ05B153PFclean
                                                  #*
                                                  akPuSoftDropZ05B153PFmatch
                                                  #*
                                                  #akPuSoftDropZ05B153PFmatchGroomed
                                                  *
                                                  akPuSoftDropZ05B153PFparton
                                                  *
                                                  akPuSoftDropZ05B153PFcorr
                                                  *
                                                  #akPuSoftDropZ05B153PFJetID
                                                  #*
                                                  akPuSoftDropZ05B153PFPatJetFlavourIdLegacy
                                                  #*
			                          #akPuSoftDropZ05B153PFPatJetFlavourId  # Use legacy algo till PU implemented
                                                  *
                                                  akPuSoftDropZ05B153PFJetTracksAssociatorAtVertex
                                                  *
                                                  akPuSoftDropZ05B153PFJetBtagging
                                                  *
                                                  akPuSoftDropZ05B153PFNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akPuSoftDropZ05B153PFpatJetsWithBtagging
                                                  *
                                                  akPuSoftDropZ05B153PFJetAnalyzer
                                                  )

akPuSoftDropZ05B153PFJetSequence_data = cms.Sequence(akPuSoftDropZ05B153PFcorr
                                                    *
                                                    #akPuSoftDropZ05B153PFJetID
                                                    #*
                                                    akPuSoftDropZ05B153PFJetTracksAssociatorAtVertex
                                                    *
                                                    akPuSoftDropZ05B153PFJetBtagging
                                                    *
                                                    akPuSoftDropZ05B153PFNjettiness 
                                                    *
                                                    akPuSoftDropZ05B153PFpatJetsWithBtagging
                                                    *
                                                    akPuSoftDropZ05B153PFJetAnalyzer
                                                    )

akPuSoftDropZ05B153PFJetSequence_jec = cms.Sequence(akPuSoftDropZ05B153PFJetSequence_mc)
akPuSoftDropZ05B153PFJetSequence_mb = cms.Sequence(akPuSoftDropZ05B153PFJetSequence_mc)

akPuSoftDropZ05B153PFJetSequence = cms.Sequence(akPuSoftDropZ05B153PFJetSequence_jec)
akPuSoftDropZ05B153PFJetAnalyzer.genPtMin = cms.untracked.double(1)
akPuSoftDropZ05B153PFJetAnalyzer.jetPtMin = cms.double(1)
akPuSoftDropZ05B153PFpatJetsWithBtagging.userData.userFloats.src += ['akPuSoftDropZ05B153PFJets:sym']
akPuSoftDropZ05B153PFpatJetsWithBtagging.userData.userInts.src += ['akPuSoftDropZ05B153PFJets:droppedBranches']
