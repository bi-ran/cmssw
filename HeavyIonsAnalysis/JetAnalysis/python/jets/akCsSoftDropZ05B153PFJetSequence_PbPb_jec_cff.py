

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akCsSoftDropZ05B153PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akCsSoftDropZ05B153PFJets"),
    matched = cms.InputTag("ak3HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.3
    )

akCsSoftDropZ05B153PFmatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B153HiSignalGenJets"),
    matched = cms.InputTag("ak3HiSignalGenJets"),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = 0.3
    )

akCsSoftDropZ05B153PFparton = patJetPartonMatch.clone(src = cms.InputTag("akCsSoftDropZ05B153PFJets")
                                                        )

akCsSoftDropZ05B153PFcorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akCsSoftDropZ05B153PFJets"),
    payload = "AK3PF_offline"
    )

akCsSoftDropZ05B153PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akCsSoftDropZ05B153CaloJets'))

#akCsSoftDropZ05B153PFclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak3HiSignalGenJets'))

akCsSoftDropZ05B153PFbTagger = bTaggers("akCsSoftDropZ05B153PF",0.3)

#create objects locally since they dont load properly otherwise
#akCsSoftDropZ05B153PFmatch = akCsSoftDropZ05B153PFbTagger.match
akCsSoftDropZ05B153PFparton = patJetPartonMatch.clone(src = cms.InputTag("akCsSoftDropZ05B153PFJets"), matched = cms.InputTag("hiSignalGenParticles"))
akCsSoftDropZ05B153PFPatJetFlavourAssociationLegacy = akCsSoftDropZ05B153PFbTagger.PatJetFlavourAssociationLegacy
akCsSoftDropZ05B153PFPatJetPartons = akCsSoftDropZ05B153PFbTagger.PatJetPartons
akCsSoftDropZ05B153PFJetTracksAssociatorAtVertex = akCsSoftDropZ05B153PFbTagger.JetTracksAssociatorAtVertex
akCsSoftDropZ05B153PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akCsSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags = akCsSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akCsSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags = akCsSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akCsSoftDropZ05B153PFCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B153PFbTagger.CombinedSecondaryVertexBJetTags
akCsSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B153PFbTagger.CombinedSecondaryVertexV2BJetTags
akCsSoftDropZ05B153PFJetBProbabilityBJetTags = akCsSoftDropZ05B153PFbTagger.JetBProbabilityBJetTags
akCsSoftDropZ05B153PFSoftPFMuonByPtBJetTags = akCsSoftDropZ05B153PFbTagger.SoftPFMuonByPtBJetTags
akCsSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags = akCsSoftDropZ05B153PFbTagger.SoftPFMuonByIP3dBJetTags
akCsSoftDropZ05B153PFTrackCountingHighEffBJetTags = akCsSoftDropZ05B153PFbTagger.TrackCountingHighEffBJetTags
akCsSoftDropZ05B153PFTrackCountingHighPurBJetTags = akCsSoftDropZ05B153PFbTagger.TrackCountingHighPurBJetTags
akCsSoftDropZ05B153PFPatJetPartonAssociationLegacy = akCsSoftDropZ05B153PFbTagger.PatJetPartonAssociationLegacy

akCsSoftDropZ05B153PFImpactParameterTagInfos = akCsSoftDropZ05B153PFbTagger.ImpactParameterTagInfos
akCsSoftDropZ05B153PFImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akCsSoftDropZ05B153PFJetProbabilityBJetTags = akCsSoftDropZ05B153PFbTagger.JetProbabilityBJetTags

akCsSoftDropZ05B153PFSecondaryVertexTagInfos = akCsSoftDropZ05B153PFbTagger.SecondaryVertexTagInfos
akCsSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags = akCsSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akCsSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags = akCsSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akCsSoftDropZ05B153PFCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B153PFbTagger.CombinedSecondaryVertexBJetTags
akCsSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B153PFbTagger.CombinedSecondaryVertexV2BJetTags

akCsSoftDropZ05B153PFSecondaryVertexNegativeTagInfos = akCsSoftDropZ05B153PFbTagger.SecondaryVertexNegativeTagInfos
akCsSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighEffBJetTags = akCsSoftDropZ05B153PFbTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akCsSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighPurBJetTags = akCsSoftDropZ05B153PFbTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akCsSoftDropZ05B153PFNegativeCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B153PFbTagger.NegativeCombinedSecondaryVertexBJetTags
akCsSoftDropZ05B153PFPositiveCombinedSecondaryVertexBJetTags = akCsSoftDropZ05B153PFbTagger.PositiveCombinedSecondaryVertexBJetTags
akCsSoftDropZ05B153PFNegativeCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B153PFbTagger.NegativeCombinedSecondaryVertexV2BJetTags
akCsSoftDropZ05B153PFPositiveCombinedSecondaryVertexV2BJetTags = akCsSoftDropZ05B153PFbTagger.PositiveCombinedSecondaryVertexV2BJetTags

akCsSoftDropZ05B153PFSoftPFMuonsTagInfos = akCsSoftDropZ05B153PFbTagger.SoftPFMuonsTagInfos
akCsSoftDropZ05B153PFSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akCsSoftDropZ05B153PFSoftPFMuonBJetTags = akCsSoftDropZ05B153PFbTagger.SoftPFMuonBJetTags
akCsSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags = akCsSoftDropZ05B153PFbTagger.SoftPFMuonByIP3dBJetTags
akCsSoftDropZ05B153PFSoftPFMuonByPtBJetTags = akCsSoftDropZ05B153PFbTagger.SoftPFMuonByPtBJetTags
akCsSoftDropZ05B153PFNegativeSoftPFMuonByPtBJetTags = akCsSoftDropZ05B153PFbTagger.NegativeSoftPFMuonByPtBJetTags
akCsSoftDropZ05B153PFPositiveSoftPFMuonByPtBJetTags = akCsSoftDropZ05B153PFbTagger.PositiveSoftPFMuonByPtBJetTags
akCsSoftDropZ05B153PFPatJetFlavourIdLegacy = cms.Sequence(akCsSoftDropZ05B153PFPatJetPartonAssociationLegacy*akCsSoftDropZ05B153PFPatJetFlavourAssociationLegacy)
#Not working with our PU sub, but keep it here for reference
#akCsSoftDropZ05B153PFPatJetFlavourAssociation = akCsSoftDropZ05B153PFbTagger.PatJetFlavourAssociation
#akCsSoftDropZ05B153PFPatJetFlavourId = cms.Sequence(akCsSoftDropZ05B153PFPatJetPartons*akCsSoftDropZ05B153PFPatJetFlavourAssociation)

akCsSoftDropZ05B153PFJetBtaggingIP       = cms.Sequence(akCsSoftDropZ05B153PFImpactParameterTagInfos *
            (akCsSoftDropZ05B153PFTrackCountingHighEffBJetTags +
             akCsSoftDropZ05B153PFTrackCountingHighPurBJetTags +
             akCsSoftDropZ05B153PFJetProbabilityBJetTags +
             akCsSoftDropZ05B153PFJetBProbabilityBJetTags 
            )
            )

akCsSoftDropZ05B153PFJetBtaggingSV = cms.Sequence(akCsSoftDropZ05B153PFImpactParameterTagInfos
            *
            akCsSoftDropZ05B153PFSecondaryVertexTagInfos
            * (akCsSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags+
                akCsSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags+
                akCsSoftDropZ05B153PFCombinedSecondaryVertexBJetTags+
                akCsSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags
              )
            )

akCsSoftDropZ05B153PFJetBtaggingNegSV = cms.Sequence(akCsSoftDropZ05B153PFImpactParameterTagInfos
            *
            akCsSoftDropZ05B153PFSecondaryVertexNegativeTagInfos
            * (akCsSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighEffBJetTags+
                akCsSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighPurBJetTags+
                akCsSoftDropZ05B153PFNegativeCombinedSecondaryVertexBJetTags+
                akCsSoftDropZ05B153PFPositiveCombinedSecondaryVertexBJetTags+
                akCsSoftDropZ05B153PFNegativeCombinedSecondaryVertexV2BJetTags+
                akCsSoftDropZ05B153PFPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akCsSoftDropZ05B153PFJetBtaggingMu = cms.Sequence(akCsSoftDropZ05B153PFSoftPFMuonsTagInfos * (akCsSoftDropZ05B153PFSoftPFMuonBJetTags
                +
                akCsSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags
                +
                akCsSoftDropZ05B153PFSoftPFMuonByPtBJetTags
                +
                akCsSoftDropZ05B153PFNegativeSoftPFMuonByPtBJetTags
                +
                akCsSoftDropZ05B153PFPositiveSoftPFMuonByPtBJetTags
              )
            )

akCsSoftDropZ05B153PFJetBtagging = cms.Sequence(akCsSoftDropZ05B153PFJetBtaggingIP
            *akCsSoftDropZ05B153PFJetBtaggingSV
            *akCsSoftDropZ05B153PFJetBtaggingNegSV
#            *akCsSoftDropZ05B153PFJetBtaggingMu
            )

akCsSoftDropZ05B153PFpatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akCsSoftDropZ05B153PFJets"),
        genJetMatch          = cms.InputTag("akCsSoftDropZ05B153PFmatch"),
        genPartonMatch       = cms.InputTag("akCsSoftDropZ05B153PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akCsSoftDropZ05B153PFcorr")),
        JetPartonMapSource   = cms.InputTag("akCsSoftDropZ05B153PFPatJetFlavourAssociationLegacy"),
	JetFlavourInfoSource   = cms.InputTag("akCsSoftDropZ05B153PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akCsSoftDropZ05B153PFJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = True,
        discriminatorSources = cms.VInputTag(cms.InputTag("akCsSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akCsSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akCsSoftDropZ05B153PFCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akCsSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akCsSoftDropZ05B153PFJetBProbabilityBJetTags"),
            cms.InputTag("akCsSoftDropZ05B153PFJetProbabilityBJetTags"),
            #cms.InputTag("akCsSoftDropZ05B153PFSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akCsSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akCsSoftDropZ05B153PFTrackCountingHighEffBJetTags"),
            cms.InputTag("akCsSoftDropZ05B153PFTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akCsSoftDropZ05B153PFJetID"),
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

akCsSoftDropZ05B153PFNjettiness = Njettiness.clone(
		    src = cms.InputTag("akCsSoftDropZ05B153PFJets"),
           	    R0  = cms.double( 0.3)
)
akCsSoftDropZ05B153PFpatJetsWithBtagging.userData.userFloats.src += ['akCsSoftDropZ05B153PFNjettiness:tau1','akCsSoftDropZ05B153PFNjettiness:tau2','akCsSoftDropZ05B153PFNjettiness:tau3']

akCsSoftDropZ05B153PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akCsSoftDropZ05B153PFpatJetsWithBtagging"),
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
                                                             bTagJetName = cms.untracked.string("akCsSoftDropZ05B153PF"),
                                                             jetName = cms.untracked.string("akCsSoftDropZ05B153PF"),
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

akCsSoftDropZ05B153PFJetSequence_mc = cms.Sequence(
                                                  #akCsSoftDropZ05B153PFclean
                                                  #*
                                                  akCsSoftDropZ05B153PFmatch
                                                  #*
                                                  #akCsSoftDropZ05B153PFmatchGroomed
                                                  *
                                                  akCsSoftDropZ05B153PFparton
                                                  *
                                                  akCsSoftDropZ05B153PFcorr
                                                  *
                                                  #akCsSoftDropZ05B153PFJetID
                                                  #*
                                                  akCsSoftDropZ05B153PFPatJetFlavourIdLegacy
                                                  #*
			                          #akCsSoftDropZ05B153PFPatJetFlavourId  # Use legacy algo till PU implemented
                                                  *
                                                  akCsSoftDropZ05B153PFJetTracksAssociatorAtVertex
                                                  *
                                                  akCsSoftDropZ05B153PFJetBtagging
                                                  *
                                                  akCsSoftDropZ05B153PFNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akCsSoftDropZ05B153PFpatJetsWithBtagging
                                                  *
                                                  akCsSoftDropZ05B153PFJetAnalyzer
                                                  )

akCsSoftDropZ05B153PFJetSequence_data = cms.Sequence(akCsSoftDropZ05B153PFcorr
                                                    *
                                                    #akCsSoftDropZ05B153PFJetID
                                                    #*
                                                    akCsSoftDropZ05B153PFJetTracksAssociatorAtVertex
                                                    *
                                                    akCsSoftDropZ05B153PFJetBtagging
                                                    *
                                                    akCsSoftDropZ05B153PFNjettiness 
                                                    *
                                                    akCsSoftDropZ05B153PFpatJetsWithBtagging
                                                    *
                                                    akCsSoftDropZ05B153PFJetAnalyzer
                                                    )

akCsSoftDropZ05B153PFJetSequence_jec = cms.Sequence(akCsSoftDropZ05B153PFJetSequence_mc)
akCsSoftDropZ05B153PFJetSequence_mb = cms.Sequence(akCsSoftDropZ05B153PFJetSequence_mc)

akCsSoftDropZ05B153PFJetSequence = cms.Sequence(akCsSoftDropZ05B153PFJetSequence_jec)
akCsSoftDropZ05B153PFJetAnalyzer.genPtMin = cms.untracked.double(1)
akCsSoftDropZ05B153PFJetAnalyzer.jetPtMin = cms.double(1)
akCsSoftDropZ05B153PFpatJetsWithBtagging.userData.userFloats.src += ['akCsSoftDropZ05B153PFJets:sym']
akCsSoftDropZ05B153PFpatJetsWithBtagging.userData.userInts.src += ['akCsSoftDropZ05B153PFJets:droppedBranches']
