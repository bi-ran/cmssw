

import FWCore.ParameterSet.Config as cms
from HeavyIonsAnalysis.JetAnalysis.patHeavyIonSequences_cff import patJetGenJetMatch, patJetPartonMatch, patJetCorrFactors, patJets
from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *
from HeavyIonsAnalysis.JetAnalysis.bTaggers_cff import *
from RecoJets.JetProducers.JetIDParams_cfi import *
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

akSoftDropZ05B153PFmatch = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B153PFJets"),
    matched = cms.InputTag("ak3HiCleanedGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.3
    )

akSoftDropZ05B153PFmatchGroomed = patJetGenJetMatch.clone(
    src = cms.InputTag("akSoftDropZ05B153HiGenJets"),
    matched = cms.InputTag("ak3HiCleanedGenJets"),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = 0.3
    )

akSoftDropZ05B153PFparton = patJetPartonMatch.clone(src = cms.InputTag("akSoftDropZ05B153PFJets")
                                                        )

akSoftDropZ05B153PFcorr = patJetCorrFactors.clone(
    useNPV = cms.bool(False),
    useRho = cms.bool(False),
#    primaryVertices = cms.InputTag("hiSelectedVertex"),
    levels   = cms.vstring('L2Relative','L3Absolute'),
    src = cms.InputTag("akSoftDropZ05B153PFJets"),
    payload = "AK3PF_offline"
    )

akSoftDropZ05B153PFJetID= cms.EDProducer('JetIDProducer', JetIDParams, src = cms.InputTag('akSoftDropZ05B153CaloJets'))

#akSoftDropZ05B153PFclean   = heavyIonCleanedGenJets.clone(src = cms.InputTag('ak3HiCleanedGenJets'))

akSoftDropZ05B153PFbTagger = bTaggers("akSoftDropZ05B153PF",0.3)

#create objects locally since they dont load properly otherwise
#akSoftDropZ05B153PFmatch = akSoftDropZ05B153PFbTagger.match
akSoftDropZ05B153PFparton = patJetPartonMatch.clone(src = cms.InputTag("akSoftDropZ05B153PFJets"), matched = cms.InputTag("selectedPartons"))
akSoftDropZ05B153PFPatJetFlavourAssociationLegacy = akSoftDropZ05B153PFbTagger.PatJetFlavourAssociationLegacy
akSoftDropZ05B153PFPatJetPartons = akSoftDropZ05B153PFbTagger.PatJetPartons
akSoftDropZ05B153PFJetTracksAssociatorAtVertex = akSoftDropZ05B153PFbTagger.JetTracksAssociatorAtVertex
akSoftDropZ05B153PFJetTracksAssociatorAtVertex.tracks = cms.InputTag("highPurityTracks")
akSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akSoftDropZ05B153PFCombinedSecondaryVertexBJetTags = akSoftDropZ05B153PFbTagger.CombinedSecondaryVertexBJetTags
akSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B153PFbTagger.CombinedSecondaryVertexV2BJetTags
akSoftDropZ05B153PFJetBProbabilityBJetTags = akSoftDropZ05B153PFbTagger.JetBProbabilityBJetTags
akSoftDropZ05B153PFSoftPFMuonByPtBJetTags = akSoftDropZ05B153PFbTagger.SoftPFMuonByPtBJetTags
akSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags = akSoftDropZ05B153PFbTagger.SoftPFMuonByIP3dBJetTags
akSoftDropZ05B153PFTrackCountingHighEffBJetTags = akSoftDropZ05B153PFbTagger.TrackCountingHighEffBJetTags
akSoftDropZ05B153PFTrackCountingHighPurBJetTags = akSoftDropZ05B153PFbTagger.TrackCountingHighPurBJetTags
akSoftDropZ05B153PFPatJetPartonAssociationLegacy = akSoftDropZ05B153PFbTagger.PatJetPartonAssociationLegacy

akSoftDropZ05B153PFImpactParameterTagInfos = akSoftDropZ05B153PFbTagger.ImpactParameterTagInfos
akSoftDropZ05B153PFImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akSoftDropZ05B153PFJetProbabilityBJetTags = akSoftDropZ05B153PFbTagger.JetProbabilityBJetTags

akSoftDropZ05B153PFSecondaryVertexTagInfos = akSoftDropZ05B153PFbTagger.SecondaryVertexTagInfos
akSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighEffBJetTags
akSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B153PFbTagger.SimpleSecondaryVertexHighPurBJetTags
akSoftDropZ05B153PFCombinedSecondaryVertexBJetTags = akSoftDropZ05B153PFbTagger.CombinedSecondaryVertexBJetTags
akSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B153PFbTagger.CombinedSecondaryVertexV2BJetTags

akSoftDropZ05B153PFSecondaryVertexNegativeTagInfos = akSoftDropZ05B153PFbTagger.SecondaryVertexNegativeTagInfos
akSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighEffBJetTags = akSoftDropZ05B153PFbTagger.NegativeSimpleSecondaryVertexHighEffBJetTags
akSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighPurBJetTags = akSoftDropZ05B153PFbTagger.NegativeSimpleSecondaryVertexHighPurBJetTags
akSoftDropZ05B153PFNegativeCombinedSecondaryVertexBJetTags = akSoftDropZ05B153PFbTagger.NegativeCombinedSecondaryVertexBJetTags
akSoftDropZ05B153PFPositiveCombinedSecondaryVertexBJetTags = akSoftDropZ05B153PFbTagger.PositiveCombinedSecondaryVertexBJetTags
akSoftDropZ05B153PFNegativeCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B153PFbTagger.NegativeCombinedSecondaryVertexV2BJetTags
akSoftDropZ05B153PFPositiveCombinedSecondaryVertexV2BJetTags = akSoftDropZ05B153PFbTagger.PositiveCombinedSecondaryVertexV2BJetTags

akSoftDropZ05B153PFSoftPFMuonsTagInfos = akSoftDropZ05B153PFbTagger.SoftPFMuonsTagInfos
akSoftDropZ05B153PFSoftPFMuonsTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
akSoftDropZ05B153PFSoftPFMuonBJetTags = akSoftDropZ05B153PFbTagger.SoftPFMuonBJetTags
akSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags = akSoftDropZ05B153PFbTagger.SoftPFMuonByIP3dBJetTags
akSoftDropZ05B153PFSoftPFMuonByPtBJetTags = akSoftDropZ05B153PFbTagger.SoftPFMuonByPtBJetTags
akSoftDropZ05B153PFNegativeSoftPFMuonByPtBJetTags = akSoftDropZ05B153PFbTagger.NegativeSoftPFMuonByPtBJetTags
akSoftDropZ05B153PFPositiveSoftPFMuonByPtBJetTags = akSoftDropZ05B153PFbTagger.PositiveSoftPFMuonByPtBJetTags
akSoftDropZ05B153PFPatJetFlavourIdLegacy = cms.Sequence(akSoftDropZ05B153PFPatJetPartonAssociationLegacy*akSoftDropZ05B153PFPatJetFlavourAssociationLegacy)
#Not working with our PU sub, but keep it here for reference
#akSoftDropZ05B153PFPatJetFlavourAssociation = akSoftDropZ05B153PFbTagger.PatJetFlavourAssociation
#akSoftDropZ05B153PFPatJetFlavourId = cms.Sequence(akSoftDropZ05B153PFPatJetPartons*akSoftDropZ05B153PFPatJetFlavourAssociation)

akSoftDropZ05B153PFJetBtaggingIP       = cms.Sequence(akSoftDropZ05B153PFImpactParameterTagInfos *
            (akSoftDropZ05B153PFTrackCountingHighEffBJetTags +
             akSoftDropZ05B153PFTrackCountingHighPurBJetTags +
             akSoftDropZ05B153PFJetProbabilityBJetTags +
             akSoftDropZ05B153PFJetBProbabilityBJetTags 
            )
            )

akSoftDropZ05B153PFJetBtaggingSV = cms.Sequence(akSoftDropZ05B153PFImpactParameterTagInfos
            *
            akSoftDropZ05B153PFSecondaryVertexTagInfos
            * (akSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags+
                akSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags+
                akSoftDropZ05B153PFCombinedSecondaryVertexBJetTags+
                akSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags
              )
            )

akSoftDropZ05B153PFJetBtaggingNegSV = cms.Sequence(akSoftDropZ05B153PFImpactParameterTagInfos
            *
            akSoftDropZ05B153PFSecondaryVertexNegativeTagInfos
            * (akSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighEffBJetTags+
                akSoftDropZ05B153PFNegativeSimpleSecondaryVertexHighPurBJetTags+
                akSoftDropZ05B153PFNegativeCombinedSecondaryVertexBJetTags+
                akSoftDropZ05B153PFPositiveCombinedSecondaryVertexBJetTags+
                akSoftDropZ05B153PFNegativeCombinedSecondaryVertexV2BJetTags+
                akSoftDropZ05B153PFPositiveCombinedSecondaryVertexV2BJetTags
              )
            )

akSoftDropZ05B153PFJetBtaggingMu = cms.Sequence(akSoftDropZ05B153PFSoftPFMuonsTagInfos * (akSoftDropZ05B153PFSoftPFMuonBJetTags
                +
                akSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags
                +
                akSoftDropZ05B153PFSoftPFMuonByPtBJetTags
                +
                akSoftDropZ05B153PFNegativeSoftPFMuonByPtBJetTags
                +
                akSoftDropZ05B153PFPositiveSoftPFMuonByPtBJetTags
              )
            )

akSoftDropZ05B153PFJetBtagging = cms.Sequence(akSoftDropZ05B153PFJetBtaggingIP
            *akSoftDropZ05B153PFJetBtaggingSV
            *akSoftDropZ05B153PFJetBtaggingNegSV
#            *akSoftDropZ05B153PFJetBtaggingMu
            )

akSoftDropZ05B153PFpatJetsWithBtagging = patJets.clone(jetSource = cms.InputTag("akSoftDropZ05B153PFJets"),
        genJetMatch          = cms.InputTag("akSoftDropZ05B153PFmatch"),
        genPartonMatch       = cms.InputTag("akSoftDropZ05B153PFparton"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("akSoftDropZ05B153PFcorr")),
        JetPartonMapSource   = cms.InputTag("akSoftDropZ05B153PFPatJetFlavourAssociationLegacy"),
	JetFlavourInfoSource   = cms.InputTag("akSoftDropZ05B153PFPatJetFlavourAssociation"),
        trackAssociationSource = cms.InputTag("akSoftDropZ05B153PFJetTracksAssociatorAtVertex"),
	useLegacyJetMCFlavour = True,
        discriminatorSources = cms.VInputTag(cms.InputTag("akSoftDropZ05B153PFSimpleSecondaryVertexHighEffBJetTags"),
            cms.InputTag("akSoftDropZ05B153PFSimpleSecondaryVertexHighPurBJetTags"),
            cms.InputTag("akSoftDropZ05B153PFCombinedSecondaryVertexBJetTags"),
            cms.InputTag("akSoftDropZ05B153PFCombinedSecondaryVertexV2BJetTags"),
            cms.InputTag("akSoftDropZ05B153PFJetBProbabilityBJetTags"),
            cms.InputTag("akSoftDropZ05B153PFJetProbabilityBJetTags"),
            #cms.InputTag("akSoftDropZ05B153PFSoftPFMuonByPtBJetTags"),
            #cms.InputTag("akSoftDropZ05B153PFSoftPFMuonByIP3dBJetTags"),
            cms.InputTag("akSoftDropZ05B153PFTrackCountingHighEffBJetTags"),
            cms.InputTag("akSoftDropZ05B153PFTrackCountingHighPurBJetTags"),
            ),
        jetIDMap = cms.InputTag("akSoftDropZ05B153PFJetID"),
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

akSoftDropZ05B153PFNjettiness = Njettiness.clone(
		    src = cms.InputTag("akSoftDropZ05B153PFJets"),
           	    R0  = cms.double( 0.3)
)
akSoftDropZ05B153PFpatJetsWithBtagging.userData.userFloats.src += ['akSoftDropZ05B153PFNjettiness:tau1','akSoftDropZ05B153PFNjettiness:tau2','akSoftDropZ05B153PFNjettiness:tau3']

akSoftDropZ05B153PFJetAnalyzer = inclusiveJetAnalyzer.clone(jetTag = cms.InputTag("akSoftDropZ05B153PFpatJetsWithBtagging"),
                                                             genjetTag = 'ak3HiGenJets',
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
                                                             bTagJetName = cms.untracked.string("akSoftDropZ05B153PF"),
                                                             jetName = cms.untracked.string("akSoftDropZ05B153PF"),
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

akSoftDropZ05B153PFJetSequence_mc = cms.Sequence(
                                                  #akSoftDropZ05B153PFclean
                                                  #*
                                                  akSoftDropZ05B153PFmatch
                                                  #*
                                                  #akSoftDropZ05B153PFmatchGroomed
                                                  *
                                                  akSoftDropZ05B153PFparton
                                                  *
                                                  akSoftDropZ05B153PFcorr
                                                  *
                                                  #akSoftDropZ05B153PFJetID
                                                  #*
                                                  akSoftDropZ05B153PFPatJetFlavourIdLegacy
                                                  #*
			                          #akSoftDropZ05B153PFPatJetFlavourId  # Use legacy algo till PU implemented
                                                  *
                                                  akSoftDropZ05B153PFJetTracksAssociatorAtVertex
                                                  *
                                                  akSoftDropZ05B153PFJetBtagging
                                                  *
                                                  akSoftDropZ05B153PFNjettiness #No constituents for calo jets in pp. Must be removed for pp calo jets but I'm not sure how to do this transparently (Marta)
                                                  *
                                                  akSoftDropZ05B153PFpatJetsWithBtagging
                                                  *
                                                  akSoftDropZ05B153PFJetAnalyzer
                                                  )

akSoftDropZ05B153PFJetSequence_data = cms.Sequence(akSoftDropZ05B153PFcorr
                                                    *
                                                    #akSoftDropZ05B153PFJetID
                                                    #*
                                                    akSoftDropZ05B153PFJetTracksAssociatorAtVertex
                                                    *
                                                    akSoftDropZ05B153PFJetBtagging
                                                    *
                                                    akSoftDropZ05B153PFNjettiness 
                                                    *
                                                    akSoftDropZ05B153PFpatJetsWithBtagging
                                                    *
                                                    akSoftDropZ05B153PFJetAnalyzer
                                                    )

akSoftDropZ05B153PFJetSequence_jec = cms.Sequence(akSoftDropZ05B153PFJetSequence_mc)
akSoftDropZ05B153PFJetSequence_mb = cms.Sequence(akSoftDropZ05B153PFJetSequence_mc)

akSoftDropZ05B153PFJetSequence = cms.Sequence(akSoftDropZ05B153PFJetSequence_mb)
akSoftDropZ05B153PFpatJetsWithBtagging.userData.userFloats.src += ['akSoftDropZ05B153PFJets:sym']
akSoftDropZ05B153PFpatJetsWithBtagging.userData.userInts.src += ['akSoftDropZ05B153PFJets:droppedBranches']
