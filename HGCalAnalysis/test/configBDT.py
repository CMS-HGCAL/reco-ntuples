import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoLocalCalo.HGCalRecProducers.HGCalLocalRecoSequence_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
from FastSimulation.Event.ParticleFilter_cfi import *
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import dEdX_weights as dEdX

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.load('EgammaTools.EgammaAnalysis.HGCalElectronFilter_cfi')

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #        'root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_3_2/RelValZEE_14/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/2E2C7104-E1A6-E711-817B-0CC47A7C3612.root'
        #        '/store/relval/CMSSW_9_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v0_D17PU200-v1/00000/082041D3-1E89-E711-80AA-0242AC130002.root'
        #       'root://polgrid4.in2p3.fr//store/relval/CMSSW_9_3_2/RelValSingleElectronPt15Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/D0226E59-00A7-E711-AA65-0CC47A78A426.root'
        #'root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_3_2/RelValZEE_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/10000/5CBBD7E3-F7AA-E711-904B-E0071B7A0670.root',
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/BE500A88-5EAD-E711-8E07-0242AC110012.root',
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/78E43CDF-5EAD-E711-BA4A-0242AC110002.root',
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/04CC740B-5FAD-E711-9A4A-0242AC110008.root',
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/A8640BD9-5EAD-E711-9B9C-0242AC11000D.root',
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/7C369899-61AD-E711-BB03-0242AC110006.root',
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/A0A41AA9-88AD-E711-8744-0242AC110002.root',
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/6AB7D61B-91AD-E711-937F-0242AC11000F.root',
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/68C29CB7-8AAD-E711-A9D8-0242AC110003.root',
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/7AED1C74-6BAD-E711-ABA3-0242AC110010.root',
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/4409A3B7-8EAD-E711-800A-0242AC110009.root',


        'root://cms-xrd-global.cern.ch///store/mc/PhaseIITDRFall17DR/SingleElectronPt5_100Eta1p6_2p8/GEN-SIM-RECO/PU200FEVT_93X_upgrade2023_realistic_v2-v2/150000/9A58FEE3-CDB4-E711-8B96-02163E013C39.root',
        'root://cms-xrd-global.cern.ch///store/mc/PhaseIITDRFall17DR/SingleElectronPt5_100Eta1p6_2p8/GEN-SIM-RECO/PU200FEVT_93X_upgrade2023_realistic_v2-v2/150000/58060787-D4B4-E711-82F3-FA163EF6B191.root',
        'root://cms-xrd-global.cern.ch///store/mc/PhaseIITDRFall17DR/SingleElectronPt5_100Eta1p6_2p8/GEN-SIM-RECO/PU200FEVT_93X_upgrade2023_realistic_v2-v2/150000/02298E27-EAB4-E711-90F2-02163E016479.root',
        'root://cms-xrd-global.cern.ch///store/mc/PhaseIITDRFall17DR/SingleElectronPt5_100Eta1p6_2p8/GEN-SIM-RECO/PU200FEVT_93X_upgrade2023_realistic_v2-v2/150000/14DDB6AE-ECB4-E711-A63B-02163E0152B7.root',
        'root://cms-xrd-global.cern.ch///store/mc/PhaseIITDRFall17DR/SingleElectronPt5_100Eta1p6_2p8/GEN-SIM-RECO/PU200FEVT_93X_upgrade2023_realistic_v2-v2/150000/6CF0C949-F1B4-E711-8C8C-FA163E648B06.root',
        'root://cms-xrd-global.cern.ch///store/mc/PhaseIITDRFall17DR/SingleElectronPt5_100Eta1p6_2p8/GEN-SIM-RECO/PU200FEVT_93X_upgrade2023_realistic_v2-v2/150000/7C0D454A-E4B4-E711-8DB7-FA163E641E47.root',
        'root://cms-xrd-global.cern.ch///store/mc/PhaseIITDRFall17DR/SingleElectronPt5_100Eta1p6_2p8/GEN-SIM-RECO/PU200FEVT_93X_upgrade2023_realistic_v2-v2/150000/F401F719-F1B4-E711-86C2-FA163EFA0A74.root',
        'root://cms-xrd-global.cern.ch///store/mc/PhaseIITDRFall17DR/SingleElectronPt5_100Eta1p6_2p8/GEN-SIM-RECO/PU200FEVT_93X_upgrade2023_realistic_v2-v2/150000/5A5EE0DB-F6B4-E711-B4B1-FA163EED17DF.root',
        'root://cms-xrd-global.cern.ch///store/mc/PhaseIITDRFall17DR/SingleElectronPt5_100Eta1p6_2p8/GEN-SIM-RECO/PU200FEVT_93X_upgrade2023_realistic_v2-v2/150000/EA2981B9-03B5-E711-B8E1-FA163EE639D0.root',
        'root://cms-xrd-global.cern.ch///store/mc/PhaseIITDRFall17DR/SingleElectronPt5_100Eta1p6_2p8/GEN-SIM-RECO/PU200FEVT_93X_upgrade2023_realistic_v2-v2/150000/98B2D3CA-FAB4-E711-9219-FA163E2A6957.root'
    ),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.prunedGenParticles = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop  *", # this is the default
        "keep++ pdgId = {Z0}",

#    "drop pdgId = {Z0} & status = 2"
    )
)

process.electrons = cms.EDFilter("PdgIdCandViewSelector",
                                 src = cms.InputTag("prunedGenParticles"),
                                 pdgId = cms.vint32( 11 )
)

process.electronFilter = cms.EDFilter("CandViewCountFilter",
                                      src = cms.InputTag("electrons"),
                                      minNumber = cms.uint32(1),
)

process.ana = cms.EDAnalyzer('HGCalAnalysis_EleID',
                             detector = cms.string("all"),
                             rawRecHits = cms.bool(True),
                             readCaloParticles = cms.bool(False),
                             storeGenParticleOrigin = cms.bool(True),
                             readGenParticles = cms.bool(True),
                             GenParticles = cms.InputTag('genParticles'),
                             storeGenParticleExtrapolation = cms.bool(False),
                             storePCAvariables = cms.bool(False),
                             storeElectrons = cms.bool(True),
                             recomputePCA = cms.bool(False),
                             includeHaloPCA = cms.bool(True),
                             dEdXWeights = dEdX,
                             layerClusterPtThreshold = cms.double(-1),  # All LayerCluster belonging to a multicluster are saved; this Pt threshold applied to the others
                             TestParticleFilter = ParticleFilterBlock.ParticleFilter,
                             EERecHits = cms.InputTag('HGCalRecHit:HGCEERecHits'),
                             FHRecHits = cms.InputTag('HGCalRecHit:HGCHEFRecHits'),
                             BHRecHits = cms.InputTag('HGCalRecHit:HGCHEBRecHits'),
                             Electrons = cms.InputTag('cleanedEcalDrivenGsfElectronsFromMultiCl'),
                             PFMultiClusters = cms.InputTag('particleFlowClusterHGCalFromMultiCl'),
                             pcaRadius = cms.double(3.),
                             vertices = cms.InputTag("offlinePrimaryVertices"),
                             puSummary = cms.InputTag("addPileupInfo"),
                             beamSpot = cms.InputTag("offlineBeamSpot"),
                             barrelLowPt = cms.FileInPath('EgammaTools/EgammaAnalysis/data/EIDmva_EB_1020_oldbarreltdrDR01_BDT.weights.xml'),
                             barrelHighPt = cms.FileInPath('EgammaTools/EgammaAnalysis/data/EIDmva_EB_20_oldbarreltdrDR01_BDT.weights.xml'),
                             endcapLowPt = cms.FileInPath('EgammaTools/EgammaAnalysis/data/HGCEIDmva_1020_trackepshowernoisolonghgcaltdrV3DR01preselmatch_BDT.weights.xml'),
                             endcapHighPt = cms.FileInPath('EgammaTools/EgammaAnalysis/data/HGCEIDmva_20_trackepshowernoisolonghgcaltdrV3DR01preselmatch_BDT.weights.xml'),
                             GsfElectrons = cms.InputTag('cleanedEcalDrivenGsfElectronsFromMultiCl'),

)

process.ana.TestParticleFilter.protonEMin = cms.double(100000)
process.ana.TestParticleFilter.etaMax = cms.double(3.1)

fileFormat = 'AOD'

process.ntuplizer = cms.EDAnalyzer('Ntuplizer',
                                   EleTag      = cms.InputTag('gedGsfElectrons'),
                                   VerticesTag = cms.InputTag('offlinePrimaryVertices'),
                                   isMC = cms.bool(True),
                                   GenParticles = cms.InputTag('prunedGenParticles')
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("hgcalNtupleEleGun_n1000.root")

                                   )
#  process.p = cms.Path(process.prunedGenParticles+process.electrons+process.electronFilter+process.ana+process.ntuplizer)
#    process.p = cms.Path(process.prunedGenParticles+process.electrons+process.electronFilter+process.cleanedEcalDrivenGsfElectronsFromMultiCl+ process.ana)
process.p = cms.Path(process.cleanedEcalDrivenGsfElectronsFromMultiCl+ process.ana)
