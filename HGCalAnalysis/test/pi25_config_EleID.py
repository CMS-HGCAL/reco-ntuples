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
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
from FastSimulation.Event.ParticleFilter_cfi import *
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import dEdX_weights as dEdX

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/2637F672-C7A6-E711-B4EF-0025905A612A.root',
        '/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/3C00B396-CBA6-E711-95E1-0025905A612A.root',
        '/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/60FC86A0-CDA6-E711-960D-0025905B856E.root',
        '/store/relval/CMSSW_9_3_2/RelValSinglePiPt25Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/E0E29FFF-C6A6-E711-93A0-003048FFCC16.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    #skipEvents = cms.untracked.uint32(497)
)

process.ana = cms.EDAnalyzer('HGCalAnalysis_EleID',
                             detector = cms.string("all"),
                             rawRecHits = cms.bool(False),
                             readCaloParticles = cms.bool(False),
                             storeGenParticleOrigin = cms.bool(False),
                             storeGenParticleExtrapolation = cms.bool(True),
                             storeElectrons = cms.bool(True),
                             storePCAvariables = cms.bool(True),
                             recomputePCA = cms.bool(True),
                             includeHaloPCA = cms.bool(True),
                             dEdXWeights = dEdX,
                             layerClusterPtThreshold = cms.double(-1),  # All LayerCluster belonging to a multicluster are saved; this Pt threshold applied to the others
                             TestParticleFilter = ParticleFilterBlock.ParticleFilter,
                             EERecHits = cms.InputTag('HGCalRecHit:HGCEERecHits'),
                             FHRecHits = cms.InputTag('HGCalRecHit:HGCHEFRecHits'),
                             BHRecHits = cms.InputTag('HGCalRecHit:HGCHEBRecHits')

)

process.ana.TestParticleFilter.protonEMin = cms.double(100000)
process.ana.TestParticleFilter.etaMax = cms.double(3.1)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("ntuples/eleIDntuple_pi25_n10000.root")
                                   )

#reRunClustering = True
reRunClustering = False

if reRunClustering:
    #process.hgcalLayerClusters.minClusters = cms.uint32(0)
    #process.hgcalLayerClusters.realSpaceCone = cms.bool(True)
    #process.hgcalLayerClusters.multiclusterRadius = cms.double(2.)  # in cm if realSpaceCone is true
    #process.hgcalLayerClusters.dependSensor = cms.bool(True)
    #process.hgcalLayerClusters.ecut = cms.double(3.)  # multiple of sigma noise if dependSensor is true
    #process.hgcalLayerClusters.kappa = cms.double(9.)  # multiple of sigma noise if dependSensor is true
    #process.hgcalLayerClusters.deltac = cms.vdouble(2.,3.,5.) #specify delta c for each subdetector separately
    #process.hgcalLayerClusters.deltac = cms.vdouble(2.3,2.3,5.) #specify delta c for each subdetector separately
    process.p = cms.Path(process.hgcalLayerClusters+process.ana)
else:
    process.p = cms.Path(process.ana)
