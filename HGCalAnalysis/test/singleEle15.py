import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoLocalCalo.HGCalRecProducers.HGCalLocalRecoSequence_cff')

from FastSimulation.Event.ParticleFilter_cfi import *
#process.load("RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://xrootd-cms.infn.it//store/relval/CMSSW_8_2_0_patch1/RelValDoubleElectronPt15Eta17_27/GEN-SIM-RECO/90X_upgrade2023_realistic_v1_2023D4-v1/10000/0E600BFE-A3E2-E611-8D74-0CC47A4C8E14.root',
        'root://xrootd-cms.infn.it//store/relval/CMSSW_8_2_0_patch1/RelValDoubleElectronPt15Eta17_27/GEN-SIM-RECO/90X_upgrade2023_realistic_v1_2023D4-v1/10000/20D4831D-A4E2-E611-8E46-0CC47A4C8F30.root',
        'root://xrootd-cms.infn.it//store/relval/CMSSW_8_2_0_patch1/RelValDoubleElectronPt15Eta17_27/GEN-SIM-RECO/90X_upgrade2023_realistic_v1_2023D4-v1/10000/5E8E5EFB-A2E2-E611-9508-0CC47A4D76AC.root',
        'root://xrootd-cms.infn.it//store/relval/CMSSW_8_2_0_patch1/RelValDoubleElectronPt15Eta17_27/GEN-SIM-RECO/90X_upgrade2023_realistic_v1_2023D4-v1/10000/68478B06-A3E2-E611-8626-0025905A48EC.root'
    ),
)



process.ana = cms.EDAnalyzer('HGCalAnalysis',
                             detector = cms.string("all"),
                             rawRecHits = cms.bool(True),
                             readOfficialReco = cms.bool(True),
                             readCaloParticles = cms.bool(False),
                             layerClusterPtThreshold = cms.double(0.01),  # All LayerCluster belonging to a multicluster are saved; this Pt threshold applied to the others
                             TestParticleFilter = ParticleFilterBlock.ParticleFilter
                             )

#Quite important
process.ana.TestParticleFilter.protonEMin = cms.double(100000)
process.ana.TestParticleFilter.etaMax = cms.double(3.1)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("singleElePt15-v3.root")

                                   )
#process.imagingClusterHGCal.ecut = cms.double(0.01)
#process.imagingClusterHGCal.eventsToDisplay = cms.untracked.uint32(2)
process.hgcalLayerClusters.minClusters = cms.uint32(0)
process.hgcalLayerClusters.realSpaceCone = cms.bool(True)
process.hgcalLayerClusters.multiclusterRadius = cms.double(2.)
process.hgcalLayerClusters.dependSensor = cms.bool(True)
process.hgcalLayerClusters.ecut = cms.double(3.)
process.hgcalLayerClusters.kappa = cms.double(9.)

process.p = cms.Path(process.HGCalLocalRecoSequence+process.ana)
#process.p = cms.Path( process.ana)
