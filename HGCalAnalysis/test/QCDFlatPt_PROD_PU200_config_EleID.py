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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'root://cms-xrd-global.cern.ch//store/relval/CMSSW_9_3_2/RelValZEE_14/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/2E2C7104-E1A6-E711-817B-0CC47A7C3612.root'
#        '/store/relval/CMSSW_9_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v0_D17PU200-v1/00000/082041D3-1E89-E711-80AA-0242AC130002.root'
#       'root://polgrid4.in2p3.fr//store/relval/CMSSW_9_3_2/RelValSingleElectronPt15Eta1p7_2p7/GEN-SIM-RECO/93X_upgrade2023_realistic_v2_2023D17noPU-v1/10000/D0226E59-00A7-E711-AA65-0CC47A78A426.root'
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/ACFA776A-5CAE-E711-A251-0025905C2CA4.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/AA832E35-77AE-E711-80F9-0025904C5DD8.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/622A6148-73AE-E711-AEB1-0025904CDDEE.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/54CC2B0C-AFAF-E711-B23D-0025905C54BA.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/F0417CC9-B2AF-E711-883A-0025904C68D8.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/70273C7A-B4AF-E711-9DD2-0025904C68DC.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/48E59B54-B2AF-E711-AC7E-0025905D1E00.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/9A9D8351-B4AF-E711-A7CF-0025904C66A4.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/34F77174-B4AF-E711-A390-0025905D1D7A.root',
'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/00000/2C0CCB5F-9AAF-E711-B301-0CC47AF9B496.root',
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.prunedGenParticles = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "keep  pt > 1.5 & abs(eta) < 5"
    )
)

process.ana = cms.EDAnalyzer('HGCalAnalysis_EleID',
                             detector = cms.string("all"),
                             rawRecHits = cms.bool(True),
                             readCaloParticles = cms.bool(False),
                             storeGenParticleOrigin = cms.bool(True),
                             readGenParticles = cms.bool(True),
                             GenParticles = cms.InputTag('prunedGenParticles'),
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
                             BHRecHits = cms.InputTag('HGCalRecHit:HGCHEBRecHits')
)

process.ana.TestParticleFilter.protonEMin = cms.double(100000)
process.ana.TestParticleFilter.etaMax = cms.double(3.1)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("hgcalNtuple.root")

                                   )

reRunClustering = False

if reRunClustering:
    #process.hgcalLayerClusters.minClusters = cms.uint32(0)
    #process.hgcalLayerClusters.realSpaceCone = cms.bool(True)
    #process.hgcalLayerClusters.multiclusterRadius = cms.double(2.)  # in cm if realSpaceCone is true
    #process.hgcalLayerClusters.dependSensor = cms.bool(True)
    #process.hgcalLayerClusters.ecut = cms.double(3.)  # multiple of sigma noise if dependSensor is true
    #process.hgcalLayerClusters.kappa = cms.double(9.)  # multiple of sigma noise if dependSensor is true
    #process.hgcalLayerClusters.deltac = cms.vdouble(2.,3.,5.) #specify delta c for each subdetector separately
    process.p = cms.Path(process.hgcalLayerClusters+process.ana)
else:
    process.p = cms.Path(process.prunedGenParticles+process.ana)
