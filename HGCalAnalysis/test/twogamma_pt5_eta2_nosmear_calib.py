import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")
process.load('Configuration.Geometry.GeometryExtended2023D3Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('RecoLocalCalo.HGCalRecHitDump.imagingClusterHGCal_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/partGun_clange_PDGid211_nPart1_Pt20_20160901/RECO/partGun_PDGid211_x40_Pt20.0To20.0_RECO_1.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.ana = cms.EDAnalyzer('HGCalAnalysis',
                             detector = cms.string("both"),
                             depthClusteringCone = cms.double(0.015),
                             rawRecHits = cms.bool(True)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("twogamma_pt5_eta2_nosmear_calib_ntuple.root")

                                   )
process.imagingClusterHGCal.ecut = cms.double(0.01)
process.imagingClusterHGCal.eventsToDisplay = cms.untracked.uint32(2)

process.p = cms.Path(process.imagingClusterHGCal+process.ana)
