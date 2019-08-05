// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Math/interface/deltaPhi.h"
// track data formats
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/GeometrySurface/interface/PlaneBuilder.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FastSimulation/CaloGeometryTools/interface/Transform3DPJ.h"
#include "FastSimulation/Event/interface/FSimEvent.h"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "FastSimulation/Event/interface/FSimVertex.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeGeometry/interface/MagVolumeOutsideValidity.h"
#include "TrackPropagation/RungeKutta/interface/defaultRKPropagator.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TH1F.h"
#include "TTree.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace HGCal_helpers {

class coordinates {
 public:
  coordinates() : x(0), y(0), z(0), eta(100), phi(0) {}
  float x, y, z, eta, phi;
  inline math::XYZTLorentzVectorD toVector() { return math::XYZTLorentzVectorD(x, y, z, 0); }
};
class simpleTrackPropagator {
 public:
  simpleTrackPropagator(MagneticField const *f)
      : field_(f), prod_(field_, alongMomentum, 5.e-5), absz_target_(0) {
    ROOT::Math::SMatrixIdentity id;
    AlgebraicSymMatrix55 C(id);
    C *= 0.001;
    err_ = CurvilinearTrajectoryError(C);
  }
  void setPropagationTargetZ(const float &z);

  bool propagate(const double px, const double py, const double pz, const double x, const double y,
                 const double z, const float charge, coordinates &coords) const;

  bool propagate(const math::XYZTLorentzVectorD &momentum, const math::XYZTLorentzVectorD &position,
                 const float charge, coordinates &coords) const;

 private:
  simpleTrackPropagator() : field_(0), prod_(field_, alongMomentum, 5.e-5), absz_target_(0) {}
  const RKPropagatorInS &RKProp() const { return prod_.propagator; }
  Plane::PlanePointer targetPlaneForward_, targetPlaneBackward_;
  MagneticField const *field_;
  CurvilinearTrajectoryError err_;
  defaultRKPropagator::Product prod_;
  float absz_target_;
};

void simpleTrackPropagator::setPropagationTargetZ(const float &z) {
  targetPlaneForward_ = Plane::build(Plane::PositionType(0, 0, std::abs(z)), Plane::RotationType());
  targetPlaneBackward_ =
      Plane::build(Plane::PositionType(0, 0, -std::abs(z)), Plane::RotationType());
  absz_target_ = std::abs(z);
}
bool simpleTrackPropagator::propagate(const double px, const double py, const double pz,
                                      const double x, const double y, const double z,
                                      const float charge, coordinates &output) const {
  output = coordinates();

  typedef TrajectoryStateOnSurface TSOS;
  GlobalPoint startingPosition(x, y, z);
  GlobalVector startingMomentum(px, py, pz);
  Plane::PlanePointer startingPlane =
      Plane::build(Plane::PositionType(x, y, z), Plane::RotationType());
  TSOS startingStateP(
      GlobalTrajectoryParameters(startingPosition, startingMomentum, charge, field_), err_,
      *startingPlane);

  TSOS trackStateP;
  if (pz > 0) {
    trackStateP = RKProp().propagate(startingStateP, *targetPlaneForward_);
  } else {
    trackStateP = RKProp().propagate(startingStateP, *targetPlaneBackward_);
  }
  if (trackStateP.isValid()) {
    output.x = trackStateP.globalPosition().x();
    output.y = trackStateP.globalPosition().y();
    output.z = trackStateP.globalPosition().z();
    output.phi = trackStateP.globalPosition().phi();
    output.eta = trackStateP.globalPosition().eta();
    return true;
  }
  return false;
}

bool simpleTrackPropagator::propagate(const math::XYZTLorentzVectorD &momentum,
                                      const math::XYZTLorentzVectorD &position, const float charge,
                                      coordinates &output) const {
  return propagate(momentum.px(), momentum.py(), momentum.pz(), position.x(), position.y(),
                   position.z(), charge, output);
}

}  // HGCal_helpers

class HGCalAnalysis : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
 public:
  //
  // constructors and destructor
  //
  typedef ROOT::Math::Transform3DPJ::Point Point;

  // approximative geometrical values
  static constexpr float hgcalOuterRadius_ = 160.;
  static constexpr float hgcalInnerRadius_ = 25.;

  HGCalAnalysis();
  explicit HGCalAnalysis(const edm::ParameterSet &);
  ~HGCalAnalysis();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  virtual void beginRun(edm::Run const &iEvent, edm::EventSetup const &) override;
  virtual void endRun(edm::Run const &iEvent, edm::EventSetup const &) override;

 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endJob() override;
  virtual int fillLayerCluster(const edm::Ptr<reco::CaloCluster> &layerCluster,
                               const bool &fillRecHits, const int &multiClusterIndex = -1);
  virtual void fillRecHit(const DetId &detid, const float &fraction, const unsigned int &layer,
                          const int &cluster_index_ = -1);

  void clearVariables();

  void retrieveLayerPositions(const edm::EventSetup &, unsigned layers);

  // ---------parameters ----------------------------
  bool readCaloParticles_;
  bool readGen_;
  bool storeMoreGenInfo_;
  bool storeGenParticleExtrapolation_;
  bool storeElectrons_;
  bool storePFCandidates_;
  bool storeGunParticles_;
  double layerClusterPtThreshold_;
  double propagationPtThreshold_;
  std::string detector_;
  std::string inputTag_HGCalMultiCluster_;
  bool rawRecHits_;
  bool verbose_;

  // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> recHitsEE_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsFH_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsBH_;
  edm::EDGetTokenT<reco::CaloClusterCollection> clusters_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticles_;
  edm::EDGetTokenT<std::vector<SimCluster>> simClusters_;
  edm::EDGetTokenT<std::vector<reco::PFCluster>> pfClusters_;
  edm::EDGetTokenT<std::vector<reco::PFCluster>> pfClustersFromMultiCl_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticles_;
  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster>> multiClusters_;
  edm::EDGetTokenT<std::vector<SimTrack>> simTracks_;
  edm::EDGetTokenT<std::vector<SimVertex>> simVertices_;
  edm::EDGetTokenT<edm::HepMCProduct> hev_;
  edm::EDGetTokenT<std::vector<reco::Track>> tracks_;
  edm::EDGetTokenT<std::vector<reco::GsfElectron>> electrons_;
  edm::EDGetTokenT<edm::ValueMap<reco::CaloClusterPtr>> electrons_ValueMapClusters_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_;
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandidates_;

  TTree *t_;

  ////////////////////
  // event
  //
  edm::RunNumber_t ev_run_;
  edm::LuminosityBlockNumber_t ev_lumi_;
  edm::EventNumber_t ev_event_;
  float vtx_x_;
  float vtx_y_;
  float vtx_z_;

  ////////////////////
  // GenParticles
  //
  std::vector<float> genpart_eta_;
  std::vector<float> genpart_phi_;
  std::vector<float> genpart_pt_;
  std::vector<float> genpart_energy_;
  std::vector<float> genpart_dvx_;
  std::vector<float> genpart_dvy_;
  std::vector<float> genpart_dvz_;
  std::vector<float> genpart_ovx_;
  std::vector<float> genpart_ovy_;
  std::vector<float> genpart_ovz_;
  std::vector<float> genpart_exx_;
  std::vector<float> genpart_exy_;
  std::vector<int> genpart_mother_;
  std::vector<float> genpart_exphi_;
  std::vector<float> genpart_exeta_;
  std::vector<float> genpart_fbrem_;
  std::vector<int> genpart_pid_;
  std::vector<int> genpart_gen_;
  std::vector<int> genpart_reachedEE_;
  std::vector<bool> genpart_fromBeamPipe_;
  std::vector<std::vector<float>> genpart_posx_;
  std::vector<std::vector<float>> genpart_posy_;
  std::vector<std::vector<float>> genpart_posz_;

  ////////////////////
  // reco::GenParticles
  //
  std::vector<float> gen_eta_;
  std::vector<float> gen_phi_;
  std::vector<float> gen_pt_;
  std::vector<float> gen_energy_;
  std::vector<int> gen_charge_;
  std::vector<int> gen_pdgid_;
  std::vector<int> gen_status_;
  std::vector<std::vector<int>> gen_daughters_;

  ////////////////////
  // RecHits
  // associated to layer clusters
  std::vector<float> rechit_eta_;
  std::vector<float> rechit_phi_;
  std::vector<float> rechit_pt_;
  std::vector<float> rechit_energy_;
  std::vector<float> rechit_x_;
  std::vector<float> rechit_y_;
  std::vector<float> rechit_z_;
  std::vector<float> rechit_time_;
  std::vector<float> rechit_thickness_;
  std::vector<int> rechit_layer_;
  std::vector<int> rechit_wafer_u_;
  std::vector<int> rechit_wafer_v_;
  std::vector<int> rechit_cell_u_;
  std::vector<int> rechit_cell_v_;
  std::vector<unsigned int> rechit_detid_;
  std::vector<bool> rechit_isHalf_;
  std::vector<int> rechit_flags_;
  std::vector<int> rechit_cluster2d_;
  std::vector<float> rechit_radius_;

  ////////////////////
  // layer clusters
  //
  std::vector<float> cluster2d_eta_;
  std::vector<float> cluster2d_phi_;
  std::vector<float> cluster2d_pt_;
  std::vector<float> cluster2d_energy_;
  std::vector<float> cluster2d_x_;
  std::vector<float> cluster2d_y_;
  std::vector<float> cluster2d_z_;
  std::vector<int> cluster2d_layer_;
  std::vector<int> cluster2d_nhitCore_;
  std::vector<int> cluster2d_nhitAll_;
  std::vector<int> cluster2d_multicluster_;
  std::vector<std::vector<unsigned int>> cluster2d_rechits_;
  std::vector<int> cluster2d_rechitSeed_;

  ////////////////////
  // multi clusters
  //
  std::vector<float> multiclus_eta_;
  std::vector<float> multiclus_phi_;
  std::vector<float> multiclus_pt_;
  std::vector<float> multiclus_energy_;
  std::vector<float> multiclus_z_;
  std::vector<float> multiclus_slopeX_;
  std::vector<float> multiclus_slopeY_;
  std::vector<std::vector<unsigned int>> multiclus_cluster2d_;
  std::vector<int> multiclus_cl2dSeed_;
  std::vector<int> multiclus_firstLay_;
  std::vector<int> multiclus_lastLay_;
  std::vector<int> multiclus_NLay_;

  ////////////////////
  // sim clusters
  //
  std::vector<float> simcluster_eta_;
  std::vector<float> simcluster_phi_;
  std::vector<float> simcluster_pt_;
  std::vector<float> simcluster_energy_;
  std::vector<float> simcluster_simEnergy_;
  std::vector<std::vector<uint32_t>> simcluster_hits_;
  std::vector<std::vector<int>> simcluster_hits_indices_;
  std::vector<std::vector<float>> simcluster_fractions_;
  std::vector<std::vector<unsigned int>> simcluster_layers_;
  std::vector<std::vector<int>> simcluster_wafers_u_;
  std::vector<std::vector<int>> simcluster_wafers_v_;
  std::vector<std::vector<int>> simcluster_cells_u_;
  std::vector<std::vector<int>> simcluster_cells_v_;

  ////////////////////
  // PF clusters
  //
  std::vector<float> pfcluster_eta_;
  std::vector<float> pfcluster_phi_;
  std::vector<float> pfcluster_pt_;
  std::vector<float> pfcluster_energy_;
  std::vector<float> pfcluster_correctedEnergy_;
  std::vector<std::vector<uint32_t>> pfcluster_hits_;
  std::vector<std::vector<float>> pfcluster_fractions_;

  ////////////////////
  // PF clusters From MultiClusters
  //
  std::vector<ROOT::Math::XYZPoint> pfclusterFromMultiCl_pos_;
  std::vector<float> pfclusterFromMultiCl_eta_;
  std::vector<float> pfclusterFromMultiCl_phi_;
  std::vector<float> pfclusterFromMultiCl_pt_;
  std::vector<float> pfclusterFromMultiCl_energy_;
  std::vector<float> pfclusterFromMultiCl_energyEE_;
  std::vector<float> pfclusterFromMultiCl_energyFH_;
  std::vector<float> pfclusterFromMultiCl_energyBH_;
  std::vector<float> pfclusterFromMultiCl_correctedEnergy_;
  std::vector<std::vector<uint32_t>> pfclusterFromMultiCl_hits_;
  std::vector<std::vector<float>> pfclusterFromMultiCl_fractions_;
  std::vector<std::vector<uint32_t>> pfclusterFromMultiCl_rechits_;

  ////////////////////
  // Ecal Driven GsfElectrons From MultiClusters
  //
  std::vector<float> ecalDrivenGsfele_charge_;
  std::vector<float> ecalDrivenGsfele_eta_;
  std::vector<float> ecalDrivenGsfele_phi_;
  std::vector<float> ecalDrivenGsfele_pt_;
  std::vector<ROOT::Math::XYZPoint> ecalDrivenGsfele_scpos_;
  std::vector<float> ecalDrivenGsfele_sceta_;
  std::vector<float> ecalDrivenGsfele_scphi_;
  std::vector<uint32_t> ecalDrivenGsfele_seedlayer_;
  std::vector<ROOT::Math::XYZPoint> ecalDrivenGsfele_seedpos_;
  std::vector<float> ecalDrivenGsfele_seedeta_;
  std::vector<float> ecalDrivenGsfele_seedphi_;
  std::vector<float> ecalDrivenGsfele_seedenergy_;
  std::vector<float> ecalDrivenGsfele_energy_;
  std::vector<float> ecalDrivenGsfele_energyEE_;
  std::vector<float> ecalDrivenGsfele_energyFH_;
  std::vector<float> ecalDrivenGsfele_energyBH_;
  std::vector<float> ecalDrivenGsfele_isEB_;
  std::vector<float> ecalDrivenGsfele_hoe_;
  std::vector<float> ecalDrivenGsfele_numClinSC_;
  std::vector<float> ecalDrivenGsfele_track_dxy_;
  std::vector<float> ecalDrivenGsfele_track_dz_;
  std::vector<float> ecalDrivenGsfele_track_simdxy_;
  std::vector<float> ecalDrivenGsfele_track_simdz_;
  std::vector<float> ecalDrivenGsfele_deltaEtaSuperClusterTrackAtVtx_;
  std::vector<float> ecalDrivenGsfele_deltaPhiSuperClusterTrackAtVtx_;
  std::vector<float> ecalDrivenGsfele_deltaEtaEleClusterTrackAtCalo_;
  std::vector<float> ecalDrivenGsfele_deltaPhiEleClusterTrackAtCalo_;
  std::vector<float> ecalDrivenGsfele_deltaEtaSeedClusterTrackAtCalo_;
  std::vector<float> ecalDrivenGsfele_deltaPhiSeedClusterTrackAtCalo_;
  std::vector<float> ecalDrivenGsfele_eSuperClusterOverP_;
  std::vector<float> ecalDrivenGsfele_eSeedClusterOverP_;
  std::vector<float> ecalDrivenGsfele_eSeedClusterOverPout_;
  std::vector<float> ecalDrivenGsfele_eEleClusterOverPout_;
  std::vector<std::vector<uint32_t>>
      ecalDrivenGsfele_pfClusterIndex_;  // the second index runs through the corresponding
                                         // PFClustersHGCalFromMultiClusters

  ////////////////////
  // calo particles
  //
  std::vector<float> calopart_eta_;
  std::vector<float> calopart_phi_;
  std::vector<float> calopart_pt_;
  std::vector<float> calopart_energy_;
  std::vector<float> calopart_simEnergy_;
  std::vector<std::vector<uint32_t>> calopart_simClusterIndex_;

  ////////////////////
  // high purity tracks
  //
  std::vector<float> track_eta_;
  std::vector<float> track_phi_;
  std::vector<float> track_pt_;
  std::vector<float> track_energy_;
  std::vector<int> track_charge_;
  std::vector<std::vector<float>> track_posx_;
  std::vector<std::vector<float>> track_posy_;
  std::vector<std::vector<float>> track_posz_;


  ////////////////////
  // PF candidates
  //
  std::vector<float> pfcandidate_eta_;
  std::vector<float> pfcandidate_phi_;
  std::vector<float> pfcandidate_pt_;
  std::vector<float> pfcandidate_energy_;
  std::vector<int> pfcandidate_pdgid_;

  ////////////////////
  // gun particles per vertex
  //
  std::vector<std::vector<int> > gunparticle_id_;
  std::vector<std::vector<float> > gunparticle_energy_;
  std::vector<std::vector<float> > gunparticle_pt_;
  std::vector<std::vector<float> > gunparticle_eta_;
  std::vector<std::vector<float> > gunparticle_phi_;


  ////////////////////
  // helper classes
  //
  unsigned int cluster_index_;
  unsigned int rechit_index_;
  std::map<DetId, const HGCRecHit *> hitmap_;
  std::map<DetId, unsigned int> detIdToRecHitIndexMap_;
  float vz_;  // primary vertex z position
  // to keep track of the 2d clusters stored within the loop on multiclusters
  std::set<edm::Ptr<reco::BasicCluster>> storedLayerClusters_;
  // to keep track of the RecHits stored within the cluster loops
  std::set<DetId> storedRecHits_;
  int algo_;
  HGCalDepthPreClusterer pre_;
  hgcal::RecHitTools recHitTools_;

  // -------convenient tool to deal with simulated tracks
  FSimEvent *mySimEvent_;
  edm::ParameterSet particleFilter_;
  std::vector<float> layerPositions_;


  // and also the magnetic field
  MagneticField const *aField_;

};

HGCalAnalysis::HGCalAnalysis() { ; }

HGCalAnalysis::HGCalAnalysis(const edm::ParameterSet &iConfig)
    : readCaloParticles_(iConfig.getParameter<bool>("readCaloParticles")),
      readGen_(iConfig.getParameter<bool>("readGenParticles")),
      storeMoreGenInfo_(iConfig.getParameter<bool>("storeGenParticleOrigin")),
      storeGenParticleExtrapolation_(iConfig.getParameter<bool>("storeGenParticleExtrapolation")),
      storeElectrons_(iConfig.getParameter<bool>("storeElectrons")),
      storePFCandidates_(iConfig.getParameter<bool>("storePFCandidates")),
      storeGunParticles_(iConfig.getParameter<bool>("storeGunParticles")),
      layerClusterPtThreshold_(iConfig.getParameter<double>("layerClusterPtThreshold")),
      propagationPtThreshold_(iConfig.getUntrackedParameter<double>("propagationPtThreshold", 3.0)),
      detector_(iConfig.getParameter<std::string>("detector")),
      inputTag_HGCalMultiCluster_(iConfig.getParameter<std::string>("inputTag_HGCalMultiCluster")),
      rawRecHits_(iConfig.getParameter<bool>("rawRecHits")),
      verbose_(iConfig.getParameter<bool>("verbose")),
      particleFilter_(iConfig.getParameter<edm::ParameterSet>("TestParticleFilter")) {
  // now do what ever initialization is needed
  mySimEvent_ = new FSimEvent(particleFilter_);

  if (detector_ == "all") {
    recHitsEE_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCEERecHits"));
    recHitsFH_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
    recHitsBH_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
    algo_ = 1;
  } else if (detector_ == "EM") {
    recHitsEE_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCEERecHits"));
    algo_ = 2;
  } else if (detector_ == "HAD") {
    recHitsFH_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
    recHitsBH_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
    algo_ = 3;
  }
  clusters_ = consumes<reco::CaloClusterCollection>(edm::InputTag("hgcalLayerClusters"));
  simClusters_ = consumes<std::vector<SimCluster>>(edm::InputTag("mix", "MergedCaloTruth"));
  hev_ = consumes<edm::HepMCProduct>(edm::InputTag("generatorSmeared"));

  simTracks_ = consumes<std::vector<SimTrack>>(edm::InputTag("g4SimHits"));
  simVertices_ = consumes<std::vector<SimVertex>>(edm::InputTag("g4SimHits"));

  if (readCaloParticles_) {
    caloParticles_ = consumes<std::vector<CaloParticle>>(edm::InputTag("mix", "MergedCaloTruth"));
  }

  if (readGen_) {
    genParticles_ = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"));
  }

  if (storePFCandidates_) {
    pfCandidates_ = consumes<std::vector<reco::PFCandidate>>(edm::InputTag("particleFlow"));
  }

  pfClusters_ = consumes<std::vector<reco::PFCluster>>(edm::InputTag("particleFlowClusterHGCal"));
  pfClustersFromMultiCl_ =
      consumes<std::vector<reco::PFCluster>>(edm::InputTag("particleFlowClusterHGCalFromMultiCl"));
  multiClusters_ =
      consumes<std::vector<reco::HGCalMultiCluster>>(edm::InputTag(inputTag_HGCalMultiCluster_));
  tracks_ = consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks"));
  electrons_ =
      consumes<std::vector<reco::GsfElectron>>(edm::InputTag("ecalDrivenGsfElectronsFromMultiCl"));
  electrons_ValueMapClusters_ = consumes<edm::ValueMap<reco::CaloClusterPtr>>(
      edm::InputTag("particleFlowSuperClusterHGCalFromMultiCl:PFClusterAssociationEBEE"));
  vertices_ = consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"));

  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs;
  fs->make<TH1F>("total", "total", 100, 0, 5.);

  t_ = fs->make<TTree>("hgc", "hgc");

  // event info
  t_->Branch("event", &ev_event_);
  t_->Branch("lumi", &ev_lumi_);
  t_->Branch("run", &ev_run_);
  t_->Branch("vtx_x", &vtx_x_);
  t_->Branch("vtx_y", &vtx_y_);
  t_->Branch("vtx_z", &vtx_z_);

  t_->Branch("genpart_eta", &genpart_eta_);
  t_->Branch("genpart_phi", &genpart_phi_);
  t_->Branch("genpart_pt", &genpart_pt_);
  t_->Branch("genpart_energy", &genpart_energy_);
  t_->Branch("genpart_dvx", &genpart_dvx_);
  t_->Branch("genpart_dvy", &genpart_dvy_);
  t_->Branch("genpart_dvz", &genpart_dvz_);
  if (storeMoreGenInfo_) {
    t_->Branch("genpart_ovx", &genpart_ovx_);
    t_->Branch("genpart_ovy", &genpart_ovy_);
    t_->Branch("genpart_ovz", &genpart_ovz_);
    t_->Branch("genpart_mother", &genpart_mother_);
  }
  if (storeGenParticleExtrapolation_) {
    t_->Branch("genpart_exphi", &genpart_exphi_);
    t_->Branch("genpart_exeta", &genpart_exeta_);
    t_->Branch("genpart_exx", &genpart_exx_);
    t_->Branch("genpart_exy", &genpart_exy_);
  }
  t_->Branch("genpart_fbrem", &genpart_fbrem_);
  t_->Branch("genpart_pid", &genpart_pid_);
  t_->Branch("genpart_gen", &genpart_gen_);
  t_->Branch("genpart_reachedEE", &genpart_reachedEE_);
  t_->Branch("genpart_fromBeamPipe", &genpart_fromBeamPipe_);
  t_->Branch("genpart_posx", &genpart_posx_);
  t_->Branch("genpart_posy", &genpart_posy_);
  t_->Branch("genpart_posz", &genpart_posz_);


  if (readGen_) {
    t_->Branch("gen_eta", &gen_eta_);
    t_->Branch("gen_phi", &gen_phi_);
    t_->Branch("gen_pt", &gen_pt_);
    t_->Branch("gen_energy", &gen_energy_);
    t_->Branch("gen_charge", &gen_charge_);
    t_->Branch("gen_pdgid", &gen_pdgid_);
    t_->Branch("gen_status", &gen_status_);
    t_->Branch("gen_daughters", &gen_daughters_);
  }

  ////////////////////
  // RecHits
  // associated to layer clusters
  t_->Branch("rechit_eta", &rechit_eta_);
  t_->Branch("rechit_phi", &rechit_phi_);
  t_->Branch("rechit_pt", &rechit_pt_);
  t_->Branch("rechit_energy", &rechit_energy_);
  t_->Branch("rechit_x", &rechit_x_);
  t_->Branch("rechit_y", &rechit_y_);
  t_->Branch("rechit_z", &rechit_z_);
  t_->Branch("rechit_time", &rechit_time_);
  t_->Branch("rechit_thickness", &rechit_thickness_);
  t_->Branch("rechit_layer", &rechit_layer_);
  t_->Branch("rechit_wafer_u", &rechit_wafer_u_);
  t_->Branch("rechit_wafer_v", &rechit_wafer_v_);
  t_->Branch("rechit_cell_u", &rechit_cell_u_);
  t_->Branch("rechit_cell_v", &rechit_cell_v_);
  t_->Branch("rechit_detid", &rechit_detid_);
  t_->Branch("rechit_isHalf", &rechit_isHalf_);
  t_->Branch("rechit_flags", &rechit_flags_);
  t_->Branch("rechit_cluster2d", &rechit_cluster2d_);
  t_->Branch("rechit_radius", &rechit_radius_);

  ////////////////////
  // layer clusters
  //
  t_->Branch("cluster2d_eta", &cluster2d_eta_);
  t_->Branch("cluster2d_phi", &cluster2d_phi_);
  t_->Branch("cluster2d_pt", &cluster2d_pt_);
  t_->Branch("cluster2d_energy", &cluster2d_energy_);
  t_->Branch("cluster2d_x", &cluster2d_x_);
  t_->Branch("cluster2d_y", &cluster2d_y_);
  t_->Branch("cluster2d_z", &cluster2d_z_);
  t_->Branch("cluster2d_layer", &cluster2d_layer_);
  t_->Branch("cluster2d_nhitCore", &cluster2d_nhitCore_);
  t_->Branch("cluster2d_nhitAll", &cluster2d_nhitAll_);
  t_->Branch("cluster2d_multicluster", &cluster2d_multicluster_);
  t_->Branch("cluster2d_rechits", &cluster2d_rechits_);
  t_->Branch("cluster2d_rechitSeed", &cluster2d_rechitSeed_);

  ////////////////////
  // multi clusters
  //
  t_->Branch("multiclus_eta", &multiclus_eta_);
  t_->Branch("multiclus_phi", &multiclus_phi_);
  t_->Branch("multiclus_pt", &multiclus_pt_);
  t_->Branch("multiclus_energy", &multiclus_energy_);
  t_->Branch("multiclus_z", &multiclus_z_);
  t_->Branch("multiclus_slopeX", &multiclus_slopeX_);
  t_->Branch("multiclus_slopeY", &multiclus_slopeY_);
  t_->Branch("multiclus_cluster2d", &multiclus_cluster2d_);
  t_->Branch("multiclus_cl2dSeed", &multiclus_cl2dSeed_);
  t_->Branch("multiclus_firstLay", &multiclus_firstLay_);
  t_->Branch("multiclus_lastLay", &multiclus_lastLay_);
  t_->Branch("multiclus_NLay", &multiclus_NLay_);

  ////////////////////
  // sim clusters
  //
  t_->Branch("simcluster_eta", &simcluster_eta_);
  t_->Branch("simcluster_phi", &simcluster_phi_);
  t_->Branch("simcluster_pt", &simcluster_pt_);
  t_->Branch("simcluster_energy", &simcluster_energy_);
  t_->Branch("simcluster_simEnergy", &simcluster_simEnergy_);
  t_->Branch("simcluster_hits", &simcluster_hits_);
  t_->Branch("simcluster_hits_indices", &simcluster_hits_indices_);
  t_->Branch("simcluster_fractions", &simcluster_fractions_);
  t_->Branch("simcluster_layers", &simcluster_layers_);
  t_->Branch("simcluster_wafers_u", &simcluster_wafers_u_);
  t_->Branch("simcluster_wafers_v", &simcluster_wafers_v_);
  t_->Branch("simcluster_cells_u", &simcluster_cells_u_);
  t_->Branch("simcluster_cells_v", &simcluster_cells_v_);

  ////////////////////
  // PF clusters
  //
  t_->Branch("pfcluster_eta", &pfcluster_eta_);
  t_->Branch("pfcluster_phi", &pfcluster_phi_);
  t_->Branch("pfcluster_pt", &pfcluster_pt_);
  t_->Branch("pfcluster_energy", &pfcluster_energy_);
  t_->Branch("pfcluster_correctedEnergy", &pfcluster_correctedEnergy_);
  t_->Branch("pfcluster_hits", &pfcluster_hits_);
  t_->Branch("pfcluster_fractions", &pfcluster_fractions_);

  ////////////////////
  // PF clusters From MultiClusters
  //
  t_->Branch("pfclusterFromMultiCl_pos", &pfclusterFromMultiCl_pos_);
  t_->Branch("pfclusterFromMultiCl_eta", &pfclusterFromMultiCl_eta_);
  t_->Branch("pfclusterFromMultiCl_phi", &pfclusterFromMultiCl_phi_);
  t_->Branch("pfclusterFromMultiCl_pt", &pfclusterFromMultiCl_pt_);
  t_->Branch("pfclusterFromMultiCl_energy", &pfclusterFromMultiCl_energy_);
  t_->Branch("pfclusterFromMultiCl_energyEE", &pfclusterFromMultiCl_energyEE_);
  t_->Branch("pfclusterFromMultiCl_energyFH", &pfclusterFromMultiCl_energyFH_);
  t_->Branch("pfclusterFromMultiCl_energyBH", &pfclusterFromMultiCl_energyBH_);
  t_->Branch("pfclusterFromMultiCl_correctedEnergy", &pfclusterFromMultiCl_correctedEnergy_);
  t_->Branch("pfclusterFromMultiCl_hits", &pfclusterFromMultiCl_hits_);
  t_->Branch("pfclusterFromMultiCl_fractions", &pfclusterFromMultiCl_fractions_);
  t_->Branch("pfclusterFromMultiCl_rechits", &pfclusterFromMultiCl_rechits_);

  ////////////////////
  // Ecal Driven Gsf Electrons From MultiClusters
  //
  if (storeElectrons_) {
    t_->Branch("ecalDrivenGsfele_charge", &ecalDrivenGsfele_charge_);
    t_->Branch("ecalDrivenGsfele_eta", &ecalDrivenGsfele_eta_);
    t_->Branch("ecalDrivenGsfele_phi", &ecalDrivenGsfele_phi_);
    t_->Branch("ecalDrivenGsfele_pt", &ecalDrivenGsfele_pt_);
    t_->Branch("ecalDrivenGsfele_scpos", &ecalDrivenGsfele_scpos_);
    t_->Branch("ecalDrivenGsfele_sceta", &ecalDrivenGsfele_sceta_);
    t_->Branch("ecalDrivenGsfele_scphi", &ecalDrivenGsfele_scphi_);
    t_->Branch("ecalDrivenGsfele_seedlayer", &ecalDrivenGsfele_seedlayer_);
    t_->Branch("ecalDrivenGsfele_seedpos", &ecalDrivenGsfele_seedpos_);
    t_->Branch("ecalDrivenGsfele_seedeta", &ecalDrivenGsfele_seedeta_);
    t_->Branch("ecalDrivenGsfele_seedphi", &ecalDrivenGsfele_seedphi_);
    t_->Branch("ecalDrivenGsfele_seedenergy", &ecalDrivenGsfele_seedenergy_);
    t_->Branch("ecalDrivenGsfele_energy", &ecalDrivenGsfele_energy_);
    t_->Branch("ecalDrivenGsfele_energyEE", &ecalDrivenGsfele_energyEE_);
    t_->Branch("ecalDrivenGsfele_energyFH", &ecalDrivenGsfele_energyFH_);
    t_->Branch("ecalDrivenGsfele_energyBH", &ecalDrivenGsfele_energyBH_);
    t_->Branch("ecalDrivenGsfele_isEB", &ecalDrivenGsfele_isEB_);
    t_->Branch("ecalDrivenGsfele_hoe", &ecalDrivenGsfele_hoe_);
    t_->Branch("ecalDrivenGsfele_numClinSC", &ecalDrivenGsfele_numClinSC_);
    t_->Branch("ecalDrivenGsfele_track_dxy", &ecalDrivenGsfele_track_dxy_);
    t_->Branch("ecalDrivenGsfele_track_dz", &ecalDrivenGsfele_track_dz_);
    t_->Branch("ecalDrivenGsfele_track_simdxy", &ecalDrivenGsfele_track_simdxy_);
    t_->Branch("ecalDrivenGsfele_track_simdz", &ecalDrivenGsfele_track_simdz_);
    t_->Branch("ecalDrivenGsfele_deltaEtaSuperClusterTrackAtVtx",
               &ecalDrivenGsfele_deltaEtaSuperClusterTrackAtVtx_);
    t_->Branch("ecalDrivenGsfele_deltaPhiSuperClusterTrackAtVtx",
               &ecalDrivenGsfele_deltaPhiSuperClusterTrackAtVtx_);
    t_->Branch("ecalDrivenGsfele_deltaEtaEleClusterTrackAtCalo",
               &ecalDrivenGsfele_deltaEtaEleClusterTrackAtCalo_);
    t_->Branch("ecalDrivenGsfele_deltaPhiEleClusterTrackAtCalo",
               &ecalDrivenGsfele_deltaPhiEleClusterTrackAtCalo_);
    t_->Branch("ecalDrivenGsfele_deltaEtaSeedClusterTrackAtCalo",
               &ecalDrivenGsfele_deltaEtaSeedClusterTrackAtCalo_);
    t_->Branch("ecalDrivenGsfele_deltaPhiSeedClusterTrackAtCalo",
               &ecalDrivenGsfele_deltaPhiSeedClusterTrackAtCalo_);
    t_->Branch("ecalDrivenGsfele_eSuperClusterOverP", &ecalDrivenGsfele_eSuperClusterOverP_);
    t_->Branch("ecalDrivenGsfele_eSeedClusterOverP", &ecalDrivenGsfele_eSeedClusterOverP_);
    t_->Branch("ecalDrivenGsfele_eSeedClusterOverPout", &ecalDrivenGsfele_eSeedClusterOverPout_);
    t_->Branch("ecalDrivenGsfele_eEleClusterOverPout", &ecalDrivenGsfele_eEleClusterOverPout_);
    t_->Branch("ecalDrivenGsfele_pfClusterIndex", &ecalDrivenGsfele_pfClusterIndex_);
  }

  ////////////////////
  // calo particles
  //
  t_->Branch("calopart_eta", &calopart_eta_);
  t_->Branch("calopart_phi", &calopart_phi_);
  t_->Branch("calopart_pt", &calopart_pt_);
  t_->Branch("calopart_energy", &calopart_energy_);
  t_->Branch("calopart_simEnergy", &calopart_simEnergy_);
  t_->Branch("calopart_simClusterIndex", &calopart_simClusterIndex_);

  ////////////////////
  // high purity trackstatep
  //
  t_->Branch("track_eta", &track_eta_);
  t_->Branch("track_phi", &track_phi_);
  t_->Branch("track_pt", &track_pt_);
  t_->Branch("track_energy", &track_energy_);
  t_->Branch("track_charge", &track_charge_);
  t_->Branch("track_posx", &track_posx_);
  t_->Branch("track_posy", &track_posy_);
  t_->Branch("track_posz", &track_posz_);

  ////////////////////
  // PF candidates
  //
  if (storePFCandidates_) {
    t_->Branch("pfcandidate_eta", &pfcandidate_eta_);
    t_->Branch("pfcandidate_phi", &pfcandidate_phi_);
    t_->Branch("pfcandidate_pt", &pfcandidate_pt_);
    t_->Branch("pfcandidate_energy", &pfcandidate_energy_);
    t_->Branch("pfcandidate_pdgid", &pfcandidate_pdgid_);
  }

  ////////////////////
  // gun particles
  //
  if (storeGunParticles_) {
    t_->Branch("gunparticle_id", &gunparticle_id_);
    t_->Branch("gunparticle_energy", &gunparticle_energy_);
    t_->Branch("gunparticle_pt", &gunparticle_pt_);
    t_->Branch("gunparticle_eta", &gunparticle_eta_);
    t_->Branch("gunparticle_phi", &gunparticle_phi_);
  }

}
HGCalAnalysis::~HGCalAnalysis() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//
void HGCalAnalysis::clearVariables() {
  ev_run_ = 0;
  ev_lumi_ = 0;
  ev_event_ = 0;
  vtx_x_ = 0;
  vtx_y_ = 0;
  vtx_z_ = 0;

  ////////////////////
  // GenParticles
  //
  genpart_eta_.clear();
  genpart_phi_.clear();
  genpart_pt_.clear();
  genpart_energy_.clear();
  genpart_dvx_.clear();
  genpart_dvy_.clear();
  genpart_dvz_.clear();
  genpart_ovx_.clear();
  genpart_ovy_.clear();
  genpart_ovz_.clear();
  genpart_exx_.clear();
  genpart_exy_.clear();
  genpart_mother_.clear();
  genpart_exphi_.clear();
  genpart_exeta_.clear();
  genpart_fbrem_.clear();
  genpart_pid_.clear();
  genpart_gen_.clear();
  genpart_reachedEE_.clear();
  genpart_fromBeamPipe_.clear();
  genpart_posx_.clear();
  genpart_posy_.clear();
  genpart_posz_.clear();

  ////////////////////
  // reco::GenParticles
  //
  gen_eta_.clear();
  gen_phi_.clear();
  gen_pt_.clear();
  gen_energy_.clear();
  gen_charge_.clear();
  gen_pdgid_.clear();
  gen_status_.clear();
  gen_daughters_.clear();

  ////////////////////
  // RecHits
  // associated to layer clusters
  rechit_eta_.clear();
  rechit_phi_.clear();
  rechit_pt_.clear();
  rechit_energy_.clear();
  rechit_x_.clear();
  rechit_y_.clear();
  rechit_z_.clear();
  rechit_time_.clear();
  rechit_thickness_.clear();
  rechit_layer_.clear();
  rechit_wafer_u_.clear();
  rechit_wafer_v_.clear();
  rechit_cell_u_.clear();
  rechit_cell_v_.clear();
  rechit_detid_.clear();
  rechit_isHalf_.clear();
  rechit_flags_.clear();
  rechit_cluster2d_.clear();
  rechit_radius_.clear();
  detIdToRecHitIndexMap_.clear();

  ////////////////////
  // layer clusters
  //
  cluster2d_eta_.clear();
  cluster2d_phi_.clear();
  cluster2d_pt_.clear();
  cluster2d_energy_.clear();
  cluster2d_x_.clear();
  cluster2d_y_.clear();
  cluster2d_z_.clear();
  cluster2d_layer_.clear();
  cluster2d_nhitCore_.clear();
  cluster2d_nhitAll_.clear();
  cluster2d_multicluster_.clear();
  cluster2d_rechits_.clear();
  cluster2d_rechitSeed_.clear();

  ////////////////////
  // multi clusters
  //
  multiclus_eta_.clear();
  multiclus_phi_.clear();
  multiclus_pt_.clear();
  multiclus_energy_.clear();
  multiclus_z_.clear();
  multiclus_slopeX_.clear();
  multiclus_slopeY_.clear();
  multiclus_cluster2d_.clear();
  multiclus_cl2dSeed_.clear();
  multiclus_firstLay_.clear();
  multiclus_lastLay_.clear();
  multiclus_NLay_.clear();

  ////////////////////
  // sim clusters
  //
  simcluster_eta_.clear();
  simcluster_phi_.clear();
  simcluster_pt_.clear();
  simcluster_energy_.clear();
  simcluster_simEnergy_.clear();
  simcluster_hits_.clear();
  simcluster_hits_indices_.clear();
  simcluster_fractions_.clear();
  simcluster_layers_.clear();
  simcluster_wafers_u_.clear();
  simcluster_wafers_v_.clear();
  simcluster_cells_u_.clear();
  simcluster_cells_v_.clear();

  ////////////////////
  // PF clusters
  //
  pfcluster_eta_.clear();
  pfcluster_phi_.clear();
  pfcluster_pt_.clear();
  pfcluster_energy_.clear();
  pfcluster_correctedEnergy_.clear();
  pfcluster_hits_.clear();
  pfcluster_fractions_.clear();

  ////////////////////
  // PF clusters From MultiClusters
  //
  pfclusterFromMultiCl_pos_.clear();
  pfclusterFromMultiCl_eta_.clear();
  pfclusterFromMultiCl_phi_.clear();
  pfclusterFromMultiCl_pt_.clear();
  pfclusterFromMultiCl_energy_.clear();
  pfclusterFromMultiCl_energyEE_.clear();
  pfclusterFromMultiCl_energyFH_.clear();
  pfclusterFromMultiCl_energyBH_.clear();
  pfclusterFromMultiCl_correctedEnergy_.clear();
  pfclusterFromMultiCl_hits_.clear();
  pfclusterFromMultiCl_fractions_.clear();
  pfclusterFromMultiCl_rechits_.clear();

  ////////////////////
  //  Ecal Driven Gsf Electrons From MultiClusters
  //
  ecalDrivenGsfele_charge_.clear();
  ecalDrivenGsfele_eta_.clear();
  ecalDrivenGsfele_phi_.clear();
  ecalDrivenGsfele_pt_.clear();
  ecalDrivenGsfele_scpos_.clear();
  ecalDrivenGsfele_sceta_.clear();
  ecalDrivenGsfele_scphi_.clear();
  ecalDrivenGsfele_seedlayer_.clear();
  ecalDrivenGsfele_seedpos_.clear();
  ecalDrivenGsfele_seedeta_.clear();
  ecalDrivenGsfele_seedphi_.clear();
  ecalDrivenGsfele_seedenergy_.clear();
  ecalDrivenGsfele_energy_.clear();
  ecalDrivenGsfele_energyEE_.clear();
  ecalDrivenGsfele_energyFH_.clear();
  ecalDrivenGsfele_energyBH_.clear();
  ecalDrivenGsfele_isEB_.clear();
  ecalDrivenGsfele_hoe_.clear();
  ecalDrivenGsfele_numClinSC_.clear();
  ecalDrivenGsfele_track_dxy_.clear();
  ecalDrivenGsfele_track_dz_.clear();
  ecalDrivenGsfele_track_simdxy_.clear();
  ecalDrivenGsfele_track_simdz_.clear();
  ecalDrivenGsfele_deltaEtaSuperClusterTrackAtVtx_.clear();
  ecalDrivenGsfele_deltaPhiSuperClusterTrackAtVtx_.clear();
  ecalDrivenGsfele_deltaEtaEleClusterTrackAtCalo_.clear();
  ecalDrivenGsfele_deltaPhiEleClusterTrackAtCalo_.clear();
  ecalDrivenGsfele_deltaEtaSeedClusterTrackAtCalo_.clear();
  ecalDrivenGsfele_deltaPhiSeedClusterTrackAtCalo_.clear();
  ecalDrivenGsfele_eSuperClusterOverP_.clear();
  ecalDrivenGsfele_eSeedClusterOverP_.clear();
  ecalDrivenGsfele_eSeedClusterOverPout_.clear();
  ecalDrivenGsfele_eEleClusterOverPout_.clear();
  ecalDrivenGsfele_pfClusterIndex_.clear();

  ////////////////////
  // calo particles
  //
  calopart_eta_.clear();
  calopart_phi_.clear();
  calopart_pt_.clear();
  calopart_energy_.clear();
  calopart_simEnergy_.clear();
  calopart_simClusterIndex_.clear();

  ////////////////////
  // high purity tracks
  //
  track_eta_.clear();
  track_phi_.clear();
  track_pt_.clear();
  track_energy_.clear();
  track_charge_.clear();
  track_posx_.clear();
  track_posy_.clear();
  track_posz_.clear();

  ////////////////////
  // PF candidates
  //
  pfcandidate_eta_.clear();
  pfcandidate_phi_.clear();
  pfcandidate_pt_.clear();
  pfcandidate_energy_.clear();
  pfcandidate_pdgid_.clear();

  ////////////////////
  // gun particles
  //
  gunparticle_id_.clear();
  gunparticle_energy_.clear();
  gunparticle_pt_.clear();
  gunparticle_eta_.clear();
  gunparticle_phi_.clear();
}

void HGCalAnalysis::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using namespace edm;

  clearVariables();

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleFH;
  Handle<HGCRecHitCollection> recHitHandleBH;
  Handle<reco::CaloClusterCollection> clusterHandle;

  iEvent.getByToken(clusters_, clusterHandle);

  Handle<edm::HepMCProduct> hevH;
  Handle<std::vector<SimTrack>> simTracksHandle;
  Handle<std::vector<SimVertex>> simVerticesHandle;

  iEvent.getByToken(hev_, hevH);

  iEvent.getByToken(simTracks_, simTracksHandle);
  iEvent.getByToken(simVertices_, simVerticesHandle);
  mySimEvent_->fill(*simTracksHandle, *simVerticesHandle);

  Handle<std::vector<SimCluster>> simClusterHandle;
  Handle<std::vector<reco::PFCluster>> pfClusterHandle;
  Handle<std::vector<reco::PFCluster>> pfClusterFromMultiClHandle;
  Handle<std::vector<CaloParticle>> caloParticleHandle;
  Handle<std::vector<reco::Track>> trackHandle;

  const std::vector<SimCluster> *simClusters = 0;
  iEvent.getByToken(simClusters_, simClusterHandle);
  simClusters = &(*simClusterHandle);

  iEvent.getByToken(pfClusters_, pfClusterHandle);
  const std::vector<reco::PFCluster> &pfClusters = *pfClusterHandle;

  iEvent.getByToken(pfClustersFromMultiCl_, pfClusterFromMultiClHandle);
  const std::vector<reco::PFCluster> &pfClustersFromMultiCl = *pfClusterFromMultiClHandle;

  const std::vector<CaloParticle> *caloParticles;
  if (readCaloParticles_) {
    iEvent.getByToken(caloParticles_, caloParticleHandle);
    caloParticles = &(*caloParticleHandle);
  }

  iEvent.getByToken(tracks_, trackHandle);
  const std::vector<reco::Track> &tracks = *trackHandle;

  Handle<std::vector<reco::GsfElectron>> eleHandle;
  iEvent.getByToken(electrons_, eleHandle);
  const std::vector<reco::GsfElectron> &electrons = *eleHandle;

  edm::Handle<edm::ValueMap<reco::CaloClusterPtr>> electrons_ValueMapClustersHandle;
  iEvent.getByToken(electrons_ValueMapClusters_, electrons_ValueMapClustersHandle);
  auto const &electrons_ValueMapClusters = *electrons_ValueMapClustersHandle;

  Handle<std::vector<reco::HGCalMultiCluster>> multiClusterHandle;
  iEvent.getByToken(multiClusters_, multiClusterHandle);
  const std::vector<reco::HGCalMultiCluster> &multiClusters = *multiClusterHandle;

  Handle<std::vector<reco::Vertex>> verticesHandle;
  iEvent.getByToken(vertices_, verticesHandle);
  auto const &vertices = *verticesHandle;

  HepMC::GenVertex *primaryVertex = *(hevH)->GetEvent()->vertices_begin();
  float vx_ = primaryVertex->position().x() / 10.;  // to put in official units
  float vy_ = primaryVertex->position().y() / 10.;
  vz_ = primaryVertex->position().z() / 10.;
  Point sim_pv(vx_, vy_, vz_);
  // std::cout << "start the fun" << std::endl;

  // fill the gunparticles per vertex
  if (storeGunParticles_) {
    HepMC::GenEvent::vertex_const_iterator vertex_it;
    for (vertex_it = hevH->GetEvent()->vertices_begin();
         vertex_it != hevH->GetEvent()->vertices_end(); vertex_it++) {
      std::vector<int> gunparticle_id;
      std::vector<float> gunparticle_energy;
      std::vector<float> gunparticle_pt;
      std::vector<float> gunparticle_eta;
      std::vector<float> gunparticle_phi;

      HepMC::GenVertex::particles_out_const_iterator particle_it;
      for (particle_it = (*vertex_it)->particles_out_const_begin();
           particle_it != (*vertex_it)->particles_out_const_end(); particle_it++) {
        gunparticle_id.push_back((*particle_it)->pdg_id());
        gunparticle_energy.push_back((*particle_it)->momentum().e());
        gunparticle_pt.push_back((*particle_it)->momentum().perp());
        gunparticle_eta.push_back((*particle_it)->momentum().eta());
        gunparticle_phi.push_back((*particle_it)->momentum().phi());
      }

      gunparticle_id_.push_back(gunparticle_id);
      gunparticle_energy_.push_back(gunparticle_energy);
      gunparticle_pt_.push_back(gunparticle_pt);
      gunparticle_eta_.push_back(gunparticle_eta);
      gunparticle_phi_.push_back(gunparticle_phi);
    }
  }

  HGCal_helpers::simpleTrackPropagator toHGCalPropagator(aField_);
  toHGCalPropagator.setPropagationTargetZ(layerPositions_[0]);
  std::vector<FSimTrack *> allselectedgentracks;
  unsigned int npart = mySimEvent_->nTracks();
  for (unsigned int i = 0; i < npart; ++i) {
    std::vector<float> xp, yp, zp;
    FSimTrack &myTrack(mySimEvent_->track(i));
    math::XYZTLorentzVectorD vtx(0, 0, 0, 0);

    int reachedEE = 0;  // compute the extrapolations for the particles reaching EE
                        // and for the gen particles
    double fbrem = -1;

    if (std::abs(myTrack.vertex().position().z()) >= layerPositions_[0]) continue;

    unsigned nlayers =   recHitTools_.lastLayerFH();
    if (myTrack.noEndVertex())  // || myTrack.genpartIndex()>=0)
    {
      HGCal_helpers::coordinates propcoords;
      bool reachesHGCal = toHGCalPropagator.propagate(
          myTrack.momentum(), myTrack.vertex().position(), myTrack.charge(), propcoords);
      vtx = propcoords.toVector();

      if (reachesHGCal && vtx.Rho() < hgcalOuterRadius_ && vtx.Rho() > hgcalInnerRadius_) {
        reachedEE = 2;
        double dpt = 0;

        for (int i = 0; i < myTrack.nDaughters(); ++i) dpt += myTrack.daughter(i).momentum().pt();
        if (abs(myTrack.type()) == 11) fbrem = dpt / myTrack.momentum().pt();
      } else if (reachesHGCal && vtx.Rho() > hgcalOuterRadius_)
        reachedEE = 1;

      HGCal_helpers::simpleTrackPropagator indiv_particleProp(aField_);
      for (unsigned il = 0; il < nlayers; ++il) {
        const float charge = myTrack.charge();
        indiv_particleProp.setPropagationTargetZ(layerPositions_[il]);
        HGCal_helpers::coordinates propCoords;
        indiv_particleProp.propagate(myTrack.momentum(), myTrack.vertex().position(), charge,
                                     propCoords);

        xp.push_back(propCoords.x);
        yp.push_back(propCoords.y);
        zp.push_back(propCoords.z);
      }
    } else {
      vtx = myTrack.endVertex().position();
    }
    auto orig_vtx = myTrack.vertex().position();

    allselectedgentracks.push_back(&mySimEvent_->track(i));
    // fill branches
    genpart_eta_.push_back(myTrack.momentum().eta());
    genpart_phi_.push_back(myTrack.momentum().phi());
    genpart_pt_.push_back(myTrack.momentum().pt());
    genpart_energy_.push_back(myTrack.momentum().energy());
    genpart_dvx_.push_back(vtx.x());
    genpart_dvy_.push_back(vtx.y());
    genpart_dvz_.push_back(vtx.z());

    genpart_ovx_.push_back(orig_vtx.x());
    genpart_ovy_.push_back(orig_vtx.y());
    genpart_ovz_.push_back(orig_vtx.z());

    HGCal_helpers::coordinates hitsHGCal;
    toHGCalPropagator.propagate(myTrack.momentum(), orig_vtx, myTrack.charge(), hitsHGCal);

    genpart_exphi_.push_back(hitsHGCal.phi);
    genpart_exeta_.push_back(hitsHGCal.eta);
    genpart_exx_.push_back(hitsHGCal.x);
    genpart_exy_.push_back(hitsHGCal.y);

    genpart_fbrem_.push_back(fbrem);
    genpart_pid_.push_back(myTrack.type());
    genpart_gen_.push_back(myTrack.genpartIndex());
    genpart_reachedEE_.push_back(reachedEE);
    genpart_fromBeamPipe_.push_back(true);

    genpart_posx_.push_back(xp);
    genpart_posy_.push_back(yp);
    genpart_posz_.push_back(zp);
  }

  // associate gen particles to mothers
  genpart_mother_.resize(genpart_posz_.size(), -1);
  for (size_t i = 0; i < allselectedgentracks.size(); i++) {
    const auto tracki = allselectedgentracks.at(i);

    for (size_t j = i + 1; j < allselectedgentracks.size(); j++) {
      const auto trackj = allselectedgentracks.at(j);

      if (!tracki->noMother()) {
        if (&tracki->mother() == trackj) genpart_mother_.at(i) = j;
      }
      if (!trackj->noMother()) {
        if (&trackj->mother() == tracki) genpart_mother_.at(j) = i;
      }
    }
  }

  // make a map detid-rechit
  hitmap_.clear();
  switch (algo_) {
    case 1: {
      iEvent.getByToken(recHitsEE_, recHitHandleEE);
      iEvent.getByToken(recHitsFH_, recHitHandleFH);
      iEvent.getByToken(recHitsBH_, recHitHandleBH);
      const auto &rechitsEE = *recHitHandleEE;
      const auto &rechitsFH = *recHitHandleFH;
      const auto &rechitsBH = *recHitHandleBH;
      for (unsigned int i = 0; i < rechitsEE.size(); ++i) {
        hitmap_[rechitsEE[i].detid()] = &rechitsEE[i];
      }
      for (unsigned int i = 0; i < rechitsFH.size(); ++i) {
        hitmap_[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for (unsigned int i = 0; i < rechitsBH.size(); ++i) {
        hitmap_[rechitsBH[i].detid()] = &rechitsBH[i];
      }
      break;
    }
    case 2: {
      iEvent.getByToken(recHitsEE_, recHitHandleEE);
      const HGCRecHitCollection &rechitsEE = *recHitHandleEE;
      for (unsigned int i = 0; i < rechitsEE.size(); i++) {
        hitmap_[rechitsEE[i].detid()] = &rechitsEE[i];
      }
      break;
    }
    case 3: {
      iEvent.getByToken(recHitsFH_, recHitHandleFH);
      iEvent.getByToken(recHitsBH_, recHitHandleBH);
      const auto &rechitsFH = *recHitHandleFH;
      const auto &rechitsBH = *recHitHandleBH;
      for (unsigned int i = 0; i < rechitsFH.size(); i++) {
        hitmap_[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for (unsigned int i = 0; i < rechitsBH.size(); i++) {
        hitmap_[rechitsBH[i].detid()] = &rechitsBH[i];
      }
      break;
    }
    default:
      break;
  }

  const reco::CaloClusterCollection &clusters = *clusterHandle;
  unsigned int nclus = clusters.size();
  cluster_index_ = 0;
  rechit_index_ = 0;
  storedLayerClusters_.clear();
  storedRecHits_.clear();
  for (unsigned int i = 0; i < multiClusters.size(); i++) {
    int cl2dSeed = 0;
    std::set<int> layers;
    std::vector<unsigned int> cl2dIndices;

    for (reco::HGCalMultiCluster::component_iterator it = multiClusters[i].begin();
         it != multiClusters[i].end(); it++) {
      if ((*it)->energy() > (*(it + cl2dSeed))->energy()) cl2dSeed = it - multiClusters[i].begin();
      cl2dIndices.push_back(cluster_index_);
      int layer = fillLayerCluster(*it, true, i);
      layers.insert(layer);
    }  // end of loop on layer clusters

    double pt = multiClusters[i].energy() / cosh(multiClusters[i].eta());

    multiclus_eta_.push_back(multiClusters[i].eta());
    multiclus_phi_.push_back(multiClusters[i].phi());
    multiclus_pt_.push_back(pt);
    multiclus_energy_.push_back(multiClusters[i].energy());
    multiclus_z_.push_back(multiClusters[i].z());
    multiclus_slopeX_.push_back(multiClusters[i].x());
    multiclus_slopeY_.push_back(multiClusters[i].y());
    multiclus_cluster2d_.push_back(cl2dIndices);
    multiclus_cl2dSeed_.push_back(cl2dSeed);
    multiclus_firstLay_.push_back(*layers.begin());
    multiclus_lastLay_.push_back(*layers.rbegin());
    multiclus_NLay_.push_back(layers.size());
  }  // end of loop on multiclusters

  // Fills the additional 2d layers
  for (unsigned ic = 0; ic < nclus; ++ic) {
    edm::Ptr<reco::BasicCluster> clusterPtr(clusterHandle, ic);
    if (storedLayerClusters_.find(clusterPtr) == storedLayerClusters_.end()) {
      double pt = clusterPtr->energy() / cosh(clusterPtr->eta());
      if (pt > layerClusterPtThreshold_) {
        fillLayerCluster(clusterPtr, rawRecHits_);
      }
    }
  }

  // Fill remaining RecHits
  if (rawRecHits_) {
    if (algo_ < 3) {
      const HGCRecHitCollection &rechitsEE = *recHitHandleEE;
      // loop over EE RecHits
      for (HGCRecHitCollection::const_iterator it_hit = rechitsEE.begin(); it_hit < rechitsEE.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->detid();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);

        if (storedRecHits_.find(detid) == storedRecHits_.end()) {
          fillRecHit(detid, -1, layer);
        }
      }
    }
    if (algo_ != 2) {
      const HGCRecHitCollection &rechitsFH = *recHitHandleFH;
      // loop over FH RecHits
      for (HGCRecHitCollection::const_iterator it_hit = rechitsFH.begin(); it_hit < rechitsFH.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->detid();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);

        if (storedRecHits_.find(detid) == storedRecHits_.end()) {
          fillRecHit(detid, -1, layer);
        }
      }
      const HGCRecHitCollection &rechitsBH = *recHitHandleBH;
      // loop over BH RecHits
      for (HGCRecHitCollection::const_iterator it_hit = rechitsBH.begin(); it_hit < rechitsBH.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->detid();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);

        if (storedRecHits_.find(detid) == storedRecHits_.end()) {
          fillRecHit(detid, -1, layer);
        }
      }
    }
  }

  if (readGen_) {
    Handle<std::vector<reco::GenParticle>> genParticlesHandle;
    iEvent.getByToken(genParticles_, genParticlesHandle);
    for (std::vector<reco::GenParticle>::const_iterator it_p = genParticlesHandle->begin();
         it_p != genParticlesHandle->end(); ++it_p) {
      gen_eta_.push_back(it_p->eta());
      gen_phi_.push_back(it_p->phi());
      gen_pt_.push_back(it_p->pt());
      gen_energy_.push_back(it_p->energy());
      gen_charge_.push_back(it_p->charge());
      gen_pdgid_.push_back(it_p->pdgId());
      gen_status_.push_back(it_p->status());
      std::vector<int> daughters(it_p->daughterRefVector().size(), 0);
      for (unsigned j = 0; j < it_p->daughterRefVector().size(); ++j) {
        daughters[j] = static_cast<int>(it_p->daughterRefVector().at(j).key());
      }
      gen_daughters_.push_back(daughters);
    }
  }

  // loop over simClusters
  for (std::vector<SimCluster>::const_iterator it_simClus = simClusters->begin();
       it_simClus != simClusters->end(); ++it_simClus) {
    const std::vector<std::pair<uint32_t, float>> hits_and_fractions =
        it_simClus->hits_and_fractions();
    std::vector<uint32_t> hits;
    std::vector<int> hits_indices;
    std::vector<float> fractions;
    std::vector<unsigned int> layers;
    std::vector<int> wafers_u;
    std::vector<int> wafers_v;
    std::vector<int> cells_u;
    std::vector<int> cells_v;
    for (std::vector<std::pair<uint32_t, float>>::const_iterator it_haf =
             hits_and_fractions.begin();
         it_haf != hits_and_fractions.end(); ++it_haf) {
      auto index_iterator = detIdToRecHitIndexMap_.find( it_haf->first);
      int index_simcluster = -1;
      if(index_iterator != detIdToRecHitIndexMap_.end())
        index_simcluster = index_iterator->second;
      hits_indices.push_back(index_simcluster);
      hits.push_back(it_haf->first);
      fractions.push_back(it_haf->second);
      layers.push_back(recHitTools_.getLayerWithOffset(it_haf->first));
      if (DetId(it_haf->first).det() == DetId::Forward || DetId(it_haf->first).det() == DetId::HGCalEE || DetId(it_haf->first).det() == DetId::HGCalHSi) {
        std::pair<int, int> this_wafer = recHitTools_.getWafer(it_haf->first);
        std::pair<int, int> this_cell = recHitTools_.getCell(it_haf->first);
        wafers_u.push_back(this_wafer.first);
        wafers_v.push_back(this_wafer.second);
        cells_u.push_back(this_cell.first);
        cells_v.push_back(this_cell.second);
      } else {
        wafers_u.push_back(std::numeric_limits<unsigned int>::max());
        wafers_v.push_back(std::numeric_limits<unsigned int>::max());
        cells_u.push_back(std::numeric_limits<unsigned int>::max());
        cells_v.push_back(std::numeric_limits<unsigned int>::max());
      }
    }

    simcluster_eta_.push_back(it_simClus->eta());
    simcluster_phi_.push_back(it_simClus->phi());
    simcluster_pt_.push_back(it_simClus->pt());
    simcluster_energy_.push_back(it_simClus->energy());
    simcluster_simEnergy_.push_back(it_simClus->simEnergy());
    simcluster_hits_.push_back(hits);
    simcluster_hits_indices_.push_back(hits_indices);
    simcluster_fractions_.push_back(fractions);
    simcluster_layers_.push_back(layers);
    simcluster_wafers_u_.push_back(wafers_u);
    simcluster_wafers_v_.push_back(wafers_v);
    simcluster_cells_u_.push_back(cells_u);
    simcluster_cells_v_.push_back(cells_v);

  }  // end loop over simClusters

  // loop over pfClusters
  for (std::vector<reco::PFCluster>::const_iterator it_pfClus = pfClusters.begin();
       it_pfClus != pfClusters.end(); ++it_pfClus) {
    const std::vector<std::pair<DetId, float>> hits_and_fractions = it_pfClus->hitsAndFractions();
    std::vector<uint32_t> hits;
    std::vector<float> fractions;
    for (std::vector<std::pair<DetId, float>>::const_iterator it_haf = hits_and_fractions.begin();
         it_haf != hits_and_fractions.end(); ++it_haf) {
      hits.push_back(it_haf->first);
      fractions.push_back(it_haf->second);
    }
    pfcluster_eta_.push_back(it_pfClus->eta());
    pfcluster_phi_.push_back(it_pfClus->phi());
    pfcluster_pt_.push_back(it_pfClus->pt());
    pfcluster_energy_.push_back(it_pfClus->energy());
    pfcluster_correctedEnergy_.push_back(it_pfClus->correctedEnergy());
    pfcluster_hits_.push_back(hits);
    pfcluster_fractions_.push_back(fractions);

  }  // end loop over pfClusters

  // loop over pfClusters From MultiClusters (python label particleFlowClusterHGCalFromMultiCl)
  for (std::vector<reco::PFCluster>::const_iterator it_pfClus = pfClustersFromMultiCl.begin();
       it_pfClus != pfClustersFromMultiCl.end(); ++it_pfClus) {
    const std::vector<std::pair<DetId, float>> hits_and_fractions = it_pfClus->hitsAndFractions();
    std::vector<uint32_t> hits;
    std::vector<float> fractions;
    float energyEE = 0.;
    float energyFH = 0.;
    float energyBH = 0.;
    for (std::vector<std::pair<DetId, float>>::const_iterator it_haf = hits_and_fractions.begin();
         it_haf != hits_and_fractions.end(); ++it_haf) {
      hits.push_back(it_haf->first);
      fractions.push_back(it_haf->second);
      HGCalDetId detid = HGCalDetId(it_haf->first);
      if (detid.subdetId() == HGCEE || detid.det() == DetId::HGCalEE)
        energyEE += hitmap_[it_haf->first]->energy() * it_haf->second;
      else if (detid.subdetId() == HGCHEF || (detid.det() == DetId::HGCalHSi))
        energyFH += hitmap_[it_haf->first]->energy() * it_haf->second;
      else if (detid.subdetId() == BHM || detid.det() == DetId::HGCalHSc)
        energyBH += hitmap_[it_haf->first]->energy() * it_haf->second;
      else
        assert(false);
    }
    pfclusterFromMultiCl_pos_.push_back(it_pfClus->position());
    pfclusterFromMultiCl_eta_.push_back(it_pfClus->eta());
    pfclusterFromMultiCl_phi_.push_back(it_pfClus->phi());
    pfclusterFromMultiCl_pt_.push_back(it_pfClus->pt());
    pfclusterFromMultiCl_energy_.push_back(it_pfClus->energy());
    pfclusterFromMultiCl_energyEE_.push_back(energyEE);
    pfclusterFromMultiCl_energyFH_.push_back(energyFH);
    pfclusterFromMultiCl_energyBH_.push_back(energyBH);
    pfclusterFromMultiCl_correctedEnergy_.push_back(it_pfClus->correctedEnergy());
    pfclusterFromMultiCl_hits_.push_back(hits);
    pfclusterFromMultiCl_fractions_.push_back(fractions);
    std::vector<uint32_t> rechits_indices;
    for (auto const i : hits) {
      assert(detIdToRecHitIndexMap_.find(i) != detIdToRecHitIndexMap_.end());
      rechits_indices.push_back(detIdToRecHitIndexMap_[i]);
    }
    pfclusterFromMultiCl_rechits_.push_back(rechits_indices);
  }  // end loop over pfClusters From MultiClusters

  // Loop over Ecal Driven Gsf Electrons From MultiClusters
  if (storeElectrons_) {
    for (auto const &ele : electrons) {
      std::vector<uint32_t> pfclustersIndex;
      auto const &sc = ele.superCluster();
      float hoe = 0.;
      float energyEE = 0.;
      float energyFH = 0.;
      float energyBH = 0.;
      for (reco::CaloCluster_iterator cl = sc->clustersBegin(); cl != sc->clustersEnd(); ++cl) {
        if ((*cl)->seed().det() == DetId::Forward || (*cl)->seed().det() == DetId::HGCalEE || (*cl)->seed().det() == DetId::HGCalHSi) {
          if (false)
            std::cout << "SuperCluster Key: " << sc.key() << " own CaloCluster Key: " << cl->key();
          if (electrons_ValueMapClusters.contains(cl->id())) {
            auto pfClusterKey = electrons_ValueMapClusters[*cl].key();
            pfclustersIndex.push_back(pfClusterKey);
            // Redefine HoE for the HGCAL case
            energyEE += pfclusterFromMultiCl_energyEE_[pfClusterKey];
            energyFH += pfclusterFromMultiCl_energyFH_[pfClusterKey];
            energyBH += pfclusterFromMultiCl_energyBH_[pfClusterKey];
            if (false) {
              std::cout << " PFCluster key: " << electrons_ValueMapClusters[*cl].key() << std::endl;
              std::cout << " PFCluster looks like: " << (*electrons_ValueMapClusters[*cl])
                        << std::endl;
              std::cout << " Own CaloCluster looks like: " << *(*cl) << std::endl;
              std::cout << " CastToPFCluster looks like: "
                        << *(dynamic_cast<const reco::PFCluster *>(&*electrons_ValueMapClusters[*cl]))
                        << std::endl;
              for (auto const &pfrh :
                   dynamic_cast<const reco::PFCluster *>(&*electrons_ValueMapClusters[*cl])
                       ->recHitFractions()) {
                std::cout << " PFRecHit key: " << pfrh.recHitRef().key() << std::endl;
                if (pfrh.recHitRef().isAvailable()) std::cout << pfrh << std::endl;
              }
            }  // end of DEBUG section
            assert(pfClusterKey <= pfclusterFromMultiCl_eta_.size());
          }
          hoe = (energyFH + energyBH) / (energyEE + energyFH + energyBH);
        }  // is within HGCAL
      }    // End of loop over clusters within the SC
      ecalDrivenGsfele_charge_.push_back(ele.charge());
      ecalDrivenGsfele_eta_.push_back(ele.eta());
      ecalDrivenGsfele_phi_.push_back(ele.phi());
      ecalDrivenGsfele_pt_.push_back(ele.pt());
      ecalDrivenGsfele_scpos_.push_back(ele.superCluster()->position());
      ecalDrivenGsfele_sceta_.push_back(ele.superCluster()->eta());
      ecalDrivenGsfele_scphi_.push_back(ele.superCluster()->phi());
      ecalDrivenGsfele_seedlayer_.push_back(
          recHitTools_.getLayerWithOffset(ele.superCluster()->seed()->seed()));
      ecalDrivenGsfele_seedpos_.push_back(ele.superCluster()->seed()->position());
      ecalDrivenGsfele_seedeta_.push_back(ele.superCluster()->seed()->eta());
      ecalDrivenGsfele_seedphi_.push_back(ele.superCluster()->seed()->phi());
      ecalDrivenGsfele_seedenergy_.push_back(ele.superCluster()->seed()->energy());
      ecalDrivenGsfele_energy_.push_back(ele.energy());
      ecalDrivenGsfele_energyEE_.push_back(energyEE);
      ecalDrivenGsfele_energyFH_.push_back(energyFH);
      ecalDrivenGsfele_energyBH_.push_back(energyBH);
      ecalDrivenGsfele_isEB_.push_back(ele.isEB());
      ecalDrivenGsfele_hoe_.push_back(hoe);
      ecalDrivenGsfele_numClinSC_.push_back(sc->clusters().size());
      ecalDrivenGsfele_track_dxy_.push_back(ele.gsfTrack()->dxy(vertices[0].position()));
      ecalDrivenGsfele_track_dz_.push_back(ele.gsfTrack()->dz(vertices[0].position()));
      ecalDrivenGsfele_track_simdxy_.push_back(ele.gsfTrack()->dxy(sim_pv));
      ecalDrivenGsfele_track_simdz_.push_back(ele.gsfTrack()->dz(sim_pv));
      ecalDrivenGsfele_deltaEtaSuperClusterTrackAtVtx_.push_back(
          ele.deltaEtaSuperClusterTrackAtVtx());
      ecalDrivenGsfele_deltaPhiSuperClusterTrackAtVtx_.push_back(
          ele.deltaPhiSuperClusterTrackAtVtx());
      ecalDrivenGsfele_deltaEtaEleClusterTrackAtCalo_.push_back(ele.deltaEtaEleClusterTrackAtCalo());
      ecalDrivenGsfele_deltaPhiEleClusterTrackAtCalo_.push_back(ele.deltaPhiEleClusterTrackAtCalo());
      ecalDrivenGsfele_deltaEtaSeedClusterTrackAtCalo_.push_back(
          ele.deltaEtaSeedClusterTrackAtCalo());
      ecalDrivenGsfele_deltaPhiSeedClusterTrackAtCalo_.push_back(
          ele.deltaPhiSeedClusterTrackAtCalo());
      ecalDrivenGsfele_eSuperClusterOverP_.push_back(ele.eSuperClusterOverP());
      ecalDrivenGsfele_eSeedClusterOverP_.push_back(ele.eSeedClusterOverP());
      ecalDrivenGsfele_eSeedClusterOverPout_.push_back(ele.eSeedClusterOverPout());
      ecalDrivenGsfele_eEleClusterOverPout_.push_back(ele.eEleClusterOverPout());
      ecalDrivenGsfele_pfClusterIndex_.push_back(pfclustersIndex);
    }  // End of loop over electrons
  }

  // loop over caloParticles
  if (readCaloParticles_) {
    for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles->begin();
         it_caloPart != caloParticles->end(); ++it_caloPart) {
      const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();
      std::vector<uint32_t> simClusterIndex;
      for (CaloParticle::sc_iterator it_sc = simClusterRefVector.begin();
           it_sc != simClusterRefVector.end(); ++it_sc) {
        simClusterIndex.push_back((*it_sc).key());
      }

      calopart_eta_.push_back(it_caloPart->eta());
      calopart_phi_.push_back(it_caloPart->phi());
      calopart_pt_.push_back(it_caloPart->pt());
      calopart_energy_.push_back(it_caloPart->energy());
      calopart_simEnergy_.push_back(it_caloPart->simEnergy());
      calopart_simClusterIndex_.push_back(simClusterIndex);

    }  // end loop over caloParticles
  }

  // loop over tracks

  //  random = new RandomEngineAndDistribution(iEvent.streamID());

  // prepare for RK propagation
  defaultRKPropagator::Product prod(aField_, alongMomentum, 5.e-5);
  auto &RKProp = prod.propagator;

  for (std::vector<reco::Track>::const_iterator it_track = tracks.begin(); it_track != tracks.end();
       ++it_track) {
    if (!it_track->quality(reco::Track::highPurity)) continue;

    double energy = it_track->pt() * cosh(it_track->eta());

    // save info about reconstructed tracks propoagation to hgcal layers (ony
    // for pt>propagationPtThreshold_ tracks)
    std::vector<float> xp, yp, zp;

    if (it_track->pt() >= propagationPtThreshold_) {
      // Define error matrix
      ROOT::Math::SMatrixIdentity id;
      AlgebraicSymMatrix55 C(id);
      C *= 0.01;
      CurvilinearTrajectoryError err(C);
      typedef TrajectoryStateOnSurface TSOS;

      GlobalPoint startingPosition(it_track->vx(), it_track->vy(), it_track->vz());
      GlobalVector startingMomentum(it_track->px(), it_track->py(), it_track->pz());

      Plane::PlanePointer startingPlane =
          Plane::build(Plane::PositionType(it_track->vx(), it_track->vy(), it_track->vz()),
                       Plane::RotationType());

      TSOS startingStateP(GlobalTrajectoryParameters(startingPosition, startingMomentum,
                                                     it_track->charge(), aField_),
                          err, *startingPlane);

      for (unsigned il = 0; il < layerPositions_.size(); ++il) {
        float xp_curr = 0;
        float yp_curr = 0;
        float zp_curr = 0;

        for (int zside = -1; zside <= 1; zside += 2) {
          // clearly try both sides
          Plane::PlanePointer endPlane = Plane::build(
              Plane::PositionType(0, 0, zside * layerPositions_[il]), Plane::RotationType());
          try {
            /*
            std::cout << "Trying from " <<
            " layer " << il <<
            " starting point " << startingStateP.globalPosition() <<
            std::endl;
            */
            TSOS trackStateP = RKProp.propagate(startingStateP, *endPlane);
            if (trackStateP.isValid()) {
              xp_curr = trackStateP.globalPosition().x();
              yp_curr = trackStateP.globalPosition().y();
              zp_curr = trackStateP.globalPosition().z();

              // std::cout << "Succesfully finished Positive track propagation
              // -------------- with RK: " << trackStateP.globalPosition() <<
              // std::endl;
            }
          } catch (...) {
            std::cout << "MagVolumeOutsideValidity not properly caught!! Lost "
                         "this track "
                      << std::endl;
          }
        }
        xp.push_back(xp_curr);
        yp.push_back(yp_curr);
        zp.push_back(zp_curr);
      }  // closes loop on layers
    }    // closes conditions pt>3

    // save info in tree
    track_pt_.push_back(it_track->pt());
    track_eta_.push_back(it_track->eta());
    track_phi_.push_back(it_track->phi());
    track_energy_.push_back(energy);
    track_charge_.push_back(it_track->charge());
    track_posx_.push_back(xp);
    track_posy_.push_back(yp);
    track_posz_.push_back(zp);

  }  // end loop over tracks

  if (storePFCandidates_) {
    Handle<std::vector<reco::PFCandidate>> pfCandidatesHandle;
    iEvent.getByToken(pfCandidates_, pfCandidatesHandle);
    for (std::vector<reco::PFCandidate>::const_iterator it_p = pfCandidatesHandle->begin();
         it_p != pfCandidatesHandle->end(); ++it_p) {
      pfcandidate_eta_.push_back(it_p->eta());
      pfcandidate_phi_.push_back(it_p->phi());
      pfcandidate_pt_.push_back(it_p->pt());
      pfcandidate_energy_.push_back(it_p->energy());
      pfcandidate_pdgid_.push_back(it_p->pdgId());
    }
  }

  ev_event_ = iEvent.id().event();
  ev_lumi_ = iEvent.id().luminosityBlock();
  ev_run_ = iEvent.id().run();

  vtx_x_ = vx_;
  vtx_y_ = vy_;
  vtx_z_ = vz_;

  t_->Fill();
}

void HGCalAnalysis::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
  edm::ESHandle<HepPDT::ParticleDataTable> pdt;
  es.getData(pdt);
  mySimEvent_->initializePdt(&(*pdt));

  recHitTools_.getEventSetup(es);
  retrieveLayerPositions(es, recHitTools_.lastLayerBH());

  edm::ESHandle<MagneticField> magfield;
  es.get<IdealMagneticFieldRecord>().get(magfield);

  aField_ = &(*magfield);
}

void HGCalAnalysis::endRun(edm::Run const &iEvent, edm::EventSetup const &) {}

void HGCalAnalysis::beginJob() { ; }

// ------------ method called once each job just after ending the event loop
// ------------
void HGCalAnalysis::endJob() {}

// ------------ method to be called once
// --------------------------------------------------

void HGCalAnalysis::retrieveLayerPositions(const edm::EventSetup &es, unsigned layers) {
  recHitTools_.getEventSetup(es);

  for (unsigned ilayer = 1; ilayer <= layers; ++ilayer) {
    const GlobalPoint pos = recHitTools_.getPositionLayer(ilayer);
    layerPositions_.push_back(pos.z());
  }
}

int HGCalAnalysis::fillLayerCluster(const edm::Ptr<reco::CaloCluster> &layerCluster,
                                    const bool &fillRecHits, const int &multiClusterIndex) {
  // std::cout << "in fillLayerCluster" << std::endl;
  const std::vector<std::pair<DetId, float>> &hf = layerCluster->hitsAndFractions();
  std::vector<unsigned int> rhIndices;
  int ncoreHit = 0;
  int layer = 0;
  int rhSeed = 0;
  if (!fillRecHits) {
    rhSeed = -1;
  }
  float maxEnergy = -1.;
  unsigned hfsize = hf.size();

  for (unsigned int j = 0; j < hfsize; j++) {
    // here we loop over detid/fraction pairs
    float fraction = hf[j].second;
    const DetId rh_detid = hf[j].first;
    layer = recHitTools_.getLayerWithOffset(rh_detid);
    const HGCRecHit *hit = hitmap_[rh_detid];
    ncoreHit += int(fraction);

    if (fillRecHits) {
      if (storedRecHits_.find(rh_detid) == storedRecHits_.end()) {
        // std::cout << "in fillLayerCluster: RecHit not yet filled" <<
        // std::endl;
        // std::cout << "in fillLayerCluster: hit energy: " << hit->energy() <<
        // std::endl;
        // std::cout << "in fillLayerCluster: first hit energy: " << maxEnergy
        // << std::endl;
        if (hit->energy() > maxEnergy) {
          rhSeed = rechit_index_;
          maxEnergy = hit->energy();
        }
        rhIndices.push_back(rechit_index_);
        fillRecHit(rh_detid, fraction, layer, cluster_index_);
      } else {
        // need to see what to do about existing rechits in case of sharing
        if (verbose_) {
          std::cout << "RecHit already filled for different layer cluster: " << int(rh_detid)
                   << std::endl;
        }
      }
    }
  }

  double pt = layerCluster->energy() / cosh(layerCluster->eta());

  cluster2d_eta_.push_back(layerCluster->eta());
  cluster2d_phi_.push_back(layerCluster->phi());
  cluster2d_pt_.push_back(pt);
  cluster2d_energy_.push_back(layerCluster->energy());
  cluster2d_x_.push_back(layerCluster->x());
  cluster2d_y_.push_back(layerCluster->y());
  cluster2d_z_.push_back(layerCluster->z());
  cluster2d_layer_.push_back(layer);
  cluster2d_nhitCore_.push_back(ncoreHit);
  cluster2d_nhitAll_.push_back(hf.size());
  cluster2d_multicluster_.push_back(multiClusterIndex);
  cluster2d_rechitSeed_.push_back(rhSeed);
  cluster2d_rechits_.push_back(rhIndices);

  storedLayerClusters_.insert(layerCluster);
  ++cluster_index_;
  return layer;
}

void HGCalAnalysis::fillRecHit(const DetId &detid, const float &fraction, const unsigned int &layer,
                               const int &cluster_index_) {
  // std::cout << "in fillRecHit" << std::endl;
  int flags = 0x0;
  if (fraction > 0. && fraction < 1.)
    flags = 0x1;
  else if (fraction < 0.)
    flags = 0x3;
  else if (fraction == 0.)
    flags = 0x2;
  const HGCRecHit *hit = hitmap_[detid];

  const GlobalPoint position = recHitTools_.getPosition(detid);
  std::pair<int, int> wafer;
  std::pair<int, int> cell;
  if (detid.det() == DetId::Forward || detid.det() == DetId::HGCalEE || detid.det() == DetId::HGCalHSi) {
    wafer = recHitTools_.getWafer(detid);
    cell = recHitTools_.getCell(detid);
  }
  else {
    wafer = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());
    cell = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());
  }
  const double cellThickness =
      ((detid.det() == DetId::Forward || detid.det() == DetId::HGCalEE || detid.det() == DetId::HGCalHSi) ? recHitTools_.getSiThickness(detid)
                                            : std::numeric_limits<std::float_t>::max());
  const bool isHalfCell = recHitTools_.isHalfCell(detid);
  const double eta = recHitTools_.getEta(position);
  const double phi = recHitTools_.getPhi(position);
  const double pt = recHitTools_.getPt(position, hit->energy());
  const double radius =
      ((detid.det() == DetId::Forward || detid.det() == DetId::HGCalEE || detid.det() == DetId::HGCalHSi) ? recHitTools_.getRadiusToSide(detid) : -1.);

  // fill the vectors
  rechit_eta_.push_back(eta);
  rechit_phi_.push_back(phi);
  rechit_pt_.push_back(pt);
  rechit_energy_.push_back(hit->energy());
  rechit_layer_.push_back(layer);
  rechit_wafer_u_.push_back(wafer.first);
  rechit_wafer_v_.push_back(wafer.second);
  rechit_cell_u_.push_back(cell.first);
  rechit_cell_v_.push_back(cell.second);
  rechit_detid_.push_back(detid);
  rechit_x_.push_back(position.x());
  rechit_y_.push_back(position.y());
  rechit_z_.push_back(position.z());
  rechit_time_.push_back(hit->time());
  rechit_thickness_.push_back(cellThickness);
  rechit_isHalf_.push_back(isHalfCell);
  rechit_flags_.push_back(flags);
  rechit_cluster2d_.push_back(cluster_index_);
  rechit_radius_.push_back(radius);

  storedRecHits_.insert(detid);
  detIdToRecHitIndexMap_[detid] = rechit_index_;
  ++rechit_index_;
}


// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------

void HGCalAnalysis::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

/*
   Surface::RotationType HGCalAnalysis::rotation( const GlobalVector& zDir)
   const
   {
   GlobalVector zAxis = zDir.unit();
   GlobalVector yAxis( zAxis.y(), -zAxis.x(), 0);
   GlobalVector xAxis = yAxis.cross( zAxis);
   return Surface::RotationType( xAxis, yAxis, zAxis);
   }
 */

// define this as a plug-in
DEFINE_FWK_MODULE(HGCalAnalysis);
