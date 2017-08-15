// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

// track data formats
#include "DataFormats/TrackReco/interface/Track.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "DataFormats/GeometrySurface/interface/PlaneBuilder.h"

#include "FastSimulation/Event/interface/FSimEvent.h"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "FastSimulation/Event/interface/FSimVertex.h"
#include "FastSimulation/Particle/interface/ParticleTable.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackPropagation/RungeKutta/interface/defaultRKPropagator.h"
#include "MagneticField/VolumeGeometry/interface/MagVolumeOutsideValidity.h"


#include "TTree.h"
#include "TH1F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include <string>
#include <map>
#include <vector>
#include <set>

#include "TTree.h"

class HGCalAnalysis : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
HGCalAnalysis();
explicit HGCalAnalysis(const edm::ParameterSet&);
~HGCalAnalysis();

static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
virtual void beginRun(edm::Run const& iEvent, edm::EventSetup const&) override;
virtual void endRun(edm::Run const& iEvent, edm::EventSetup const&) override;

private:
virtual void beginJob() override;
virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
virtual void endJob() override;
virtual void fillLayerCluster(const edm::Ptr<reco::CaloCluster>& layerCluster, const bool& fillRecHits, const int& multiClusterIndex = -1);
virtual void fillRecHit(const DetId& detid, const float& fraction, const unsigned int& layer, const int& cluster_index = -1);

void clearVariables();

void retrieveLayerPositions(const edm::EventSetup&, unsigned layers);

// ---------parameters ----------------------------
bool readOfficialReco;
bool readCaloParticles;
double layerClusterPtThreshold;
double propagationPtThreshold;
std::string detector;
bool rawRecHits;
bool readGen;

// ----------member data ---------------------------

edm::EDGetTokenT<HGCRecHitCollection> _recHitsEE;
edm::EDGetTokenT<HGCRecHitCollection> _recHitsFH;
edm::EDGetTokenT<HGCRecHitCollection> _recHitsBH;
edm::EDGetTokenT<reco::CaloClusterCollection> _clusters;
edm::EDGetTokenT<std::vector<TrackingVertex> > _vtx;
edm::EDGetTokenT<std::vector<TrackingParticle> > _part;
edm::EDGetTokenT<std::vector<SimCluster> > _simClusters;
edm::EDGetTokenT<std::vector<reco::PFCluster> > _pfClusters;
edm::EDGetTokenT<std::vector<CaloParticle> > _caloParticles;
edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster> > _multiClusters;
edm::EDGetTokenT<std::vector<SimTrack> > _simTracks;
edm::EDGetTokenT<std::vector<SimVertex> > _simVertices;
edm::EDGetTokenT<edm::HepMCProduct> _hev;
edm::EDGetTokenT<std::vector<reco::Track> > _tracks;
edm::EDGetTokenT<std::vector<reco::GenParticle> > _genParticles;

TTree                     *t;

////////////////////
// event
//
edm::RunNumber_t ev_run;
edm::LuminosityBlockNumber_t ev_lumi;
edm::EventNumber_t ev_event;
float vtx_x;
float vtx_y;
float vtx_z;

////////////////////
// GenParticles
//
std::vector<float> genpart_eta;
std::vector<float> genpart_phi;
std::vector<float> genpart_pt;
std::vector<float> genpart_energy;
std::vector<float> genpart_dvx;
std::vector<float> genpart_dvy;
std::vector<float> genpart_dvz;
std::vector<float> genpart_fbrem;
std::vector<int> genpart_pid;
std::vector<int> genpart_gen;
std::vector<bool> genpart_reachedEE;
std::vector<bool> genpart_fromBeamPipe;
std::vector<std::vector<float> > genpart_posx;
std::vector<std::vector<float> > genpart_posy;
std::vector<std::vector<float> > genpart_posz;

////////////////////
// GenParticles
//
std::vector<float> gen_eta;
std::vector<float> gen_phi;
std::vector<float> gen_pt;
std::vector<float> gen_energy;
std::vector<int> gen_pdgid;
std::vector<int> gen_status;


////////////////////
// RecHits
// associated to layer clusters
std::vector<float> rechit_eta;
std::vector<float> rechit_phi;
std::vector<float> rechit_pt;
std::vector<float> rechit_energy;
std::vector<float> rechit_x;
std::vector<float> rechit_y;
std::vector<float> rechit_z;
std::vector<float> rechit_time;
std::vector<float> rechit_thickness;
std::vector<int> rechit_layer;
std::vector<int> rechit_wafer;
std::vector<int> rechit_cell;
std::vector<unsigned int> rechit_detid;
std::vector<bool> rechit_isHalf;
std::vector<int> rechit_flags;
std::vector<int> rechit_cluster2d;

////////////////////
// layer clusters
//
std::vector<float> cluster2d_eta;
std::vector<float> cluster2d_phi;
std::vector<float> cluster2d_pt;
std::vector<float> cluster2d_energy;
std::vector<float> cluster2d_x;
std::vector<float> cluster2d_y;
std::vector<float> cluster2d_z;
std::vector<int> cluster2d_layer;
std::vector<int> cluster2d_nhitCore;
std::vector<int> cluster2d_nhitAll;
std::vector<int> cluster2d_multicluster;
std::vector<std::vector<unsigned int> > cluster2d_rechits;
std::vector<int> cluster2d_rechitSeed;

////////////////////
// multi clusters
//
std::vector<float> multiclus_eta;
std::vector<float> multiclus_phi;
std::vector<float> multiclus_pt;
std::vector<float> multiclus_energy;
std::vector<float> multiclus_z;
std::vector<float> multiclus_slopeX;
std::vector<float> multiclus_slopeY;
std::vector<std::vector<unsigned int> > multiclus_cluster2d;
std::vector<int> multiclus_cl2dSeed;

////////////////////
// sim clusters
//
std::vector<float> simcluster_eta;
std::vector<float> simcluster_phi;
std::vector<float> simcluster_pt;
std::vector<float> simcluster_energy;
std::vector<float> simcluster_simEnergy;
std::vector<std::vector<uint32_t> > simcluster_hits;
std::vector<std::vector<float> > simcluster_fractions;
std::vector<std::vector<unsigned int> > simcluster_layers;
std::vector<std::vector<unsigned int> > simcluster_wafers;
std::vector<std::vector<unsigned int> > simcluster_cells;


////////////////////
// PF clusters
//
std::vector<float> pfcluster_eta;
std::vector<float> pfcluster_phi;
std::vector<float> pfcluster_pt;
std::vector<float> pfcluster_energy;


////////////////////
// calo particles
//
std::vector<float> calopart_eta;
std::vector<float> calopart_phi;
std::vector<float> calopart_pt;
std::vector<float> calopart_energy;
std::vector<float> calopart_simEnergy;
std::vector<std::vector<uint32_t> > calopart_simClusterIndex;


////////////////////
// high purity tracks
//
std::vector<float> track_eta;
std::vector<float> track_phi;
std::vector<float> track_pt;
std::vector<float> track_energy;
std::vector<int> track_charge;
std::vector<std::vector<float> > track_posx;
std::vector<std::vector<float> > track_posy;
std::vector<std::vector<float> > track_posz;

////////////////////
// helper classes
//
unsigned int cluster_index;
unsigned int rechit_index;
std::map<DetId,const HGCRecHit*> hitmap;
float vz; // primary vertex z position
// to keep track of the 2d clusters stored within the loop on multiclusters
std::set<edm::Ptr<reco::BasicCluster> > storedLayerClusters;
// to keep track of the RecHits stored within the cluster loops
std::set<DetId> storedRecHits;
int algo;
HGCalDepthPreClusterer pre;
hgcal::RecHitTools recHitTools;

// -------convenient tool to deal with simulated tracks
FSimEvent * mySimEvent;
edm::ParameterSet particleFilter;
std::vector <float> layerPositions;

// and also the magnetic field
MagneticField const * aField;
};

HGCalAnalysis::HGCalAnalysis() {
	;
}

HGCalAnalysis::HGCalAnalysis(const edm::ParameterSet& iConfig) :
	readOfficialReco(iConfig.getParameter<bool>("readOfficialReco")),
	readCaloParticles(iConfig.getParameter<bool>("readCaloParticles")),
	layerClusterPtThreshold(iConfig.getParameter<double>("layerClusterPtThreshold")),
	propagationPtThreshold(iConfig.getUntrackedParameter<double>("propagationPtThreshold",3.0)),
	detector(iConfig.getParameter<std::string >("detector")),
	rawRecHits(iConfig.getParameter<bool>("rawRecHits")),
	readGen(iConfig.getParameter<bool>("readGenParticles")),
	particleFilter(iConfig.getParameter<edm::ParameterSet>("TestParticleFilter"))
{
	// now do what ever initialization is needed
	mySimEvent = new FSimEvent(particleFilter);

	if(detector=="all") {
		_recHitsEE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"));
		_recHitsFH = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"));
		_recHitsBH = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEBRecHits"));
		algo = 1;
	}else if(detector=="EM") {
		_recHitsEE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"));
		algo = 2;
	}else if(detector=="HAD") {
		_recHitsFH = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"));
		_recHitsBH = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEBRecHits"));
		algo = 3;
	}
	_clusters = consumes<reco::CaloClusterCollection>(edm::InputTag("hgcalLayerClusters"));
	_simClusters = consumes<std::vector<SimCluster> >(edm::InputTag("mix","MergedCaloTruth"));
	_hev = consumes<edm::HepMCProduct>(edm::InputTag("generatorSmeared") );
	if(!readOfficialReco) {
		_vtx = consumes<std::vector<TrackingVertex> >(edm::InputTag("mix","MergedTrackTruth"));
		_part = consumes<std::vector<TrackingParticle> >(edm::InputTag("mix","MergedTrackTruth"));
	}
	else {
		_simTracks = consumes<std::vector<SimTrack> >(edm::InputTag("g4SimHits"));
		_simVertices = consumes<std::vector<SimVertex> >(edm::InputTag("g4SimHits"));
	}
	if (readCaloParticles) {
		_caloParticles = consumes<std::vector<CaloParticle> >(edm::InputTag("mix","MergedCaloTruth"));
	}
  if (readGen) {
    _genParticles = consumes<std::vector<reco::GenParticle> >(edm::InputTag("genParticles"));
  }
	_pfClusters = consumes<std::vector<reco::PFCluster> >(edm::InputTag("particleFlowClusterHGCal"));
	_multiClusters = consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("hgcalLayerClusters"));
	_tracks = consumes<std::vector<reco::Track> >(edm::InputTag("generalTracks"));


	usesResource(TFileService::kSharedResource);
	edm::Service<TFileService> fs;
	fs->make<TH1F>("total", "total", 100, 0, 5.);

	t = fs->make<TTree>("hgc","hgc");

	// event info
	t->Branch("event", &ev_event);
	t->Branch("lumi", &ev_lumi);
	t->Branch("run", &ev_run);
	t->Branch("vtx_x", &vtx_x);
	t->Branch("vtx_y", &vtx_y);
	t->Branch("vtx_z", &vtx_z);

	t->Branch("genpart_eta", &genpart_eta);
	t->Branch("genpart_phi", &genpart_phi);
	t->Branch("genpart_pt", &genpart_pt);
	t->Branch("genpart_energy", &genpart_energy);
	t->Branch("genpart_dvx", &genpart_dvx);
	t->Branch("genpart_dvy", &genpart_dvy);
	t->Branch("genpart_dvz", &genpart_dvz);
	t->Branch("genpart_fbrem", &genpart_fbrem);
	t->Branch("genpart_pid", &genpart_pid);
	t->Branch("genpart_gen", &genpart_gen);
	t->Branch("genpart_reachedEE", &genpart_reachedEE);
	t->Branch("genpart_fromBeamPipe", &genpart_fromBeamPipe);
	t->Branch("genpart_posx", &genpart_posx);
	t->Branch("genpart_posy", &genpart_posy);
	t->Branch("genpart_posz", &genpart_posz);


  //////////////////
	// reco::GenParticles
	t->Branch("gen_eta", &gen_eta);
	t->Branch("gen_phi", &gen_phi);
	t->Branch("gen_pt", &gen_pt);
	t->Branch("gen_energy", &gen_energy);
	t->Branch("gen_pdgid", &gen_pdgid);
	t->Branch("gen_status", &gen_status);

  //////////////////
	// RecHits
	// associated to layer clusters
	t->Branch("rechit_eta", &rechit_eta);
	t->Branch("rechit_phi", &rechit_phi);
	t->Branch("rechit_pt", &rechit_pt);
	t->Branch("rechit_energy", &rechit_energy);
	t->Branch("rechit_x", &rechit_x);
	t->Branch("rechit_y", &rechit_y);
	t->Branch("rechit_z", &rechit_z);
	t->Branch("rechit_time", &rechit_time);
	t->Branch("rechit_thickness", &rechit_thickness);
	t->Branch("rechit_layer", &rechit_layer);
	t->Branch("rechit_wafer", &rechit_wafer);
	t->Branch("rechit_cell", &rechit_cell);
	t->Branch("rechit_detid", &rechit_detid);
	t->Branch("rechit_isHalf", &rechit_isHalf);
	t->Branch("rechit_flags", &rechit_flags);
	t->Branch("rechit_cluster2d", &rechit_cluster2d);

	////////////////////
	// layer clusters
	//
	t->Branch("cluster2d_eta", &cluster2d_eta);
	t->Branch("cluster2d_phi", &cluster2d_phi);
	t->Branch("cluster2d_pt", &cluster2d_pt);
	t->Branch("cluster2d_energy", &cluster2d_energy);
	t->Branch("cluster2d_x", &cluster2d_x);
	t->Branch("cluster2d_y", &cluster2d_y);
	t->Branch("cluster2d_z", &cluster2d_z);
	t->Branch("cluster2d_layer", &cluster2d_layer);
	t->Branch("cluster2d_nhitCore", &cluster2d_nhitCore);
	t->Branch("cluster2d_nhitAll", &cluster2d_nhitAll);
	t->Branch("cluster2d_multicluster", &cluster2d_multicluster);
	t->Branch("cluster2d_rechits", &cluster2d_rechits);
	t->Branch("cluster2d_rechitSeed", &cluster2d_rechitSeed);

	////////////////////
	// multi clusters
	//
	t->Branch("multiclus_eta", &multiclus_eta);
	t->Branch("multiclus_phi", &multiclus_phi);
	t->Branch("multiclus_pt", &multiclus_pt);
	t->Branch("multiclus_energy", &multiclus_energy);
	t->Branch("multiclus_z", &multiclus_z);
	t->Branch("multiclus_slopeX", &multiclus_slopeX);
	t->Branch("multiclus_slopeY", &multiclus_slopeY);
	t->Branch("multiclus_cluster2d", &multiclus_cluster2d);
	t->Branch("multiclus_cl2dSeed", &multiclus_cl2dSeed);

	////////////////////
	// sim clusters
	//
	t->Branch("simcluster_eta", &simcluster_eta);
	t->Branch("simcluster_phi", &simcluster_phi);
	t->Branch("simcluster_pt", &simcluster_pt);
	t->Branch("simcluster_energy", &simcluster_energy);
	t->Branch("simcluster_simEnergy", &simcluster_simEnergy);
	t->Branch("simcluster_hits", &simcluster_hits);
	t->Branch("simcluster_fractions", &simcluster_fractions);
	t->Branch("simcluster_layers", &simcluster_layers);
	t->Branch("simcluster_wafers", &simcluster_wafers);
	t->Branch("simcluster_cells", &simcluster_cells);


	////////////////////
	// PF clusters
	//
	t->Branch("pfcluster_eta", &pfcluster_eta);
	t->Branch("pfcluster_phi", &pfcluster_phi);
	t->Branch("pfcluster_pt", &pfcluster_pt);
	t->Branch("pfcluster_energy", &pfcluster_energy);


	////////////////////
	// calo particles
	//
	t->Branch("calopart_eta", &calopart_eta);
	t->Branch("calopart_phi", &calopart_phi);
	t->Branch("calopart_pt", &calopart_pt);
	t->Branch("calopart_energy", &calopart_energy);
	t->Branch("calopart_simEnergy", &calopart_simEnergy);
	t->Branch("calopart_simClusterIndex", &calopart_simClusterIndex);


	////////////////////
	// high purity tracks
	//
	t->Branch("track_eta", &track_eta);
	t->Branch("track_phi", &track_phi);
	t->Branch("track_pt", &track_pt);
	t->Branch("track_energy", &track_energy);
	t->Branch("track_charge", &track_charge);
	t->Branch("track_posx", &track_posx);
	t->Branch("track_posy", &track_posy);
	t->Branch("track_posz", &track_posz);

}
HGCalAnalysis::~HGCalAnalysis()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

//
// member functions
//
void HGCalAnalysis::clearVariables() {

	ev_run = 0;
	ev_lumi = 0;
	ev_event = 0;
	vtx_x = 0;
	vtx_y = 0;
	vtx_z = 0;

	////////////////////
	// GenParticles
	//
	genpart_eta.clear();
	genpart_phi.clear();
	genpart_pt.clear();
	genpart_energy.clear();
	genpart_dvx.clear();
	genpart_dvy.clear();
	genpart_dvz.clear();
	genpart_fbrem.clear();
	genpart_pid.clear();
	genpart_gen.clear();
	genpart_reachedEE.clear();
	genpart_fromBeamPipe.clear();
	genpart_posx.clear();
	genpart_posy.clear();
	genpart_posz.clear();

	////////////////////
	// reco::GenParticles
	//
  gen_eta.clear();
	gen_phi.clear();
	gen_pt.clear();
	gen_energy.clear();
	gen_pdgid.clear();
	gen_status.clear();

	////////////////////
	// RecHits
	// associated to layer clusters
	rechit_eta.clear();
	rechit_phi.clear();
	rechit_pt.clear();
	rechit_energy.clear();
	rechit_x.clear();
	rechit_y.clear();
	rechit_z.clear();
	rechit_time.clear();
	rechit_thickness.clear();
	rechit_layer.clear();
	rechit_wafer.clear();
	rechit_cell.clear();
	rechit_detid.clear();
	rechit_isHalf.clear();
	rechit_flags.clear();
	rechit_cluster2d.clear();

	////////////////////
	// layer clusters
	//
	cluster2d_eta.clear();
	cluster2d_phi.clear();
	cluster2d_pt.clear();
	cluster2d_energy.clear();
	cluster2d_x.clear();
	cluster2d_y.clear();
	cluster2d_z.clear();
	cluster2d_layer.clear();
	cluster2d_nhitCore.clear();
	cluster2d_nhitAll.clear();
	cluster2d_multicluster.clear();
	cluster2d_rechits.clear();
	cluster2d_rechitSeed.clear();

	////////////////////
	// multi clusters
	//
	multiclus_eta.clear();
	multiclus_phi.clear();
	multiclus_pt.clear();
	multiclus_energy.clear();
	multiclus_z.clear();
	multiclus_slopeX.clear();
	multiclus_slopeY.clear();
	multiclus_cluster2d.clear();
	multiclus_cl2dSeed.clear();

	////////////////////
	// sim clusters
	//
	simcluster_eta.clear();
	simcluster_phi.clear();
	simcluster_pt.clear();
	simcluster_energy.clear();
	simcluster_simEnergy.clear();
	simcluster_hits.clear();
	simcluster_fractions.clear();
	simcluster_layers.clear();
	simcluster_wafers.clear();
	simcluster_cells.clear();


	////////////////////
	// PF clusters
	//
	pfcluster_eta.clear();
	pfcluster_phi.clear();
	pfcluster_pt.clear();
	pfcluster_energy.clear();


	////////////////////
	// calo particles
	//
	calopart_eta.clear();
	calopart_phi.clear();
	calopart_pt.clear();
	calopart_energy.clear();
	calopart_simEnergy.clear();
	calopart_simClusterIndex.clear();


	////////////////////
	// high purity tracks
	//
	track_eta.clear();
	track_phi.clear();
	track_pt.clear();
	track_energy.clear();
	track_charge.clear();
	track_posx.clear();
	track_posy.clear();
	track_posz.clear();

}


void
HGCalAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	clearVariables();
	// std::cout << "after clearVariables" << std::endl;

	ParticleTable::Sentry ptable(mySimEvent->theTable());
	recHitTools.getEventSetup(iSetup);

	Handle<HGCRecHitCollection> recHitHandleEE;
	Handle<HGCRecHitCollection> recHitHandleFH;
	Handle<HGCRecHitCollection> recHitHandleBH;
	Handle<reco::CaloClusterCollection> clusterHandle;

	iEvent.getByToken(_clusters,clusterHandle);
	Handle<std::vector<TrackingVertex> > vtxHandle;
	Handle<std::vector<TrackingParticle> > partHandle;
	const std::vector<TrackingParticle> * part;

	Handle<edm::HepMCProduct> hevH;
	Handle<std::vector<SimTrack> >simTracksHandle;
	Handle<std::vector<SimVertex> >simVerticesHandle;

	iEvent.getByToken(_hev,hevH);
	if(!readOfficialReco) {
		iEvent.getByToken(_vtx,vtxHandle);
		iEvent.getByToken(_part,partHandle);
		part = &(*partHandle);
	} else // use SimTracks and HepMCProduct
	{
		iEvent.getByToken(_simTracks,simTracksHandle);
		iEvent.getByToken(_simVertices,simVerticesHandle);
		mySimEvent->fill(*simTracksHandle,*simVerticesHandle);

	}


	Handle<std::vector<SimCluster> > simClusterHandle;
	Handle<std::vector<reco::PFCluster> > pfClusterHandle;
	Handle<std::vector<CaloParticle> > caloParticleHandle;
	Handle<std::vector<reco::Track> > trackHandle;

	const std::vector<SimCluster> * simClusters = 0;
	iEvent.getByToken(_simClusters, simClusterHandle);
	simClusters = &(*simClusterHandle);

	iEvent.getByToken(_pfClusters, pfClusterHandle);
	const std::vector<reco::PFCluster>& pfClusters = *pfClusterHandle;

	const std::vector<CaloParticle>* caloParticles;
	if(readCaloParticles) {
		iEvent.getByToken(_caloParticles, caloParticleHandle);
		caloParticles = &(*caloParticleHandle);
	}

	iEvent.getByToken(_tracks, trackHandle);
	const std::vector<reco::Track>& tracks = *trackHandle;

	Handle<std::vector<reco::HGCalMultiCluster> > multiClusterHandle;
	iEvent.getByToken(_multiClusters, multiClusterHandle);
	const std::vector<reco::HGCalMultiCluster>& multiClusters = *multiClusterHandle;
	HepMC::GenVertex * primaryVertex = *(hevH)->GetEvent()->vertices_begin();
	float vx = primaryVertex->position().x()/10.; // to put in official units
	float vy = primaryVertex->position().y()/10.;
	vz = primaryVertex->position().z()/10.;

	// std::cout << "start the fun" << std::endl;

  Handle<std::vector<reco::GenParticle>> genParticlesHandle;
  if (readGen) {
    iEvent.getByToken(_genParticles, genParticlesHandle);
  }


	if( !readOfficialReco) {
		unsigned int npart = part->size();
		for(unsigned int i=0; i<npart; ++i) {

			// event=0 is the hard scattering (all PU have event()>=1)
			// bunchCrossing == 0 intime, buncCrossing!=0 offtime, standard generation has [-12,+3]
			if((*part)[i].eventId().event() ==0 and (*part)[i].eventId().bunchCrossing()==0) {

				// look for the generator particle, set to -1 for Geant produced paticles
				int tp_genpart=-1;
				if (!(*part)[i].genParticles().empty()) tp_genpart=(*part)[i].genParticle_begin().key();

				// default values for decay position is outside detector, i.e. ~stable
				int reachedEE=1;
				double fbrem=-1;
				float dvx=999.;
				float dvy=999.;
				float dvz=999.;
				if((*part)[i].decayVertices().size()>=1) { //they can be delta rays, in this case you have multiple decay verices
					dvx=(*part)[i].decayVertices()[0]->position().x();
					dvy=(*part)[i].decayVertices()[0]->position().y();
					dvz=(*part)[i].decayVertices()[0]->position().z();
					if ((*part)[i].decayVertices()[0]->inVolume()) reachedEE=0; //if it decays inside the tracker volume
				}

				// look for origin of particle inside the beampipe so that we can assess if a track can potentially give a reco::Track
				float r_origin=(*part)[i].parentVertex()->position().Pt();
				bool fromBeamPipe=true;
				if (r_origin>2.0) fromBeamPipe=false;
				genpart_eta.push_back((*part)[i].eta());
				genpart_phi.push_back((*part)[i].phi());
				genpart_pt.push_back((*part)[i].pt());
				genpart_energy.push_back((*part)[i].energy());
				genpart_dvx.push_back(dvx);
				genpart_dvy.push_back(dvy);
				genpart_dvz.push_back(dvz);
				genpart_fbrem.push_back(fbrem);
				genpart_pid.push_back((*part)[i].pdgId());
				genpart_gen.push_back(tp_genpart);
				genpart_reachedEE.push_back(reachedEE);
				genpart_fromBeamPipe.push_back(fromBeamPipe);
			}
		}
	} else
	{
		unsigned int npart = mySimEvent->nTracks();
		for (unsigned int i=0; i<npart; ++i) {
			std::vector<float> xp,yp,zp;
			FSimTrack &myTrack(mySimEvent->track(i));
			math::XYZTLorentzVectorD vtx(0,0,0,0);

			int reachedEE=0; // compute the extrapolations for the particles reaching EE and for the gen particles
			double fbrem=-1;
			if (myTrack.noEndVertex() || myTrack.genpartIndex()>=0)
			{

				RawParticle part(myTrack.momentum(),myTrack.vertex().position());
				part.setID(myTrack.type());
				BaseParticlePropagator myPropag(part,160,layerPositions[0],3.8);
				myPropag.propagate();
				unsigned result=myPropag.getSuccess();
				vtx=myPropag.propagated().vertex();
				unsigned nlayers=40;

				if (myTrack.noEndVertex()) {
					if (result==2 && vtx.Rho()> 25) {
						reachedEE=2;
						double dpt=0;

						for(int i=0; i<myTrack.nDaughters(); ++i)
							dpt+=myTrack.daughter(i).momentum().pt();
						if(abs(myTrack.type())==11)
							fbrem = dpt/myTrack.momentum().pt();
					}
					if (result==1) reachedEE=1;
				}
				else // in that case we propagate only to the first layers
					nlayers=1;

				for(unsigned il=0; il<nlayers; ++il) {
					myPropag.setPropagationConditions(140,layerPositions[il],false);
					if(il>0) // set PID 22 for a straight-line extrapolation after the 1st layer
						myPropag.setID(22);
					myPropag.propagate();
					RawParticle propParticle=myPropag.propagated();
					xp.push_back(propParticle.vertex().x());
					yp.push_back(propParticle.vertex().y());
					zp.push_back(propParticle.vertex().z());
				}
			}
			else
			{
				vtx = myTrack.endVertex().position();
			}

			// fill branches
			genpart_eta.push_back(myTrack.momentum().eta());
			genpart_phi.push_back(myTrack.momentum().phi());
			genpart_pt.push_back(myTrack.momentum().pt());
			genpart_energy.push_back(myTrack.momentum().energy());
			genpart_dvx.push_back(vtx.x());
			genpart_dvy.push_back(vtx.y());
			genpart_dvz.push_back(vtx.z());
			genpart_fbrem.push_back(fbrem);
			genpart_pid.push_back(myTrack.type());
			genpart_gen.push_back(myTrack.genpartIndex());
			genpart_reachedEE.push_back(reachedEE);
			genpart_fromBeamPipe.push_back(true);
		}
	}
	//make a map detid-rechit
	hitmap.clear();
	switch(algo) {
	case 1:
	{
		iEvent.getByToken(_recHitsEE,recHitHandleEE);
		iEvent.getByToken(_recHitsFH,recHitHandleFH);
		iEvent.getByToken(_recHitsBH,recHitHandleBH);
		const auto& rechitsEE = *recHitHandleEE;
		const auto& rechitsFH = *recHitHandleFH;
		const auto& rechitsBH = *recHitHandleBH;
		for(unsigned int i = 0; i < rechitsEE.size(); ++i) {
			hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
		}
		for(unsigned int i = 0; i < rechitsFH.size(); ++i) {
			hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
		}
		for(unsigned int i = 0; i < rechitsBH.size(); ++i) {
			hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
		}
		break;
	}
	case 2:
	{
		iEvent.getByToken(_recHitsEE,recHitHandleEE);
		const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
		for(unsigned int i = 0; i < rechitsEE.size(); i++) {
			hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
		}
		break;
	}
	case 3:
	{
		iEvent.getByToken(_recHitsFH,recHitHandleFH);
		iEvent.getByToken(_recHitsBH,recHitHandleBH);
		const auto& rechitsFH = *recHitHandleFH;
		const auto& rechitsBH = *recHitHandleBH;
		for(unsigned int i = 0; i < rechitsFH.size(); i++) {
			hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
		}
		for(unsigned int i = 0; i < rechitsBH.size(); i++) {
			hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
		}
		break;
	}
	default:
		break;
	}


	const reco::CaloClusterCollection &clusters = *clusterHandle;
	unsigned int nclus = clusters.size();
	cluster_index = 0;
	rechit_index = 0;
	storedLayerClusters.clear();
	storedRecHits.clear();
	for(unsigned int i = 0; i < multiClusters.size(); i++) {

		int cl2dSeed = 0;
		std::vector<unsigned int> cl2dIndices;

		for(reco::HGCalMultiCluster::component_iterator it = multiClusters[i].begin();
		    it!=multiClusters[i].end(); it++) {

			if((*it)->energy() > (*(it+cl2dSeed))->energy()) cl2dSeed = it - multiClusters[i].begin();
			cl2dIndices.push_back(cluster_index);

			fillLayerCluster(*it, true, i);

		}

		double pt = multiClusters[i].energy() / cosh(multiClusters[i].eta());

		multiclus_eta.push_back(multiClusters[i].eta());
		multiclus_phi.push_back(multiClusters[i].phi());
		multiclus_pt.push_back(pt);
		multiclus_energy.push_back(multiClusters[i].energy());
		multiclus_z.push_back(multiClusters[i].z());
		multiclus_slopeX.push_back(multiClusters[i].x());
		multiclus_slopeY.push_back(multiClusters[i].y());
		multiclus_cluster2d.push_back(cl2dIndices);
		multiclus_cl2dSeed.push_back(cl2dSeed);

	}

	// Fills the additional 2d layers
	for(unsigned ic=0; ic < nclus; ++ic )  {
		edm::Ptr<reco::BasicCluster> clusterPtr(clusterHandle,ic);
		if(storedLayerClusters.find(clusterPtr)==storedLayerClusters.end() ) {

			double pt = clusterPtr->energy() / cosh(clusterPtr->eta());
			if(pt>layerClusterPtThreshold) {

				fillLayerCluster(clusterPtr, rawRecHits);

			}
		}
	}

	// Fill remaining RecHits
	if (rawRecHits) {
		if (algo < 3) {
			const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
			// loop over EE RecHits
			for (HGCRecHitCollection::const_iterator it_hit = rechitsEE.begin(); it_hit < rechitsEE.end(); ++it_hit) {

				const HGCalDetId detid = it_hit->detid();
				unsigned int layer = recHitTools.getLayerWithOffset(detid);

				if(storedRecHits.find(detid)==storedRecHits.end() ) {
					fillRecHit(detid, -1, layer);
				}

			}
		}
		if (algo != 2) {
			const HGCRecHitCollection& rechitsFH = *recHitHandleFH;
			// loop over FH RecHits
			for (HGCRecHitCollection::const_iterator it_hit = rechitsFH.begin(); it_hit < rechitsFH.end(); ++it_hit) {

				const HGCalDetId detid = it_hit->detid();
				unsigned int layer = recHitTools.getLayerWithOffset(detid);

				if(storedRecHits.find(detid)==storedRecHits.end() ) {
					fillRecHit(detid, -1, layer);
				}

			}
			const HGCRecHitCollection& rechitsBH = *recHitHandleBH;
			// loop over BH RecHits
			for (HGCRecHitCollection::const_iterator it_hit = rechitsBH.begin(); it_hit < rechitsBH.end(); ++it_hit) {

				const HGCalDetId detid = it_hit->detid();
				unsigned int layer = recHitTools.getLayerWithOffset(detid);

				if(storedRecHits.find(detid)==storedRecHits.end() ) {
					fillRecHit(detid, -1, layer);
				}

			}
		}
	}

	// loop over simClusters
	for (std::vector<SimCluster>::const_iterator it_simClus = simClusters->begin(); it_simClus != simClusters->end(); ++it_simClus) {
		const std::vector<std::pair<uint32_t,float> > hits_and_fractions = it_simClus->hits_and_fractions();
		std::vector<uint32_t> hits;
		std::vector<float> fractions;
		std::vector<unsigned int> layers;
		std::vector<unsigned int> wafers;
		std::vector<unsigned int> cells;
		for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
			hits.push_back(it_haf->first);
			fractions.push_back(it_haf->second);
			layers.push_back(recHitTools.getLayerWithOffset(it_haf->first));
			if( DetId::Forward == DetId(it_haf->first).det() ) {
				wafers.push_back(recHitTools.getWafer(it_haf->first));
				cells.push_back(recHitTools.getCell(it_haf->first));
			} else {
				wafers.push_back(std::numeric_limits<unsigned int>::max());
				cells.push_back(std::numeric_limits<unsigned int>::max());
			}
		}

		simcluster_eta.push_back(it_simClus->eta());
		simcluster_phi.push_back(it_simClus->phi());
		simcluster_pt.push_back(it_simClus->pt());
		simcluster_energy.push_back(it_simClus->energy());
		simcluster_simEnergy.push_back(it_simClus->simEnergy());
		simcluster_hits.push_back(hits);
		simcluster_fractions.push_back(fractions);
		simcluster_layers.push_back(layers);
		simcluster_wafers.push_back(wafers);
		simcluster_cells.push_back(cells);

	} // end loop over simClusters

	// loop over pfClusters
	for (std::vector<reco::PFCluster>::const_iterator it_pfClus = pfClusters.begin(); it_pfClus != pfClusters.end(); ++it_pfClus) {

		pfcluster_eta.push_back(it_pfClus->eta());
		pfcluster_phi.push_back(it_pfClus->phi());
		pfcluster_pt.push_back(it_pfClus->pt());
		pfcluster_energy.push_back(it_pfClus->energy());

	} // end loop over pfClusters

	// loop over caloParticles
	if (readCaloParticles) {
		for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles->begin(); it_caloPart != caloParticles->end(); ++it_caloPart) {
			const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();
			std::vector<uint32_t> simClusterIndex;
			for (CaloParticle::sc_iterator it_sc = simClusterRefVector.begin(); it_sc != simClusterRefVector.end(); ++it_sc) {
				simClusterIndex.push_back((*it_sc).key());
			}

			calopart_eta.push_back(it_caloPart->eta());
			calopart_phi.push_back(it_caloPart->phi());
			calopart_pt.push_back(it_caloPart->pt());
			calopart_energy.push_back(it_caloPart->energy());
			calopart_simEnergy.push_back(it_caloPart->simEnergy());
			calopart_simClusterIndex.push_back(simClusterIndex);

		} // end loop over caloParticles
	}

	// loop over tracks

	//  random = new RandomEngineAndDistribution(iEvent.streamID());

	// prepare for RK propagation
	defaultRKPropagator::Product prod( aField, alongMomentum, 5.e-5);
	auto & RKProp = prod.propagator;


	for (std::vector<reco::Track>::const_iterator it_track = tracks.begin(); it_track != tracks.end(); ++it_track) {

		if (!it_track->quality(reco::Track::highPurity)) continue;

		double energy = it_track->pt() * cosh(it_track->eta());

		//save info about reconstructed tracks propoagation to hgcal layers (ony for pt>propagationPtThreshold tracks)
		std::vector<float> xp,yp,zp;

		if (it_track->pt()>= propagationPtThreshold) {

			// Define error matrix
			ROOT::Math::SMatrixIdentity id;
			AlgebraicSymMatrix55 C(id);
			C *= 0.01;
			CurvilinearTrajectoryError err(C);
			typedef TrajectoryStateOnSurface TSOS;

			GlobalPoint startingPosition(it_track->vx(),it_track->vy(),it_track->vz());
			GlobalVector startingMomentum(it_track->px(),it_track->py(),it_track->pz());

			Plane::PlanePointer startingPlane = Plane::build( Plane::PositionType (it_track->vx(),it_track->vy(),it_track->vz()), Plane::RotationType () );

			TSOS startingStateP(GlobalTrajectoryParameters(startingPosition,startingMomentum, it_track->charge(), aField), err, *startingPlane);

			for(unsigned il=0; il<layerPositions.size(); ++il) {
				float xp_curr=0;
				float yp_curr=0;
				float zp_curr=0;

				for (int zside = -1; zside <=1; zside+=2)
				{
					// clearly try both sides
					Plane::PlanePointer endPlane = Plane::build( Plane::PositionType (0,0,zside*layerPositions[il]), Plane::RotationType());
					try{
						/*
						   std::cout << "Trying from " <<
						   " layer " << il <<
						   " starting point " << startingStateP.globalPosition() <<
						   std::endl;
						 */
						TSOS trackStateP = RKProp.propagate( startingStateP, *endPlane);
						if (trackStateP.isValid()) {
							xp_curr=trackStateP.globalPosition().x();
							yp_curr=trackStateP.globalPosition().y();
							zp_curr=trackStateP.globalPosition().z();

							// std::cout << "Succesfully finished Positive track propagation  -------------- with RK: " << trackStateP.globalPosition() << std::endl;
						}
					}
					catch (...) {
						std::cout << "MagVolumeOutsideValidity not properly caught!! Lost this track " << std::endl;
					}
				}
				xp.push_back(xp_curr);
				yp.push_back(yp_curr);
				zp.push_back(zp_curr);
			} // closes loop on layers
		} // closes conditions pt>3

		// save info in tree
		track_pt.push_back(it_track->pt());
		track_eta.push_back(it_track->eta());
		track_phi.push_back(it_track->phi());
		track_energy.push_back(energy);
		track_charge.push_back(it_track->charge());
		track_posx.push_back(xp);
		track_posy.push_back(yp);
		track_posz.push_back(zp);

	} // end loop over tracks

  for (std::vector<reco::GenParticle>::const_iterator it_p = genParticlesHandle->begin(); it_p != genParticlesHandle->end(); ++it_p) {
    gen_eta.push_back(it_p->eta());
    gen_phi.push_back(it_p->phi());
    gen_pt.push_back(it_p->pt());
    gen_energy.push_back(it_p->energy());
    gen_pdgid.push_back(it_p->pdgId());
    gen_status.push_back(it_p->status());
  }

	ev_event = iEvent.id().event();
	ev_lumi = iEvent.id().luminosityBlock();
	ev_run = iEvent.id().run();

	vtx_x = vx;
	vtx_y = vy;
	vtx_z = vz;

	t->Fill();
}

void HGCalAnalysis::beginRun(edm::Run const& iEvent, edm::EventSetup const& es) {

	edm::ESHandle < HepPDT::ParticleDataTable > pdt;
	es.getData(pdt);
	mySimEvent->initializePdt(&(*pdt));

	retrieveLayerPositions(es,52);

	edm::ESHandle<MagneticField> magfield;
	es.get<IdealMagneticFieldRecord>().get(magfield);

	aField=&(*magfield);

}

void HGCalAnalysis::endRun(edm::Run const& iEvent, edm::EventSetup const&) {
}

void
HGCalAnalysis::beginJob()
{
	;
}

// ------------ method called once each job just after ending the event loop  ------------
void
HGCalAnalysis::endJob()
{
}

// ------------ method to be called once --------------------------------------------------


void HGCalAnalysis::retrieveLayerPositions(const edm::EventSetup& es, unsigned layers)
{
	recHitTools.getEventSetup(es);

	DetId id;
	for(unsigned ilayer=1; ilayer<=layers; ++ilayer) {
		if (ilayer<=28) id=HGCalDetId(ForwardSubdetector::HGCEE,1,ilayer,1,50,1);
		if (ilayer>28 && ilayer<=40) id=HGCalDetId(ForwardSubdetector::HGCHEF,1,ilayer-28,1,50,1);
		if (ilayer>40) id=HcalDetId(HcalSubdetector::HcalEndcap, 50, 100, ilayer-40);
		const GlobalPoint pos = recHitTools.getPosition(id);
		// std::cout << "GEOM " ;
		// std::cout << " layer " << ilayer << " " << pos.z() << std::endl;
		layerPositions.push_back(pos.z());
	}
}


void HGCalAnalysis::fillLayerCluster(const edm::Ptr<reco::CaloCluster>& layerCluster, const bool& fillRecHits, const int& multiClusterIndex) {

	// std::cout << "in fillLayerCluster" << std::endl;
	const std::vector< std::pair<DetId, float> > &hf = layerCluster->hitsAndFractions();
	std::vector<unsigned int> rhIndices;
	int ncoreHit = 0;
	int layer = 0;
	int rhSeed = 0;
	if (!fillRecHits) {
		rhSeed = -1;
	}
	float maxEnergy = -1.;

	for(unsigned int j = 0; j < hf.size(); j++) {
		//here we loop over detid/fraction pairs
		float fraction = hf[j].second;
		const DetId rh_detid = hf[j].first;
		layer = recHitTools.getLayerWithOffset(rh_detid);
		const HGCRecHit *hit = hitmap[rh_detid];
		ncoreHit += int(fraction);

		if (fillRecHits) {
			if(storedRecHits.find(rh_detid)==storedRecHits.end() ) {
				// std::cout << "in fillLayerCluster: RecHit not yet filled" << std::endl;
				// std::cout << "in fillLayerCluster: hit energy: " << hit->energy() << std::endl;
				// std::cout << "in fillLayerCluster: first hit energy: " << maxEnergy << std::endl;
				if(hit->energy() > maxEnergy) {
					rhSeed = rechit_index;
					maxEnergy = hit->energy();
				}
				rhIndices.push_back(rechit_index);
				fillRecHit(rh_detid, fraction, layer, cluster_index);
			}
			else {
				// need to see what to do about existing rechits in case of sharing
				std::cout << "RecHit already filled for different layer cluster: " << int(rh_detid) << std::endl;
			}
		}
	}

	double pt = layerCluster->energy() / cosh(layerCluster->eta());

	cluster2d_eta.push_back(layerCluster->eta());
	cluster2d_phi.push_back(layerCluster->phi());
	cluster2d_pt.push_back(pt);
	cluster2d_energy.push_back(layerCluster->energy());
	cluster2d_x.push_back(layerCluster->x());
	cluster2d_y.push_back(layerCluster->y());
	cluster2d_z.push_back(layerCluster->z());
	cluster2d_layer.push_back(layer);
	cluster2d_nhitCore.push_back(ncoreHit);
	cluster2d_nhitAll.push_back(hf.size());
	cluster2d_multicluster.push_back(multiClusterIndex);
	cluster2d_rechitSeed.push_back(rhSeed);
	cluster2d_rechits.push_back(rhIndices);

	storedLayerClusters.insert(layerCluster);
	++cluster_index;

}


void HGCalAnalysis::fillRecHit(const DetId& detid, const float& fraction, const unsigned int& layer, const int& cluster_index) {

	// std::cout << "in fillRecHit" << std::endl;
	int flags = 0x0;
	if (fraction>0. && fraction<1.)
		flags = 0x1;
	else if(fraction<0.)
		flags = 0x3;
	else if(fraction==0.)
		flags = 0x2;
	const HGCRecHit *hit = hitmap[detid];

	const GlobalPoint position = recHitTools.getPosition(detid);
	const unsigned int wafer = ( DetId::Forward == DetId(detid).det() ? recHitTools.getWafer(detid) : std::numeric_limits<unsigned int>::max() );
	const unsigned int cell  = ( DetId::Forward == DetId(detid).det() ? recHitTools.getCell(detid) :  std::numeric_limits<unsigned int>::max() );
	const double cellThickness = ( DetId::Forward == DetId(detid).det() ? recHitTools.getSiThickness(detid) : std::numeric_limits<std::float_t>::max() );
	const bool isHalfCell = recHitTools.isHalfCell(detid);
	const double eta = recHitTools.getEta(position, vz);
	const double phi = recHitTools.getPhi(position);
	const double pt = recHitTools.getPt(position, hit->energy(), vz);


	// fill the vectors
	rechit_eta.push_back(eta);
	rechit_phi.push_back(phi);
	rechit_pt.push_back(pt);
	rechit_energy.push_back(hit->energy());
	rechit_layer.push_back(layer);
	rechit_wafer.push_back(wafer);
	rechit_cell.push_back(cell);
	rechit_detid.push_back(detid);
	rechit_x.push_back(position.x());
	rechit_y.push_back(position.y());
	rechit_z.push_back(position.z());
	rechit_time.push_back(hit->time());
	rechit_thickness.push_back(cellThickness);
	rechit_isHalf.push_back(isHalfCell);
	rechit_flags.push_back(flags);
	rechit_cluster2d.push_back(cluster_index);

	storedRecHits.insert(detid);
	++rechit_index;

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void
HGCalAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

/*
   Surface::RotationType HGCalAnalysis::rotation( const GlobalVector& zDir) const
   {
   GlobalVector zAxis = zDir.unit();
   GlobalVector yAxis( zAxis.y(), -zAxis.x(), 0);
   GlobalVector xAxis = yAxis.cross( zAxis);
   return Surface::RotationType( xAxis, yAxis, zAxis);
   }
 */

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalAnalysis);
