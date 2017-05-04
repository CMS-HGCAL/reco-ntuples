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
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
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

#include "RecoNtuples/HGCalAnalysis/interface/AEvent.h"
#include "RecoNtuples/HGCalAnalysis/interface/AObData.h"

#include <string>
#include <map>
#include <vector>
#include <set>

class HGCalAnalysis : public  edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  HGCalAnalysis();
  explicit HGCalAnalysis(const edm::ParameterSet&);
  ~HGCalAnalysis();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  virtual void beginRun(edm::Run const& iEvent, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const& iEvent, edm::EventSetup const&) override ;

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  void retrieveLayerPositions(const edm::EventSetup&, unsigned layers) ;

  //  Surface::RotationType rotation( const GlobalVector& zDir) const;
  
  // ---------parameters ----------------------------
  bool readOfficialReco;
  bool readCaloParticles;
  double layerClusterPtThreshold;
  double propagationPtThreshold;
  std::string                detector;
  bool                       rawRecHits;

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

  TTree                     *tree;
  AEvent                    *event;
  AGenPartCollection        *agpc;
  ARecHitCollection         *arhc;
  ARecHitCollection         *arhc_raw;
  ACluster2dCollection      *acdc;
  AMultiClusterCollection   *amcc;
  ASimClusterCollection     *ascc;
  APFClusterCollection      *apfcc;
  ACaloParticleCollection   *acpc;
  ATrackCollection          *atrc;
  int                        algo;
  HGCalDepthPreClusterer     pre;
  hgcal::RecHitTools         recHitTools;

  // -------convenient tool to deal with simulated tracks
  FSimEvent * mySimEvent;
  edm::ParameterSet particleFilter;
  std::vector <float> layerPositions;

  // and also the magnetic field
  MagneticField const * aField;
};

HGCalAnalysis::HGCalAnalysis() {;}

HGCalAnalysis::HGCalAnalysis(const edm::ParameterSet& iConfig) :
  readOfficialReco(iConfig.getParameter<bool>("readOfficialReco")),
  readCaloParticles(iConfig.getParameter<bool>("readCaloParticles")),
  layerClusterPtThreshold(iConfig.getParameter<double>("layerClusterPtThreshold")),
  propagationPtThreshold(iConfig.getUntrackedParameter<double>("propagationPtThreshold",3.0)),
  detector(iConfig.getParameter<std::string >("detector")),
  rawRecHits(iConfig.getParameter<bool>("rawRecHits")),
  particleFilter(iConfig.getParameter<edm::ParameterSet>("TestParticleFilter"))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
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
  _pfClusters = consumes<std::vector<reco::PFCluster> >(edm::InputTag("particleFlowClusterHGCal"));
  _multiClusters = consumes<std::vector<reco::HGCalMultiCluster> >(edm::InputTag("hgcalLayerClusters"));
  _tracks = consumes<std::vector<reco::Track> >(edm::InputTag("generalTracks"));



  edm::Service<TFileService> fs;
  fs->make<TH1F>("total", "total", 100, 0, 5.);
  tree = new TTree("hgc","Analysis");
  agpc = new AGenPartCollection();
  arhc = new ARecHitCollection();
  if (rawRecHits)
    arhc_raw = new ARecHitCollection();
  acdc = new ACluster2dCollection();
  amcc = new AMultiClusterCollection();
  ascc = new ASimClusterCollection();
  apfcc = new APFClusterCollection();
  acpc = new ACaloParticleCollection();
  atrc = new ATrackCollection();
  event = new AEvent();
  tree->Branch("event","AEvent",&event,16000,99);
  tree->Branch("particles","AGenPartCollection",&agpc,16000,0);
  tree->Branch("rechits","ARecHitCollection",&arhc,16000,0);
  if (rawRecHits)
    tree->Branch("rechits_raw","ARecHitCollection",&arhc_raw,16000,0);
  tree->Branch("cluster2d","ACluster2dCollection",&acdc,16000,0);
  tree->Branch("multicluster","AMultiClusterCollection",&amcc,16000,0);
  tree->Branch("simcluster","ASimClusterCollection",&ascc,16000,0);
  tree->Branch("pfcluster","APFClusterCollection",&apfcc,16000,0);
  tree->Branch("caloparticles","ACaloParticleCollection",&acpc,16000,0);
  tree->Branch("tracks","ATrackCollection", &atrc, 16000, 0);

  
}
HGCalAnalysis::~HGCalAnalysis()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void

HGCalAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  agpc->clear();
  arhc->clear();
  arhc_raw->clear();
  acdc->clear();
  amcc->clear();
  ascc->clear();
  apfcc->clear();
  acpc->clear();
  atrc->clear();

  ParticleTable::Sentry ptable(mySimEvent->theTable());
  recHitTools.getEventSetup(iSetup);

  unsigned int npart = 0;
  unsigned int nhit  = 0;
  unsigned int nhit_raw = 0;
  unsigned int nclus = 0;
  unsigned int nmclus = 0;
  unsigned int nsimclus = 0;
  unsigned int npfclus = 0;
  unsigned int ncalopart = 0;
  unsigned int ntracks = 0;

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleFH;
  Handle<HGCRecHitCollection> recHitHandleBH;
  Handle<reco::CaloClusterCollection> clusterHandle;

  iEvent.getByToken(_clusters,clusterHandle);
  Handle<std::vector<TrackingVertex> > vtxHandle;
  Handle<std::vector<TrackingParticle> > partHandle;
  //  const std::vector<TrackingVertex>* vtxs;
  const std::vector<TrackingParticle> * part;

  Handle<edm::HepMCProduct> hevH;
  Handle<std::vector<SimTrack> >simTracksHandle;
  Handle<std::vector<SimVertex> >simVerticesHandle;

  iEvent.getByToken(_hev,hevH);
  if(!readOfficialReco) {
    iEvent.getByToken(_vtx,vtxHandle);
    iEvent.getByToken(_part,partHandle);
    //    vtxs = &(*vtxHandle);
    part = &(*partHandle);
  } else  // use SimTracks and HepMCProduct
    {
      iEvent.getByToken(_simTracks,simTracksHandle);
      iEvent.getByToken(_simVertices,simVerticesHandle);
      //      std::cout << " Filling FSimEvent " << simTracksHandle->size() << " " << simVerticesHandle->size() << std::endl;
      //      for(unsigned i=0; i<simTracksHandle->size();++i)
      //	std::cout << "i " << (*simTracksHandle)[i].type() << std::endl;
      mySimEvent->fill(*simTracksHandle,*simVerticesHandle);

    }


  Handle<std::vector<SimCluster> > simClusterHandle;
  Handle<std::vector<reco::PFCluster> > pfClusterHandle;
  Handle<std::vector<CaloParticle> > caloParticleHandle;
  Handle<std::vector<reco::Track> > trackHandle;

  const std::vector<SimCluster> * simClusters = 0 ;
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
  float vz = primaryVertex->position().z()/10.;

  if( !readOfficialReco) {
    npart = part->size();
    for(unsigned int i=0;i<npart;++i){

	// event=0 is the hard scattering (all PU have event()>=1)
	// bunchCrossing == 0 intime, buncCrossing!=0 offtime, standard generation has [-12,+3]
        if((*part)[i].eventId().event() ==0 and (*part)[i].eventId().bunchCrossing()==0){

	  //look for the generator particle, set to -1 for Geant produced paticles
	  int tp_genpart=-1;
	  if (! (*part)[i].genParticles().empty()) tp_genpart=(*part)[i].genParticle_begin().key();

	  // default values for decay position is outside detector, i.e. ~stable 
	  int reachedEE=1;
	  float dvx=999.;
	  float dvy=999.;
	  float dvz=999.;
	  if((*part)[i].decayVertices().size()>=1){ //they can be delta rays, in this case you have multiple decay verices
	    dvx=(*part)[i].decayVertices()[0]->position().x();
	    dvy=(*part)[i].decayVertices()[0]->position().y();
	    dvz=(*part)[i].decayVertices()[0]->position().z();
	    if ((*part)[i].decayVertices()[0]->inVolume()) reachedEE=0; //if it decays inside the tracker volume
	  }
	  
	  //look for origin of particle inside the beampipe so that we can assess if a track can potentially give a reco::Track
	  float r_origin=(*part)[i].parentVertex()->position().Pt();
	  bool fromBeamPipe=true;
	  if (r_origin>2.0) fromBeamPipe=false;
	  agpc->push_back(AGenPart((*part)[i].eta(),(*part)[i].phi(),(*part)[i].pt(),(*part)[i].energy(),dvx,dvy,dvz,-1.,(*part)[i].pdgId(),tp_genpart,reachedEE,fromBeamPipe));
	}
    }
  } else
    {
      npart = mySimEvent->nTracks();
      for (unsigned int i=0;i<npart ; ++i) {
	std::vector<float> xp,yp,zp;
	FSimTrack &myTrack(mySimEvent->track(i));
	math::XYZTLorentzVectorD vtx(0,0,0,0);

	int reachedEE=0; // compute the extrapolations for the particles reaching EE and for the gen particles
	double fbrem=-1.;
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
		
		for(int i=0;i<myTrack.nDaughters();++i)
		  dpt+=myTrack.daughter(i).momentum().pt();
		if(abs(myTrack.type())==11)
		  fbrem = dpt/myTrack.momentum().pt();
	      }
	      if (result==1) reachedEE=1;
	    }
	    else  // in that case we propagate only to the first layers
	      nlayers=1;
	    
	    for(unsigned il=0;il<nlayers;++il) {
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

	AGenPart part(myTrack.momentum().eta(),myTrack.momentum().phi(),myTrack.momentum().pt(),myTrack.momentum().energy(),vtx.x(),vtx.y(),vtx.z(),fbrem,myTrack.type(),myTrack.genpartIndex(),reachedEE);
	part.setExtrapolations(xp,yp,zp);
	agpc->push_back(part);
      }
    }
  //make a map detid-rechit
  std::map<DetId,const HGCRecHit*> hitmap;
  switch(algo){
  case 1:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      iEvent.getByToken(_recHitsFH,recHitHandleFH);
      iEvent.getByToken(_recHitsBH,recHitHandleBH);
      const auto& rechitsEE = *recHitHandleEE;
      const auto& rechitsFH = *recHitHandleFH;
      const auto& rechitsBH = *recHitHandleBH;
      for(unsigned int i = 0; i < rechitsEE.size(); ++i){
	hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
      }
      for(unsigned int i = 0; i < rechitsFH.size(); ++i){
	hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for(unsigned int i = 0; i < rechitsBH.size(); ++i){
	hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
      }
      break;
    }
  case 2:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
      for(unsigned int i = 0; i < rechitsEE.size(); i++){
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
      for(unsigned int i = 0; i < rechitsFH.size(); i++){
	hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for(unsigned int i = 0; i < rechitsBH.size(); i++){
	hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
      }
      break;
    }
  default:
    break;
  }

  // dump raw RecHits
  if (rawRecHits) {
     if (algo < 3) {
        const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
        // loop over EE RecHits
        for (HGCRecHitCollection::const_iterator it_hit = rechitsEE.begin(); it_hit < rechitsEE.end(); ++it_hit) {
            int clusterIndex = -1;
            int flags = 0x0;

	    const HGCalDetId detid = it_hit->detid();
	    unsigned int layer = recHitTools.getLayerWithOffset(detid);
            const GlobalPoint position = recHitTools.getPosition(it_hit->detid());
            const unsigned int wafer = recHitTools.getWafer(detid);
            const unsigned int cell  = recHitTools.getCell(detid);
            const double cellThickness = recHitTools.getSiThickness(detid);
            const bool isHalfCell = recHitTools.isHalfCell(detid);
            const double eta = recHitTools.getEta(position, vz);
            const double phi = recHitTools.getPhi(position);
            const double pt = recHitTools.getPt(position, it_hit->energy(), vz);

        	++nhit_raw;

            arhc_raw->push_back(ARecHit(layer, wafer, cell, detid,
                        position.x(), position.y(), position.z(),
                        eta,phi,pt,
                        it_hit->energy(), it_hit->time(), cellThickness,
                        isHalfCell, flags, clusterIndex));

        }
     }
     if (algo != 2) {
         const HGCRecHitCollection& rechitsFH = *recHitHandleFH;
         // loop over EE RecHits
         for (HGCRecHitCollection::const_iterator it_hit = rechitsFH.begin(); it_hit < rechitsFH.end(); ++it_hit) {
             int clusterIndex = -1;
             int flags = 0x0;

         	 const HGCalDetId detid = it_hit->detid();
         	 unsigned int layer = recHitTools.getLayerWithOffset(detid);
             const GlobalPoint position = recHitTools.getPosition(it_hit->detid());
             const unsigned int wafer = recHitTools.getWafer(detid);
             const unsigned int cell  = recHitTools.getCell(detid);
             const double cellThickness = recHitTools.getSiThickness(detid);
             const bool isHalfCell = recHitTools.isHalfCell(detid);
             const double eta = recHitTools.getEta(position, vz);
             const double phi = recHitTools.getPhi(position);
             const double pt = recHitTools.getPt(position, it_hit->energy(), vz);

         	++nhit_raw;

             arhc_raw->push_back(ARecHit(layer, wafer, cell, detid,
                         position.x(), position.y(), position.z(),
                         eta,phi,pt,
                         it_hit->energy(), it_hit->time(), cellThickness,
                         isHalfCell, flags, clusterIndex));

         }
	 const HGCRecHitCollection& rechitsBH = *recHitHandleBH;
         // loop over EE RecHits
         for (HGCRecHitCollection::const_iterator it_hit = rechitsBH.begin(); it_hit < rechitsBH.end(); ++it_hit) {
             int clusterIndex = -1;
             int flags = 0x0;

	     const HcalDetId detid = it_hit->detid();
	     unsigned int layer = recHitTools.getLayerWithOffset(detid);
             const GlobalPoint position = recHitTools.getPosition(it_hit->detid());
             const unsigned int wafer = std::numeric_limits<unsigned int>::max();//recHitTools.getWafer(detid);
             const unsigned int cell  = std::numeric_limits<unsigned int>::max();//recHitTools.getCell(detid);
             const double cellThickness = std::numeric_limits<std::float_t>::max();//recHitTools.getSiThickness(detid);
             const bool isHalfCell = recHitTools.isHalfCell(detid);
             const double eta = recHitTools.getEta(position, vz);
             const double phi = recHitTools.getPhi(position);
             const double pt = recHitTools.getPt(position, it_hit->energy(), vz);

         	++nhit_raw;

             arhc_raw->push_back(ARecHit(layer, wafer, cell, detid,
                         position.x(), position.y(), position.z(),
                         eta,phi,pt,
                         it_hit->energy(), it_hit->time(), cellThickness,
                         isHalfCell, flags, clusterIndex));

         }
     }
  }

  const reco::CaloClusterCollection &clusters = *clusterHandle;
  nclus = clusters.size();
  nmclus = multiClusters.size();
  unsigned int cluster_index = 0;
  // to keep track of the 2d clusters stored within the loop on multiclusters
  std::set<edm::Ptr<reco::BasicCluster> > storedLayerClusters;
  for(unsigned int i = 0; i < multiClusters.size(); i++){
      int cl2dSeed = 0;
    for(reco::HGCalMultiCluster::component_iterator it = multiClusters[i].begin();
	it!=multiClusters[i].end(); it++){
      if((*it)->energy() > (*(it+cl2dSeed))->energy()) cl2dSeed = it - multiClusters[i].begin();

      const std::vector< std::pair<DetId, float> > &hf = (*it)->hitsAndFractions();
      int ncoreHit = 0;
      int layer = 0;
      int rhSeed = 0;

      for(unsigned int j = 0; j < hf.size(); j++){
	//here we loop over detid/fraction pairs
	ncoreHit += int(hf[j].second);
	int flags = 0x0;
	if (hf[j].second>0. && hf[j].second<1.)
	  flags = 0x1;
	else if(hf[j].second<0.)
	  flags = 0x3;
	else if(hf[j].second==0.)
	  flags = 0x2;
	const HGCRecHit *hit = hitmap[hf[j].first];
	layer = recHitTools.getLayerWithOffset(hf[j].first);

	const GlobalPoint position = recHitTools.getPosition(hf[j].first);
	const unsigned int detid = hf[j].first;
	const unsigned int wafer = ( DetId::Forward == DetId(hf[j].first).det() ? recHitTools.getWafer(hf[j].first) : std::numeric_limits<unsigned int>::max() );
	const unsigned int cell  = ( DetId::Forward == DetId(hf[j].first).det() ? recHitTools.getCell(hf[j].first) :  std::numeric_limits<unsigned int>::max() );
	const double cellThickness = ( DetId::Forward == DetId(hf[j].first).det() ? recHitTools.getSiThickness(hf[j].first) : std::numeric_limits<std::float_t>::max() );
	const bool isHalfCell = recHitTools.isHalfCell(hf[j].first);
	const double eta = recHitTools.getEta(position, vz);
	const double phi = recHitTools.getPhi(position);
	const double pt = recHitTools.getPt(position, hit->energy(), vz);

	if(hit->energy() > hitmap[hf[rhSeed].first]->energy()) rhSeed = j;
	++nhit;

	arhc->push_back(ARecHit(layer, wafer, cell, detid,
				position.x(), position.y(), position.z(),
				eta,phi,pt,
				hit->energy(), hit->time(), cellThickness,
				isHalfCell, flags, cluster_index));

      }
      double pt = (*it)->energy() / cosh((*it)->eta());
      acdc->push_back(ACluster2d((*it)->x(),(*it)->y(),(*it)->z(),
				 (*it)->eta(),(*it)->phi(),pt,(*it)->energy(),
				 layer, ncoreHit,hf.size(),i, rhSeed));
      storedLayerClusters.insert((*it));
      cluster_index++;
    }

    double pt = multiClusters[i].energy() / cosh(multiClusters[i].eta());
    amcc->push_back(AMultiCluster(multiClusters[i].eta(),
				  multiClusters[i].phi(),
                                  pt,
				  multiClusters[i].z(),
				  multiClusters[i].x(),
				  multiClusters[i].y(),
                                  multiClusters[i].energy(),
				  multiClusters[i].size(),
				  cl2dSeed));

  }

  // Fills the additional 2d layers
  for(unsigned ic=0 ; ic < nclus ; ++ic )  {
    edm::Ptr<reco::BasicCluster> clusterPtr(clusterHandle,ic);
    if(storedLayerClusters.find(clusterPtr)==storedLayerClusters.end() ) {
      double pt = clusterPtr->energy() / cosh(clusterPtr->eta());
      if(pt>layerClusterPtThreshold) {
	int layer = recHitTools.getLayerWithOffset(clusterPtr->hitsAndFractions()[0].first);
	acdc->push_back(ACluster2d(clusterPtr->x(),clusterPtr->y(),clusterPtr->z(),
				   clusterPtr->eta(),clusterPtr->phi(),pt,clusterPtr->energy(),
				   layer, 0, (int)clusterPtr->hitsAndFractions().size(),-1, -1));
      }
    }
  }

  // loop over simClusters
  for (std::vector<SimCluster>::const_iterator it_simClus = simClusters->begin(); it_simClus != simClusters->end(); ++it_simClus) {
    ++nsimclus;
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
    ascc->push_back(ASimCluster(it_simClus->pt(),
  			  it_simClus->eta(),
  			  it_simClus->phi(),
  			  it_simClus->energy(),
  			  it_simClus->simEnergy(),
  			  hits,
  			  fractions,
  			  layers,
  			  wafers,
  			  cells));
  } // end loop over simClusters

  // loop over pfClusters
  for (std::vector<reco::PFCluster>::const_iterator it_pfClus = pfClusters.begin(); it_pfClus != pfClusters.end(); ++it_pfClus) {
    ++npfclus;
    apfcc->push_back(APFCluster(it_pfClus->pt(),
                  it_pfClus->eta(),
                  it_pfClus->phi(),
                  it_pfClus->energy()));

  } // end loop over pfClusters

  // loop over caloParticles
  if (readCaloParticles) {
    for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles->begin(); it_caloPart != caloParticles->end(); ++it_caloPart) {
      ++ncalopart;
      const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();
      std::vector<uint32_t> simClusterIndex;
      for (CaloParticle::sc_iterator it_sc = simClusterRefVector.begin(); it_sc != simClusterRefVector.end(); ++it_sc) {
        simClusterIndex.push_back((*it_sc).key());
      }
      acpc->push_back(ACaloParticle(it_caloPart->pt(),
				    it_caloPart->eta(),
				    it_caloPart->phi(),
				    it_caloPart->energy(),
				    it_caloPart->simEnergy(),
				    simClusterIndex));

    } // end loop over caloParticles
  }

  // loop over tracks

  //  random = new RandomEngineAndDistribution(iEvent.streamID());

  // prepare for RK propagation
  defaultRKPropagator::Product  prod( aField, alongMomentum, 5.e-5); 
  auto & RKProp = prod.propagator;
          

  for (std::vector<reco::Track>::const_iterator it_track = tracks.begin(); it_track != tracks.end(); ++it_track) {
    
    if (! it_track->quality(reco::Track::highPurity)) continue;

    ++ntracks;
    double energy = it_track->pt() * cosh(it_track->eta());
    
    //    atrc->push_back(ATrack(it_track->pt(),
    //			     it_track->eta(),
    //			     it_track->phi(),
    //			     energy));
    

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
      
      for(unsigned il=0;il<layerPositions.size();++il) {
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
	    catch (...){
	      std::cout << "MagVolumeOutsideValidity not properly caught!! Lost this track " << std::endl;
	    }
	  }
	xp.push_back(xp_curr);
	yp.push_back(yp_curr);
	zp.push_back(zp_curr);
      } // closes loop on layers
    } // closes conditions pt>3

    // save info in tree
    ATrack thistrack(it_track->pt(),it_track->eta(),it_track->phi(),energy,it_track->charge());
    thistrack.setExtrapolations(xp,yp,zp);
    atrc->push_back(thistrack);

  } // end loop over tracks

  event->set(iEvent.run(),iEvent.id().event(),npart,nhit,nhit_raw,nclus,nmclus,
         nsimclus, npfclus, ncalopart, ntracks,
	     vx,vy,vz);
  tree->Fill();
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

void HGCalAnalysis::endRun(edm::Run const& iEvent, edm::EventSetup const&) {}

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
  for(unsigned ilayer=1;ilayer<=layers;++ilayer) {
    if (ilayer<=28) id=HGCalDetId(ForwardSubdetector::HGCEE,1,ilayer,1,2,1);
    if (ilayer>28 && ilayer<=40) id=HGCalDetId(ForwardSubdetector::HGCHEF,1,ilayer-28,1,2,1);
    if (ilayer>40) id=HcalDetId(HcalSubdetector::HcalEndcap, 50, 100, ilayer-40);
    const GlobalPoint pos = recHitTools.getPosition(id);
    // std::cout << "GEOM " ;
    // std::cout << " layer " << ilayer << " " << pos.z() << std::endl;
    layerPositions.push_back(pos.z());
  }
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
