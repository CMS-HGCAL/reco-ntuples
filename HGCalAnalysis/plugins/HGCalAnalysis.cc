// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"



#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"



#include "TTree.h"
#include "TH1F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"



#include "RecoLocalCalo/HGCalRecHitDump/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecHitDump/interface/HGCalMultiCluster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "RecoNtuples/HGCalAnalysis/interface/AEvent.h"
#include "RecoNtuples/HGCalAnalysis/interface/AObData.h"

#include <string>
#include <map>

class HGCalAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  explicit HGCalAnalysis(const edm::ParameterSet&);
  ~HGCalAnalysis();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

 // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> _recHitsEE;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsHE;
  edm::EDGetTokenT<reco::CaloClusterCollection> _clusters;
  edm::EDGetTokenT<std::vector<TrackingVertex> > _vtx;
  edm::EDGetTokenT<std::vector<TrackingParticle> > _part;
  edm::EDGetTokenT<std::vector<SimCluster> > _simClusters;
  edm::EDGetTokenT<std::vector<reco::PFCluster> > _pfClusters;
  edm::EDGetTokenT<std::vector<CaloParticle> > _caloParticles;

  TTree                     *tree;
  AEvent                    *event;
  AGenPartCollection        *agpc;
  ARecHitCollection         *arhc;
  ACluster2dCollection      *acdc;
  AMultiClusterCollection   *amcc;
  ASimClusterCollection     *ascc;
  APFClusterCollection      *apfcc;
  ACaloParticleCollection   *acpc;
  std::string                detector;
  int                        algo;
  HGCalDepthPreClusterer     pre;
  hgcal::RecHitTools         recHitTools;
};



HGCalAnalysis::HGCalAnalysis(const edm::ParameterSet& iConfig) :
  detector(iConfig.getParameter<std::string >("detector")),
  pre(iConfig.getParameter<double>("depthClusteringCone"))
{
  //now do what ever initialization is needed
  usesResource("TFileService");


  if(detector=="both"){
    _recHitsEE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"));
    _recHitsHE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"));
    _clusters = consumes<reco::CaloClusterCollection>(edm::InputTag("imagingClusterHGCal"));
    algo = 1;
  }else if(detector=="EE"){
    _recHitsEE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"));
    _clusters = consumes<reco::CaloClusterCollection>(edm::InputTag("imagingClusterHGCal"));
    algo = 2;
  }else{
    _recHitsHE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"));
    _clusters = consumes<reco::CaloClusterCollection>(edm::InputTag("imagingClusterHGCal"));
    algo = 3;
  }
  _vtx = consumes<std::vector<TrackingVertex> >(edm::InputTag("mix","MergedTrackTruth"));
  _part = consumes<std::vector<TrackingParticle> >(edm::InputTag("mix","MergedTrackTruth"));
  _simClusters = consumes<std::vector<SimCluster> >(edm::InputTag("mix","MergedCaloTruth"));
  _pfClusters = consumes<std::vector<reco::PFCluster> >(edm::InputTag("particleFlowClusterHGCal"));
  _caloParticles = consumes<std::vector<CaloParticle> >(edm::InputTag("mix","MergedCaloTruth"));

  edm::Service<TFileService> fs;
  fs->make<TH1F>("totale", "totale", 100, 0, 5.);
  tree = new TTree("hgc","Analysis");
  agpc = new AGenPartCollection();
  arhc = new ARecHitCollection();
  acdc = new ACluster2dCollection();
  amcc = new AMultiClusterCollection();
  ascc = new ASimClusterCollection();
  apfcc = new APFClusterCollection();
  acpc = new ACaloParticleCollection();
  event = new AEvent();
  tree->Branch("event","AEvent",&event,16000,99);
  tree->Branch("particles","AGenPartCollection",&agpc,16000,0);
  tree->Branch("rechits","ARecHitCollection",&arhc,16000,0);
  tree->Branch("cluster2d","ACluster2dCollection",&acdc,16000,0);
  tree->Branch("multicluster","AMultiClusterCollection",&amcc,16000,0);
  tree->Branch("simcluster","ASimClusterCollection",&ascc,16000,0);
  tree->Branch("pfcluster","APFClusterCollection",&apfcc,16000,0);
  tree->Branch("caloparticles","ACaloParticleCollection",&acpc,16000,0);
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
  acdc->clear();
  amcc->clear();
  ascc->clear();
  apfcc->clear();
  acpc->clear();


  recHitTools.getEventSetup(iSetup);

  int npart = 0;
  int nhit  = 0;
  int nclus = 0;
  int nmclus = 0;
  int nsimclus = 0;
  int npfclus = 0;
  int ncalopart = 0;

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleHE;
  Handle<reco::CaloClusterCollection> clusterHandle;

  std::vector<HGCalMultiCluster> multiClusters;
  iEvent.getByToken(_clusters,clusterHandle);
  Handle<std::vector<TrackingVertex> > vtxHandle;
  Handle<std::vector<TrackingParticle> > partHandle;
  iEvent.getByToken(_vtx,vtxHandle);
  iEvent.getByToken(_part,partHandle);
  const std::vector<TrackingVertex>& vtxs = *vtxHandle;
  const std::vector<TrackingParticle>& part = *partHandle;

  Handle<std::vector<SimCluster> > simClusterHandle;
  Handle<std::vector<reco::PFCluster> > pfClusterHandle;
  Handle<std::vector<CaloParticle> > caloParticleHandle;

  iEvent.getByToken(_simClusters, simClusterHandle);
  iEvent.getByToken(_pfClusters, pfClusterHandle);
  iEvent.getByToken(_caloParticles, caloParticleHandle);

  const std::vector<SimCluster>& simClusters = *simClusterHandle;
  const std::vector<reco::PFCluster>& pfClusters = *pfClusterHandle;
  const std::vector<CaloParticle>& caloParticles = *caloParticleHandle;

  float vx = 0.;
  float vy = 0.;
  float vz = 0.;
  if(vtxs.size()!=0){
    vx = vtxs[0].position().x();
    vy = vtxs[0].position().y();
    vz = vtxs[0].position().z();
  }
  // TODO: should fall back to beam spot if no vertex
  npart = part.size();
  for(unsigned int i=0;i<part.size();++i){
    if(part[i].parentVertex()->nGenVertices()>0){
      float dvx=0.;
      float dvy=0.;
      float dvz=0.;
      if(part[i].decayVertices().size()==1){
	 dvx=part[i].decayVertices()[0]->position().x();
	 dvy=part[i].decayVertices()[0]->position().y();
	 dvz=part[i].decayVertices()[0]->position().z();
      }
      agpc->push_back(AGenPart(part[i].eta(),part[i].phi(),part[i].pt(),part[i].energy(),dvx,dvy,dvz,part[i].pdgId()));
    }
  }
  //make a map detid-rechit
  std::map<HGCalDetId,const HGCRecHit*> hitmap;
  switch(algo){
  case 1:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      iEvent.getByToken(_recHitsHE,recHitHandleHE);
      const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
      const HGCRecHitCollection& rechitsHE = *recHitHandleHE;
      for(unsigned int i = 0; i < rechitsEE.size(); i++){
	hitmap[HGCalDetId(rechitsEE[i].detid())] = &rechitsEE[i];
      }
      for(unsigned int i = 0; i < rechitsHE.size(); i++){
	hitmap[HGCalDetId(rechitsHE[i].detid())] = &rechitsHE[i];
      }
      break;
    }
  case 2:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
      for(unsigned int i = 0; i < rechitsEE.size(); i++){
	hitmap[HGCalDetId(rechitsEE[i].detid())] = &rechitsEE[i];
      }
      break;
    }
  case 3:
    {
      iEvent.getByToken(_recHitsHE,recHitHandleHE);
      const HGCRecHitCollection& rechitsHE = *recHitHandleHE;
      for(unsigned int i = 0; i < rechitsHE.size(); i++){
	hitmap[HGCalDetId(rechitsHE[i].detid())] = &rechitsHE[i];
      }
      break;
    }
  default:
    break;
  }




  const reco::CaloClusterCollection &clusters = *clusterHandle;
  nclus = clusters.size();
  multiClusters = pre.makePreClusters(clusters);
  nmclus = multiClusters.size();
  unsigned int cluster_index = 0;
  for(unsigned int i = 0; i < multiClusters.size(); i++){
      int cl2dSeed = 0;
    for(HGCalMultiCluster::component_iterator it = multiClusters[i].begin();
	it!=multiClusters[i].end(); it++){
      if(it->energy() > (it+cl2dSeed)->energy()) cl2dSeed = it - multiClusters[i].begin();

      const std::vector< std::pair<DetId, float> > &hf = it->hitsAndFractions();
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
	layer = HGCalDetId(hf[j].first).layer();

    const GlobalPoint position = recHitTools.getPosition(hf[j].first);
    const unsigned int wafer = recHitTools.getWafer(hf[j].first);
    const unsigned int cell  = recHitTools.getCell(hf[j].first);
    const double cellThickness = recHitTools.getSiThickness(hf[j].first);
    const bool isHalfCell = recHitTools.isHalfCell(hf[j].first);
    const double eta = recHitTools.getEta(position, vz);
    const double phi = recHitTools.getPhi(position);
    const double pt = recHitTools.getPt(position, hit->energy(), vz);

	if(hit->energy() > hitmap[hf[rhSeed].first]->energy()) rhSeed = j;

	if(layer < 29){
	  nhit++;
	}
    arhc->push_back(ARecHit(layer, wafer, cell,
                position.x(), position.y(), position.z(),
                eta,phi,pt,
                hit->energy(), hit->time(), cellThickness,
                isHalfCell, flags, cluster_index));

      }
      double pt = it->energy() / cosh(it->eta());
      acdc->push_back(ACluster2d(it->x(),it->y(),it->z(),
				 it->eta(),it->phi(),pt,it->energy(),
				 layer, ncoreHit,hf.size(),i, rhSeed));
      cluster_index++;
    }

    double pt = multiClusters[i].total_uncalibrated_energy() / cosh(multiClusters[i].simple_eta(vz));
    amcc->push_back(AMultiCluster(multiClusters[i].simple_eta(vz),
				  multiClusters[i].simple_phi(),
                  pt,
				  multiClusters[i].simple_z(vz),
				  multiClusters[i].simple_slope_x(vz),
				  multiClusters[i].simple_slope_y(vz),
				  multiClusters[i].total_uncalibrated_energy(),
				  multiClusters[i].size(),
				  cl2dSeed));
  }

  // loop over simClusters
  for (std::vector<SimCluster>::const_iterator it_simClus = simClusters.begin(); it_simClus != simClusters.end(); ++it_simClus) {
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
        layers.push_back(HGCalDetId(it_haf->first).layer());
        wafers.push_back(recHitTools.getWafer(it_haf->first));
        cells.push_back(recHitTools.getCell(it_haf->first));
    }
    ascc->push_back(ASimCluster(it_simClus->pt(),
				  it_simClus->eta(),
				  it_simClus->phi(),
				  it_simClus->energy(),
				  it_simClus->simEnergy(),
                  it_simClus->numberOfSimHits(),
                  it_simClus->numberOfRecHits(),
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
  for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles.begin(); it_caloPart != caloParticles.end(); ++it_caloPart) {
    ++ncalopart;
    const std::vector<std::pair<uint32_t,float> > hits_and_fractions = it_caloPart->hits_and_fractions();
    std::vector<uint32_t> hits;
    std::vector<float> fractions;
    for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
        hits.push_back(it_haf->first);
        fractions.push_back(it_haf->second);
    }
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
                  it_caloPart->numberOfSimHits(),
                  it_caloPart->numberOfRecHits(),
				  hits,
				  fractions,
                  simClusterIndex));

  } // end loop over caloParticles


  event->set(iEvent.run(),iEvent.id().event(),npart,nhit,nclus,nmclus,
         nsimclus, npfclus, ncalopart,
	     vx,vy,vz);
  tree->Fill();
}

void
HGCalAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
HGCalAnalysis::endJob()
{
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

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalAnalysis);
