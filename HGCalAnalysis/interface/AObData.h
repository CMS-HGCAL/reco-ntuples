#ifndef aobdata_h
#define aobdata_h

#include "TObject.h"

class AGenPart : public TObject
{
public:

  AGenPart(): eta(-1000.),phi(-1000.),pt(-1000.),energy(-1000.),dvx(0.),dvy(0.),dvz(0.),pid(0)
  {
  }
  AGenPart(float i_eta, float i_phi, float i_pt, float i_energy,
	   float i_dvx, float i_dvy,float i_dvz, int i_pid) :
    eta(i_eta),phi(i_phi),pt(i_pt),energy(i_energy),dvx(i_dvx),dvy(i_dvy),dvz(i_dvz),pid(i_pid)
  {
  }

  float eta;
  float phi;
  float pt;
  float energy;
  float dvx;
  float dvy;
  float dvz;
  int pid;

  ClassDef(AGenPart,1)
};


class ARecHit : public TObject
{
public:
  ARecHit() :  layer(0),wafer(0),cell(0),detid(0),
                x(0.),y(0.),z(0.),
                eta(-1000.),phi(-1000.),pt(-1000.),
                energy(-1000.),time(-1),thickness(0.),
                isHalf(false), flags(0), cluster2d(-1)
  {
  }
 ARecHit(int i_layer, int i_wafer, int i_cell, unsigned int i_detid,
	 float i_x, float i_y, float i_z,
     float i_eta, float i_phi, float i_pt,
	 float i_energy, float i_time, float i_thickness,
	 bool i_isHalf, int i_flags, int i_cluster2d):
    layer(i_layer),wafer(i_wafer),cell(i_cell),detid(i_detid),
    x(i_x),y(i_y),z(i_z),
    eta(i_eta),phi(i_phi),pt(i_pt),
    energy(i_energy),time(i_time),thickness(i_thickness),
    isHalf(i_isHalf),
    flags(i_flags),
    cluster2d(i_cluster2d)
    {
  }


  virtual ~ARecHit() { }

  int   layer,wafer,cell;
  unsigned int detid;
  float x,y,z;
  float eta,phi,pt;
  float energy,time,thickness;
  bool  isHalf;
  int   flags;
  int   cluster2d;

  ClassDef(ARecHit,1)
};


class ACluster2d : public TObject
{
public:

  ACluster2d() :  x(0.),y(0.),z(0.),eta(-1000.),phi(-1000.),pt(-1000.),
    energy(-1000.),layer(0),nhitCore(0),nhitAll(0), multicluster(-1), rechitSeed(-1)
  {
  }
  ACluster2d(  float i_x,
	       float i_y,
	       float i_z,
	       float i_eta,
	       float i_phi,
           float i_pt,
	       float i_energy,
	       int i_layer,
	       int i_nhitCore,
	       int i_nhitAll,
	       int i_multicluster,
	       int i_rechitSeed) :
  x(i_x),y(i_y),z(i_z),eta(i_eta),phi(i_phi),pt(i_pt),
    energy(i_energy),layer(i_layer),
    nhitCore(i_nhitCore),nhitAll(i_nhitAll),multicluster(i_multicluster), rechitSeed(i_rechitSeed)
  {
  }


  float x;
  float y;
  float z;
  float eta;
  float phi;
  float pt;
  float energy;
  int layer;
  int nhitCore;
  int nhitAll;
  int multicluster;
  int rechitSeed;

  ClassDef(ACluster2d,1)
};


class AMultiCluster : public TObject
{
public:

 AMultiCluster() :  eta(-1000.),phi(-1000.),pt(-1000.), z(0),
    slopeX(0.), slopeY(0.),
    energy(0.), nclus(0), cl2dSeed(-1)
  {
  }
  AMultiCluster(  float i_eta,
		  float i_phi,
          float i_pt,
		  float i_z,
		  float i_slopeX,
		  float i_slopeY,
		  float i_energy,
		  int i_nclus,
		  int i_cl2dSeed) :
  eta(i_eta),phi(i_phi),pt(i_pt),z(i_z),
    slopeX(i_slopeX), slopeY(i_slopeY),
    energy(i_energy),nclus(i_nclus), cl2dSeed(i_cl2dSeed)
  {
  }


  float eta;
  float phi;
  float pt;
  float z;
  float slopeX;
  float slopeY;
  float energy;
  int nclus;
  int cl2dSeed;


  ClassDef(AMultiCluster,1)
};

class ASimCluster : public TObject
{
public:
  ASimCluster() :  pt(0), eta(0), phi(0), energy(0),
    simEnergy(0), hits(0), fractions(0), layers(0),
    wafers(0), cells(0)
  {
  }
 ASimCluster(float i_pt, float i_eta, float i_phi, float i_energy,
     float i_simEnergy,
     std::vector<uint32_t> i_hits, std::vector<float> i_fractions,
     std::vector<unsigned int> i_layers, std::vector<unsigned int> i_wafers, std::vector<unsigned int> i_cells):
    pt(i_pt), eta(i_eta), phi(i_phi), energy(i_energy),
    simEnergy(i_simEnergy)
    {
    hits = i_hits;
    fractions = i_fractions;
    layers = i_layers;
    wafers = i_wafers;
    cells = i_cells;
  }


  virtual ~ASimCluster() { }

  float pt, eta, phi, energy, simEnergy;
  std::vector<uint32_t> hits;
  std::vector<float> fractions;
  std::vector<unsigned int> layers;
  std::vector<unsigned int> wafers;
  std::vector<unsigned int> cells;


  ClassDef(ASimCluster,1)
};


class APFCluster : public TObject
{
public:
  APFCluster() :  pt(0), eta(0), phi(0), energy(0)
    // simEnergy(0), hits(0), fractions(0)
  {
  }
 APFCluster(float i_pt, float i_eta, float i_phi, float i_energy
    //  float i_simEnergy, int i_numberOfSimHits, int i_numberOfRecHits,
    //  std::vector<uint32_t> i_hits, std::vector<float> i_fractions
 ):
    pt(i_pt), eta(i_eta), phi(i_phi), energy(i_energy)
    // simEnergy(i_simEnergy), numberOfSimHits(i_numberOfSimHits), numberOfRecHits(i_numberOfRecHits)
    {
    // hits = i_hits;
    // fractions = i_fractions;
  }


  virtual ~APFCluster() { }

  float pt, eta, phi, energy; //, simEnergy;
  // int numberOfSimHits, numberOfRecHits;
  // std::vector<uint32_t> hits;
  // std::vector<float> fractions;


  ClassDef(APFCluster,1)
};


class ACaloParticle : public TObject
{
public:
  ACaloParticle() :  pt(0), eta(0), phi(0), energy(0),
    simEnergy(0), simClusterIndex(0)
  {
  }
 ACaloParticle(float i_pt, float i_eta, float i_phi, float i_energy,
     float i_simEnergy,
     std::vector<uint32_t> i_simClusterIndex):
    pt(i_pt), eta(i_eta), phi(i_phi), energy(i_energy),
    simEnergy(i_simEnergy)
    {
    simClusterIndex = i_simClusterIndex;
  }


  virtual ~ACaloParticle() { }

  float pt, eta, phi, energy, simEnergy;
  std::vector<uint32_t> simClusterIndex;


  ClassDef(ACaloParticle,1)
};


typedef std::vector<AGenPart> AGenPartCollection;
typedef std::vector<ARecHit> ARecHitCollection;
typedef std::vector<ACluster2d> ACluster2dCollection;
typedef std::vector<AMultiCluster> AMultiClusterCollection;
typedef std::vector<ASimCluster> ASimClusterCollection;
typedef std::vector<APFCluster> APFClusterCollection;
typedef std::vector<ACaloParticle> ACaloParticleCollection;

#endif
