#ifndef aobdata_h
#define aobdata_h

#include "TObject.h"

class AGenPart : public TObject
{
public: 

  AGenPart(): eta(-1000.),phi(-1000.),pt(-1000.),dvx(0.),dvy(0.),dvz(0.),pid(0)
  {
  }  
  AGenPart(float i_eta, float i_phi, float i_pt, 
	   float i_dvx, float i_dvy,float i_dvz, int i_pid) :
    eta(i_eta),phi(i_phi),pt(i_pt),dvx(i_dvx),dvy(i_dvy),dvz(i_dvz),pid(i_pid)
  {
  }

  float eta;
  float phi;
  float pt;
  float dvx;
  float dvy;
  float dvz;
  int pid;

  ClassDef(AGenPart,1)
};


class ARecHit : public TObject
{
public:
  ARecHit() :  layer(0),x(0.),y(0.),z(0.),
	       energy(0),time(-1),thickness(0.),isHalf(false), flags(0), cluster2d(-1)
  { 
  }
  ARecHit(int i_layer, float i_x, float i_y, float i_z,  
	  float i_energy, float i_time, float i_thickness,
	  bool i_isHalf, int i_flags, int i_cluster2d):  layer(i_layer),x(i_x),y(i_y),z(i_z),
							 energy(i_energy),time(i_time),
							 thickness(i_thickness),
							 isHalf(i_isHalf),
							 flags(i_flags),
							 cluster2d(i_cluster2d)
  {
  }


  virtual ~ARecHit() { }
  
  int   layer;
  float x,y,z;
  float eta,phi;
  float energy,time,thickness;
  bool  isHalf;
  int   flags;
  int   cluster2d;

  ClassDef(ARecHit,1)
};


class ACluster2d : public TObject
{
public:
  
  ACluster2d() :  x(0.),y(0.),z(0.),eta(-1000.),phi(-1000.),
		  energy(0.),layer(0),nhitCore(0),nhitAll(0), multicluster(-1)
  { 
  }
  ACluster2d(  float i_x,
	       float i_y,
	       float i_z,
	       float i_eta,
	       float i_phi,
	       float i_energy,
	       int i_layer,
	       int i_nhitCore,
	       int i_nhitAll,
	       int i_multicluster) :  
    x(i_x),y(i_y),z(i_z),eta(i_eta),phi(i_phi),
    energy(i_energy),layer(i_layer),
    nhitCore(i_nhitCore),nhitAll(i_nhitAll),multicluster(i_multicluster)
  {
  }


  float x;
  float y;
  float z;
  float eta;
  float phi;
  float energy;
  int layer;
  int nhitCore;
  int nhitAll;
  int multicluster;
  
  ClassDef(ACluster2d,1)
};


class AMultiCluster : public TObject
{
public:
  
  AMultiCluster() :  eta(-1000.),phi(-1000.),
		     energy(0.), nclus(0)
  { 
  }
  AMultiCluster(  float i_eta,
		  float i_phi,
		  float i_energy,
		  int i_nclus) :  
    eta(i_eta),phi(i_phi),
    energy(i_energy),nclus(i_nclus)
  {
  }


  float eta;
  float phi;
  float energy;
  int nclus;

  
  ClassDef(AMultiCluster,1)
};


typedef std::vector<AGenPart> AGenPartCollection;
typedef std::vector<ARecHit> ARecHitCollection;
typedef std::vector<ACluster2d> ACluster2dCollection;
typedef std::vector<AMultiCluster> AMultiClusterCollection;

#endif
