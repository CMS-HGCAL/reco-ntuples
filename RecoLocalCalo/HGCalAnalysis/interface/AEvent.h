#ifndef aevent_h
#define aevent_h

#include "TObject.h"

#include <vector>

class AGenPart;
class ARecHit;
class ACluster2d;
class AMultiCluster;
class AEvent : public TObject
{
public:

  AEvent(): run(0), evn(0), ngen(0), nhit(0), nclus2d(0),nclus3d(0),
	    vx(0.), vy(0.), vz(0.)
  {
  }
  void set(int i_run, int i_event, 
	   int i_ngen, int i_nhit, int i_nclus, int i_nmclus,
	   float i_vx,float i_vy,float i_vz)
  {
    run =   i_run;
    evn = i_event; 
    ngen = i_ngen;
    nhit = i_nhit;
    nclus2d = i_nclus;
    nclus3d = i_nmclus;
    vx = i_vx; 
    vy = i_vy; 
    vz = i_vz;
  
  }

  unsigned int run;
  unsigned int evn;
  int   ngen;
  int   nhit;
  int   nclus2d;
  int   nclus3d;
  float vx;
  float vy;
  float vz;


  ClassDef(AEvent,1)
};
#endif
