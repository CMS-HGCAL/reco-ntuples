#ifndef aevent_h
#define aevent_h

#include "TObject.h"

#include <vector>

class AGenPart;
class ARecHit;
class ACluster2d;
class AMultiCluster;
class ASimCluster;
class ACaloParticle;
class AEvent : public TObject
{
public:

  AEvent(): run(0), evn(0), ngen(0), nhit(0), nhit_raw(0), nclus2d(0),nclus3d(0),
        nsimclus(0), npfclus(0), ncalopart(0),
	    vx(0.), vy(0.), vz(0.)
  {
  }
  void set(int i_run, int i_event,
	   int i_ngen, int i_nhit, int i_nhit_raw, int i_nclus, int i_nmclus,
       int i_nsimclus, int i_npfclus, int i_ncalopart,
	   float i_vx,float i_vy,float i_vz)
  {
    run =   i_run;
    evn = i_event;
    ngen = i_ngen;
    nhit = i_nhit;
    nhit_raw = i_nhit_raw;
    nclus2d = i_nclus;
    nclus3d = i_nmclus;
    nsimclus = i_nsimclus;
    npfclus = i_npfclus;
    ncalopart = i_ncalopart;
    vx = i_vx;
    vy = i_vy;
    vz = i_vz;

  }

  unsigned int run;
  unsigned int evn;
  int   ngen;
  int   nhit;
  int   nhit_raw;
  int   nclus2d;
  int   nclus3d;
  int   nsimclus;
  int   npfclus;
  int   ncalopart;
  float vx;
  float vy;
  float vz;


  ClassDef(AEvent,1)
};
#endif
