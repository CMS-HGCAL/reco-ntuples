#ifndef arechit_h
#define arechit_h

#include "TObject.h"


#include <vector>

class ARecHit : public TObject
{
public:
  ARecHit() :  layer(0),x(0.),y(0.),z(0.),
	       energy(0),time(-1),thickness(0.),isHalf(false), cluster2d(-1)
  { 
  }
  ARecHit(int i_layer, float i_x, float i_y, float i_z,  
	  float i_energy, float i_time, float i_thickness,
	  bool i_isHalf, int i_cluster2d):  layer(i_layer),x(i_x),y(i_y),z(i_z),
					    energy(i_energy),time(i_time),
					    thickness(i_thickness),
					    isHalf(i_isHalf), 
					    cluster2d(i_cluster2d)
  {
  }


  virtual ~ARecHit() { }
  
  int   layer;
  float x,y,z;
  float eta,phi;
  float energy,time,thickness;
  bool  isHalf;
  int   cluster2d;

  ClassDef(ARecHit,1)
};

typedef std::vector<ARecHit> ARecHitCollection;

#endif
