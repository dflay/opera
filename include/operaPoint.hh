#ifndef OPERA_POINT_HH
#define OPERA_POINT_HH

// a data struct to store opera data points (from track files) 

#include <cstdlib>  

namespace opera { 
   typedef struct point { 
      double x;      // x position  
      double y;      // y position
      double z;      // z position
      double vx;     // speed along x
      double vy;     // speed along y 
      double vz;     // speed along z
      int ID;        // track number 

      point():
	 x(0),y(0),z(0),vx(0),vy(0),vz(0),ID(-1) 
      { } 

   } point_t;
} // ::opera

#endif 
