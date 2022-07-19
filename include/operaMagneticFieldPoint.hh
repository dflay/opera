#ifndef OPERA_MAGNETIC_FIELD_POINT_HH
#define OPERA_MAGNETIC_FIELD_POINT_HH

// convenient struct for handling Opera data 

namespace opera {
   enum bUnits { kTesla, kGauss };

   typedef struct magneticFieldPt {
      double x;    // x coordinate
      double y;    // y coordinate 
      double z;    // z coordinate
      double bx;   // bx value  
      double by;   // by value 
      double bz;   // bz value
      double bmod; // |B| = sqrt( bx*bx + by*by + bz*bz )  

     // constructor 
     magneticFieldPt(): 
	x(0),y(0),z(0),bx(0),by(0),bz(0),bmod(0)
     {}

   } magneticFieldPt_t;

   // sorting functions 
   bool compareByX(const magneticFieldPt &a,const magneticFieldPt &b){
      return (a.x<b.x);
   }
   bool compareByY(const magneticFieldPt &a,const magneticFieldPt &b){
      return (a.y<b.y);
   }
   bool compareByZ(const magneticFieldPt &a,const magneticFieldPt &b){
      return (a.z<b.z);
   }

}

#endif 
