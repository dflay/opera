#ifndef OPERA_MAGNETIC_FIELD_POINT_HH
#define OPERA_MAGNETIC_FIELD_POINT_HH

// convenient struct for handling Opera data 

namespace opera {
   enum bUnits { kTesla, kGauss };

   typedef struct magneticFieldPt {
      double x;
      double y;
      double z;
      double bx;
      double by;
      double bz;
      double bmod;

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
