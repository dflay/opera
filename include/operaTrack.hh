#ifndef OPERA_TRACK_HH
#define OPERA_TRACK_HH

// a track class for Opera data 

#include <cstdlib>
#include <vector>

// #include "operaPoint.hh"

// namespace opera {
//    typedef struct track { 
//       std::vector<point_t> p;  // x,y,z,vx,vy,vz
//       int ID;                  // track ID 
// 
//       // constructor 
//       track():
// 	 ID(-1)
//       {}
// 
//    } track_t; 
// } // ::opera

namespace opera { 

   class Track { 

      private: 
	 std::vector<double> fx,fy,fz;     //  x,  y,  z coordinates 
	 std::vector<double> fvx,fvy,fvz;  // vx, vy, vz components 
	 int fID;                          // track ID 

      public: 
	 Track(int id=-1); 
	 ~Track();

         void ClearData();

	 void PushBackCoordinateComp(char axis,double q);  
	 void PushBackCoordinate(double x,double y,double z);  

	 void SetCoordinateComp(char axis,int i,double q); 
	 void SetCoordinate(int i,double x,double y,double z);  
 
	 void PushBackVelocityComp(char axis,double v);  
	 void PushBackVelocity(double vx,double vy,double vz);

	 void SetVelocityComp(char axis,int i,double v);  
	 void SetVelocity(int i,double vx,double vy,double vz);

         int GetNumPoints(); 

         double GetX(int i)  const { return fx[i];  }  
         double GetY(int i)  const { return fy[i];  }  
         double GetZ(int i)  const { return fz[i];  } 
 
         double GetVX(int i) const { return fvx[i]; }  
         double GetVY(int i) const { return fvy[i]; }  
         double GetVZ(int i) const { return fvz[i]; }  
   }; 

} // ::opera 

#endif 
