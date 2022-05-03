#include "../include/operaTrack.hh"
//______________________________________________________________________________
namespace opera {
   //______________________________________________________________________________
   Track::Track(int id){
      fID = id;
   }
   //______________________________________________________________________________
   Track::~Track(){
      ClearData();
   }
   //______________________________________________________________________________
   void Track::ClearData(){
      fx.clear();  fy.clear();  fz.clear();
      fvx.clear(); fvy.clear(); fvz.clear();
      fID = -1;
   }
   //______________________________________________________________________________
   int Track::GetNumPoints(){
      int NX  = fx.size();
      int NY  = fy.size();
      int NZ  = fz.size();
      int NVX = fvx.size();
      int NVY = fvy.size();
      int NVZ = fvz.size();
  
      // check coordinates 
      int NP = 0; 
      if( NX==NY   && NX==NZ  && NY==NZ ) NP = NX;  
      
      // check velocity  
      int NV = 0; 
      if( NVX==NVY   && NVX==NVZ  && NVY==NVZ ) NV = NVX; 

      char msg[200]; 

      // check coordinates against velocity 
      if(NP==NV){
	 // check successful, do nothing 
      }else{
	 std::cout << "[Track::GetNumPoints]: ERROR! Number of points for (x,y,z) and (vx,vy,vz) don't match! " << std::endl;
	 sprintf(msg,"                       NX = %d, NY = %d, NZ = %d, NVX = %d, NVY = %d, NVZ = %d",NX,NY,NZ,NVX,NVY,NVZ); 
         std::cout << msg << std::endl;
      } 

      return NP; 
   }
   //______________________________________________________________________________
   void Track::PushBackCoordinateComp(char axis,double q){
      if(axis=='x') fx.push_back(q); 
      if(axis=='y') fy.push_back(q); 
      if(axis=='z') fz.push_back(q);
   }
   //______________________________________________________________________________
   void Track::SetCoordinateComp(char axis,int i,double q){
      int NX = fx.size();
      int NY = fy.size();
      int NZ = fz.size();
      if(axis=='x'){
	 if(i<NX) fx[i] = q;
      } 
      if(axis=='y'){
	 if(i<NY) fy[i] = q;
      } 
      if(axis=='z'){
	 if(i<NZ) fz[i] = q;
      }
   }
   //______________________________________________________________________________
   void Track::PushBackCoordinate(double x,double y,double z){
      fx.push_back(x); fy.push_back(y); fz.push_back(z); 
   }
   //______________________________________________________________________________
   void Track::SetCoordinate(int i,double x,double y,double z){
      int NX = fx.size();
      int NY = fy.size();
      int NZ = fz.size();
      if(i<NX) fx[i] = x;
      if(i<NY) fy[i] = y;
      if(i<NZ) fz[i] = z;
   }
   //______________________________________________________________________________
   void Track::PushBackVelocityComp(char axis,double q){
      if(axis=='x') fvx.push_back(q); 
      if(axis=='y') fvy.push_back(q); 
      if(axis=='z') fvz.push_back(q); 
   }
   //______________________________________________________________________________
   void Track::SetVelocityComp(char axis,int i,double q){
      int NX = fvx.size();
      int NY = fvy.size();
      int NZ = fvz.size();
      if(axis=='x'){
	 if(i<NX) fvx[i] = q;
      } 
      if(axis=='y'){
	 if(i<NY) fvy[i] = q;
      } 
      if(axis=='z'){
	 if(i<NZ) fvz[i] = q;
      } 
   }
   //______________________________________________________________________________
   void Track::PushBackVelocity(double x,double y,double z){
      fvx.push_back(x); fvy.push_back(y); fvz.push_back(z); 
   }
   //______________________________________________________________________________
   void Track::SetVelocity(int i,double x,double y,double z){
      int NX = fvx.size();
      int NY = fvy.size();
      int NZ = fvz.size();
      if(i<NX) fvx[i] = x;
      if(i<NY) fvy[i] = y;
      if(i<NZ) fvz[i] = z;
   }
} // ::opera 
