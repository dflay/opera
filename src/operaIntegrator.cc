#include "../include/operaIntegrator.hh"
//______________________________________________________________________________
namespace opera { 
   //______________________________________________________________________________
   Integrator::Integrator(){

   }
   //______________________________________________________________________________
   Integrator::~Integrator(){
      ClearData();
   }
   //______________________________________________________________________________
   void Integrator::SetData(util_df::CSVManager *data,double x0,double y0){
      // set data vectors
      std::vector<magneticFieldPt_t> bf;
      opera::ConsolidateData(data,bf);

      // sort the opera data
      std::sort(bf.begin(),bf.end(),opera::compareByZ);
      int NN = bf.size();
      for(int i=0;i<NN;i++){
	 // cut on data for x = x0, y = y0
	 if( abs(bf[i].x-x0)<1E-3 && abs(bf[i].y-y0)<1E-3 ){
	    fz.push_back(bf[i].z);
	    fBx.push_back(bf[i].bx);
	    fBy.push_back(bf[i].by);
	    fBz.push_back(bf[i].bz);
	 }
      }
   }
   //______________________________________________________________________________
   void Integrator::ClearData(){
      fBx.clear();
      fBy.clear();
      fBz.clear();
      fz.clear();
   }
   //______________________________________________________________________________
   double Integrator::Integrate(std::string Bq,double zMin,double zMax){
      // integrate field component Bq (q = x,y,z) over the z axis
      // FIXME: The function call isn't working! 
      double ans=0;
      if(Bq.compare("BX")==0) ans = util_df::Math::SimpsonIntegral(void(*Integrator::fieldFunc_x)(const double),zMin,zMax); 
      if(Bq.compare("BY")==0) ans = util_df::Math::SimpsonIntegral(void(*Integrator::fieldFunc_y)(const double),zMin,zMax); 
      if(Bq.compare("BZ")==0) ans = util_df::Math::SimpsonIntegral(void(*Integrator::fieldFunc_z)(const double),zMin,zMax); 
      return ans;
   }
   //______________________________________________________________________________
   double Integrator::fieldFunc_x(const double z){
      // Bx as a function of z
      int lo=0,hi=0;
      util_df::Algorithm::BinarySearch<double>(fz,z,lo,hi);
      double arg = 0.5*(fBx[lo] + fBx[hi]);
      return arg;
   }
   //______________________________________________________________________________
   double Integrator::fieldFunc_y(const double z){
      // By as a function of z
      int lo=0,hi=0;
      util_df::Algorithm::BinarySearch<double>(fz,z,lo,hi);
      double arg = 0.5*(fBy[lo] + fBy[hi]);
      return arg;
   }
   //______________________________________________________________________________
   double Integrator::fieldFunc_z(const double z){
      // Bz as a function of z
      int lo=0,hi=0;
      util_df::Algorithm::BinarySearch<double>(fz,z,lo,hi);
      double arg = 0.5*(fBz[lo] + fBz[hi]);
      return arg;
   }

} // ::opera 
