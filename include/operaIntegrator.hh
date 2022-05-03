#ifndef OPERA_INTEGRATOR_HH
#define OPERA_INTEGRATOR_HH

// a class to simplify/manage integrating the magnetic field produced by Opera

#include "CSVManager.hh"
#include "UtilDFAlgorithm.hh"
#include "UtilDFMath.hh"

namespace opera {

   class Integrator { 

      private: 
	 std::vector<double> fz;
	 std::vector<double> fBx;
	 std::vector<double> fBy;
	 std::vector<double> fBz;

	 double fieldFunc_x(const double); 
	 double fieldFunc_y(const double); 
	 double fieldFunc_z(const double); 

      public: 
	 Integrator();
	 ~Integrator();

	 void SetData(util_df::CSVManager *data,double x0,double y0); 
	 void ClearData(); 

	 double Integrate(std::string Bq,double zMin,double zMax); 

   }; 

} // ::opera

#endif 
