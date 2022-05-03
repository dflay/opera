#ifndef OPERA_PARAMETERS_HH
#define OPERA_PARAMETERS_HH

// opera data parameters
// these are found in *.res files

#include <cstdlib>
#include <string>
#include <vector>  

namespace opera {
   typedef struct parameters {
      std::string symmetry;            // symmetry  
      int NActiveElem;                 // number of active elements         
      int NNodes;                      // number of nodes
      int NEqns;                       // number of equations 
      int NNonZeros;                   // number of non-zeros 
      std::vector<double> I;           // current densities  
      std::vector<double> drive;       // drive scale factors
      std::vector<std::string> iLabel; // labels for currents
      std::vector<std::string> dLabel; // labels for drive scale factors 

      // constructor 
      parameters():
	 symmetry("UNKNOWN"),NActiveElem(0),NNodes(0),NEqns(0),NNonZeros(0)
      {} 

   } parameters_t;
} // ::opera 

#endif 
