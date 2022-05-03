#ifndef OPERA_MAGNET_COIL_HH
#define OPERA_MAGNET_COIL_HH

// a data struct for an SBS corrector magnet coil in Opera  

namespace opera {
   typedef struct magnetCoil { 
      std::string name;         
      double current;            // [A]
      double currentDensity;     // [A/cm^2]  
      double xSize;              // size in x (horizontal) [cm]
      double ySize;              // size in y (vertical) [cm] 
      double driveScale;         // drive scale (from OPERA); this is applied when calculating current  
      int nTurns;                // number of turns
      int nCoils;                // number of coils 
      // constructor
      magnetCoil():
	 name("NONE"),current(0),currentDensity(0),xSize(0),ySize(0),driveScale(0),nTurns(0),nCoils(0)
      {} 
   } magnetCoil_t;
} // ::opera 

#endif 
