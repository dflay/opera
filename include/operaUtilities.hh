// namespace for Opera utility functions 

#include <cstdlib> 
#include <iostream>
#include <vector>

#include "TGraph.h"

#include "UtilDFFunc.hh"
#include "UtilDFAlgorithm.hh"
#include "CSVManager.hh"

#include "operaPoint.hh"
#include "operaMagneticFieldPoint.hh"
#include "operaParameters.hh"
#include "operaMagnetCoil.hh"

namespace opera {
   void PrintParameters(parameters_t data);
   void PrintCoil(magnetCoil_t mc); 
   void CalculateCurrents(parameters_t pars,std::vector<magnetCoil_t> &data,std::string fileName); 
   void CalculateCorrectorCurrents(parameters_t pars,std::vector<std::string> &label,
                                   std::vector<double> &Icm2,std::vector<double> &I); 
   void CalculateCorrectorCurrents(std::vector<std::string> label,std::vector<double> Icm2,std::vector<double> &I);

   void GetVectors(util_df::CSVManager *data,std::string yAxis,double x0,double y0,int units,
                   std::vector<double> &Z,std::vector<double> &B); 
   void GetVectors_Interp(util_df::CSVManager *data,std::string yAxis,double x0,double y0,int units,
                   std::vector<double> &Z,std::vector<double> &B); 

   int ReadResFile(const char *inpath,parameters_t &data,bool isDebug=false); 
   int ReadTrackFile(const char *inpath,std::vector<point_t> &track,bool isDebug=false); 
   int GetBMod(util_df::CSVManager *data,std::vector<double> &bmod); 
   int ConsolidateData(util_df::CSVManager *in,std::vector<magneticFieldPt_t> &data); 

   double GetSBSAngle(std::string fileName); 

   std::string getResNameFromTrackName(std::string fileName);
   std::string getResNameFromTrackName_E(std::string fileName,double &E,int eIndexOffset=0);

   TGraph *GetTGraph(util_df::CSVManager *data,std::string yAxis,double x0,double y0,int units); 
   TGraph *GetTGraph_xySlice(util_df::CSVManager *data,std::string yAxis,double x0,double y0,int units);
   TGraph *GetTGraph(std::vector<point_t> data,std::string xAxis,std::string yAxis); 
} // ::opera 
