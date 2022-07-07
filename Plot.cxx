// Read a single field map and plot x, y, z components 
// User inputs:
// - (x,y) trajectory to plot as a function of z  
// - Range of z to integrate  

#include <cstdlib> 
#include <iostream>

#include "TGraph2D.h"
#include "TArc.h"
#include "TStopwatch.h"

#include "CSVManager.hh"
#include "JSONManager.hh"
#include "UtilDFAlgorithm.hh"
#include "UtilDFImport.hh"
#include "UtilDFMath.hh"
#include "UtilDFGraph.hh"

#include "./include/operaMagneticFieldPoint.hh"
#include "./include/operaParameters.hh"
#include "./src/operaUtilities.cc"
// #include "./src/operaIntegrator.cc"

int StoreToGlobalVectors(util_df::CSVManager *data,double x0,double y0); 

double IntegrateField_xySlice(util_df::CSVManager *data,std::string yAxis,double x0,double y0,int units); 
double fieldFunc_x(const double z);
double fieldFunc_y(const double z);
double fieldFunc_z(const double z);
double circFunc(double *x,double *p); 

// set to true if the inital map is in (mm,T) and want (cm,G) 
bool conv_mm_to_cm = false;
bool conv_T_to_G   = false;

// for field integral 
std::vector<double> gL,gbx,gby,gbz; 
std::vector<double> zz,bbx,bby,bbz; 

int Plot(){

   // subdirectory prefix in the data directory
   char prefix[200]; 
   sprintf(prefix,"fin");   // "final" study files => correct corrector coil geometry, final currents/drives

   gStyle->SetPalette(kRainBow);
   // gStyle->SetPalette(kLightTemperature);
   // gStyle->SetPalette(kTemperatureMap);

   util_df::JSONManager *jpars = new util_df::JSONManager(); 
   jpars->ReadFile("./input/json/plot.json"); 

   std::string fileName = jpars->GetValueFromKey_str("file"); 

   // read parameters 
   char inpath_pars[200];
   sprintf(inpath_pars,"./data/%s/%s.res",prefix,fileName.c_str());
 
   opera::parameters_t data; 
   int rc = opera::ReadResFile(inpath_pars,data);
   if(rc==0){
      opera::PrintParameters(data);  
   }else{
      return 1;
   } 

   // the last row is actually (x,y,zmin,zmax)
   double x0   = jpars->GetValueFromSubKey<double>("integral","x0" );
   double y0   = jpars->GetValueFromSubKey<double>("integral","y0" ); 
   double zMin = jpars->GetValueFromSubKey<double>("integral","min");
   double zMax = jpars->GetValueFromSubKey<double>("integral","max"); 

   // get currents from res file 
   std::vector<double> Icm2,I; 
   std::vector<std::string> ll; 
   opera::CalculateCorrectorCurrents(data,ll,Icm2,I);  

   int NN = ll.size();
   std::vector<double> USL(NN),USR(NN),DSL(NN),DSR(NN); 
   std::vector<double> USRi,USLi,DSRi,DSLi; 
   for(int i=0;i<NN;i++){
      if(ll[i].compare("USR")==0){
	 USR[0]  = Icm2[i]; 
	 USRi.push_back(I[i]);
      } 
      if(ll[i].compare("USL")==0){
	 USL[0] = Icm2[i];
	 USLi.push_back(I[i]);
      } 
      if(ll[i].compare("DSR")==0){
	 DSR[0] = Icm2[i]; 
	 DSRi.push_back(I[i]);
      } 
      if(ll[i].compare("DSL")==0){
	 DSL[0] = Icm2[i];
	 DSLi.push_back(I[i]); 
      }  
   } 

   std::cout << "Corrector current densities: " << std::endl;
   std::cout << Form("USL = %.3lf A/cm^2, USR = %.3lf A/cm^2, DSL = %.3lf A/cm^2, DSR = %.3lf A/cm^2",
                     USL[0],USR[0],DSL[0],DSR[0]) << std::endl;
   std::cout << "Corrector currents: " << std::endl;
   std::cout << Form("USL = %.2lf A, USR = %.2lf A, DSL = %.2lf A, DSR = %.2lf A",
                     USLi[0],USRi[0],DSLi[0],DSRi[0]) << std::endl;

   // data file header 
   std::string header = "X,Y,Z,BX,BY,BZ"; 

   char table_path[200],inpath2[200],inpath3[200]; 
   sprintf(table_path,"./data/%s/%s.table",prefix,fileName.c_str()); 

   int NSkip = 8; // for "raw" opera headers  
   // int NSkip = 11; // for g4sbs headers
 
   // field map 
   util_df::CSVManager *f1 = new util_df::CSVManager("tsv"); 
   rc = f1->ReadFile(table_path,false,NSkip);
   if(rc!=0) return 1;
   f1->SetHeader(header);

   int units = opera::kGauss;

   std::string xAxis      = "Z"; 
   std::string yAxisUnits = "Gauss";  

   double yMin = -1E+4; 
   double yMax =  2E+4;

   if(units==opera::kTesla){
      yAxisUnits = "T";
      yMin /= 1E+4; 
      yMax /= 1E+4; 
   } 

   double LENGTH = zMax-zMin;

   double bInt_simp[3];    
   StoreToGlobalVectors(f1,x0,y0);
   bInt_simp[0] = util_df::Math::SimpsonIntegral(&fieldFunc_x,zMin,zMax); 
   bInt_simp[1] = util_df::Math::SimpsonIntegral(&fieldFunc_y,zMin,zMax); 
   bInt_simp[2] = util_df::Math::SimpsonIntegral(&fieldFunc_z,zMin,zMax);

   // opera::Integrator *myInt = new opera::Integrator(); 
   // myInt->SetData(f1,x0,y0); 
   // bInt_simp[0] = myInt->Integrate("BX",zMin,zMax); 
   // bInt_simp[1] = myInt->Integrate("BY",zMin,zMax); 
   // bInt_simp[2] = myInt->Integrate("BZ",zMin,zMax); 

   char axis[3] = {'x','y','z'};

   std::cout << Form("Integration for (x,y) = (%.2lf,%.2lf) cm, z = %.2lf to %.2lf",x0,y0,zMin,zMax) << std::endl;
   for(int i=0;i<3;i++){
      std::cout << Form("B%c*dz = %.3lf %s-cm",axis[i],bInt_simp[i],yAxisUnits.c_str()) << std::endl;
   }

   TStopwatch *watch = new TStopwatch(); 

   // get plots
   watch->Start(); 
   TGraph *gBx = opera::GetTGraph(f1,"BX",x0,y0,units);
   double t_x = watch->RealTime(); 

   watch->Start();
   TGraph *gBy = opera::GetTGraph(f1,"BY",x0,y0,units);
   double t_y = watch->RealTime();
 
   watch->Start();
   TGraph *gBz = opera::GetTGraph(f1,"BZ",x0,y0,units);
   double t_z = watch->RealTime();

   std::cout << Form("Time: %.1lf s, %.1lf s, %.1lf s",t_x,t_y,t_z) << std::endl;

   util_df::Graph::SetParameters(gBx,20,kBlack  ,0.5,2); 
   util_df::Graph::SetParameters(gBy,20,kGreen+2,0.5,2); 
   util_df::Graph::SetParameters(gBz,20,kBlue   ,0.5,2); 

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gBx,"l"); 
   mg->Add(gBy,"l"); 
   // mg->Add(gBz,"l"); 

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(gBx,"B_{x}","l");  
   L->AddEntry(gBy,"B_{y}","l");  
   // L->AddEntry(gBz,"B_{z}","l");  

   TString Title      = Form("(x,y) = (%.2lf,%.2lf) cm, #int_{%.1lf}^{%.1lf} Bx dz = %.1lf %s-cm, #int_{%.1lf}^{%.1lf} By dz = %.1lf %s-cm",
                             x0,y0,zMin,zMax,bInt_simp[0],yAxisUnits.c_str(),zMin,zMax,bInt_simp[1],yAxisUnits.c_str());
   TString xAxisTitle = Form("%s [cm]",xAxis.c_str());
   TString yAxisTitle = Form("B_{q} [%s]",yAxisUnits.c_str());

   TCanvas *c1 = new TCanvas("c1","Opera Field Map",1000,800);  

   c1->cd();   
   mg->Draw("a"); 
   util_df::Graph::SetLabels(mg,Title,xAxisTitle,yAxisTitle);
   mg->Draw("a");
   L->Draw("same"); 
   c1->Update(); 

   return 0;
}
//______________________________________________________________________________
double IntegrateField_xySlice(util_df::CSVManager *data,std::string yAxis,double x0,double y0,int units){
   // compute the field integral of a field component yAxis along the z axis
   std::vector<double> x,y,z; 
   data->GetColumn_byName<double>("X",x);   
   data->GetColumn_byName<double>("Y",y);   
   data->GetColumn_byName<double>("Z",z);  

   std::vector<double> b;
   if(yAxis.compare("bmod")==0){
      opera::GetBMod(data,b);
   }else{
      data->GetColumn_byName<double>(yAxis,b); 
   }

   // B units; input data are in Gauss  
   std::string unitName = "Gauss";
   const int N = x.size();
   if(units==opera::kTesla){
      for(int i=0;i<N;i++) b[i] /= 1E+4; 
      unitName = "Tesla";  
   }

   if(conv_T_to_G){
      for(int i=0;i<N;i++) b[i] *= 1E+4;  
   }

   double sf=1;
   if(conv_mm_to_cm) sf = 1E-1;  

   // choose values for a given (x0,y0)
   double step=0;
   double z_prev=z[0];  
   std::vector<double> L,B;  
   for(int i=0;i<N;i++){
      if( abs(x[i]-x0/sf)<1E-3 && abs(y[i]-y0/sf)<1E-3 ){
	 L.push_back(z[i]); 
	 B.push_back(b[i]); 
	 step = z[i]-z_prev;
      }
      z_prev = z[i]; 
   }

   // compute integral; units are [cm*B-field], with B-field in Gauss or Tesla 
   double sum=0,stepSum=0,bSum=0;
   const int NN = L.size();
   for(int i=0;i<NN;i++){
      sum     += step*B[i];
      stepSum += step; 
      bSum    += B[i]; 
   }
   sum += 0; // get rid of -0 result   
   // std::cout << Form("[IntegrateField_xySlice]: Step size = %.3lf cm, step sum = %.3lf cm, B sum = %.3lf %s",
   //                   step,stepSum,bSum,unitName.c_str()) << std::endl;
   return sum; 
}
//______________________________________________________________________________
int StoreToGlobalVectors(util_df::CSVManager *data,double x0,double y0){

   // clear existing data 
   zz.clear(); 
   bbx.clear(); 
   bby.clear(); 
   bbz.clear();
  
   // double sf=1;
   // if(conv_mm_to_cm) sf = 1E-1;  

   // // consolidate into a magneticFieldPt struct  
   // std::vector<opera::magneticFieldPt_t> bf; 
   // opera::ConsolidateData(data,bf); 
   // // sort the opera data
   // std::sort(bf.begin(),bf.end(),opera::compareByZ); 

   // int NN = bf.size();
   // for(int i=0;i<NN;i++){
   //    // cut on data for x = x0, y = y0
   //    if( abs(bf[i].x-x0/sf)<1E-3 && abs(bf[i].y-y0/sf)<1E-3 ){
   //       zz.push_back(bf[i].z); 
   //       bbx.push_back(bf[i].bx); 
   //       bby.push_back(bf[i].by); 
   //       bbz.push_back(bf[i].bz); 
   //    }
   // }

   int units = opera::kGauss;

   std::vector<double> z,bx,by,bz;
   opera::GetVectors_Interp(data,"BX",x0,y0,units,z,bx);
   z.clear();
   opera::GetVectors_Interp(data,"BY",x0,y0,units,z,by);
   z.clear();
   opera::GetVectors_Interp(data,"BZ",x0,y0,units,z,bz);

   opera::magneticFieldPt_t apt;
   std::vector<opera::magneticFieldPt_t> bf; 
   const int N = z.size();
   for(int i=0;i<N;i++){
      apt.z  = z[i]; 
      apt.bx = bx[i]; apt.by = by[i]; apt.bz = bz[i];
      bf.push_back(apt);
   }

   // sort by z! (not really needed) 
   std::sort(bf.begin(),bf.end(),opera::compareByZ);  

   for(int i=0;i<N;i++){
      zz.push_back(bf[i].z);
      bbx.push_back(bf[i].bx);
      bby.push_back(bf[i].by);
      bbz.push_back(bf[i].bz);
   }

   return 0;
}
//______________________________________________________________________________
double fieldFunc_x(const double z){  
   int lo=0,hi=0;
   util_df::Algorithm::BinarySearch<double>(zz,z,lo,hi); 
   double arg = 0.5*(bbx[lo] + bbx[hi]);
   return arg; 
}
//______________________________________________________________________________
double fieldFunc_y(const double z){
   int lo=0,hi=0;
   util_df::Algorithm::BinarySearch<double>(zz,z,lo,hi); 
   double arg = 0.5*(bby[lo] + bby[hi]);
   return arg; 
}
//______________________________________________________________________________
double fieldFunc_z(const double z){
   int lo=0,hi=0;
   util_df::Algorithm::BinarySearch<double>(zz,z,lo,hi); 
   double arg = 0.5*(bbz[lo] + bbz[hi]);
   return arg; 
}
