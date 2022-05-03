// Read a single field map and plot x, y, z components
// for different trajectories in (x,y).  Integrate over z for each (x,y) trajectory 
// and show total field integral over some (x,y) surface  
// User inputs:
// - (x,y) trajectory to plot as a function of z  
// - Range of z to integrate  

#include <cstdlib> 
#include <iostream>

#include "TGraph2D.h"
#include "TArc.h"

#include "CSVManager.hh"
#include "UtilDFAlgorithm.hh"
#include "UtilDFImport.hh"
#include "UtilDFMath.hh"
#include "UtilDFGraph.hh"

#include "./include/operaMagneticFieldPoint.hh"
#include "./include/operaParameters.hh"
#include "./src/operaUtilities.cc"
// #include "./src/operaIntegrator.cc"

int StoreToGlobalVectors(util_df::CSVManager *data,double x0,double y0,int units); 

double IntegrateField_xySlice(util_df::CSVManager *data,std::string yAxis,double x0,double y0,int units); 
double fieldFunc_x(const double z);
double fieldFunc_y(const double z);
double fieldFunc_z(const double z);

// set to true if the inital map is in (mm,T) and want (cm,G) 
bool conv_mm_to_cm = false;
bool conv_T_to_G   = false;

// for field integral 
std::vector<double> gL,gbx,gby,gbz; 
std::vector<double> zz,bbx,bby,bbz; 

int xyInt(){

   gStyle->SetPalette(kRainBow);
   // gStyle->SetPalette(kLightTemperature);
   // gStyle->SetPalette(kTemperatureMap);

   util_df::CSVManager *pars = new util_df::CSVManager("csv"); 
   pars->ReadFile("./input/xy-int.csv",true);

   std::vector<std::string> fileName; 
   pars->GetColumn_byName_str("name",fileName);

   // only first row has the file name 
   char inpath_pars[200];
   sprintf(inpath_pars,"./data/fin/%s.res",fileName[0].c_str());

   // now get bounds on x, y trajectories
   std::vector<int> np;
   pars->GetColumn_byIndex<int>(2,np); 

   int nx = np[1];   
   int ny = np[2];  
 
   opera::parameters_t data; 
   int rc = opera::ReadResFile(inpath_pars,data);
   if(rc==0){
      std::cout << "Opening file: " << inpath_pars << std::endl;
      opera::PrintParameters(data);  
   }else{
      return 1;
   } 

   // get min and max of x and y for trajectories 
   std::vector<double> min,max; 
   pars->GetColumn_byIndex<double>(3,min); 
   pars->GetColumn_byIndex<double>(4,max); 
    
   double xMin = min[1]; 
   double yMin = min[2]; 
 
   double xMax = max[1]; 
   double yMax = max[2]; 

   // the last row is actually (zmin,zmax)
   int N = min.size(); 
   double zMin = min[N-1];
   double zMax = max[N-1]; 

   // determine corrector current densities and currents 
   std::vector<std::string> ll; 
   std::vector<double> Icm2,I; 
   opera::CalculateCorrectorCurrents(data,ll,Icm2,I);
   // current densities [A/cm2] 
   double USR = Icm2[0];  
   double USL = Icm2[1];  
   double DSR = Icm2[2];  
   double DSL = Icm2[3];
   // currents [A]  
   double iUSR = Icm2[0];  
   double iUSL = Icm2[1];  
   double iDSR = Icm2[2];  
   double iDSL = Icm2[3];

   std::cout << "Corrector current densities: " << std::endl;
   std::cout << Form("USL = %.3lf A/cm^2, USR = %.3lf A/cm^2, DSL = %.3lf A/cm^2, DSR = %.3lf A/cm^2",
                     USL,USR,DSL,DSR) << std::endl;
   std::cout << "Corrector currents: " << std::endl;
   std::cout << Form("USL = %.2lf A, USR = %.2lf A, DSL = %.2lf A, DSR = %.2lf A",
                     iUSL,iUSR,iDSL,iDSR) << std::endl;

   // data file header 
   std::string header = "X,Y,Z,BX,BY,BZ"; 

   char inpath1[200],inpath2[200],inpath3[200]; 
   sprintf(inpath1,"./data/fin/%s.table",fileName[0].c_str()); 

   int NSkip = 8; // for "raw" opera headers  
   // int NSkip = 11; // for g4sbs headers
 
   // field map 
   util_df::CSVManager *f1 = new util_df::CSVManager("tsv"); 
   rc = f1->ReadFile(inpath1,false,NSkip);
   if(rc!=0) return 1;
   f1->SetHeader(header);

   int units = opera::kGauss;

   std::string xAxis      = "Z"; 
   std::string yAxisUnits = "Gauss";  

   // double yMin = -1E+4; 
   // double yMax =  2E+4;

   if(units==opera::kTesla){
      yAxisUnits = "T";
      // yMin /= 1E+4; 
      // yMax /= 1E+4; 
   } 

   double xStep = (xMax-xMin)/( (double)nx ); 
   double yStep = (yMax-yMin)/( (double)ny ); 
   
   int NL = nx*ny;
   TGraph2D *g2D = new TGraph2D(NL);

   // // artificially make the plot area larger than the data
   // const int NPTS = 50;
   // double xMIN = -2.0;
   // double xMAX =  2.0;
   // double yMIN = -2.0;
   // double yMAX =  2.0;
   // double xs = (xMAX-xMIN)/( (double)NPTS );
   // double ys = (yMAX-yMIN)/( (double)NPTS );
   // TGraph2D *g2D = new TGraph2D(NPTS);

   // int k=0;
   // double xx=0,yy=0;
   // for(int i=0;i<NPTS;i++){
   //    xx = xMIN + ( (double)i )*xs;
   //    for(int j=0;j<NPTS;j++){
   //       yy = yMIN + ( (double)j )*ys;
   //       g2D->SetPoint(k,xx,yy,0);
   //       k++;
   //    }
   // }

   int k=0; 
   double arg_x=0,arg_y=0,arg_z=0,arg_xy=0,xi=0,yi=0;
   for(int i=0;i<=nx;i++){
      xi = xMin + ( (double)i )*xStep; 
      for(int j=0;j<=ny;j++){
	 yi = yMin + ( (double)j )*yStep; 
	 StoreToGlobalVectors(f1,xi,yi,units);
	 arg_x = util_df::Math::SimpsonIntegral(&fieldFunc_x,zMin,zMax); 
	 arg_y = util_df::Math::SimpsonIntegral(&fieldFunc_y,zMin,zMax); 
	 arg_z = util_df::Math::SimpsonIntegral(&fieldFunc_z,zMin,zMax);
         arg_xy = arg_x + arg_y;  
	 std::cout << "--------------------------" << std::endl;
	 std::cout << Form("x = %.2lf, y = %.2lf",xi,yi) << std::endl;
         std::cout << Form("  int Bx*dz = %.3lf %s-cm"     ,arg_x ,yAxisUnits.c_str()) << std::endl; 
         std::cout << Form("  int By*dz = %.3lf %s-cm"     ,arg_y ,yAxisUnits.c_str()) << std::endl; 
         std::cout << Form("  int Bz*dz = %.3lf %s-cm"     ,arg_z ,yAxisUnits.c_str()) << std::endl; 
         std::cout << Form("  int (Bx+By)*dz = %.3lf %s-cm",arg_xy,yAxisUnits.c_str()) << std::endl;
	 // fill graph
	 g2D->SetPoint(k,xi,yi,arg_xy);
	 k++;
      }
   }
 
   TArc *circle = new TArc(0,0,0.25);
   circle->SetLineColor(kBlack);
   circle->SetLineWidth(2);
   circle->SetFillStyle(0);

   TCanvas *c1 = new TCanvas("c1","Opera Field Map",1000,800);  

   c1->cd();   
   g2D->Draw("colz");
   g2D->GetXaxis()->SetTitle("x [cm]");  
   g2D->GetXaxis()->CenterTitle();  
   g2D->GetYaxis()->SetTitle("y [cm]");  
   g2D->GetYaxis()->CenterTitle();  
   g2D->Draw("colz");
   // circle->Draw("same"); 
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
int StoreToGlobalVectors(util_df::CSVManager *data,double x0,double y0,int units){
   // clear existing data 
   zz.clear(); 
   bbx.clear(); 
   bby.clear(); 
   bbz.clear(); 

   // // store the data for B(x,y) = B(x0,y0)
   // // coordinates x0,y0 may not exist in intial file; use TGraph to do the interpolation 
   // TGraph *gbx = opera::GetTGraph(data,"BX",x0,y0,units);  
   // TGraph *gby = opera::GetTGraph(data,"BY",x0,y0,units);  
   // TGraph *gbz = opera::GetTGraph(data,"BZ",x0,y0,units); 

   // const int N = gbx->GetN();
   // double *z  = gbx->GetX();  
   // double *bx = gbx->GetY();  
   // double *by = gby->GetY();  
   // double *bz = gbz->GetY(); 

   // for(int i=0;i<N;i++){
   //    zz.push_back(z[i]); 
   //    bbx.push_back(bx[i]); 
   //    bby.push_back(by[i]); 
   //    bbz.push_back(bz[i]); 
   // }

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
