// Read multiple field maps and plot them 
// - Mainly to compare different corrector current settings  

#include <cstdlib> 
#include <iostream>

#include "TGraph2D.h"

#include "JSONManager.hh"
#include "CSVManager.hh"
#include "UtilDFAlgorithm.hh"
#include "UtilDFImport.hh"
#include "UtilDFMath.hh"
#include "UtilDFGraph.hh"

#include "./include/operaMagneticFieldPoint.hh"
#include "./include/operaParameters.hh"
#include "./src/operaUtilities.cc"

int StoreToGlobalVectors(util_df::CSVManager *data,double x0,double y0); 

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
// std::vector<opera::magneticFieldPt_t> opera;  

int CompareCC(){

   // subdirectory prefix in the data directory
   char prefix[200];
   sprintf(prefix,"fin");   // "final" study files => correct corrector coil geometry, final currents/drives

   // read in parameters
   util_df::JSONManager *jpars = new util_df::JSONManager(); 
   jpars->ReadFile("./input/json/cc-compare.json");

   // integration trajectory 
   double x0 = jpars->GetValueFromSubKey<double>("integral","x0"); 
   double y0 = jpars->GetValueFromSubKey<double>("integral","y0"); 

   // integration bounds on z
   double min = jpars->GetValueFromSubKey<double>("integral","min");
   double max = jpars->GetValueFromSubKey<double>("integral","max");
   double LENGTH = max - min;
 
   // plot range 
   double xMin = min; 
   double xMax = max; 

   // file names
   std::vector<std::string> fileName; 
   jpars->GetVectorFromKey_str("files",fileName);

   // config file 
   std::string confName = jpars->GetValueFromKey_str("config");
   std::cout << "Reading configuration: " << confName << std::endl;
   util_df::JSONManager *cf = new util_df::JSONManager();
   cf->ReadFile("./input/json/gen-conf.json");

   // SBS angle 
   double sbsAngle = cf->GetValueFromSubKey<double>(confName,"sbs-angle");

   // currents 
   std::vector<double> USL,USR,DSL,DSR; 

   std::string header = "X,Y,Z,BX,BY,BZ"; 

   int NSkip = 8; // for "raw" opera headers  
   // int NSkip = 11; // for g4sbs headers

   int rc=0;
   util_df::CSVManager *ff = new util_df::CSVManager("tsv"); 

   char inpath[200],inpath_pars[200]; 

   opera::parameters_t sim_pars; 

   const int NF = fileName.size(); 
   double bInt_xs[NF],bInt_ys[NF],bInt_zs[NF]; // simpson integrals 

   TGraph **gBx = new TGraph*[NF]; 
   TGraph **gBy = new TGraph*[NF];
   TGraph **gBz = new TGraph*[NF];

   TMultiGraph *mgx  = new TMultiGraph();
   TMultiGraph *mgy  = new TMultiGraph(); 
   TMultiGraph *mgz  = new TMultiGraph();

   double I_corr=0;

   TLegend *LX = new TLegend(0.6,0.6,0.8,0.8);
   TLegend *LY = new TLegend(0.6,0.6,0.8,0.8);
   TLegend *LZ = new TLegend(0.6,0.6,0.8,0.8);

   int color[10] = {kBlack ,kRed  ,kBlue  ,kGreen+1,kOrange+1,
                    kViolet,kRed+2,kBlue+2,kGreen+3,kOrange+3}; 

   int c=0;
   int lineWidth=2;
   int NC=0,ND=0;
   int units = opera::kGauss; 
   std::string xAxis      = "Z"; 
   std::string yAxisUnits = "Gauss";  
 
   TString labelName[NF]; 

   std::vector<std::string> label,ll; 
   std::vector<double> Icm2,I;
   std::vector<double> USLi,USRi,DSLi,DSRi; 
 
   double cc=0;

   bool isDebug = true;

   // read all files 
   for(int i=0;i<NF;i++){
      // get data table (field map) 
      sprintf(inpath,"./data/%s/%s.table",prefix,fileName[i].c_str()); 
      std::cout << "Reading file: " << inpath << "..." << std::endl;
      rc = ff->ReadFile(inpath,false,NSkip);
      if(rc!=0) return 1;
      ff->SetHeader(header);
      // get simulation parameters
      sprintf(inpath_pars,"./data/%s/%s.res",prefix,fileName[i].c_str()); 
      std::cout << "Reading file: " << inpath_pars << "..." << std::endl;
      rc = opera::ReadResFile(inpath_pars,sim_pars,isDebug);
      if(rc==0){
	 opera::PrintParameters(sim_pars);
      }else{
	 return 1;
      }
      // determine corrector current densities and currents 
      opera::CalculateCorrectorCurrents(sim_pars,ll,Icm2,I);  
      // current densities [A/cm2]  
      USR.push_back(Icm2[0]); USL.push_back(Icm2[1]); 
      DSR.push_back(Icm2[2]); DSL.push_back(Icm2[3]);
      // currents [A]  
      USRi.push_back(I[0]); USLi.push_back(I[1]); 
      DSRi.push_back(I[2]); DSLi.push_back(I[3]); 
      // set label 
      // labelName[i] = Form("US-L = %.0lf, US-R = %.0lf, DS-L = %.0lf, DS-R = %.0lf",USL[i],USR[i],DSL[i],DSR[i]);
      labelName[i] = Form("US-L = %.2lf A, US-R = %.2lf A, DS-L = %.2lf A, DS-R = %.2lf A",USLi[i],USRi[i],DSLi[i],DSRi[i]);
      // do integrals
      StoreToGlobalVectors(ff,x0,y0);
      bInt_xs[i] = util_df::Math::SimpsonIntegral(&fieldFunc_x,min,max); 
      bInt_ys[i] = util_df::Math::SimpsonIntegral(&fieldFunc_y,min,max); 
      bInt_zs[i] = util_df::Math::SimpsonIntegral(&fieldFunc_z,min,max);
      // get graphs  
      gBx[i] = opera::GetTGraph(ff,"BX",x0,y0,units); 
      gBy[i] = opera::GetTGraph(ff,"BY",x0,y0,units); 
      gBz[i] = opera::GetTGraph(ff,"BZ",x0,y0,units); 
      // set colors, markers
      c = i+1; 
      if(i<10) c = color[i];  
      util_df::Graph::SetParameters(gBx[i],20,c,0.5,lineWidth);
      util_df::Graph::SetParameters(gBy[i],20,c,0.5,lineWidth);
      util_df::Graph::SetParameters(gBz[i],20,c,0.5,lineWidth);
      // add to multigraphs 
      mgx->Add(gBx[i],"l");  
      mgy->Add(gBy[i],"l");  
      mgz->Add(gBz[i],"l");  
      // add labels 
      LX->AddEntry(gBx[i],labelName[i],"l"); 
      LY->AddEntry(gBy[i],labelName[i],"l"); 
      LZ->AddEntry(gBz[i],labelName[i],"l");
      // set up for next file 
      ff->ClearData();
      sim_pars.I.clear(); 
      sim_pars.iLabel.clear(); 
      sim_pars.drive.clear(); 
      sim_pars.dLabel.clear(); 
      ll.clear();
      Icm2.clear();
      I.clear();
   }

   double yMin = -1E+4; 
   double yMax =  2E+4;

   if(units==opera::kTesla){
      yAxisUnits = "T";
      yMin /= 1E+4; 
      yMax /= 1E+4; 
   } 

   TString Title      = Form("SBS Field Map for GMn-13 (Q^{2} = 12 GeV^{2})");
   TString xAxisTitle = Form("%s [cm]",xAxis.c_str());
   TString yAxisTitle = Form("B Field Component [%s]",yAxisUnits.c_str());
   
   // double bInt_sum[N];
   // for(int i=0;i<N;i++) bInt_sum[i]  = bInt_x[i]  + bInt_y[i]  + bInt_z[i];  
   double bInt_sums[NF];
   for(int i=0;i<NF;i++) bInt_sums[i] = bInt_xs[i] + bInt_ys[i] + bInt_zs[i];  

   std::cout << Form("Field integrals at (x,y) = (%.2lf,%.2lf) cm",x0,y0) << std::endl;
   for(int i=0;i<NF;i++){
      std::cout << "----------------" << std::endl;
      std::cout << Form("Corrector currents [A/cm^2]: US-L = %.1lf, US-R = %.1lf, DS-L = %.1lf, DS-R = %.1lf",
	                USL[i],USR[i],DSL[i],DSR[i]) << std::endl;
      std::cout << Form("Corrector currents [A]: US-L = %.1lf, US-R = %.1lf, DS-L = %.1lf, DS-R = %.1lf",
	                USLi[i],USRi[i],DSLi[i],DSRi[i]) << std::endl;
      std::cout << Form("   Bx*dL: simpson = %.3lf %s-cm",bInt_xs[i],yAxisUnits.c_str())     << std::endl;
      std::cout << Form("   By*dL: simpson = %.3lf %s-cm",bInt_ys[i],yAxisUnits.c_str())     << std::endl;
      std::cout << Form("   Bz*dL: simpson = %.3lf %s-cm",bInt_zs[i],yAxisUnits.c_str())     << std::endl;
      std::cout << Form("   total: simpson = %.3lf %s-cm",bInt_sums[i],yAxisUnits.c_str()) << std::endl;
   }
   
   std::cout << "**** For cut and paste use: " << std::endl;
   for(int i=0;i<NF;i++){
      std::cout << Form("%.1lf,%.1lf,%.1lf,%.1lf",bInt_xs[i],bInt_ys[i],bInt_zs[i],bInt_sums[i]) << std::endl;
   }

   TString Title0       = Form("SBS Angle = %.1lf#circ (x = %.1lf cm, y = %.1lf cm)",sbsAngle,x0,y0);
   TString yAxisTitleB  = Form("Field Strength [%s]",yAxisUnits.c_str()); 
   TString yAxisTitleBX = Form("B_{x} [%s]",yAxisUnits.c_str()); 
   TString yAxisTitleBY = Form("B_{y} [%s]",yAxisUnits.c_str()); 
   TString yAxisTitleBZ = Form("B_{z} [%s]",yAxisUnits.c_str()); 

   TCanvas *c1 = new TCanvas("c1","Opera Field Map Bx and By (xy slices)",1000,800);  
   c1->Divide(1,2); 

   c1->cd(1);   
   mgx->Draw("a"); 
   util_df::Graph::SetLabels(mgx,Title0,xAxisTitle,yAxisTitleBX);
   util_df::Graph::SetLabelSizes(mgx,0.05,0.06);
   mgx->GetXaxis()->SetLimits(xMin,xMax); 
   mgx->Draw("a");
   LX->Draw("same"); 
   c1->Update(); 

   c1->cd(2);   
   mgy->Draw("a"); 
   util_df::Graph::SetLabels(mgy,"",xAxisTitle,yAxisTitleBY);
   util_df::Graph::SetLabelSizes(mgy,0.05,0.06);
   mgy->GetXaxis()->SetLimits(xMin,xMax); 
   mgy->Draw("a");
   c1->Update(); 

   TCanvas *c2 = new TCanvas("c2","Bz Component",1000,600); 

   c2->cd();
   mgz->Draw("a"); 
   util_df::Graph::SetLabels(mgz,"",xAxisTitle,yAxisTitleBZ);
   // util_df::Graph::SetLabelSizes(mgz,0.05,0.06);
   mgz->GetXaxis()->SetLimits(xMin,xMax); 
   mgz->Draw("a");
   LZ->Draw("same");
   c2->Update();

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
