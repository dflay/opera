// compare trajectories from different files 

#include <cstdlib> 
#include <iostream>

#include "CSVManager.hh"
#include "JSONManager.hh"
#include "UtilDFGraph.hh"
#include "UtilDFFunc.hh"

#include "./include/operaParameters.hh"
#include "./src/operaUtilities.cc"

int CompareTraj(){

   bool useBeamE = false; // use alternate trajectory with specific beam energy tag 

   // read in parameters
   util_df::JSONManager *jpars = new util_df::JSONManager();
   jpars->ReadFile("./input/json/compare-traj.json");

   // file names
   std::vector<std::string> fileName;
   jpars->GetVectorFromKey_str("files",fileName);

   // kinematics 
   double sbsAngle = jpars->GetValueFromSubKey<double>("config","sbs-angle");
   double Ebeam    = jpars->GetValueFromSubKey<double>("config","beam-energy"); 
   double Ibeam    = jpars->GetValueFromSubKey<double>("config","beam-current"); 

   std::vector<double> USL,USR,DSL,DSR; 

   const int N = fileName.size();

   // vector of graphs 
   TGraphAsymmErrors **gBand = new TGraphAsymmErrors*[N]; 
 
   // for the track files
   util_df::CSVManager *data = new util_df::CSVManager("csv"); 

   char inpath[200]; 
   std::vector<double> z,min,max;
  
   int color[5] = {kBlack,kRed+2,kBlue+2,kOrange,kViolet}; 

   TString label; 
   TLegend *L = new TLegend(0.6,0.6,0.8,0.8);

   int rc=0;
   char inpath_pars[200];
   opera::parameters_t sim_pars;

   double I_corr=0,cc=0;
   std::vector<double> Icm2,I; 
   std::vector<double> USLi,USRi,DSLi,DSRi; 
   std::vector<std::string> ll; 

   // int M=0;
   // std::vector<std::string> sv;
   // std::string suffix; 
   std::string resName;

   for(int i=0;i<N;i++){
      // read in data
      sprintf(inpath,"./output/%s_band.csv",fileName[i].c_str()); 
      std::cout << "Reading file: " << inpath << std::endl; 
      rc = data->ReadFile(inpath,true);
      data->GetColumn_byName<double>("Z(cm)"   ,z  );  
      data->GetColumn_byName<double>("minR(cm)",min);  
      data->GetColumn_byName<double>("maxR(cm)",max);
      if(rc!=0){
	 // file read failed
	 return 1; 
      } 
      // get simulation parameters
      // build the filename
      if(useBeamE){
	 resName = opera::getResNameFromTrackName_E(fileName[i],Ebeam,1);
         std::cout << Form("Beam energy = %.2lf GeV",Ebeam) << std::endl; 
      }else{
	 // resName = opera::getResNameFromTrackName(fileName[i]); 
	 resName = fileName[i] + ".res";  
      }
      sprintf(inpath_pars,"./data/fin/%s",resName.c_str());
      // sprintf(inpath_pars,"./data/fin/%s.res",fileName[i].c_str());
      std::cout << "Reading file: " << inpath_pars << "..." << std::endl;
      rc = opera::ReadResFile(inpath_pars,sim_pars);
      if(rc==0){
         opera::PrintParameters(sim_pars);
      }else{
         return 1;
      }
      // get currents 
      opera::CalculateCorrectorCurrents(sim_pars,ll,Icm2,I);
      // parse into vectors  
      USR.push_back(Icm2[0]); USL.push_back(Icm2[1]);
      DSR.push_back(Icm2[2]); DSL.push_back(Icm2[3]);        
      USRi.push_back(I[0]); USLi.push_back(I[1]);
      DSRi.push_back(I[2]); DSLi.push_back(I[3]);        
      // create band 
      gBand[i] = util_df::Graph::GetBand(z,min,max,color[i],1001,0.20);
      // add to legend
      // label = Form("USL = %.0lf, USR = %.0lf, DSL = %.0lf, DSR = %.0lf",USL[i],USR[i],DSL[i],DSR[i]);  // A/cm^2
      label = Form("USL = %.2lf A, USR = %.2lf A, DSL = %.2lf A, DSR = %.2lf A",USLi[i],USRi[i],DSLi[i],DSRi[i]); // A 
      L->AddEntry(gBand[i],label,"f");
      // clean up for next file 
      z.clear(); 
      min.clear(); 
      max.clear(); 
      ll.clear();
      I.clear();
      Icm2.clear();
      sim_pars.I.clear();
      sim_pars.iLabel.clear();
      sim_pars.drive.clear();
      sim_pars.dLabel.clear();
      data->ClearData();
   }

   double xMin = -1300; 
   double xMax = 1; 
   double yMin = 0; 
   double yMax = 1; 

   TString Title      = Form("I = %.0lf #muA, E = %.2lf GeV, #theta_{SBS} = %.1lf#circ",Ibeam,Ebeam,sbsAngle);
   TString xAxisTitle = Form("z [cm]");
   TString yAxisTitle = Form("r [cm]");

   TCanvas *c1 = new TCanvas("c1","Trajectory Ranges",1000,600);

   c1->cd();
   for(int i=0;i<N;i++){
      if(i==0){
	 gBand[i]->Draw("a e4");
         gBand[i]->SetTitle(Title); 
         gBand[i]->GetXaxis()->SetLimits(xMin,xMax); 
         gBand[i]->GetYaxis()->SetRangeUser(yMin,yMax);  
         gBand[i]->GetXaxis()->SetTitle(xAxisTitle); 
         gBand[i]->GetXaxis()->CenterTitle(); 
         gBand[i]->GetYaxis()->SetTitle(yAxisTitle); 
         gBand[i]->GetYaxis()->CenterTitle(); 
      }else{
	 gBand[i]->Draw("e4 same");
      } 
      c1->Update();
   }
   L->Draw("same"); 
   c1->Update();

   return 0;
}
