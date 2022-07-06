// plot arbitrary BH curves

#include <cstdlib> 
#include <iostream>

#include "TStyle.h"
#include "TPad.h"

#include "UtilDFUnits.hh"
#include "UtilDFConstants.hh"
#include "UtilDFGraph.hh"
#include "CSVManager.hh"

#include "./src/NonLinearBHCurve.cc"

TGraph *GetTGraph(NonLinearBHCurve *bh,std::string xAxis,std::string yAxis); 

int BH(){

   std::string units = "CGS"; 
   std::string steel = "annealed"; // "cold-roll"; 

   NonLinearBHCurve *bh_1030_old = new NonLinearBHCurve("1030_old");
   bh_1030_old->LoadData("./data/bh/tenthirty.bh","tsv",true,0);

   NonLinearBHCurve *bh_1030_new = new NonLinearBHCurve("1030");
   bh_1030_new->LoadData("./data/bh/tenthirty.bh","tsv",true,0);

   NonLinearBHCurve *bh_1006_old = new NonLinearBHCurve("1006_old");
   bh_1006_old->LoadData("./data/bh/default.bh","tsv",true,0);

   NonLinearBHCurve *bh_1008_new = new NonLinearBHCurve("1008");
   bh_1008_new->LoadData("./data/bh/1008.bh","tsv",true,kConvertMKStoCGS);

   std::string xAxis = "H";
   std::string yAxis = "B"; 

   int lineWidth = 4; 
   double markerSize = 1.0; 

   TGraph *g1030_old = GetTGraph(bh_1030_old,xAxis,yAxis);
   util_df::Graph::SetParameters(g1030_old,21,kBlack,markerSize,lineWidth);

   TGraph *g1030_new = GetTGraph(bh_1030_new,xAxis,yAxis);
   util_df::Graph::SetParameters(g1030_new,20,kRed,markerSize,lineWidth);
   g1030_new->SetLineStyle(2);

   TGraph *g1006_old = GetTGraph(bh_1006_old,xAxis,yAxis);
   util_df::Graph::SetParameters(g1006_old,21,kBlue,markerSize,lineWidth);

   TGraph *g1008_new = GetTGraph(bh_1008_new,xAxis,yAxis);
   util_df::Graph::SetParameters(g1008_new,20,kCyan+2,markerSize,lineWidth);
   g1008_new->SetLineStyle(2);

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(g1030_old,"c");  
   mg->Add(g1030_new,"c");  
   mg->Add(g1006_old,"c");  
   mg->Add(g1008_new,"c");  

   TLegend *LL = new TLegend(0.6,0.6,0.8,0.8);
   LL->AddEntry(g1030_old,"Beamline Shielding: 1030 (Old build)","l"); 
   LL->AddEntry(g1030_new,"Beamline Shielding: 1030 (New build)","l"); 
   LL->AddEntry(g1006_old,"SBS and Corrector Yokes: Default (Old build)"       ,"l"); 
   LL->AddEntry(g1008_new,"SBS and Corrector Yokes: 1008 (New build)"          ,"l"); 

   TString Title      = Form("B vs H for Various Materials");
   TString xAxisTitle = Form("H (A/m)");
   TString yAxisTitle = Form("B (T)");

   if(units.compare("CGS")==0){
      xAxisTitle = Form("H (Oerstead)");
      yAxisTitle = Form("B (Gauss)");
   }

   double xMin = 1E-2; 
   double xMax = 5E+3; 

   TCanvas *c1 = new TCanvas("c1","B vs H",1000,600);
 
   c1->cd(); 
   gPad->SetLogx();  
   mg->Draw("a"); 
   util_df::Graph::SetLabels(mg,Title,xAxisTitle,yAxisTitle); 
   mg->GetXaxis()->SetLimits(xMin,xMax);  
   mg->Draw("a");
   LL->Draw("same");
   c1->Update(); 

   return 0;
}
//______________________________________________________________________________
TGraph *GetTGraph(NonLinearBHCurve *bh,std::string xAxis,std::string yAxis){

   std::vector<double> B,H; 
   bh->GetData(B,H);  

   double boh=0;
   std::vector<double> x,y; 
   const int N = B.size();
   for(int i=0;i<N;i++){
      boh = 0;
      if(H[i]!=0) boh = B[i]/H[i];
      // x axis 
      if(xAxis.compare("B")==0)   x.push_back(B[i]); 
      if(xAxis.compare("H")==0)   x.push_back(H[i]); 
      if(xAxis.compare("BoH")==0) x.push_back(boh); 
      // y axis 
      if(yAxis.compare("B")==0)   y.push_back(B[i]); 
      if(yAxis.compare("H")==0)   y.push_back(H[i]); 
      if(yAxis.compare("BoH")==0) y.push_back(boh); 
   }

   TGraph *g = util_df::Graph::GetTGraph(x,y);
   return g;
}
