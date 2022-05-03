// optimize corrector currents based on OPERA simulations 
// to minimize the field integral in the xy plane 

#include <cstdlib> 

#include "CSVManager.hh"
#include "UtilDFGraph.hh"

#include "./include/operaPoint.hh"
#include "./include/operaParameters.hh"

#include "./src/operaUtilities.cc"
#include "./src/operaTrack.cc"

int Optimize(){

   gStyle->SetPalette(kRainBow);
   // gStyle->SetPalette(kLightTemperature);
   // gStyle->SetPalette(kTemperatureMap);

   util_df::CSVManager *data = new util_df::CSVManager(); 
   data->ReadFile("./input/opt-gmn_bl03_233.csv",true); 
   // data->Print();  

   util_df::CSVManager *data2 = new util_df::CSVManager();
   data2->ReadFile("./input/gmn-field-integrals.csv",true);
   // data2->Print();

   TGraph *gf = util_df::Graph::GetTGraph(data2,"SBS-angle(deg)","Bxdz");
   util_df::Graph::SetParameters(gf,20,kBlack); 

   std::vector<double> USL,USR,DSL,DSR,xInt,yInt; 
   data->GetColumn_byName<double>("USL" ,USL);   
   data->GetColumn_byName<double>("USR" ,USR);   
   data->GetColumn_byName<double>("DSL" ,DSL);   
   data->GetColumn_byName<double>("DSR" ,DSR);   
   data->GetColumn_byName<double>("xInt",xInt);  
   data->GetColumn_byName<double>("yInt",yInt);  

   const int N = USL.size();
   TGraph2D *g2D = new TGraph2D(N);

   double arg=0,arg_us=0,arg_ds=0,arg_xy=0;
   std::vector<double> US,DS,It,xyInt; 
   for(int i=0;i<N;i++){
      // arg = fabs(USL[i]) + fabs(USR[i]) + fabs(DSL[i]) + fabs(DSR[i]);
      arg = DSL[i]; // fabs(USL[i]) + fabs(USR[i]);
      arg_us = 1E+3*(USL[i] + USR[i]); 
      arg_ds = 1E+3*(DSL[i] + DSR[i]); 
      arg_xy = xInt[i] + yInt[i];
      xyInt.push_back(arg_xy);
      US.push_back(arg_us); 
      DS.push_back(arg_ds);
      It.push_back(arg);
      std::cout << Form("i = %d, US = %.3lf, DS = %.3lf, int = %.3E",i,arg_us,arg_ds,arg_xy) << std::endl;
      g2D->SetPoint(i,arg_us,arg_ds,arg_xy); 
   }

   TGraph *g = util_df::Graph::GetTGraph(It,xyInt);
   util_df::Graph::SetParameters(g,20,kBlack,1,1); 

   TCanvas *c1 = new TCanvas("c1","Test",1000,600); 

   c1->cd(); 
   g->Draw("ap");
   util_df::Graph::SetLabels(g,"Bx Integral vs DSL Corrector Current","I [A/cm^{2}]","#int (Bx+By) dz [G-cm]"); 
   g->Draw("ap");
   c1->Update();
 
   // TCanvas *c2 = new TCanvas("c2","Corrector Currents and Field Integral",1000,600); 

   // TString xAxisTitle = Form("Total US Corrector Current [A/cm^{2}]"); 
   // TString yAxisTitle = Form("Total DS Corrector Current [A/cm^{2}]"); 

   // c2->cd(); 
   // g2D->Draw("colz");
   // g2D->GetXaxis()->SetTitle(xAxisTitle); 
   // g2D->GetXaxis()->CenterTitle(); 
   // g2D->GetYaxis()->SetTitle(yAxisTitle); 
   // g2D->GetYaxis()->CenterTitle(); 
   // g2D->Draw("colz");
   // c2->Update(); 

   TCanvas *c3 = new TCanvas("c3","Bx integral for all settings",600,600); 

   c3->cd();
   gf->Draw("ap"); 
   util_df::Graph::SetLabels(gf,"Bx Integral for SBS Angle Settings (Correctors Off)","SBS Angle (deg)","#int Bx dz (G-cm)");
   c3->Update();  

   return 0;
}
