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

   int N = 45; 
   double bMin = 0; 
   double bMax = 2.5; // in Tesla 

   // annealed 1030 steel 
   double k1_an = 50; 
   double k2_an = 1.371; 
   double k3_an = 645.3; 

   // cold rolled 1030 steel 
   double k1_cr = 40; 
   double k2_cr = 1.416; 
   double k3_cr = 1212; 

   double k1=0,k2=0,k3=0;
   if(steel.compare("annealed")==0){
      k1 = k1_an; k2 = k2_an; k3 = k3_an;
   }else{
      k1 = k1_cr; k2 = k2_cr; k3 = k3_cr;
   }

   int unts = util_df::Units::kMKS;
   if(units.compare("CGS")==0) unts = util_df::Units::kCGS;

   NonLinearBHCurve *bh = new NonLinearBHCurve(steel.c_str(),k1,k2,k3,unts);
   bh->Calculate(N,bMin,bMax);
   bh->WriteFile("tenthirty.bh");  

   std::vector<double> B,H; 
   bh->GetData(B,H); 

   NonLinearBHCurve *bh_cr = new NonLinearBHCurve("1030_cold-rolled",k1_cr,k2_cr,k3_cr,unts);
   bh_cr->Calculate(N,bMin,bMax); 
   bh_cr->WriteFile("tenthirty_cr.bh"); 
   
   std::vector<double> B_cr,H_cr,BoH_cr; 
   bh_cr->GetData(B_cr,H_cr);  

   NonLinearBHCurve *bh_an = new NonLinearBHCurve("1030_annealed",k1_an,k2_an,k3_an,unts);
   bh_an->Calculate(N,bMin,bMax); 
   bh_an->WriteFile("tenthirty_an.bh"); 
   
   std::vector<double> B_an,H_an; 
   bh_an->GetData(B_an,H_an); 

   // compare to GEP12long BH curve 
   util_df::CSVManager *gep = new util_df::CSVManager("tsv");
   gep->ReadFile("./data/bh/tenthirty_gep.bh",true);
   std::vector<double> bb,hh; 
   gep->GetColumn_byIndex<double>(0,bb);   
   gep->GetColumn_byIndex<double>(1,hh);   
 
   NonLinearBHCurve *bhGEp = new NonLinearBHCurve(); 
   bhGEp->SetData(bb,hh);  

   // mild steel 
   util_df::CSVManager *mild = new util_df::CSVManager("tsv");
   mild->ReadFile("./data/bh/mildaverage.bh",true);
   std::vector<double> bbm,hhm; 
   mild->GetColumn_byIndex<double>(0,bbm);   
   mild->GetColumn_byIndex<double>(1,hhm);  
 
   NonLinearBHCurve *bhMil = new NonLinearBHCurve(); 
   bhMil->SetData(bbm,hhm);  

   // default steel 
   util_df::CSVManager *def = new util_df::CSVManager("tsv");
   def->ReadFile("./data/bh/default.bh",true);
   std::vector<double> bbd,hhd; 
   def->GetColumn_byIndex<double>(0,bbd);   
   def->GetColumn_byIndex<double>(1,hhd);  
   
   NonLinearBHCurve *bhDef = new NonLinearBHCurve(); 
   bhDef->SetData(bbd,hhd);  

   std::string alloy = "1008";
   char inpath_test[200]; 
   sprintf(inpath_test,"./data/bh/%s.bh",alloy.c_str()); 
   NonLinearBHCurve *bhTest = new NonLinearBHCurve(alloy.c_str(),0,0,0,unts);
   bhTest->LoadData(inpath_test,"tsv",true,kConvertMKStoCGS);
   bhTest->Print();
   bhTest->WriteFile("tenzeroeight.bh");

   std::string xAxis = "H";
   std::string yAxis = "B"; 

   TGraph *g = GetTGraph(bh,xAxis,yAxis);
   util_df::Graph::SetParameters(g,20,kBlack);
   
   TGraph *gGEP = GetTGraph(bhGEp,xAxis,yAxis);
   util_df::Graph::SetParameters(gGEP,20,kBlue); 

   TGraph *gDEF = GetTGraph(bhDef,xAxis,yAxis); 
   util_df::Graph::SetParameters(gDEF,20,kViolet+2);

   TGraph *gMIL = GetTGraph(bhMil,xAxis,yAxis); 
   util_df::Graph::SetParameters(gMIL,20,kGreen+2);

   TGraph *gTEST = GetTGraph(bhTest,xAxis,yAxis); 
   util_df::Graph::SetParameters(gTEST,20,kMagenta);

   TGraph *gAN = util_df::Graph::GetTGraph(H_an,B_an); 
   TGraph *gCR = util_df::Graph::GetTGraph(H_cr,B_cr);
   util_df::Graph::SetParameters(gAN,21,kBlack);
   util_df::Graph::SetParameters(gCR,20,kRed);

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(g   ,"lp");  
   mg->Add(gGEP,"lp");  
   mg->Add(gCR ,"lp"); 
   mg->Add(gDEF,"lp"); 
   mg->Add(gMIL,"lp");
   mg->Add(gTEST,"lp"); 

   TLegend *LL = new TLegend(0.6,0.6,0.8,0.8);
   LL->AddEntry(g    ,"1030 (DF calc, annealed)"   ,"p"); 
   LL->AddEntry(gCR  ,"1030 (DF calc, cold-rolled)","p"); 
   LL->AddEntry(gGEP ,"GEp Build (Opera)"          ,"p"); 
   LL->AddEntry(gDEF ,"Default (Opera)"            ,"p"); 
   LL->AddEntry(gMIL ,"Mild Avg (Opera)"           ,"p"); 
   LL->AddEntry(gTEST,alloy.c_str()                ,"p"); 

   TMultiGraph *mgc = new TMultiGraph();
   mgc->Add(gAN,"lp"); 
   mgc->Add(gCR,"lp"); 

   TLegend *L = new TLegend(0.6,0.6,0.8,0.8); 
   L->AddEntry(gAN,"Annealed"   ,"p"); 
   L->AddEntry(gCR,"Cold Rolled","p"); 

   TString TitleC     = Form("B vs H for 1030 Steel");
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

   TCanvas *c2 = new TCanvas("c2","B vs H (annealed vs cold rolled)",1000,600);
 
   c2->cd();
   gPad->SetLogx();  
   mgc->Draw("a"); 
   util_df::Graph::SetLabels(mgc,TitleC,xAxisTitle,yAxisTitle); 
   mgc->GetXaxis()->SetLimits(xMin,xMax); 
   mgc->Draw("a");
   L->Draw("same"); 
   c2->Update();

   delete bh; 
   delete bh_cr; 
   delete bh_an; 

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
