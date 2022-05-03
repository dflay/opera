// Test integration 

#include <cstdlib>
#include <iostream>

#include "TGraph2D.h"

#include "CSVManager.hh"
#include "UtilDFMath.hh"
#include "UtilDFGraph.hh"

double lineFunc(const double x); 
double sinFunc(const double x); 

int TestInt(){
  
   double clkFreq = 10E+6;  
   double freq    = 100E+3;
   double omega   = 2.*TMath::Pi()*freq; 

   const int NPTS = 200; 
   double zMax    = 1000.; 
   double zMin    = 0; 
   double step    = (zMax-zMin)/( (double)NPTS );
   double A       = 1;
   double z=0,arg=0,sum=0;
   std::vector<double> Z,Y;   
   for(int i=0;i<NPTS;i++){
      z    = (zMin + ( (double)i )*step)/clkFreq;
      arg  = A*TMath::Sin(omega*z); 
      sum += arg*step; 
      Z.push_back(z);  
      Y.push_back(arg);  
   }
   sum += 0;  

   double zz    = Z[NPTS-1]; 
   double myInt = util_df::Math::SimpsonIntegral(&sinFunc,0,zz); 
   std::cout << Form("sum = %.3E, simpson = %.3E",sum,myInt) << std::endl; 

   double myInt_line = util_df::Math::SimpsonIntegral(&lineFunc,0,1); 
   std::cout << Form("Integral of line: simpson = %.3E",myInt_line) << std::endl; 

   TGraph *g = util_df::Graph::GetTGraph(Z,Y); 
   util_df::Graph::SetParameters(g,20,kBlack,0.5,2); 
 
   TCanvas *c1 = new TCanvas("c1","Test Integral",1000,500);

   c1->cd();
   g->Draw("alp");
   c1->Update(); 

   return 0;
}
//______________________________________________________________________________
double lineFunc(const double x){
   // test function f(x) = 2x + 3 
   double f = 2.*x + 3; 
   return f;  
}
//______________________________________________________________________________
double sinFunc(const double x){
   double A     = 1; 
   double freq  = 20E+3; 
   double omega = 2.*TMath::Pi()*freq;
   double f     = A*TMath::Sin(omega*x);
   return f;  
}
