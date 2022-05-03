#include <cstdlib> 

#include "CSVManager.hh"
#include "UtilDFGraph.hh"

#include "./include/operaPoint.hh"
#include "./include/operaParameters.hh"
#include "./include/operaMagnetCoil.hh"

#include "./src/operaUtilities.cc"
#include "./src/operaTrack.cc"

int Test(){

   // opera::point_t data;
   // std::cout << Form("ID = %d, x = %.3E, y = %.3E, z = %.3E",data.ID,data.x,data.y,data.z) << std::endl; 

   // opera::Track *aTrack = new opera::Track(1);
   // aTrack->PushBackCoordinate(0,1,2); 
   // delete aTrack;  

   char inpath[200];
   sprintf(inpath,"./data/fin/GMN-BL03-SBS-38-4_final_correctors-1000_drives-1_coil-update_lgv2.res");
   // sprintf(inpath,"./data/fin/GMN-13_final_correctors-1000_drives-4_coil-update.res");
 
   opera::parameters_t pars;
   int rc = opera::ReadResFile(inpath,pars); 
   if(rc==0){
      opera::PrintParameters(pars);
   }else{
      std::cout << "Error!" << std::endl;
      return rc;
   }

   // compute currents 
   std::string inpath_c = "./input/coil-dimensions.csv"; 
   std::vector<opera::magnetCoil_t> MC; 
   opera::CalculateCurrents(pars,MC,inpath_c);

   const int NC = MC.size(); 
   for(int i=0;i<NC;i++) opera::PrintCoil(MC[i]);   

   // old way 
   // determine corrector current densities and currents 
   std::vector<std::string> ll; 
   std::vector<double> Icm2,I;
   std::vector<double> usl,usr,dsl,dsr; 
   std::vector<double> usli,usri,dsli,dsri; 
   opera::CalculateCorrectorCurrents(pars,ll,Icm2,I);
   // current densities [A/cm2]  
   usr.push_back(Icm2[0]); usl.push_back(Icm2[1]);
   dsr.push_back(Icm2[2]); dsl.push_back(Icm2[3]);
   // currents [A]  
   usri.push_back(I[0]); usli.push_back(I[1]);
   dsri.push_back(I[2]); dsli.push_back(I[3]);

   std::cout << "OLD WAY: " << std::endl;
   const int NN = usr.size();
   for(int i=0;i<NN;i++){
      std::cout << Form("US-L = %.2lf A, US-R = %.2lf A, DS-L = %.2lf A, DS-R = %.2lf A",usli[i],usri[i],dsli[i],dsri[i]) << std::endl;
   }

   util_df::CSVManager *data = new util_df::CSVManager(); 
   data->ReadFile("./input/opt-gmn_lgv2.csv",true); 
   data->Print();  

   std::vector<double> USL,USR,DSL,DSR,xInt,yInt; 
   data->GetColumn_byName<double>("USL" ,USL);   
   data->GetColumn_byName<double>("USR" ,USR);   
   data->GetColumn_byName<double>("DSL" ,DSL);   
   data->GetColumn_byName<double>("DSR" ,DSR);   
   data->GetColumn_byName<double>("xInt",xInt);  
   data->GetColumn_byName<double>("yInt",yInt);  

   double arg=0;
   const int N = USL.size();
   std::vector<double> It,xyInt; 
   for(int i=0;i<N;i++){
      // arg = fabs(USL[i]) + fabs(USR[i]) + fabs(DSL[i]) + fabs(DSR[i]);
      arg = DSL[i]; // fabs(USL[i]) + fabs(USR[i]);
      xyInt.push_back( xInt[i]+yInt[i] );
      It.push_back(arg); 
   }

   TGraph *g = util_df::Graph::GetTGraph(It,xyInt);
   util_df::Graph::SetParameters(g,20,kBlack,1,1); 

   TCanvas *c1 = new TCanvas("c1","Test",1000,600); 

   c1->cd(); 
   g->Draw("ap");
   util_df::Graph::SetLabels(g,"Bx Integral vs DSL Corrector Current","I [A/cm^{2}]","#int (Bx+By) dz [G-cm]"); 
   g->Draw("ap");
   c1->Update(); 

   return rc;
}
