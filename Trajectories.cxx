// plot Opera-generated trajectory files (*.lp) 

#include <cstdlib>
#include <iostream>

#include "UtilDFGraph.hh"
#include "UtilDFFunc.hh"
#include "CSVManager.hh"

#include "./include/operaPoint.hh"
#include "./src/operaUtilities.cc"

bool gIsDebug = false; 

int getTrack(int track,std::vector<opera::point_t> data,std::vector<opera::point_t> &out); 

void findExtrema(std::vector<opera::point_t> track,std::vector<double> &z,std::vector<double> &minR,std::vector<double> &maxR); 

int Trajectories(){

   TString Title      = Form("Particle Trajectories");    // GMn 
   TString xAxisTitle = Form("z [cm]");
   TString yAxisTitle = Form("r [cm]");

   double xMin = -1500; 
   double xMax = 0; 

   double yMin = 0; 
   double yMax = 1; 

   int c=0; 
   int lineWidth = 2;
   double markerSize = 0.5; 
  
   TString label;  
   char inpath[200],outpath[200],fileName[200];
   sprintf(fileName,"GMN-BL03-SBS-22-5_tracks_correctors-zeroed_coil-update_lgv2_new_E-400");
   sprintf(inpath  ,"./data/tracks/%s.lp" ,fileName);
   sprintf(outpath ,"./output/%s_band.csv",fileName);
 
   std::vector<opera::point_t> tracks; 
   int rc = opera::ReadTrackFile(inpath,tracks); 
   if(rc!=0) return 1; 

   const int NN = tracks.size();
   int numTracks = tracks[0].ID;
   for(int i=0;i<NN;i++) if(tracks[i].ID!=numTracks) numTracks++; 
   std::cout << "Found " << numTracks << " tracks in file " << inpath << std::endl;  

   const int NTracks = numTracks;
   TGraph **g      = new TGraph*[NTracks]; 
   TMultiGraph *mg = new TMultiGraph();
   TLegend *L      = new TLegend(0.6,0.6,0.8,0.8);

   int color[10] = {kBlack ,kRed  ,kBlue  ,kGreen+1,kOrange+1,
                    kViolet,kRed+2,kBlue+2,kGreen+3,kOrange+3};

   vector<double> minR,maxR,Z;

   std::vector<opera::point_t> aTrack;  

   for(int i=0;i<NTracks;i++){
      std::cout << Form("Parsing track %02d",i+1) << std::endl;
      // set color 
      c = i+1; 
      if(i<10) c = color[i];
      // get the track 
      rc = getTrack(i+1,tracks,aTrack);
      if(rc!=0){
	 std::cout << "Error! No points found! " << std::endl;
	 exit(1);
      }
      // create graph 
      g[i] = opera::GetTGraph(aTrack,"Z","R"); 
      util_df::Graph::SetParameters(g[i],20,c,markerSize,lineWidth); 
      mg->Add(g[i],"l"); 
      // add to legend 
      label = Form("Track %02d",i+1);
      L->AddEntry(g[i],label,"l"); 
      // find extrema 
      findExtrema(aTrack,Z,minR,maxR);
      // reset for next file
      aTrack.clear(); 
   }

   // print band to file
   const int NROW = Z.size();
   const int NCOL = 3; 
   util_df::CSVManager *data = new util_df::CSVManager();
   data->InitTable(NROW,NCOL);
   data->SetHeader("Z(cm),minR(cm),maxR(cm)");
   data->SetColumn<double>(0,Z   ); 
   data->SetColumn<double>(1,minR); 
   data->SetColumn<double>(2,maxR);
   data->WriteFile(outpath); 

   TGraph *gMin = util_df::Graph::GetTGraph(Z,minR);  
   TGraph *gMax = util_df::Graph::GetTGraph(Z,maxR);
 
   util_df::Graph::SetParameters(gMin,20,kViolet+1,0.5,lineWidth+1);  
   util_df::Graph::SetParameters(gMax,20,kViolet+1,0.5,lineWidth+1);  

   TGraphAsymmErrors *gBand = util_df::Graph::GetBand(Z,minR,maxR,kViolet+1);

   mg->Add(gMin,"l"); 
   mg->Add(gMax,"l");

   TCanvas *c1 = new TCanvas("c1","Trajectories",1000,500); 

   c1->cd();
   mg->Draw("a"); 
   util_df::Graph::SetLabels(mg,Title,xAxisTitle,yAxisTitle); 
   mg->GetXaxis()->SetLimits(xMin,xMax);  
   mg->GetYaxis()->SetRangeUser(yMin,yMax);  
   mg->Draw("a");
   gBand->Draw("e4"); 
   L->Draw("same"); 
   c1->Update(); 

   return 0;
}
//______________________________________________________________________________
int getTrack(int trackNo,std::vector<opera::point_t> data,std::vector<opera::point_t> &out){
   // get track info for a single track ID 
   const int N = data.size();
   for(int i=0;i<N;i++){
      if(data[i].ID==trackNo){
	 out.push_back(data[i]); 
      }
   }

   int NN = out.size(); 
   if(NN==0){
      return 1;
   }else{
      if(gIsDebug) std::cout << "Found " << NN << " points" << std::endl;
   }
 
   return 0;
}
//______________________________________________________________________________
void findExtrema(std::vector<opera::point_t> track,std::vector<double> &z,std::vector<double> &minR,std::vector<double> &maxR){
   // find min and max R as a function of z
   z.clear(); // will always be the same z values 

   const int N = track.size();
   double arg=0;
   std::vector<double> r;
   for(int i=0;i<N;i++){
      arg = TMath::Sqrt( track[i].x*track[i].x + track[i].y*track[i].y);
      r.push_back(arg);
      z.push_back(track[i].z);  
   } 

   double MIN  =  1E+7; 
   double MAX  = -1E-7; 

   if(track[0].ID==1){
      // first track, will use this to establish initial min and max r 
      for(int i=0;i<N;i++){
	 if(r[i]>MAX) MAX = r[i];
	 if(r[i]<MIN) MIN = r[i];
	 minR.push_back(MIN); 
	 maxR.push_back(MAX);
	 // reset for next z 
	 MAX = -1E-7;
	 MIN =  1E+7;  
      }
   }else{
      // min,max r vectors established already, just update
      for(int i=0;i<N;i++){
	 if(r[i]>maxR[i]) maxR[i] = r[i];
	 if(r[i]<minR[i]) minR[i] = r[i];
      }
   }
  
}
