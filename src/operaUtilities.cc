#include "../include/operaUtilities.hh"
//______________________________________________________________________________
namespace opera { 
   //______________________________________________________________________________
   void CalculateCurrents(parameters_t pars,std::vector<magnetCoil_t> &data,std::string inpath){
      // compute magnet coil currents from the parameters object

      // load coil dimensions
      util_df::CSVManager *csv = new util_df::CSVManager();
      csv->ReadFile(inpath.c_str(),true); 

      std::vector<std::string> name; 
      csv->GetColumn_byName_str("name",name);

      const int M = name.size();

      std::vector<double> X,Y,NTURNS,NCOILS; 
      csv->GetColumn_byName<double>("xlen(cm)",X     ); 
      csv->GetColumn_byName<double>("ylen(cm)",Y     ); 
      csv->GetColumn_byName<double>("nturns"  ,NTURNS); 
      csv->GetColumn_byName<double>("ncoils"  ,NCOILS);
      delete csv; 

      // get values for the racetrack and bedstead  
      int NI = pars.iLabel.size(); 
      double Icm2_rt=0,Icm2_bs=0,Icm2_cc=0;
      for(int i=0;i<NI;i++){
	 if(pars.iLabel[i].compare("arc")==0)        Icm2_rt = pars.I[i]; 
	 if(pars.iLabel[i].compare("bedstead")==0)   Icm2_bs = pars.I[i]; 
	 if(pars.iLabel[i].compare("correctors")==0) Icm2_cc = pars.I[i]; 
      }

      // get drive factors
      double sf_rt=0,sf_bs=0; 
      int ND = pars.dLabel.size(); 
      for(int i=0;i<ND;i++){
	 if(pars.dLabel[i].compare("right")==0) sf_rt = pars.drive[i]; 
	 if(pars.dLabel[i].compare("left")==0)  sf_bs = pars.drive[i]; 
       }

      // now get coil dimensions
      double x_bs=0,y_bs=0,nt_bs=0,nc_bs=0;
      double x_rt=0,y_rt=0,nt_rt=0,nc_rt=0;
      for(int i=0;i<M;i++){
	 if(name[i].compare("bedstead")==0){
	    x_bs  = X[i]; 
	    y_bs  = Y[i]; 
	    nt_bs = NTURNS[i];
	    nc_bs = NCOILS[i];
         }
         if(name[i].compare("racetrack")==0){
	    x_rt  = X[i]; 
	    y_rt  = Y[i]; 
	    nt_rt = NTURNS[i];
	    nc_rt = NCOILS[i];
         } 
      } 
 
      opera::magnetCoil_t bs;
      bs.name           = "bedstead";  
      bs.xSize          = x_bs; 
      bs.ySize          = y_bs; 
      bs.nCoils         = nc_bs; 
      bs.nTurns         = nt_bs;
      bs.driveScale     = sf_bs; 
      bs.currentDensity = Icm2_bs;  
      bs.current        = sf_bs*Icm2_bs*x_bs*y_bs/(nc_bs*nt_bs); 

      opera::magnetCoil_t rt; 
      rt.name           = "racetrack";  
      rt.xSize          = x_rt; 
      rt.ySize          = y_rt; 
      rt.nCoils         = nc_rt; 
      rt.nTurns         = nt_rt;
      rt.driveScale     = sf_rt; 
      rt.currentDensity = Icm2_rt;  
      rt.current        = sf_rt*Icm2_rt*x_rt*y_rt/(nc_rt*nt_rt); 

      // save the bedstead and racetrack 
      data.push_back(bs);
      data.push_back(rt); 
 
      // loop over all correctors 
      std::vector<std::string> clabel; 
      clabel.push_back("corr_usr"); 
      clabel.push_back("corr_usl"); 
      clabel.push_back("corr_dsr"); 
      clabel.push_back("corr_dsl");
 
      opera::magnetCoil_t mc;
      int NC = clabel.size(); 
      double Icm2=0,I=0,nt=0,nc=0,x=0,y=0,sf=0;
      for(int i=0;i<NC;i++){
         // search in the file 
	 for(int j=0;j<M;j++){
	    if(clabel[i].compare(name[j])==0){
	       x  = X[j];
	       y  = Y[j];
	       nt = NTURNS[j];
               nc = NCOILS[j];  
            }
	 }
	 // search pars for drives
         for(int j=0;j<ND;j++){
	    if(clabel[i].compare(pars.dLabel[j])==0) sf = pars.drive[j]; 
         } 
	 // now put it together 
	 Icm2 = Icm2_cc; 
	 I    = sf*Icm2*(x*y)/(nt*nc);
	 // fill data struct
	 mc.name           = clabel[i]; 
	 mc.currentDensity = Icm2; 
	 mc.current        = I;  
	 mc.xSize          = x; 
	 mc.ySize          = y;
	 mc.nTurns         = nt; 
	 mc.nCoils         = nc;
         mc.driveScale     = sf;  
	 data.push_back(mc);
      }
       
   }
   //______________________________________________________________________________
   void CalculateCorrectorCurrents(parameters_t pars,std::vector<std::string> &label,
                                   std::vector<double> &Icm2,std::vector<double> &I){
      // extract corrector currents in (A/cm2 and A) from the parameters object
      double cc=0,I_corr=0;
      std::vector<double> icm2;
      std::vector<std::string> ll;
      // determine corrector current density (A/cm^2)  
      int NC = pars.I.size();
      for(int j=0;j<NC;j++) if(pars.iLabel[j].compare("correctors")==0) I_corr = pars.I[j];
      // use drive scale factors to get the true densities
      int ND = pars.drive.size();
      // note we account for arbitrary ordering and all drive factors from the res file 
      double Icm2_usl=0,Icm2_usr=0,Icm2_dsl=0,Icm2_dsr=0;
      for(int j=0;j<ND;j++){
         cc = I_corr*pars.drive[j];
         // std::cout << I_corr << " " << pars.dLabel[j] << " " << cc << " " << pars.drive[j] << std::endl; 
         if(pars.dLabel[j].compare("corr_usr")==0) Icm2_usr = cc;
         if(pars.dLabel[j].compare("corr_usl")==0) Icm2_usl = cc;
         if(pars.dLabel[j].compare("corr_dsr")==0) Icm2_dsr = cc;
         if(pars.dLabel[j].compare("corr_dsl")==0) Icm2_dsl = cc;
      }
      // note the fixed order! 
      Icm2.push_back(Icm2_usr); label.push_back("USR");
      Icm2.push_back(Icm2_usl); label.push_back("USL");
      Icm2.push_back(Icm2_dsr); label.push_back("DSR");
      Icm2.push_back(Icm2_dsl); label.push_back("DSL");
      // now get the current in A
      CalculateCorrectorCurrents(label,Icm2,I);  
   }
   //______________________________________________________________________________
   void CalculateCorrectorCurrents(std::vector<std::string> label,std::vector<double> Icm2,std::vector<double> &I){
      // input: vector of corrector current densities (A/cm^2)  
      // output: vector of corrector currents (A)  

      // load corrector dimensions
      util_df::CSVManager *data = new util_df::CSVManager();
      data->ReadFile("./input/csv/corr-dim.csv",true); 

      std::vector<std::string> name; 
      data->GetColumn_byName_str("corr",name);

      std::vector<double> xl,yl,nturns,ncoils; 
      data->GetColumn_byName<double>("xlen(cm)",xl    ); 
      data->GetColumn_byName<double>("ylen(cm)",yl    ); 
      data->GetColumn_byName<double>("nturns"  ,nturns); 
      data->GetColumn_byName<double>("ncoils"  ,ncoils); 

      std::cout << "Calculating corrector currents: " << std::endl;

      int iusl=0,iusr=0,idsl=0,idsr=0;
      const int N = name.size(); 
      for(int i=0;i<N;i++){
	 // std::cout << name[i] << std::endl;
	 if(name[i].compare("USL")==0) iusl = i; 
	 if(name[i].compare("USR")==0) iusr = i; 
	 if(name[i].compare("DSL")==0) idsl = i; 
	 if(name[i].compare("DSR")==0) idsr = i; 
      }
     
      std::vector<int> index;
      index.push_back(iusr);  
      index.push_back(iusl);  
      index.push_back(idsr);  
      index.push_back(idsl);  
 
      int j=0;
      double area=0,nt=0,current=0;
      for(int i=0;i<N;i++){
	 j       = index[i]; 
         // std::cout << i << " " << j << " " << name[i] << " " << Icm2[i] << std::endl;
	 area    = (ncoils[j]*xl[j])*yl[j];   // WARNING: factor of ncoils on x since coils are stacked along x! 
	 nt      = ncoils[j]*nturns[j];
	 current = Icm2[j]*area/(nt);
	 I.push_back(current);
      }

      // print for confirmation 
      for(int i=0;i<N;i++){
	 std::cout << Form("  %s: rho = %.3lf A/cm^2, I = %.3lf A",name[i].c_str(),Icm2[i],I[i]) << std::endl;
      }

      delete data;
   }
   //______________________________________________________________________________
   void PrintParameters(parameters_t data){
      std::cout << "---------------- Opera Parameters ----------------" << std::endl;
      std::cout << Form("Number of active elements: %d",data.NActiveElem)      << std::endl; 
      std::cout << Form("Number of nodes:           %d",data.NNodes)           << std::endl; 
      std::cout << Form("Number of equations:       %d",data.NEqns)            << std::endl; 
      std::cout << Form("Number of non-zeros:       %d",data.NNonZeros)        << std::endl;
      std::cout << Form("Symmetry:                  %s",data.symmetry.c_str()) << std::endl;
 
      std::cout << "Current Densities [A/cm^2]: " << std::endl; 
      int NC = data.I.size();
      for(int i=0;i<NC;i++) std::cout << Form("   %s: %.4lf",data.iLabel[i].c_str(),data.I[i]) << std::endl; 

      std::cout << "Drive Scale Factors: " << std::endl; 
      int ND = data.drive.size();
      for(int i=0;i<ND;i++) std::cout << Form("   %s: %.4lf",data.dLabel[i].c_str(),data.drive[i]) << std::endl; 
      std::cout << "--------------------------------------------------" << std::endl;
   }
   //______________________________________________________________________________
   void PrintCoil(magnetCoil_t mc){
         std::cout << Form("------------ MAGNET COIL DATA ------------") << std::endl;
         std::cout << Form("name   = %s"   ,mc.name.c_str()  ) << std::endl;
         std::cout << Form("Icm2   = %.2lf",mc.currentDensity) << std::endl;
         std::cout << Form("I      = %.2lf",mc.current       ) << std::endl;
         std::cout << Form("x      = %.2lf",mc.xSize         ) << std::endl;
         std::cout << Form("y      = %.2lf",mc.ySize         ) << std::endl;
         std::cout << Form("nTurns = %d"   ,mc.nTurns        ) << std::endl;
         std::cout << Form("nCoils = %d"   ,mc.nCoils        ) << std::endl;
         std::cout << Form("drive  = %.3lf",mc.driveScale    ) << std::endl;
   }
   //______________________________________________________________________________
   int ReadResFile(const char *inpath,parameters_t &data,bool isDebug){
      // read the *.res file and gather relevant parameters
      // NOTE: labels for currents and drives are specific to SBS

      std::string aLine,junk,sym,label;
      std::vector<std::string> col,dl;  
      char phrase[200]; 

      int lineNum=0,M=0;

      bool foundConductors = false; 
      bool foundDrives     = false; 
      bool foundModelSize  = false;
      bool foundSymmetry   = false;

      // currents 
      double ii=0;
      std::vector<double> I; 

      // drive settings 
      double dd=0;
      std::vector<double> D; 
    
      // model size info 
      int NElem=0,NNode=0,NEqn=0,NNZ=0; 

      int NLines = 14; // 12; // 12 for GMn 

      char msg[200]; 
 
      std::ifstream infile;
      infile.open(inpath);
      if( infile.fail() ){
	 std::cout << "Cannot open the file: " << inpath << std::endl;
	 return 1;
      }else{
	 while( !infile.eof() ){
	    std::getline(infile,aLine);
	    lineNum++;
	    // find conductor info 
            sprintf(phrase,"Biot-Savart conductors"); 
            foundConductors = util_df::IsInCharString(aLine.c_str(),phrase); 
            if(foundConductors){
	       if(isDebug) std::cout << "Found conductors" << std::endl;
	       // read 8 lines; odd lines have current densities 
               for(int i=0;i<8;i++){
		  std::getline(infile,aLine);
                  if(i%2!=0){
		     if(isDebug) std::cout << aLine << std::endl;
		     util_df::SplitString(':',aLine,col);
                     M = col.size();
		     if(isDebug) for(int j=0;j<M;j++) std::cout << col[j] << std::endl; 
                     if(i==1) ii = std::atof(col[1].c_str()); 
                     if(i==3) ii = std::atof(col[1].c_str()); 
                     if(i==5) ii = std::atof(col[1].c_str()); 
                     if(i==7) ii = std::atof(col[1].c_str());
		     I.push_back(ii); 
		     col.clear();
                  }
		  lineNum++;
               }
	       foundConductors = false;
            }
            // find drive settings 
            sprintf(phrase,"Drive sets and functions");
            foundDrives = util_df::IsInCharString(aLine.c_str(),phrase);
            if(foundDrives){
	       if(isDebug) std::cout << "[operaUtilities::ReadResFile]: Found drives" << std::endl;
               // read NLines lines; odd lines have the scaling factor 
	       for(int i=0;i<NLines;i++){
		  std::getline(infile,aLine);
		  util_df::SplitString(':',aLine,col);
		  if(isDebug){
		     std::cout << aLine << std::endl;
                     std::cout << "Parsed: [";
		     M = col.size();
		     if(M>0){
			for(int j=0;j<M-1;j++){
			   util_df::RemoveWhiteSpace(col[j]);
			   std::cout << col[j] << ",";
                        }
			util_df::RemoveWhiteSpace(col[M-1]);
			std::cout << col[M-1];
		     }
		     std::cout << "]" << std::endl;
		     std::cout << "------" << std::endl;
		  } 
                  if(i%2!=0){
		     if(i==1)  dd = std::atof(col[2].c_str()); 
                     if(i==3)  dd = std::atof(col[2].c_str()); 
                     if(i==5)  dd = std::atof(col[2].c_str()); 
                     if(i==7)  dd = std::atof(col[2].c_str()); 
                     if(i==9)  dd = std::atof(col[2].c_str()); 
                     if(i==11) dd = std::atof(col[2].c_str()); 
                     if(i==13) dd = std::atof(col[2].c_str()); 
		     D.push_back(dd);
                  }else{
		     if(i==0)  label = col[0]; 
                     if(i==2)  label = col[0]; 
                     if(i==4)  label = col[0]; 
                     if(i==6)  label = col[0]; 
                     if(i==8)  label = col[0]; 
                     if(i==10) label = col[0];
                     if(i==12) label = col[0];
		     dl.push_back(label); 
                  }
		  col.clear();
		  lineNum++;
	       }
	       foundDrives = false;
            }
	    // find model details  
            sprintf(phrase,"Model size information"); 
            foundModelSize = util_df::IsInCharString(aLine.c_str(),phrase); 
            if(foundModelSize){
	       if(isDebug) std::cout << "[operaUtilities::ReadResFile]: Found model information" << std::endl;
	       for(int i=0;i<5;i++){
		  std::getline(infile,aLine);  
                  if(i>0){
		     util_df::SplitString(':',aLine,col);
                     if(i==1) NElem = std::atof(col[1].c_str()); 
                     if(i==2) NNode = std::atof(col[1].c_str()); 
                     if(i==3) NEqn  = std::atof(col[1].c_str()); 
                     if(i==4) NNZ   = std::atof(col[1].c_str()); 
                  } 
		  col.clear();
		  lineNum++;
	       }
            }
            // find model symmetry 
            foundSymmetry = util_df::IsInCharString(aLine.c_str(),"Symmetry"); 
            if(foundSymmetry){
	       if(isDebug) std::cout << "[operaUtilities::ReadResFile]: Found model symmetry" << std::endl;
	       std::getline(infile,aLine); 
	       sym = aLine;
	       lineNum++; 
	    }
         } 
      }
   
      std::vector<std::string> il; 
      il.push_back("correctors");
      il.push_back("bedstead");  
      il.push_back("arc");  
      il.push_back("straight");  

      // dl.push_back("right");
      // dl.push_back("corr_usr");  
      // dl.push_back("corr_usl");  
      // dl.push_back("left");
      // dl.push_back("corr_dsl");  
      // dl.push_back("corr_dsr"); 

      // standardize the drive names 
      M = dl.size(); 
      for(int i=0;i<M;i++){
	 util_df::RemoveWhiteSpace(dl[i]);
         if(dl[i].compare("RIGHT")==0  || dl[i].compare("right")==0)  dl[i] = "right";  
         if(dl[i].compare("LEFT")==0   || dl[i].compare("left")==0 )  dl[i] = "left";  
         if(dl[i].compare("COIL")==0   || dl[i].compare("coil")==0)   dl[i] = "coil";  
         if(dl[i].compare("SADDLE")==0 || dl[i].compare("saddle")==0) dl[i] = "saddle";  
	 if(dl[i].compare("COMP2R")==0 || dl[i].compare("CORUR")==0)  dl[i] = "corr_usr"; 
	 if(dl[i].compare("COMP2L")==0 || dl[i].compare("CORUL")==0)  dl[i] = "corr_usl"; 
	 if(dl[i].compare("COMP3R")==0 || dl[i].compare("CORDR")==0)  dl[i] = "corr_dsr"; 
	 if(dl[i].compare("COMP3L")==0 || dl[i].compare("CORDL")==0)  dl[i] = "corr_dsl"; 
	 if(dl[i].compare("CUR01")==0 )                               dl[i] = "cur01"; 
      } 

      // set struct info 
      data.symmetry    = sym; 
      data.NActiveElem = NElem; 
      data.NNodes      = NNode; 
      data.NEqns       = NEqn; 
      data.NNonZeros   = NNZ; 
 
      int NC = I.size();
      for(int i=0;i<NC;i++){ 
	 data.I.push_back(I[i]);     
	 data.iLabel.push_back(il[i]);    
      } 
 
      int ND = D.size();  
      for(int i=0;i<ND;i++){
	 data.drive.push_back(D[i]);     
	 data.dLabel.push_back(dl[i]);   
      } 

      if(isDebug) PrintParameters(data);  

      return 0;
   }
   //______________________________________________________________________________
   int ReadTrackFile(const char *inpath,std::vector<point_t> &track,bool isDebug){
      // read the file line by line 
      // when the phrase "Track number i" (i = 1, 2,...) is matched, 
      // read the immediately following lines into a vector  
      // output 
      // - vector of *all* tracks; vector elements are of type opera::point_t  

      std::string aLine,junk;
      std::vector<std::string> line,col,val;

      point_t aPoint;

      std::string here = "[operaUtilities::ReadTrackFile]"; 
    
      int lineNum=0,trackNum=1,M=0,NC=0,NP=0,trackCnt=0;
      char phrase[200],msg[200],ans='t';

      bool foundTrack = false;
      bool foundNStep = false;
      bool foundEntry = false;

      std::ifstream infile;
      infile.open(inpath);
      if( infile.fail() ){
	 std::cout << "Cannot open the file: " << inpath << std::endl;
	 return 1;
      }else{
         std::cout << here << ": Reading file: " << inpath << std::endl;
	 while( !infile.eof() ){
	    std::getline(infile,aLine);
	    // find number of points 
	    sprintf(phrase,"NSTEP");
	    foundNStep = util_df::IsInCharString(aLine.c_str(),phrase);
	    if(foundNStep){
	       if(isDebug) std::cout << here << ": " << aLine << std::endl;
	       util_df::SplitString_whiteSpace(aLine,col);
	       NC = col.size();
	       for(int i=0;i<NC;i++){
		  if(isDebug) std::cout << here << ": " << col[i] << std::endl;
		  foundEntry = util_df::IsInCharString(col[i].c_str(),"NSTEP");
		  if(foundEntry){
		     util_df::SplitString('=',col[i],val);
		     NP = std::atof( val[1].c_str() );
		     std::cout << here << ": Found number of points = " << NP << std::endl;
		     foundEntry = false;
		     val.clear();
		  }
	       }
	       col.clear();
	    }
	    // find track 
	    sprintf(phrase,"Track number %d",trackNum);
	    foundTrack = util_df::IsInCharString(aLine.c_str(),phrase);
	    if(foundTrack){
	       sprintf(msg,"%s: Found track %d starting on line %d",here.c_str(),trackNum,lineNum);
	       if(isDebug) std::cout << msg << std::endl;
	       // skip three lines
	       for(int i=0;i<3;i++){
		  std::getline(infile,junk);
		  if(isDebug) std::cout << here << ": junk: " << junk << std::endl;
		  lineNum++;
	       }
	       // fill the vector of structs
	       for(int i=0;i<NP;i++){
		  std::getline(infile,aLine);
		  // split to vector
		  if(isDebug) std::cout << here << ": parsing: " << aLine << std::endl;
		  util_df::SplitString_whiteSpace(aLine,col);
		  // extract data 
		  aPoint.x  = std::atof( col[0].c_str() );
		  aPoint.y  = std::atof( col[1].c_str() );
		  aPoint.z  = std::atof( col[2].c_str() );
		  aPoint.vx = std::atof( col[3].c_str() );
		  aPoint.vy = std::atof( col[4].c_str() );
		  aPoint.vz = std::atof( col[5].c_str() );
		  aPoint.ID = trackNum;
		  sprintf(msg,"%s: line %d, track %d: x = %.3E, y = %.3E, z = %.3E, vx = %.3E, vy = %.3E, vz = %.3E",
			here.c_str(),lineNum,trackNum,aPoint.x,aPoint.y,aPoint.z,aPoint.vx,aPoint.vy,aPoint.vz);
		  if(isDebug) std::cout << msg << std::endl;
		  // push back the vector
		  track.push_back(aPoint);
		  // clear for next line 
		  col.clear();
		  lineNum++;
	       }
	       // increment track number  
	       trackNum++;
	       trackCnt++;
	       sprintf(msg,"%s: Stored data for %d lines, ending on line %d",here.c_str(),NP,lineNum);
	       if(isDebug) std::cout << msg << std::endl;
	    }else{
	       // do nothing; continue to next line 
	       lineNum++;
	    }
	 }
      }

      int rc=0;
      if(trackCnt==0){
	 rc = 1;
	 std::cout << here << ": ERROR! No tracks found!" << std::endl;
      }else{
	 std::cout << here << ": Found " << trackCnt << " tracks" << std::endl;
      }

      return rc;
   }
   //______________________________________________________________________________
   int GetBMod(util_df::CSVManager *data,std::vector<double> &bmod){
      // get |B|
      bmod.clear();

      std::vector<double> bx,by,bz;
      data->GetColumn_byName<double>("BX",bx);
      data->GetColumn_byName<double>("BY",by);
      data->GetColumn_byName<double>("BZ",bz);

      const int N = bx.size();
      double arg=0;
      for(int i=0;i<N;i++){
	 arg = TMath::Sqrt( bx[i]*bx[i] + by[i]*by[i] + bz[i]*bz[i] );
	 bmod.push_back(arg);
      }
      return 0;
   }
   //______________________________________________________________________________
   int ConsolidateData(util_df::CSVManager *in,std::vector<magneticFieldPt_t> &data){
      // consolidate data into a vector of a struct 
      std::vector<double> x,y,z,bx,by,bz,bmod;
      in->GetColumn_byName<double>("X" ,x);
      in->GetColumn_byName<double>("Y" ,y);
      in->GetColumn_byName<double>("Z" ,z);
      in->GetColumn_byName<double>("BX",bx);
      in->GetColumn_byName<double>("BY",by);
      in->GetColumn_byName<double>("BZ",bz);

      const int N = x.size();
      double arg=0;
      for(int i=0;i<N;i++){
	 arg = TMath::Sqrt( bx[i]*bx[i] + by[i]*by[i] + bz[i]*bz[i] );
	 bmod.push_back(arg);
      }

      magneticFieldPt_t d;
      for(int i=0;i<N;i++){
	 d.x  = x[i];  d.y  = y[i];  d.z  = z[i];
	 d.bx = bx[i]; d.by = by[i]; d.bz = bz[i]; d.bmod = bmod[i];
	 data.push_back(d);
      }

      return 0;
   }
   //______________________________________________________________________________
   double GetSBSAngle(std::string fileName){
      // extract the SBS angle from the input file name  
      std::vector<std::string> sv;
      // split on the underscore
      util_df::SplitString('_',fileName.c_str(),sv);
      int M = sv.size();
      // get the angle; it's in the first entry of the parsed string
      std::string tag = sv[0];
      std::vector<std::string> av;
      util_df::SplitString('-',tag.c_str(),av);
      // angle is the 4th and 5th entries (so, 3 and 4 in C++ indexing) 
      double a1 = std::atof(av[3].c_str());   
      double a2 = std::atof(av[4].c_str())/10.;  // this is the decimal of the angle 
      double angle = a1 + a2;  
      return angle; 
   }
   //______________________________________________________________________________
   void GetVectors(util_df::CSVManager *data,std::string yAxis,double x0,double y0,int units,
                   std::vector<double> &Z,std::vector<double> &B){
      // Get B as a function of z for a given x0 and y0 

      std::vector<double> x,y,z,b;
      data->GetColumn_byName<double>("X",x);
      data->GetColumn_byName<double>("Y",y);
      data->GetColumn_byName<double>("Z",z);

      if(yAxis.compare("bmod")==0){
	 GetBMod(data,b);
      }else{
	 data->GetColumn_byName<double>(yAxis,b);
      }

      // B units; input data are in Gauss  
      const int N = x.size();
      if(units==opera::kTesla){
	 for(int i=0;i<N;i++) b[i] /= 1E+4;
      }

      // if(conv_T_to_G){
      //    for(int i=0;i<N;i++) b[i] *= 1E+4;
      // }

      // double sf=1;
      // if(conv_mm_to_cm) sf = 1E-1;

      // choose values to fill for graph for a given (x0,y0) 
      double xdiff=0,ydiff=0;
      std::vector<double> X,Y;
      for(int i=0;i<N;i++){
	 if( abs(x[i]-x0)<1E-3 && abs(y[i]-y0)<1E-3 ){
	    Z.push_back(z[i]);
	    B.push_back(b[i]);
	 }
      }
   }
   //______________________________________________________________________________
   void GetVectors_Interp(util_df::CSVManager *data,std::string yAxis,double x0,double y0,int units,
                          std::vector<double> &Z,std::vector<double> &B){
      // grab the B field corresponding to (x0,y0) as a function of z
      // do bilinear interpolation to estimate B(x0,y0,z)  
 
      // get all x, y 
      std::vector<double> x,y;
      data->GetColumn_byName<double>("X",x);
      data->GetColumn_byName<double>("Y",y);

      // sort and strip out repeats
      util_df::Algorithm::SortedRemoveDuplicates<double>(x);  
      util_df::Algorithm::SortedRemoveDuplicates<double>(y);  

      // find bounding values for x0 and y0  
      int ixlo=0,ixhi=0,iylo=0,iyhi=0;
      util_df::Algorithm::BinarySearch<double>(x,x0,ixlo,ixhi);
      util_df::Algorithm::BinarySearch<double>(y,y0,iylo,iyhi);

      // get bounding (x,y) values
      double xl   = x[ixlo]; 
      double xh   = x[ixhi]; 
      double yl   = y[iylo]; 
      double yh   = y[iyhi]; 
      // weight factors for interpolation 
      double xwl  = (x0-xl)/(xh-xl);
      double xwh  = (xh-x0)/(xh-xl);
      double ywl  = (y0-yl)/(yh-yl);
      double ywh  = (yh-y0)/(yh-yl);

      // construct B(xlo,ylo), B(xlo,yhi), B(xhi,ylo), B(xhi,yhi) as a function of z  
      std::vector<double> bll,blh,bhl,bhh;
      GetVectors(data,yAxis,xl,yl,units,Z,bll);
      Z.clear();
      GetVectors(data,yAxis,xl,yh,units,Z,blh);
      Z.clear();
      GetVectors(data,yAxis,xh,yl,units,Z,bhl);
      Z.clear();
      GetVectors(data,yAxis,xh,yh,units,Z,bhh);

      // construct interpolated B 
      double arg=0;
      const int N = Z.size();
      for(int i=0;i<N;i++){
         arg = ywh*(xwh*bll[i] + xwl*bhl[i]) + ywl*(xwh*blh[i] + xwl*bhh[i]); 
	 B.push_back(arg);
      }

   }
   //______________________________________________________________________________
   TGraph *GetTGraph_xySlice(util_df::CSVManager *data,std::string yAxis,double x0,double y0,int units){
      // get a TGraph of B vs z; values for (x0,y0) must exist in the CSV object 
      std::vector<double> z,b; 
      GetVectors(data,yAxis,x0,y0,units,z,b); 
      TGraph *g = util_df::Graph::GetTGraph(z,b);
      return g;
   }
   //______________________________________________________________________________
   TGraph *GetTGraph(util_df::CSVManager *data,std::string yAxis,double x0,double y0,int units){
      // get a TGraph of the B(x0,y0) as a function of z
      // do bilinear interpolation to estimate B(x0,y0,z)  
      std::vector<double> z,b; 
      GetVectors_Interp(data,yAxis,x0,y0,units,z,b);
      TGraph *g = util_df::Graph::GetTGraph(z,b);
      return g;
   }
   //______________________________________________________________________________
   TGraph *GetTGraph(std::vector<point_t> data,std::string xAxis,std::string yAxis){

      double arg=0;
      std::vector<double> x,y,r;
      const int N = data.size();
      // construct the radius if we need it 
      for(int i=0;i<N;i++){
	 arg = TMath::Sqrt( data[i].x*data[i].x + data[i].y*data[i].y );
	 r.push_back(arg);
      }

      for(int i=0;i<N;i++){
	 // x axis 
	 if( xAxis.compare("X")==0 )  x.push_back(data[i].x);
	 if( xAxis.compare("Y")==0 )  x.push_back(data[i].y);
	 if( xAxis.compare("Z")==0 )  x.push_back(data[i].z);
	 if( xAxis.compare("R")==0 )  x.push_back(r[i]);
	 if( xAxis.compare("VX")==0 ) x.push_back(data[i].vx);
	 if( xAxis.compare("VY")==0 ) x.push_back(data[i].vy);
	 if( xAxis.compare("VZ")==0 ) x.push_back(data[i].vz);
	 // y axis 
	 if( yAxis.compare("X")==0 )  y.push_back(data[i].x);
	 if( yAxis.compare("Y")==0 )  y.push_back(data[i].y);
	 if( yAxis.compare("Z")==0 )  y.push_back(data[i].z);
	 if( yAxis.compare("R")==0 )  y.push_back(r[i]);
	 if( yAxis.compare("VX")==0 ) y.push_back(data[i].vx);
	 if( yAxis.compare("VY")==0 ) y.push_back(data[i].vy);
	 if( yAxis.compare("VZ")==0 ) y.push_back(data[i].vz);
      }

      int NX = x.size();
      int NY = y.size();
      if(NX!=NY) std::cout << "[GetTGraph]: Error! No points!" << std::endl;

      TGraph *g = util_df::Graph::GetTGraph(x,y);
      return g;
   }
   //______________________________________________________________________________
   std::string getResNameFromTrackName_E(std::string fileName,double &E,int eIndexOffset){
      // build the res filename from the track name (remove E-xyz tag) 
      // and compute the beam energy based on the tag 
      // eIndexOffset => if there's trailing characters *after* the etag, need to jump back further in string
      // generally, need to replace 'track' with 'final'
      std::vector<std::string> sv;
      // std::cout << "PARSING: " << fileName << std::endl;
      // split on the underscore
      util_df::SplitString('_',fileName.c_str(),sv);
      int M = sv.size();
      std::cout << "M = " << M << std::endl;
      for(int i=0;i<M;i++) std::cout << i << " " << sv[i] << std::endl;
      // get the beam energy: split on the dash on the last term
      std::string etag = sv[M-1-eIndexOffset];
      // std::cout << "etag = " << etag << std::endl;
      std::vector<std::string> ev;
      util_df::SplitString('-',etag.c_str(),ev);
      E = std::atof(ev[1].c_str())/100.; // tag written as "xyz" => x.yz GeV, divide by 100 
      // std::cout << "Determined E = " << E << std::endl; 
      // build suffix: remove etag  
      std::vector<std::string> svf; 
      for(int i=0;i<M;i++){
	 if(sv[i].compare(etag)!=0) svf.push_back(sv[i]);  
      } 
      const int MF = svf.size();
      std::string suffix=""; 
      for(int i=2;i<MF-1;i++) suffix += svf[i] + "_";
      suffix += svf[MF-1];  
      // std::string suffix="";
      // int start = 2; 
      // int end   = start + 2; // M - 2;
      // std::cout << "start = " << start << ", end = " << end << std::endl;
      // for(int j=start;j<end;j++) suffix += sv[j] + "_";
      // suffix += sv[end];
      // char suffix[200]; 
      // sprintf(suffix,"%s_%s_%s_%s",sv[2].c_str(),sv[3].c_str(),sv[4].c_str(),sv[5].c_str()); 
      // std::string suffix = sv[2] + "_" + sv[3] + "_" + sv[4]; //  + "_" + sv[5];
      std::cout << suffix << std::endl; 
      char path[200];
      sprintf(path,"%s_final_%s.res",sv[0].c_str(),suffix.c_str());
      std::string theStr = path;
      return theStr;
   }
   //______________________________________________________________________________
   std::string getResNameFromTrackName(std::string fileName){
      // build the res filename from the track name 
      // generally, need to replace 'track' with 'final'
      std::string suffix;
      std::vector<std::string> sv;
      util_df::SplitString('_',fileName.c_str(),sv);
      int M = sv.size();
      for(int j=2;j<M-1;j++) suffix += sv[j] + "_";
      suffix += sv[M-1];
      char path[200];
      sprintf(path,"%s_final_%s.res",sv[0].c_str(),suffix.c_str());
      std::string theStr = path;
      return theStr;
   }
} // ::opera 
