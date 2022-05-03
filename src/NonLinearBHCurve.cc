#include "../include/NonLinearBHCurve.hh"
//______________________________________________________________________________
NonLinearBHCurve::NonLinearBHCurve(const char *name,double k1,double k2,double k3,int units){
   fName  = name;
   fk1    = k1; 
   fk2    = k2; 
   fk3    = k3;
   fUnits = units;  
}
//______________________________________________________________________________
NonLinearBHCurve::~NonLinearBHCurve(){
   ClearData();
}
//______________________________________________________________________________
void NonLinearBHCurve::ClearData(){
   fName  = "NONE"; 
   fk1    = 0;
   fk2    = 0;
   fk3    = 0;
   // fUnits = util_df::Units::kCGS; 
   fB.clear();
   fH.clear();
}
//______________________________________________________________________________
double NonLinearBHCurve::GetMu(double B){
   double arg = fk1*exp(fk2*B*B) + fk3;
   double f=0;
   if(arg!=0){
      f = util_df::Constants::mu_0 + 1./arg;
   }
   return f;
}
//______________________________________________________________________________
void NonLinearBHCurve::Calculate(int N,double bMin,double bMax){
   // calculate B vs H in range (bMin,bMax)

   // first clear vectors 
   fB.clear();
   fH.clear(); 

   // scale factors (for units)
   // assume mks (B = Tesla, H = A/m)   
   double sf_b = 1; 
   double sf_h = 1;

   if(fUnits==util_df::Units::kCGS){
      sf_b = util_df::Constants::Gauss;
      sf_h = util_df::Constants::Oerstead;
   } 

   double arg=0,arg_mu=0,arg_b=0,arg_h=0;
   double step = (bMax-bMin)/( (double)N );
   for(int i=0;i<N;i++){
      arg_b  = bMin + ( (double)i )*step;
      arg_mu = GetMu(arg_b); 
      arg_h  = arg_b/arg_mu; 
      arg_b *= 1./sf_b; 
      arg_h *= 1./sf_h;
      fB.push_back(arg_b);  
      fH.push_back(arg_h);  
   } 
}
//______________________________________________________________________________
void NonLinearBHCurve::GetData(std::vector<double> &b,std::vector<double> &h){
   const int N = fB.size();
   for(int i=0;i<N;i++){
      b.push_back(fB[i]); 
      h.push_back(fH[i]); 
   }
}
//______________________________________________________________________________
void NonLinearBHCurve::SetData(std::vector<double> b,std::vector<double> h){
   const int N = b.size();
   for(int i=0;i<N;i++){
      fB.push_back(b[i]); 
      fH.push_back(h[i]); 
   }
}
//______________________________________________________________________________
void NonLinearBHCurve::Print(){
   const int N = fB.size();
   for(int i=0;i<N;i++){
      std::cout << Form("H = %.3lf, B = %.3lf",fH[i],fB[i]) << std::endl;
   }
}
//______________________________________________________________________________
int NonLinearBHCurve::LoadData(const char *inpath,const char *delim,bool isHeader,int scaling){
   // load data from a file
   util_df::CSVManager *data = new util_df::CSVManager(delim);
   int rc = data->ReadFile(inpath,isHeader);
   if(rc!=0) return rc;

   std::vector<double> b,h; 
   if(isHeader){
      data->GetColumn_byName<double>("h",h); 
      data->GetColumn_byName<double>("b",b); 
   }else{
      // assume H is the first column, B is the second column
      data->GetColumn_byIndex<double>(0,h);  
      data->GetColumn_byIndex<double>(1,b);  
   }

   // apply a scale factor if desired
   double sfb=1.,sfh=1.;
   if(scaling==kConvertMKStoCGS){
      sfb = 1./util_df::Constants::Gauss;
      sfh = 1./util_df::Constants::Oerstead;
   }else if(scaling==kConvertCGStoMKS){
      sfb = util_df::Constants::Gauss;
      sfh = util_df::Constants::Oerstead;
   }
   const int N = b.size();
   for(int i=0;i<N;i++){
      b[i] *= sfb;
      h[i] *= sfh;
   }

   SetData(b,h);
 
   return 0;
}
//______________________________________________________________________________
void NonLinearBHCurve::WriteFile(const char *outpath,const char *delim){
   // write data to a tsv file 
   const int NROW = fB.size();
   const int NCOL = 2; 
   util_df::CSVManager *data = new util_df::CSVManager(delim);
   data->InitTable(NROW,NCOL); 
   data->SetColumn<double>(0,fB);  
   data->SetColumn<double>(1,fH);
   data->WriteFile(outpath); 
   delete data; 
}
