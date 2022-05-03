#ifndef NON_LINEAR_BH_CURVE_HH
#define NON_LINEAR_BH_CURVE_HH

// a class for drawing B vs H curves for non-linear materials

#include "UtilDFUnits.hh" 
#include "UtilDFConstants.hh"
#include "CSVManager.hh" 

enum sfChoice { kNoConversion=0, kConvertMKStoCGS=1, kConvertCGStoMKS=2 }; 

class NonLinearBHCurve {

   private:
      std::string fName;           // name of material 
      int fNPTS;                   // number of points in graph
      int fUnits;                  // units: CGS or MKS 
      double fBMin,fBMax;          // min and max B    
      double fk1,fk2,fk3;          // k parameters for B vs H curve  
      std::vector<double> fB,fH;

      double GetMu(double B); 

   public:
      NonLinearBHCurve(const char *name="NONE",double k1=0,double k2=0,double k3=0,int units=util_df::Units::kCGS);
      ~NonLinearBHCurve(); 

      void Print();  
      void ClearData();
      void WriteFile(const char *outpath,const char *delim="tsv"); 
      void SetParameters(double k1,double k2,double k3,int units)   { fk1 = k1; fk2 = k2; fk3 = k3; fUnits = units; }

      void SetData(std::vector<double> b,std::vector<double> h); 
 
      void Calculate(int N,double bMin,double bMax);                // compute B vs H for N points in a range B = (bMin,bMax)
      void GetData(std::vector<double> &b,std::vector<double> &h);    
      
      int LoadData(const char *inpath,const char *delim="csv",bool header=false,int scaling=kNoConversion); 
      
};

#endif 
