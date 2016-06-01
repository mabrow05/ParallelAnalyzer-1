#ifndef CALTOOLS_HH
#define CALTOOLS_HH

#include <TF1.h>
#include <vector>
#include <fstream>
#include <cstdlib>


class LinearityCurve {
 
public:

  LinearityCurve(Int_t srcPeriod);  
  ~LinearityCurve(); 
  void readLinearityCurveParams(Int_t srcPeriod);      
  Double_t applyLinCurve(Int_t pmt, Double_t x); 
  Double_t applyInverseLinCurve(Int_t pmt, Double_t x); 
  Int_t getCurrentSrcCalPeriod() { return sourceCalPeriod; }; 

private:                                                                                                                      

  Int_t sourceCalPeriod;  
  std::vector < std::vector < Double_t > > pmtParams;                                                                                        
  TF1 *linCurve;                                                                    
};        

#endif
