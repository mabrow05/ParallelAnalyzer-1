#include "calibrationTools.hh"




LinearityCurve::LinearityCurve(Int_t period) : sourceCalPeriod(period) {
  
  pmtParams.resize(8, std::vector <Double_t> (6,0.));
  linCurve = new TF1("linCurve","([0] + [1]*x + [2]*x*x)*(0.5+0.5*TMath::TanH((x-[4])/[5]))+([3]*x)*(0.5-0.5*TMath::TanH((x-[4])/[5]))", -100., 4096.0);
  readLinearityCurveParams(sourceCalPeriod);

}

LinearityCurve::~LinearityCurve() {
  delete linCurve;
}


void LinearityCurve::readLinearityCurveParams(Int_t srcPeriod) {
  // Read linearity curve                                                                                                                      
  char tempFileLinearityCurve[500];

  sprintf(tempFileLinearityCurve, "%s/lin_curves_srcCal_Period_%i.dat",getenv("LINEARITY_CURVES"),srcPeriod);
  std::cout << "Reading in linearity curve from " << tempFileLinearityCurve << ":\n";
  std::cout << "p0\tp1\tp2\tp3\tp4\tp5\n";
  
  ifstream fileLinearityCurve(tempFileLinearityCurve);
  Double_t p[6];//p0,p1,p2;                                                                                                                    
  Int_t i=0;
  while (fileLinearityCurve >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5]) {
    Int_t ii=0;
    while(ii<6) {pmtParams[i][ii] = p[ii]; ii++;}
    i++;
    std::cout << p[0] << " " << p[1] << " " << p[2]  << " " << p[3] << " " << p[4] << " " << p[5] << std::endl;
    if (fileLinearityCurve.fail()) break;
  }

  sourceCalPeriod = srcPeriod;
}

Double_t LinearityCurve::applyLinCurve(Int_t pmt, Double_t x) {
  return linCurve->EvalPar(&x,&pmtParams[pmt][0]);
}

Double_t LinearityCurve::applyInverseLinCurve(Int_t pmt, Double_t y) {
  linCurve->SetParameters(pmtParams[pmt][0], pmtParams[pmt][1], pmtParams[pmt][2], pmtParams[pmt][3], pmtParams[pmt][4], pmtParams[pmt][5]);
  return linCurve->GetX(y, -50., 4096.);
}
