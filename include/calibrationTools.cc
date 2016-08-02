#include "calibrationTools.hh"




LinearityCurve::LinearityCurve(Int_t period, bool useTanhSmear) : sourceCalPeriod(period), useTanh(useTanhSmear) {
  
  if (useTanh) {
    nParams = 6;
    linCurve = new TF1("linCurve","([0] + [1]*x + [2]*x*x)*(0.5+0.5*TMath::TanH((x-[4])/[5]))+([3]*x)*(0.5-0.5*TMath::TanH((x-[4])/[5]))", -100., 4096.0);
  }
  else {
    nParams = 3;
    linCurve = new TF1("linCurve","([0] + [1]*x + [2]*x*x)", -100., 4096.0);
  }
  
  pmtParams.resize(8, std::vector <Double_t> (nParams,0.));
  
  
  //
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

  if (useTanh) std::cout << "p0\tp1\tp2\tp3\tp4\tp5\n";
  else std::cout << "p0\tp1\tp2\n";

  ifstream fileLinearityCurve(tempFileLinearityCurve);
  //std::vector <Double_t> p(nParams);                                                           
  Double_t p;

  Int_t i = 0;
  Int_t ii=0;
  while (fileLinearityCurve >> p) {
    pmtParams[i][ii] = p;
    std::cout << p << " ";
    ii++;
    if (ii==nParams) { ii = 0; i++; std::cout << "\n"; }

    if (fileLinearityCurve.fail()) break;
  }
    
                                                       
  /*Int_t i=0;
  while (fileLinearityCurve >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5]) {
    Int_t ii=0;
    while(ii<6) {pmtParams[i][ii] = p[ii]; ii++;}
    i++;
    std::cout << p[0] << " " << p[1] << " " << p[2]  << " " << p[3] << " " << p[4] << " " << p[5] << std::endl;
    if (fileLinearityCurve.fail()) break;
  }
  */
  sourceCalPeriod = srcPeriod;
}

Double_t LinearityCurve::applyLinCurve(Int_t pmt, Double_t x) {
  return linCurve->EvalPar(&x,&pmtParams[pmt][0]);
}

Double_t LinearityCurve::applyInverseLinCurve(Int_t pmt, Double_t y) {
  if (useTanh) linCurve->SetParameters(pmtParams[pmt][0], pmtParams[pmt][1], pmtParams[pmt][2], pmtParams[pmt][3], pmtParams[pmt][4], pmtParams[pmt][5]);
  else linCurve->SetParameters(pmtParams[pmt][0], pmtParams[pmt][1], pmtParams[pmt][2]);
  return linCurve->GetX(y, -50., 4096.);
}


///////////////////////////////////////////////////////////////////////

TriggerFunctions::TriggerFunctions(Int_t run) : currentRun(run) {

  triggerFuncs.resize(8,std::vector<Double_t>(8,0.));
  readTriggerFunctions(currentRun);
  
  func = new TF1("func","([0]+[1]*TMath::Erf((x-[2])/[3]))*(0.5-.5*TMath::TanH((x-[2])/[4]))+(0.5+.5*TMath::TanH((x-[2])/[4]))*([5]+[6]*TMath::TanH((x-[2])/[7]))",-80.,4096.);

}

TriggerFunctions::~TriggerFunctions() {
  
  delete func;

}


void TriggerFunctions::readTriggerFunctions(Int_t run) {

  ifstream file(TString::Format("%s/runs/%i_thresholds.dat",getenv("TRIGGER_FUNC"),run).Data());

  std::vector <Double_t> params(8,0.);

  if (file.is_open()) {
    Int_t pmt = 0;
    while (file >> params[0] >> params[1] >> params[2] >> params[3] >> params[4] >>
	   params[5] >> params[6] >> params[7]) {
      
      for (Int_t param=0; param<8; param++) {
	triggerFuncs[pmt][param] = params[param];
	std::cout << triggerFuncs[pmt][param] << " ";
      }
      std::cout << std::endl;
      pmt++;
    }
  }
  else { 
    std::cout << "Couldn't open trigger function file for run " << run << std::endl;
    exit(0);
  }
  file.close();

  //Making sure the trigger probability is zero for WPMT2 when in the bad run range...
  if (run<=17249 || run>=16983) {
    triggerFuncs[5][0] = triggerFuncs[5][1] = triggerFuncs[5][5] = triggerFuncs[5][6] = 0.;
  }

}

bool TriggerFunctions::decideEastTrigger(std::vector<Double_t> ADC, std::vector<Double_t> rands) {

  Int_t numTriggs = 0;
  Double_t adc = 0.;
  
  adc = ADC[0];
  numTriggs = rands[0] < func->EvalPar(&adc,&triggerFuncs[0][0]) ? numTriggs+1 : numTriggs;
  adc = ADC[1];
  numTriggs = rands[1] < func->EvalPar(&adc,&triggerFuncs[1][0]) ? numTriggs+1 : numTriggs;
  
  if (numTriggs<2) {
    adc = ADC[2];
    numTriggs = rands[2] < func->EvalPar(&adc,&triggerFuncs[2][0]) ? numTriggs+1 : numTriggs;
  }
  else return true;
  
  if (numTriggs<2) {
    adc = ADC[3];
    numTriggs = rands[3] < func->EvalPar(&adc,&triggerFuncs[3][0]) ? numTriggs+1 : numTriggs;
  }
  else return true;

  return numTriggs>=2 ? true : false;

}

bool TriggerFunctions::decideWestTrigger(std::vector<Double_t> ADC, std::vector<Double_t> rands) {

  Int_t numTriggs = 0;
  Double_t adc = 0.;
  
  adc = ADC[0];
  numTriggs = rands[0] < func->EvalPar(&adc,&triggerFuncs[4][0]) ? numTriggs+1 : numTriggs;

  //std::cout << rands[0] << " " << func->EvalPar(&adc,&triggerFuncs[4][0]) << std::endl;

  adc = ADC[1];
  numTriggs = rands[1] < func->EvalPar(&adc,&triggerFuncs[5][0]) ? numTriggs+1 : numTriggs;
  
  if (numTriggs<2) {
    adc = ADC[2];
    numTriggs = rands[2] < func->EvalPar(&adc,&triggerFuncs[6][0]) ? numTriggs+1 : numTriggs;
  }
  else return true;

  if (numTriggs<2) {
    adc = ADC[3];
    numTriggs = rands[3] < func->EvalPar(&adc,&triggerFuncs[7][0]) ? numTriggs+1 : numTriggs;
  }
  else return true;

  return numTriggs>=2 ? true : false;

}
