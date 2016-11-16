/*
Class which allows for generic reading of data and construction 
of histograms for rates of each event type. This can be any type of data coming
from UCNA root files. What is stored should be histograms of event rates binned
by energy along with the run time for that particular run.

Also included is a class to handle reading in reverse calibrated 
simulation data
*/

#include "EvtRateHandler_ref.hh"
#include "MBUtils.hh"
#include <cstdio>
#include <fstream>

#include <TString.h>


int separate23(int side, double mwpcEn) {
  int type = 2;
  if (side==0)  type = ( mwpcEn>4.14 ) ? 3 : 2;  
  if (side==1)  type = ( mwpcEn>4.14 ) ? 3 : 2;
  return type;
};


EvtRateHandler::EvtRateHandler(std::vector<int> rn, bool fg, std::string anaCh, double enBinWidth, double fidCut, bool ukdata, bool unblind) : runs(rn),FG(fg),analysisChoice(anaCh),fiducialCut(fidCut),UKdata(ukdata),unblinded(unblind),pol(0) {

  numEnergyBins = (int)(1200./enBinWidth);

  pol = polarization(runs[0]); //THIS NEEDS TO BE REWRITTEN TO READ IN MASTER FILE
  
  rateEvec.resize(numEnergyBins,0.); 
  rateWvec.resize(numEnergyBins,0.);
  rateEerr.resize(numEnergyBins,0.); 
  rateWerr.resize(numEnergyBins,0.);

  refSpectraE.resize(numEnergyBins,0.);
  refSpectraW.resize(numEnergyBins,0.);
  totalRefTime = 0.;
  totalRefCountsE = 0.;
  totalRefCountsW = 0.;

  loadReferenceSpectra(); // Loading the proper reference spectra
  
  runLength.resize(runs.size(), std::vector<double> (2,0.));
  totalRunLengthE = totalRunLengthW = 0.;
  UCNMonIntegral.resize(runs.size(), 0.);
  
  //Loading information from log file

  for (unsigned int i = 0; i<runs.size(); i++) {
    std::string logFilePath = std::string(getenv("RUN_INFO_FILES"))+"runInfo_"+itos(runs[i])+".dat";
    std::ifstream logFile(logFilePath.c_str());
    std::vector <std::string> title(4);
    std::vector <double> value(4);
    for (int i=0; i<4; i++) {
      logFile >> title[i] >> value[i];
    }
    logFile.close();
    
    if (unblinded) {
      runLength[i][0] = value[2];
      runLength[i][1] = value[2];
    }
    else {
      runLength[i][0] = value[0];
      runLength[i][1] = value[1];
    }
    UCNMonIntegral[i] = value[3];

    totalRunLengthE += runLength[i][0];
    totalRunLengthW += runLength[i][1];
    
    std::cout << title[0] << "\t" << runLength[i][0] << std::endl;
    std::cout << title[1] << "\t" << runLength[i][1] << std::endl;
    std::cout << title[3] << "\t" << UCNMonIntegral[i] << std::endl;
    //std::cout << holdTitle << "\t" << runLength[0] << std::endl;
  }
};

EvtRateHandler::~EvtRateHandler() {
  
    if (hisCounts[0]) delete hisCounts[0];
    if (hisCounts[1]) delete hisCounts[1]; 

};

int EvtRateHandler::polarization(int run) {

  std::ifstream infile(TString::Format("%s/masterBetaRunList.txt",getenv("OCTET_LIST")).Data());

  int rn = 0;
  std::string type = "";
  bool runFound = false;

  while ( infile >> rn >> type ) {
    if ( rn == run ) { runFound = true; infile.close(); continue; }
  }

  if ( !runFound ) throw "You didn't pick a run in an octet...";
  //Flipper OFF is p = -1
  if ( type=="A1" || type=="A2" || type=="A10" || type=="A12" || type=="B4" || type=="B5" || type=="B7" || type=="B9" ) return -1;
  else if ( type=="B1" || type=="B2" || type=="B10" || type=="B12" || type=="A4" || type=="A5" || type=="A7" || type=="A9" ) return 1;
  
  else throw "You chose a Depol Run.. no polarization";
  
};

void EvtRateHandler::loadReferenceSpectra() {

  // NEED TO REWRITE THE REFERENCE THINGIES TO ALSO SPIT OUT TOTAL EVENTS RATHER THAN RATES
  
  TString dir = ( FG ? TString::Format("%s/reference_rates/foregroundRatesByAnaChoice/",getenv("ANALYSIS_CODE")) :
		  TString::Format("%s/reference_rates/backgroundRatesByAnaChoice/",getenv("ANALYSIS_CODE")) );
    
  TString filename = ( dir + TString("ReferenceSpectra_Octets-") +
		       ( runs[0]<20000 ? TString("0-59") : TString("60-121") ) +
		       ( pol==-1 ? TString("_sfOFF-") : TString("_sfON-") ) +
		       TString::Format("AnaCh-%s.txt",analysisChoice.c_str()) );

  std::ifstream infile(filename.Data());

  std::string label;
  infile >> label >> totalRefTime;

  double binMid=0., eastRef=0., eastErr=0., westRef=0., westErr=0.;
  int inc = 0;

  while ( infile >> binMid >> eastRef >> eastErr >> westRef >> westErr ) {

    refSpectraE[inc] = eastRef;
    refSpectraW[inc] = westRef;

    totalRefCountsE += eastRef;
    totalRefCountsW += westRef;
    
    ++inc;

  }

  std::cout << "Finished Loading Reference Rate data...\n\n";
  
};

double EvtRateHandler::referenceError(int side, int bin) {

  if ( side == 0 ) return totalRefCountsE>0. ? sqrt( refSpectraE[bin] * ( totalCountsE / totalRefCountsE ) ) : 0.;
  else if ( side == 1 ) return totalRefCountsW>0. ? sqrt( refSpectraW[bin] * ( totalCountsW / totalRefCountsW ) ) : 0.;
  else throw "BAD SIDE GIVEN TO referenceError(int side, int bin)";
};

void EvtRateHandler::CalcRates() {
  
  this->dataReader();

  Double_t binContentE = 0., binContentW = 0.;

  //NEED THE TOTAL NUMBER OF EVENTS IN EACH SIDES HISTOGRAM FOR SCALING
  for (unsigned int i=0; i<numEnergyBins; i++) {

    binContentE = hisCounts[0]->GetBinContent(i+1);
    binContentW = hisCounts[1]->GetBinContent(i+1);
    
    rateEvec[i] = binContentE / totalRunLengthE;
    rateWvec[i] = binContentW / totalRunLengthW;

    rateEerr[i] = binContentE > 25. ? sqrt(binContentE) : referenceError(0,i);
    rateWerr[i] = binContentW > 25. ? sqrt(binContentW) : referenceError(1,i);

    rateEerr[i] /= totalRunLengthE;
    rateWerr[i] /= totalRunLengthW;
	
  }
      
};


void EvtRateHandler::dataReader() {

  hisCounts[0] = new TH1D("East Event Hist","EastCounts",numEnergyBins,0.,1200.);
  hisCounts[1] = new TH1D("West Event Hist","WestCounts",numEnergyBins,0.,1200.);

  bool sep23 = false;

  // 2/3 separation
  if ( analysisChoice==std::string("C") ||
       analysisChoice==std::string("H") ||
       analysisChoice==std::string("J") ||
       analysisChoice==std::string("K") )
    sep23 = true;

  //Event types
  bool Type0=false, Type1=false, Type2=false, Type3=false;
  if ( analysisChoice=="A" || analysisChoice=="B" || analysisChoice=="C" || analysisChoice=="D" ) Type0 = true;
  if ( analysisChoice=="A" || analysisChoice=="B" || analysisChoice=="C" || analysisChoice=="D" ) Type1 = true;
  if ( analysisChoice=="A" || analysisChoice=="C" || analysisChoice=="G" || analysisChoice=="H" || analysisChoice=="J" ) Type2 = true;
  if ( analysisChoice=="A" || analysisChoice=="C" || analysisChoice=="G" || analysisChoice=="H" || analysisChoice=="K") Type3 = true;


  
  char temp[100];
  
  std::string infile;
  TFile *input;
  TTree *Tin;

  double EmwpcX=0., EmwpcY=0., WmwpcX=0., WmwpcY=0., TimeE=0., TimeW=0., Erecon=0., MWPCEnergyE=0., MWPCEnergyW=0.; //Branch Variables being read in
  float EmwpcX_f=0., EmwpcY_f=0., WmwpcX_f=0., WmwpcY_f=0., TimeE_f=0., TimeW_f=0., Erecon_f=0., MWPCEnergyE_f=0., MWPCEnergyW_f=0.; // For reading in data from MPM replays

  int PID, Side, Type;

  for ( unsigned int i=0; i<runs.size(); i++ ) {
  
    //Set branch addresses
    if (UKdata) {
      sprintf(temp,"%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),runs[i]);
      std::string infile = std::string(temp);
      input = new TFile(infile.c_str(), "READ");
      Tin = (TTree*)input->Get("pass3");
      
      Tin->SetBranchAddress("PID", &PID);
      Tin->SetBranchAddress("Type", &Type);
      Tin->SetBranchAddress("Side", &Side); 
      Tin->SetBranchAddress("Erecon",&Erecon);
      Tin->SetBranchAddress("TimeE",&TimeE);
      Tin->SetBranchAddress("TimeW",&TimeW);
      Tin->GetBranch("xE")->GetLeaf("center")->SetAddress(&EmwpcX);
      Tin->GetBranch("yE")->GetLeaf("center")->SetAddress(&EmwpcY);
      Tin->GetBranch("xW")->GetLeaf("center")->SetAddress(&WmwpcX);
      Tin->GetBranch("yW")->GetLeaf("center")->SetAddress(&WmwpcY);
      
      //NEED TO ADD IN WIRECHAMBER ENERGY FOR 2/3 SEPARATION
    }
    else {
      sprintf(temp,"%s/spec_%i.root",getenv("UK_SPEC_REPLAY"),runs[i]);
      std::string infile = std::string(temp);
      input = new TFile(infile.c_str(), "READ");
      Tin = (TTree*)input->Get("phys");
      
      Tin->SetBranchAddress("PID", &PID);
      Tin->SetBranchAddress("Type", &Type);
      Tin->SetBranchAddress("Side", &Side); 
      Tin->SetBranchAddress("Erecon",&Erecon_f);
      Tin->SetBranchAddress("TimeE",&TimeE_f);
      Tin->SetBranchAddress("TimeW",&TimeW_f);
      Tin->GetBranch("xEmpm")->GetLeaf("center")->SetAddress(&EmwpcX_f);
      Tin->GetBranch("yEmpm")->GetLeaf("center")->SetAddress(&EmwpcY_f);
      Tin->GetBranch("xWmpm")->GetLeaf("center")->SetAddress(&WmwpcX_f);
      Tin->GetBranch("yWmpm")->GetLeaf("center")->SetAddress(&WmwpcY_f);
    }
    unsigned int nevents = Tin->GetEntriesFast();
    std::cout << "Number of Events: " << nevents << std::endl;
    
    
    //Run over all events in file, fill output histogram with rate
    
    double r2E = 0.; //position of event squared
    double r2W = 0.; //position of event squared
    
    for (unsigned int i=0; i<nevents; i++) {
      
      Tin->GetEvent(i);
      
      //Cast float variables to double if they are from MPM replay
      if (!UKdata) {
	EmwpcX = (double) EmwpcX_f;
	EmwpcY = (double) EmwpcY_f;
	WmwpcX = (double) WmwpcX_f;
	WmwpcY = (double) WmwpcY_f;
	Erecon = (double) Erecon_f;
	TimeE = (double) TimeE_f;
	TimeW = (double) TimeW_f;
      }
      
      if (PID==1) {
	
	r2E=EmwpcX*EmwpcX+EmwpcY*EmwpcY;
	r2W=WmwpcX*WmwpcX+WmwpcY*WmwpcY;
	
	if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut ) ) {
	  
	  
	  if ( sep23 ) {
	    if (Erecon>0. && Type==2) {
	      
	      if (Side==0) {
		Type = separate23(Side,MWPCEnergyE);
		Side = Type==2 ? 1 : 0;
	      }
	      else if (Side==1) {
		Type = separate23(Side,MWPCEnergyW);
		Side = Type==2 ? 0 : 1;
	      }
	    }
	  }
	  
	  //Type0
	  if ( Type0 && Type==0 ) hisCounts[Side]->Fill(Erecon);
	  //Type1
	  if ( Type1 && Type==1 ) hisCounts[Side]->Fill(Erecon);
	  //Type0
	  if ( Type2 && Type==2 ) hisCounts[Side]->Fill(Erecon);
	  //Type0
	  if ( Type3 && Type==3 ) hisCounts[Side]->Fill(Erecon); 
	  
	}
      }
    }
    input->Close();
  }

  totalCountsE = (double) hisCounts[0]->Integral();
  totalCountsW = (double) hisCounts[1]->Integral();
    
  std::cout << "Beta Events: " << totalCountsE+totalCountsW << std::endl;
  
  if (input) delete input;
  
};

void SimEvtRateHandler::dataReader() {

  hisCounts[0] = new TH1D("East Event Hist","EastCounts",numEnergyBins,0.,1200.);
  hisCounts[1] = new TH1D("West Event Hist","WestCounts",numEnergyBins,0.,1200.);

  bool sep23 = false;

  // 2/3 separation
  if ( analysisChoice==std::string("C") ||
       analysisChoice==std::string("H") ||
       analysisChoice==std::string("J") ||
       analysisChoice==std::string("K") )
    sep23 = true;

  //Event types
  bool Type0=false, Type1=false, Type2=false, Type3=false;
  if ( analysisChoice=="A" || analysisChoice=="B" || analysisChoice=="C" || analysisChoice=="D" ) Type0 = true;
  if ( analysisChoice=="A" || analysisChoice=="B" || analysisChoice=="C" || analysisChoice=="D" ) Type1 = true;
  if ( analysisChoice=="A" || analysisChoice=="C" || analysisChoice=="G" || analysisChoice=="H" || analysisChoice=="J" ) Type2 = true;
  if ( analysisChoice=="A" || analysisChoice=="C" || analysisChoice=="G" || analysisChoice=="H" || analysisChoice=="K") Type3 = true;
  
  char temp[100];
  std::string infile;
  TFile *input;
  TTree *Tin;

  double TimeE=0., TimeW=0., Erecon=0., MWPCEnergyE=0., MWPCEnergyW=0.; //Branch Variables being read in
  int PID, Side, Type;

  double AsymWeight=1.;
  double mwpcPosE[3]={0.};
  double mwpcPosW[3]={0.}; //holds the position of the event in the MWPC for simulated data

  for ( unsigned int i=0; i<runs.size(); i++ ) {

    sprintf(temp,"%s/beta/revCalSim_%i_Beta.root",getenv("REVCALSIM"),runs[i]);
    infile = std::string(temp);
    input = new TFile(infile.c_str(), "READ");
    Tin = (TTree*)input->Get("revCalSim");
    //std::cout << "made it here in dataReader\n";
    //Set branch addresses
    
    Tin->SetBranchAddress("PID", &PID);
    Tin->SetBranchAddress("type", &Type);
    Tin->SetBranchAddress("side", &Side); 
    Tin->SetBranchAddress("Erecon",&Erecon);
    Tin->GetBranch("time")->GetLeaf("timeE")->SetAddress(&TimeE);
    Tin->GetBranch("time")->GetLeaf("timeW")->SetAddress(&TimeW);
    Tin->GetBranch("MWPCPos")->GetLeaf("MWPCPosE")->SetAddress(mwpcPosE);
    Tin->GetBranch("MWPCPos")->GetLeaf("MWPCPosW")->SetAddress(mwpcPosW);
    Tin->SetBranchAddress("AsymWeight",&AsymWeight);
    // ADD IN WIRECHAMBER ENERGY DEPOSITION
    
    unsigned int nevents = Tin->GetEntriesFast();
    //std::cout << nevents << std::endl;
    
    //Run over all events in file, fill output histogram with rate
    
    double r2E = 0.; //position of event squared
    double r2W = 0.;
    
    for (unsigned int i=0; i<nevents; i++) {
      
      Tin->GetEvent(i);
      
      if (PID==1) {
        
	r2W=mwpcPosW[0]*mwpcPosW[0]+mwpcPosW[1]*mwpcPosW[1];
	r2E=mwpcPosE[0]*mwpcPosE[0]+mwpcPosE[1]*mwpcPosE[1];
	
	if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut ) ) {
	  
	  
	  if ( sep23 ) {
	    if (Erecon>0. && Type==2) {
	      
	      if (Side==0) {
		Type = separate23(Side,MWPCEnergyE);
		Side = Type==2 ? 1 : 0;
	      }
	      else if (Side==1) {
		Type = separate23(Side,MWPCEnergyW);
		Side = Type==2 ? 0 : 1;
	      }
	    }
	  }
	  
	  //Type0
	  if ( Type0 && Type==0 ) hisCounts[Side]->Fill(Erecon);
	  //Type1
	  if ( Type1 && Type==1 ) hisCounts[Side]->Fill(Erecon);
	  //Type2
	  if ( Type2 && Type==2 ) hisCounts[Side]->Fill(Erecon);
	  //Type3
	  if ( Type3 && Type==3 ) hisCounts[Side]->Fill(Erecon); 
	  
	}
      }
    }
    input->Close();
  }
  totalCountsE = (double) hisCounts[0]->Integral();
  totalCountsW = (double) hisCounts[1]->Integral();
    
  std::cout << "Beta Events: " << totalCountsE+totalCountsW << std::endl;

  if (input) delete input;
  
};


BGSubtractedRate::BGSubtractedRate(std::vector <int> fgruns, std::vector <int> bgruns, std::string anaCh,
				   double enBin, double fidCut, bool ukdata, bool sim, bool unblind) : FGruns(fgruns), BGruns(bgruns), analysisChoice(anaCh), EnergyBinWidth(enBin), fiducialCut(fidCut), UKdata(ukdata), Simulation(sim), UNBLIND(unblind) {
  
  int numBins = int(1200./double(EnergyBinWidth));

  BetaRateE.resize(numBins,0.);
  BGRateE.resize(numBins,0.);
  FinalRateE.resize(numBins,0.);
  BetaRateW.resize(numBins,0.);
  BGRateW.resize(numBins,0.);
  FinalRateW.resize(numBins,0.);
  BetaRateErrorE.resize(numBins,0.);
  BGRateErrorE.resize(numBins,0.);
  FinalRateErrorE.resize(numBins,0.);
  BetaRateErrorW.resize(numBins,0.);
  BGRateErrorW.resize(numBins,0.);
  FinalRateErrorW.resize(numBins,0.);

  runLengthBeta.resize(2,0.);
  runLengthBG.resize(2,0.);
};


void BGSubtractedRate::calcBGSubtRates() {

  if (Simulation) 
    {
      SimEvtRateHandler *evt = new SimEvtRateHandler(FGruns, true, analysisChoice, EnergyBinWidth, fiducialCut, UNBLIND);
      evt->CalcRates();
      FinalRateE = evt->getRateVectors(0);
      FinalRateW = evt->getRateVectors(1);
      FinalRateErrorE = evt->getRateErrors(0);
      FinalRateErrorW = evt->getRateErrors(1);
      delete evt;
      //for (unsigned int t=0;t<FinalRateE.size();t++) {
      //for (unsigned int i=0;i<FinalRateE[0].size();i++) {
      //  FinalRateErrorE[t][i] = sqrt(FinalRateE[t][i]);
      //  FinalRateErrorW[t][i] = sqrt(FinalRateW[t][i]);
      //}
      //}
    }
  else 
    {
      LoadRatesByBin();
      std::cout << "Loaded Rates for run " << FGruns[0]  << std::endl;
      CalcFinalRate();
    }
};

std::vector<double > BGSubtractedRate::returnBGSubtRate(int side) {  

  if (side==0) return FinalRateE;
  else if (side==1) return FinalRateW;
  else throw "Invalid side when returning BG subtracted rate";
  
};

std::vector<double > BGSubtractedRate::returnBGSubtRateError(int side) {
  
  if (side==0) return FinalRateErrorE;
  else if (side==1) return FinalRateErrorW;
  else throw "Invalid side when returning BG subtracted rate error";
  
};


void BGSubtractedRate::LoadRatesByBin() {
  
  EvtRateHandler *evtBG = new EvtRateHandler(BGruns, false, analysisChoice, EnergyBinWidth, fiducialCut, UKdata, UNBLIND);
  evtBG->CalcRates();
  BGRateE = evtBG->getRateVectors(0);
  BGRateErrorE = evtBG->getRateErrors(0);
  BGRateW = evtBG->getRateVectors(1);
  BGRateErrorW = evtBG->getRateErrors(1);
  runLengthBG[0] = evtBG->returnRunLength(0);
  runLengthBG[1] = evtBG->returnRunLength(1);
  delete evtBG;
    
  EvtRateHandler *evt = new EvtRateHandler(FGruns, true, analysisChoice, EnergyBinWidth, fiducialCut, UKdata, UNBLIND);
  evt->CalcRates();
  BetaRateE = evt->getRateVectors(0);
  BetaRateErrorE = evt->getRateErrors(0);
  BetaRateW = evt->getRateVectors(1);
  BetaRateErrorW = evt->getRateErrors(1);
  runLengthBeta[0] = evt->returnRunLength(0);
  runLengthBeta[1] = evt->returnRunLength(1);
  delete evt;
  
};

std::vector<double> BGSubtractedRate::returnRunLengths(bool beta) {
  if (beta) return runLengthBeta;
  else return runLengthBG;
}

void BGSubtractedRate::CalcFinalRate()  {
  std::cout << "Calculating Final Rates for run " << FGruns[0] << std::endl;
  
  if (BetaRateE.size()==BGRateE.size()) {
    FinalRateE.resize(BetaRateE.size(),0.);
    for (unsigned int i=0; i<BetaRateE.size(); i++)
      {
	FinalRateE[i] = BetaRateE[i]-BGRateE[i];
	FinalRateErrorE[i] = sqrt(power(BetaRateErrorE[i],2)+power(BGRateErrorE[i],2));
	//std::cout << BetaRateE[i] << " " << BGRateE[i] << " " << FinalRateE[i] << std::endl;
      }
  }
  else throw "Number of energy bins do not agree between Beta and BG runs. Can't calculate final rate!";
  
  if (BetaRateW.size()==BGRateW.size()) {
    FinalRateW.resize(BetaRateW.size(),0.);
    for (unsigned int i=0; i<BetaRateW.size(); i++)
      {
	FinalRateW[i] = BetaRateW[i]-BGRateW[i];
	FinalRateErrorW[i] = sqrt(power(BetaRateErrorW[i],2)+power(BGRateErrorW[i],2));
	//std::cout << BetaRate[i] << " " << BGRate[i] << " " << FinalRate[i] << std::endl;
      }
  }
  else throw "Number of energy bins do not agree between Beta and BG runs. Can't calculate final rate!";

};


