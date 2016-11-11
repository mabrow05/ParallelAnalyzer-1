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


EvtRateHandler::EvtRateHandler(std::vector<int> rn, bool fg, std::string anaCh, double enBinWidth=10., double fidCut=100., bool ukdata, bool unblind) : runs(rn),FG(fg),analysisChoice(anaCh),fiducialCut(fidCut),UKdata(ukdata), pol(0), unblinded(ub) {

  numEnergyBins = (int)(1200./enBinWidth);

  pol = polarization(runs[0]); //THIS NEEDS TO BE REWRITTEN TO READ IN MASTER FILE
  
  rateEvec.resize(numEnergyBins,0.); 
  rateWvec.resize(numEnergyBins,0.);
  rateEerr.resize(numEnergyBins,0.); 
  rateWerr.resize(numEnergyBins,0.);

  refRate.resize(numEnergyBins,0.);
  
  runLength.resize(runs.size(), std::vector<double> (2,0.));
  UCNMonIntegral.resize(runs.size(), 0.);
  
  //Loading information from log file

  for (unsigned int i = 0; i<runs.size(); i++) {
    std::string logFilePath = std::string(getenv("RUN_INFO_FILES"))+"runInfo_"+itos(runNumber)+".dat";
    ifstream logFile(logFilePath.c_str());
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
  // NEED TO CHANGE THIS TO NOT USE DB
  std::string dbAddress = std::string(getenv("UCNADBADDRESS"));
  std::string dbname = std::string(getenv("UCNADB"));
  std::string dbUser = std::string(getenv("UCNADBUSER"));
  std::string dbPass = std::string(getenv("UCNADBPASS"));
  
  char cmd[200];
  sprintf(cmd,"SELECT flipper FROM run WHERE run_number=%i;",run);

  SQLdatabase *db = new SQLdatabase(dbname, dbAddress, dbUser, dbPass);
  db->fetchQuery(cmd);
  std::string flipperStatus = db->returnQueryEntry();
  std::cout << flipperStatus << std::endl;
  delete db;
  if (flipperStatus=="On") return 1;
  else if (flipperStatus=="Off") return -1;
  else throw "Polarization isn't applicaple or you chose a Depol Run";
  
};

void EvtRateHandler::loadReferenceRate() {

  // NEED TO REWRITE THE REFERENCE THINGIES TO ALSO SPIT OUT TOTAL EVENTS RATHER THAN RATES
  
  TString dir = ( FG ? TString::Format("%s/reference_rates/foregroundRatesByAnaChoice/",getenv("ANALYSIS_CODE")) :
		  TString::Format("%s/reference_rates/backgroundRatesByAnaChoice/",getenv("ANALYSIS_CODE")) );
    
  TString filename = ( dir + TString("ReferenceRates_Octets-") +
		       ( runs[0]<20000 ? TString("0-59") : TString("60-121") ) +
		       ( pol==-1 ? TString("_sfOFF-") : TString("_sfON-") ) +
		       TString::Format("AnaCh-%s.txt",analysisChoice.c_str()) );

  std::ifstream infile(filename.Data());

  double binMid=0., eastRef=0., eastErr=0., westRef=0., westErr=0.;
  int inc = 0;

  while ( infile >> binMid >> eastRef >> eastErr >> westRef >> westErr ) {

    refRateE[i] = eastRef;
    refRateW[i] = westRef;

  }

  std::cout << "Finished Loading Reference Rate data...\n\n";
  
};

double EvtRateHandler::referenceError(int side, int bin) {

  if ( side == 0 ) return refRateE[i]
};

void EvtRateHandler::CalcRates() {
  
  this->dataReader();

  //NEED THE TOTAL NUMBER OF EVENTS IN EACH SIDES HISTOGRAM FOR SCALING
  
  // Calculate rates from data histograms and run length
  for (unsigned int j=0; j<3; j++) { //Right now we leave type 3 empty since no separation occurs
    for (unsigned int i=0; i<numEnergyBins; i++) {
      
      rateEvec[j][i] = hisE[j]->GetBinContent(i+1) / runLength[0];
      rateWvec[j][i] = hisW[j]->GetBinContent(i+1) / runLength[1];
      //// How to handle low statistics?

      if ( setErrorFloor ) {
	rateEerr[j][i] = hisE[j]->GetBinError(i+1) > 1. ? hisE[j]->GetBinError(i+1) / runLength[0] : 1./runLength[0];
	rateWerr[j][i] = hisW[j]->GetBinError(i+1) > 1. ? hisW[j]->GetBinError(i+1) / runLength[1] : 1./runLength[1];
      }
      else {
	rateEerr[j][i] = hisE[j]->GetBinError(i+1) / runLength[0]; 
	rateWerr[j][i] = hisW[j]->GetBinError(i+1) / runLength[1]; 
      }
      
    }
  }
      
};



TH1D EvtRateHandler::getCountHist(int side, int evtType) {
  if (evtType==0 || evtType==1 || evtType==2 || evtType==3) {
    if (side==0) return *hisE[evtType];
    else if (side==1) return *hisW[evtType];
    else throw "Bad value for side thrown in getCountHist";
  }
  else throw "Bad Evt type when trying to get count hist";
};

void EvtRateHandler::dataReader() {

  hisCounts[0] = new TH1D("East Event Hist","EastCounts",numEnergyBins,0.,1200.);
  hisCounts[1] = new TH1D("West Event Hist","WestCounts",numEnergyBins,0.,1200.);

  bool separate23 = false;

  // 2/3 separation
  if ( analysisChoice==std::string("C") ||
       analysisChoice==std::string("H") ||
       analysisChoice==std::string("J") ||
       analysisChoice==std::string("K") )
    separate23 = true;

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
  
  //Set branch addresses
  if (UKdata) {
    sprintf(temp,"replay_pass3_%i.root",runNumber);
    std::string infile = inputDir+"/"+std::string(temp);
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
    sprintf(temp,"spec_%i.root",runNumber);
    std::string infile = inputDir+"/"+std::string(temp);
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
  int counter=0;
  
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
      counter++;
      
      r2E=EmwpcX*EmwpcX+EmwpcY*EmwpcY;
      r2W=WmwpcX*WmwpcX+WmwpcY*WmwpcY;
      
      if ( r2E<(fiducialCut*fiducialCut) && r2W<(fiducialCut*fiducialCut ) ) {
	
	
	if ( separate23 ) {
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
  std::cout << "Beta Events: " << counter << std::endl;
  input->Close();
  if (input) delete input;
};

void SimEvtRateHandler::dataReader() {

  hisCounts[0] = new TH1D("East Event Hist","EastCounts",numEnergyBins,0.,1200.);
  hisCounts[1] = new TH1D("West Event Hist","WestCounts",numEnergyBins,0.,1200.);

  bool separate23 = false;

  // 2/3 separation
  if ( analysisChoice==std::string("C") ||
       analysisChoice==std::string("H") ||
       analysisChoice==std::string("J") ||
       analysisChoice==std::string("K") )
    separate23 = true;

  //Event types
  bool Type0=false, Type1=false, Type2=false, Type3=false;
  if ( analysisChoice=="A" || analysisChoice=="B" || analysisChoice=="C" || analysisChoice=="D" ) Type0 = true;
  if ( analysisChoice=="A" || analysisChoice=="B" || analysisChoice=="C" || analysisChoice=="D" ) Type1 = true;
  if ( analysisChoice=="A" || analysisChoice=="C" || analysisChoice=="G" || analysisChoice=="H" || analysisChoice=="J" ) Type2 = true;
  if ( analysisChoice=="A" || analysisChoice=="C" || analysisChoice=="G" || analysisChoice=="H" || analysisChoice=="K") Type3 = true;
  
  char temp[100];

  double TimeE=0., TimeW=0., Erecon=0., MWPCEnergyE=0., MWPCEnergyW=0.; //Branch Variables being read in
  int PID, Side, Type;

  double AsymWeight=1.;
  double mwpcPosE[3]={0.};
  double mwpcPosW[3]={0.}; //holds the position of the event in the MWPC for simulated data

  sprintf(temp,"/beta/revCalSim_%i_Beta.root",runNumber);
  std::string infile = inputDir+"/"+std::string(temp);
  TFile *input = new TFile(infile.c_str(), "READ");
  TTree *Tin = (TTree*)input->Get("revCalSim");
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
	
	
	if ( separate23 ) {
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
  if (input) delete input;
  
};


BGSubtractedRate::BGSubtractedRate(int run, int bgRun, double enBin, double fidCut, bool ukdata, bool sim, bool applyAsym, bool unblind): runNumber(run), BGrunNumber(bgRun), EnergyBinWidth(enBin), fiducialCut(fidCut), UKdata(ukdata), Simulation(sim), applyAsymmetry(applyAsym), UNBLIND(unblind) {
  int numBins = int(1200./double(EnergyBinWidth));
  BetaRateE.resize(4, std::vector <double> (numBins,0.));
  BGRateE.resize(4, std::vector <double> (numBins,0.));
  FinalRateE.resize(4, std::vector <double> (numBins,0.));
  BetaRateW.resize(4, std::vector <double> (numBins,0.));
  BGRateW.resize(4, std::vector <double> (numBins,0.));
  FinalRateW.resize(4, std::vector <double> (numBins,0.));
  BetaRateErrorE.resize(4, std::vector <double> (numBins,0.));
  BGRateErrorE.resize(4, std::vector <double> (numBins,0.));
  FinalRateErrorE.resize(4, std::vector <double> (numBins,0.));
  BetaRateErrorW.resize(4, std::vector <double> (numBins,0.));
  BGRateErrorW.resize(4, std::vector <double> (numBins,0.));
  FinalRateErrorW.resize(4, std::vector <double> (numBins,0.));

  runLengthBeta.resize(2,0.);
  runLengthBG.resize(2,0.);
};


void BGSubtractedRate::calcBGSubtRates() {

  if (Simulation) 
    {
      std::string indir = std::string(getenv("REVCALSIM"));
      SimEvtRateHandler *evt = new SimEvtRateHandler(runNumber, indir, EnergyBinWidth,fiducialCut,applyAsymmetry, UNBLIND);
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
      std::cout << "Loaded Rates for run " << runNumber << std::endl;
      CalcFinalRate();
    }
};

std::vector<double > BGSubtractedRate::returnBGSubtRate(int side, int etype) {
  if (etype==0 || etype==1 || etype==2 || etype==3) {  
    if (side==0) return FinalRateE[etype];
    else if (side==1) return FinalRateW[etype];
    else throw "Invalid side when returning BG subtracted rate";
  }
  else throw "Invalid event type when returning BG subtracted rate";
};

std::vector<double > BGSubtractedRate::returnBGSubtRateError(int side, int etype) {
  if (etype==0 || etype==1 || etype==2 || etype==3) {  
    if (side==0) return FinalRateErrorE[etype];
    else if (side==1) return FinalRateErrorW[etype];
    else throw "Invalid side when returning BG subtracted rate error";
  }
  else throw "Invalid event type when returning BG subtracted rate error";
};

int BGSubtractedRate::getBackgroundRun(int run) {
  std::string dbAddress = std::string(getenv("UCNADBADDRESS"));
  std::string dbname = std::string(getenv("UCNADB"));
  std::string dbUser = std::string(getenv("UCNADBUSER"));
  std::string dbPass = std::string(getenv("UCNADBPASS"));
  
  char cmd[500];
  sprintf(cmd,"SELECT asym_oct FROM run WHERE run_number=%i;",run);

  SQLdatabase *db = new SQLdatabase(dbname, dbAddress, dbUser, dbPass);
  db->fetchQuery(cmd);
  std::string asymOct = db->returnQueryEntry();
  std::string type="";
  int ntries=0;

  if (asymOct=="A2" || asymOct=="A5" || asymOct=="B2" || asymOct=="B5") {
    while (type!="A1" && type!="A4" && type!="B1" && type!="B4")
      {
	run--;
	sprintf(cmd,"SELECT asym_oct FROM run WHERE run_number=%i;",run);
	db->fetchQuery(cmd);
	//if (db->queryReturnsTrue()) 
	type = db->returnQueryEntry();
	ntries++;
      }
  }
  else if (asymOct=="A7" || asymOct=="A10" || asymOct=="B7" || asymOct=="B10") {
    while (type!="A9" && type!="A12" && type!="B9" && type!="B12") {
      run++;
      sprintf(cmd,"SELECT asym_oct FROM run WHERE run_number=%i;",run);
      db->fetchQuery(cmd);
      //if (db->queryReturnsTrue()) 
      type = db->returnQueryEntry();
      ntries++;
    }
  }
  delete db;
  if (ntries<10) return run;
  else throw "Run Number chosen is not a Beta Decay Run and no BG run exists";
};

void BGSubtractedRate::CreateRateHistograms() {

};

void BGSubtractedRate::LoadRatesByBin() {
  std::string indir;

  if (UKdata) indir = std::string(getenv("REPLAY_PASS3"));
  else indir = std::string(getenv("UCNAOUTPUTDIR"))+"/hists";

  EvtRateHandler *evtBG = new EvtRateHandler(BGrunNumber, indir, EnergyBinWidth, fiducialCut, UKdata, UNBLIND);
  evtBG->CalcRates();
  BGRateE = evtBG->getRateVectors(0);
  BGRateErrorE = evtBG->getRateErrors(0);
  BGRateW = evtBG->getRateVectors(1);
  BGRateErrorW = evtBG->getRateErrors(1);
  runLengthBG[0] = evtBG->returnRunLength(0);
  runLengthBG[1] = evtBG->returnRunLength(1);
  delete evtBG;
    
  EvtRateHandler *evt = new EvtRateHandler(runNumber, indir, EnergyBinWidth, fiducialCut, UKdata, UNBLIND);
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
  std::cout << "Calculating Final Rates for run " << runNumber << std::endl;
  for (unsigned int evtType = 0; evtType<4; evtType++) {
    if (BetaRateE[evtType].size()==BGRateE[evtType].size()) {
      FinalRateE.resize(4,std::vector<double>(BetaRateE.size(),0.));
	for (unsigned int i=0; i<BetaRateE[evtType].size(); i++)
	  {
	    FinalRateE[evtType][i] = BetaRateE[evtType][i]-BGRateE[evtType][i];
	    FinalRateErrorE[evtType][i] = sqrt(power(BetaRateErrorE[evtType][i],2)+power(BGRateErrorE[evtType][i],2));
	    //std::cout << BetaRateE[evtType][i] << " " << BGRateE[evtType][i] << " " << FinalRateE[evtType][i] << std::endl;
	  }
      }
    else throw "Number of energy bins do not agree between Beta and BG runs. Can't calculate final rate!";

    if (BetaRateW[evtType].size()==BGRateW[evtType].size()) {
      FinalRateW.resize(4,std::vector<double>(BetaRateW.size(),0.));
	for (unsigned int i=0; i<BetaRateW[evtType].size(); i++)
	  {
	    FinalRateW[evtType][i] = BetaRateW[evtType][i]-BGRateW[evtType][i];
	    FinalRateErrorW[evtType][i] = sqrt(power(BetaRateErrorW[evtType][i],2)+power(BGRateErrorW[evtType][i],2));
	    //std::cout << BetaRate[i] << " " << BGRate[i] << " " << FinalRate[i] << std::endl;
	  }
      }
    else throw "Number of energy bins do not agree between Beta and BG runs. Can't calculate final rate!";
  }
};


