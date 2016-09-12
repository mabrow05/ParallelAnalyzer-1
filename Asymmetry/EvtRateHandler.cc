/*
Class which allows for generic reading of data and construction 
of histograms for rates of each event type. This can be any type of data coming
from UCNA root files. What is stored should be histograms of event rates binned
by energy along with the run time for that particular run.

Also included is a class to handle reading in reverse calibrated 
simulation data
*/

#include "EvtRateHandler.hh"
#include "MBUtils.hh"
#include <cstdio>
#include <fstream>

EvtRateHandler::EvtRateHandler(int run, const std::string& inDir, double enBinWidth, double fidCut, bool ukdata, bool ub) : runNumber(run),inputDir(inDir),fiducialCut(fidCut),UKdata(ukdata), pol(0), unblinded(ub) {
  rateE.resize(4,NULL);
  rateW.resize(4,NULL);

  numEnergyBins = (int)(1200./enBinWidth);

  rateEvec.resize(4, std::vector<double> (numEnergyBins,0.)); 
  rateWvec.resize(4,std::vector <double> (numEnergyBins,0.));
  rateEerr.resize(4, std::vector<double> (numEnergyBins,0.)); 
  rateWerr.resize(4,std::vector <double> (numEnergyBins,0.));

  //Loading information from log file
  std::string logFilePath = std::string(getenv("RUN_INFO_FILES"))+"runInfo_"+itos(runNumber)+".dat";
  ifstream logFile(logFilePath.c_str());
  std::vector <std::string> title(4);
   std::vector <double> value(4);
  for (int i=0; i<4; i++) {
    logFile >> title[i] >> value[i];
  }
  logFile.close();

  if (unblinded) {
    runLength[0] = value[2];
    runLength[1] = value[2];
  }
  else {
    runLength[0] = value[0];
    runLength[1] = value[1];
  }
  UCNMonIntegral = value[3];
 
  std::cout << title[0] << "\t" << runLength[0] << std::endl;
  std::cout << title[1] << "\t" << runLength[1] << std::endl;
  std::cout << title[3] << "\t" << UCNMonIntegral << std::endl;
  //std::cout << holdTitle << "\t" << runLength[0] << std::endl;
}

EvtRateHandler::~EvtRateHandler() {
  for (unsigned int i = 0; i<4; i++) {
    if (rateE[i]) delete rateE[i];
    if (rateW[i]) delete rateW[i]; 
  }
};

int EvtRateHandler::polarization(int run) {
 
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

void EvtRateHandler::CalcRates() {
  char temp1[200], temp2[200], temp3[100], temp4[100];
  
  for (unsigned int evtType = 0; evtType<4; evtType++) {
 
    sprintf(temp1,"East Type %i event rate",evtType);
    sprintf(temp2,"EastRate%i",evtType);
    sprintf(temp3,"West Type %i event rate",evtType);
    sprintf(temp4,"WestRate%i",evtType);
    
    //rateE.push_back(new TH1D(temp2,temp1,numEnergyBins,0.,1200.));
    //rateW.push_back(new TH1D(temp4,temp3,numEnergyBins,0.,1200.));
         
    rateE[evtType] = new TH1D(temp2,temp1,numEnergyBins,0.,1200.);
    rateW[evtType] = new TH1D(temp4,temp3,numEnergyBins,0.,1200.);
  }
  this->dataReader();
      
};

std::vector< std::vector<double> > EvtRateHandler::getRateVectors(int side) {

  if (side==0)
    {
      for (unsigned int j=0; j<4; j++) {
	for (unsigned int i=0; i<numEnergyBins; i++) {rateEvec[j][i]=rateE[j]->GetBinContent(i);}
      }
      return rateEvec;
    }
  else if (side==1) {
    for (unsigned int j=0; j<4; j++) {
      for (unsigned int i=0; i<numEnergyBins; i++) {rateWvec[j][i]=rateW[j]->GetBinContent(i);}
    }
    return rateWvec;
  }
  else throw "Bad value for side thrown in getRateVector";
};

std::vector< std::vector<double> > EvtRateHandler::getRateErrors(int side) {
  rateEerr.resize(4, std::vector<double> (numEnergyBins,0.)); 
  rateWerr.resize(4,std::vector <double> (numEnergyBins,0.));

  if (side==0) {
    for (unsigned int j=0; j<4; j++) {
      for (unsigned int i=0; i<numEnergyBins; i++) {
	rateEerr[j][i]=rateE[j]->GetBinContent(i)*runLength[0]>=1. ? sqrt(rateE[j]->GetBinContent(i)/runLength[0]): 0.;//1./runLength[0]; // How to handle low statistics?
      }
    }
    return rateEerr;
  }
  else if (side==1) {
    for (unsigned int j=0; j<4; j++) {
      for (unsigned int i=0; i<numEnergyBins; i++) {
	rateWerr[j][i]=rateW[j]->GetBinContent(i)*runLength[1]>=1. ? sqrt(rateW[j]->GetBinContent(i)/runLength[1]): 0.;//1./runLength[1];
      }
    }
    return rateWerr;
  }
  else throw "Bad value for side thrown in getRateErrors";
};

TH1D EvtRateHandler::getRateHist(int side, int evtType) {
  if (evtType==0 || evtType==1 || evtType==2 || evtType==3) {
    if (side==0) return *rateE[evtType];
    else if (side==1) return *rateW[evtType];
    else throw "Bad value for side thrown in getRateHist";
  }
  else throw "Bad Evt type when trying to get rate hist";
};

void EvtRateHandler::dataReader() {
  char temp[100];
  
  std::string infile;
  TFile *input;
  TTree *Tin;

  double EmwpcX=0., EmwpcY=0., WmwpcX=0., WmwpcY=0., TimeE=0., TimeW=0., Erecon=0.; //Branch Variables being read in
  float EmwpcX_f=0., EmwpcY_f=0., WmwpcX_f=0., WmwpcY_f=0., TimeE_f=0., TimeW_f=0., Erecon_f=0.; // For reading in data from MPM replays

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
  
  double EastWeight = 1./runLength[0];
  double WestWeight = 1./runLength[1];
  
  std::cout << EastWeight << std::endl;

  //Run over all events in file, fill output histogram with rate
  
  float r2 = 0.; //position of event squared
  int counter=0;
  for (unsigned int i=0; i<nevents; i++)
    {
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
      if (PID==1) 
	{
	  counter++;
	  if (Side==0)
	    {
	      r2=EmwpcX*EmwpcX+EmwpcY*EmwpcY;
	      if (r2<(fiducialCut*fiducialCut))
		{
		  rateE[Type]->Fill(Erecon, EastWeight);
		}
	    }
	  else if (Side==1)
	    {
	      r2=WmwpcX*WmwpcX+WmwpcY*WmwpcY;
	      if (r2<(fiducialCut*fiducialCut))
		{
		  rateW[Type]->Fill(Erecon, WestWeight);
		}
	    }
	}
    }
  std::cout << "Beta Events: " << counter << std::endl;
  input->Close();
  if (input) delete input;
};

void SimEvtRateHandler::dataReader() {
  char temp[100];

  double TimeE=0., TimeW=0., Erecon=0.; //Branch Variables being read in
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


  unsigned int nevents = Tin->GetEntriesFast();
  //std::cout << nevents << std::endl;

  double EastWeight = 1./runLength[0];
  double WestWeight = 1./runLength[1];

  //Run over all events in file, fill output histogram with rate
  
  double r2 = 0.; //position of event squared

  for (unsigned int i=0; i<nevents; i++)
    {
      Tin->GetEvent(i);

      if (!applyAsymmetry) AsymWeight=1.;
   
      if (Type<4 && PID==1) 
	{
	  if (Side==0)
	    {
	      r2=mwpcPosE[0]*mwpcPosE[0]+mwpcPosE[1]*mwpcPosE[1];
	      //std::cout << r2 << std::endl;
	      if (r2<(fiducialCut*fiducialCut))
		{
		  rateE[Type]->Fill(Erecon,AsymWeight*EastWeight);
		  //std::cout << EreconSim << std::endl;
		}
	    }
	  else if (Side==1)
	    {
	      r2=mwpcPosW[0]*mwpcPosW[0]+mwpcPosW[1]*mwpcPosW[1]; 
	      //std::cout << r2 << std::endl;
	      if (r2<(fiducialCut*fiducialCut))
		{
		  rateW[Type]->Fill(Erecon,AsymWeight*WestWeight);
		  //std::cout << EreconSim << std::endl;
		}
	    }
	}
    }

  input->Close();
  if (input) delete input;
  
};


BGSubtractedRate::BGSubtractedRate(int run, double enBin, double fidCut, bool ukdata, bool sim, bool applyAsym): runNumber(run), EnergyBinWidth(enBin), fiducialCut(fidCut), UKdata(ukdata), Simulation(sim), applyAsymmetry(applyAsym) {
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
      SimEvtRateHandler *evt = new SimEvtRateHandler(runNumber, indir, EnergyBinWidth,fiducialCut,applyAsymmetry);
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
  int bgRun = getBackgroundRun(runNumber);
  std::string indir;

  if (UKdata) indir = std::string(getenv("REPLAY_PASS3"));
  else indir = std::string(getenv("UCNAOUTPUTDIR"))+"/hists";

  EvtRateHandler *evtBG = new EvtRateHandler(bgRun, indir, EnergyBinWidth, fiducialCut, UKdata);
  evtBG->CalcRates();
  BGRateE = evtBG->getRateVectors(0);
  BGRateErrorE = evtBG->getRateErrors(0);
  BGRateW = evtBG->getRateVectors(1);
  BGRateErrorW = evtBG->getRateErrors(1);
  runLengthBG[0] = evtBG->returnRunLength(0);
  runLengthBG[1] = evtBG->returnRunLength(1);
  delete evtBG;
    
  EvtRateHandler *evt = new EvtRateHandler(runNumber, indir, EnergyBinWidth, fiducialCut, UKdata);
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


