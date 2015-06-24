/*
Class which allows for generic reading of data and construction 
of histograms for rates of each event type. This can be any type of data coming
from UCNA root files. What is stored should be histograms of event rates binned
by energy along with the run time for that particular run.

Also included is a class to handle reading in reverse calibrated 
simulation data
*/

#include "EvtRateHandler.hh"
#include <cstdio>

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
  if (flipperStatus=="On") return 1;
  else if (flipperStatus=="Off") return -1;
  else throw "Polarization isn't applicaple or you chose a Depol Run";
  
};

void EvtRateHandler::CalcRates(int evtType, double enBins, double fidCut) {
  char temp1[200], temp2[200];
  fiducialCut = fidCut;
  numEnergyBins = int(1200./enBins);
  
  if (evtType==0 || evtType==1 || evtType==23)
    {
      if (evtType==23) {
	sprintf(temp1,"East Type 2/3 event rate");
	sprintf(temp2,"West Type 2/3 event rate");
      }
      else {
	sprintf(temp1,"East Type %i event rate",evtType);
	sprintf(temp2,"West Type %i event rate",evtType);
      }   
      
      rateE = new TH1D("EastRate",temp1,numEnergyBins,0.,1200.);
      rateW = new TH1D("WestRate",temp2,numEnergyBins,0.,1200.);
      dataReader(evtType);
    }
  else throw "Attempt to calculate rate for non-existent event type";

};

std::vector<double> EvtRateHandler::getRateVector(int side) {
  rateEvec.resize(numEnergyBins,0);
  rateWvec.resize(numEnergyBins,0);
  if (side==0)
    {
      for (int i=0; i<numEnergyBins; i++) {rateEvec[i]=rateE->GetBinContent(i);}
      return rateEvec;
    }
  else if (side==1) {
    for (int i=0; i<numEnergyBins; i++) {rateWvec[i]=rateW->GetBinContent(i);}
    return rateWvec;
  }
  else throw "Bad value for side thrown in getRateVector";
};

TH1D EvtRateHandler::getRateHist(int side) {
  if (side==0) return *rateE;
  else if (side==1) return *rateW;
  else throw "Bad value for side thrown in getRateHist";
};

void EvtRateHandler::dataReader(int evtType) {
  char temp[100];
  sprintf(temp,"spec_%i.root",runNumber);
  std::string infile = inputDir+"/"+std::string(temp);
  TFile *input = new TFile(infile.c_str(), "READ");
  TTree *Tin = (TTree*)input->Get("phys"); //TODO: make sure this is the correct tree name
  
  //Set branch addresses
  Tin->SetBranchAddress("PID", &PID);
  Tin->SetBranchAddress("Type", &Type);
  Tin->SetBranchAddress("Side", &Side); 
  Tin->SetBranchAddress("Etrue",&Erecon);
  Tin->SetBranchAddress("TimeE",&TimeE);
  Tin->SetBranchAddress("TimeW",&TimeW);
  Tin->GetBranch("xEmpm")->GetLeaf("center")->SetAddress(&EmwpcX);
  Tin->GetBranch("yEmpm")->GetLeaf("center")->SetAddress(&EmwpcY);
  Tin->GetBranch("xWmpm")->GetLeaf("center")->SetAddress(&WmwpcX);
  Tin->GetBranch("yWmpm")->GetLeaf("center")->SetAddress(&WmwpcY);

  unsigned int nevents = Tin->GetEntriesFast();
 
  //Determine total time on each side
  Tin->GetEvent(nevents-1);
 
  float totalTimeE = TimeE;
  double EastWeight = 1./double(totalTimeE);
  float totalTimeW = TimeW;
  double WestWeight = 1./double(totalTimeW);

  //Run over all events in file, fill output histogram with rate
  
  float r2 = 0.; //position of event squared

  for (unsigned int i=0; i<nevents; i++)
    {
      Tin->GetEvent(i);
      if (Type==evtType && PID==1) 
	{
	  if (Side==0)
	    {
	      r2=EmwpcX*EmwpcX+EmwpcY*EmwpcY;
	      if (r2<(fiducialCut*fiducialCut))
		{
		  rateE->Fill(double(Erecon), EastWeight);
		}
	    }
	  else if (Side==1)
	    {
	      r2=WmwpcX*WmwpcX+WmwpcY*WmwpcY;
	      if (r2<(fiducialCut*fiducialCut))
		{
		  rateW->Fill(double(Erecon), WestWeight);
		}
	    }
	}
    }

  input->Close();
  if (input) delete input;
};

void SimEvtRateHandler::dataReader(int evtType) {
  char temp[100];
  sprintf(temp,"%i_sim.root",runNumber);
  std::string infile = inputDir+"/"+std::string(temp);
  TFile *input = new TFile(infile.c_str(), "READ");
  TTree *Tin = (TTree*)input->Get("SimOut"); //TODO: make sure this is the correct tree name
  double mwpcPos[2][3]; //holds the position of the event in the MWPC for simulated data
  double EreconSim; //This is double whereas data comes in as float
  //Set branch addresses
  Tin->GetBranch("RevCalSimData")->GetLeaf("PID")->SetAddress(&PID);
  Tin->GetBranch("RevCalSimData")->GetLeaf("type")->SetAddress(&Type);
  Tin->GetBranch("RevCalSimData")->GetLeaf("side")->SetAddress(&Side);
  Tin->GetBranch("RevCalSimData")->GetLeaf("Erecon")->SetAddress(&EreconSim);
  Tin->GetBranch("SimData")->GetLeaf("mwpcPos")->SetAddress(mwpcPos);

  unsigned int nevents = Tin->GetEntriesFast();
  //std::cout << nevents << std::endl;

  //Run over all events in file, fill output histogram with rate
  
  double r2 = 0.; //position of event squared

  for (unsigned int i=0; i<nevents; i++)
    {
      Tin->GetEvent(i);
      if (Type==evtType && PID==1) 
	{
	  if (Side==0)
	    {
	      r2=mwpcPos[0][0]*mwpcPos[0][0]+mwpcPos[0][1]*mwpcPos[0][1];
	      //std::cout << r2 << std::endl;
	      if (r2<(fiducialCut*fiducialCut))
		{
		  rateE->Fill(EreconSim);
		  //std::cout << EreconSim << std::endl;
		}
	    }
	  else if (Side==1)
	    {
	      r2=mwpcPos[1][0]*mwpcPos[1][0]+mwpcPos[1][1]*mwpcPos[1][1]; 
	      //std::cout << r2 << std::endl;
	      if (r2<(fiducialCut*fiducialCut))
		{
		  rateW->Fill(EreconSim);
		  //std::cout << EreconSim << std::endl;
		}
	    }
	}
    }

  input->Close();
  if (input) delete input;
  
};

std::vector<double > BGSubtractedRate::ReturnBGSubtRate(int side) {
  double numBins = int(1200./double(EnergyBinWidth));
  //BetaRate.resize(numBins,0.);
  //BGRate.resize(numBins,0.);
  //FinalRate.resize(numBins,0.);

  if (Simulation) //We don't save the histograms for simulations because they are just the histograms on file
    {
      std::string indir = std::string(getenv("UCNA_REVCAL_DIR"));
      SimEvtRateHandler *evt = new SimEvtRateHandler(runNumber, indir);
      evt->CalcRates(evtType,EnergyBinWidth,fiducialCut);
      FinalRate = evt->getRateVector(side);
      delete evt;
      return FinalRate; //Note that there is no background in the simulation so we just return the data counts
    }
  else 
    {
      LoadRatesByBin(side);
      CalcFinalRate();
      return FinalRate;
    }
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
  
  if (asymOct=="A2" || asymOct=="A5" || asymOct=="B2" || asymOct=="B5") return run-1;
  else if (asymOct=="A7" || asymOct=="A10" || asymOct=="B7" || asymOct=="B10") return run+2;
  else throw "Run Number chosen is not a Beta Decay Run and no BG run exists";

};

void BGSubtractedRate::CreateRateHistograms() {

};

void BGSubtractedRate::LoadRatesByBin(int side) {
  int bgRun = getBackgroundRun(runNumber);

  std::string indir = std::string(getenv("UCNAOUTPUTDIR"))+"/hists";
  EvtRateHandler *evtBG = new EvtRateHandler(bgRun, indir);
  evtBG->CalcRates(evtType,EnergyBinWidth,fiducialCut);
  BGRate = evtBG->getRateVector(side);
  delete evtBG;

  EvtRateHandler *evt = new EvtRateHandler(runNumber, indir);
  evt->CalcRates(evtType,EnergyBinWidth,fiducialCut);
  BetaRate = evt->getRateVector(side);
  delete evt;
  
};

void BGSubtractedRate::CalcFinalRate()  {

  if (BetaRate.size()==BGRate.size())
    {
      FinalRate.resize(BetaRate.size(),0.);
      for (unsigned int i=0; i<BetaRate.size(); i++)
	{
	  FinalRate[i] = BetaRate[i]-BGRate[i];
	  std::cout << BetaRate[i] << " " << BGRate[i] << " " << FinalRate[i] << std::endl;
	}
    }
  else throw "Number of energy bins do not agree between Beta and BG runs. Can't calculate final rate!";

};


