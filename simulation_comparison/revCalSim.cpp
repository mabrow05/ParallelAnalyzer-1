/* Code to take a run number, retrieve it's runperiod, and contruct the 
weighted spectra which would be seen as a reconstructed energy on one side
of the detector. Also applies the trigger functions */

#include "revCalSim.h"

#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <TRandom3.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TChain.h>
#include <TMath.h>

using namespace std;

vector < vector < Double_t > > getTriggerFunctionParams(Int_t XeRunPeriod, Int_t nParams) {
  Char_t file[500];
  sprintf(file,"%s/trigger_functions_XePeriod_%i.dat",getenv("TRIGGER_FUNC"),XeRunPeriod);
  ifstream infile(file);
  vector < vector <Double_t> > func;
  func.resize(2,vector <Double_t> (nParams,0.));
  //cout << "made it here\n";
  for (Int_t side = 0; side<2; side++) {
    Int_t param = 0;
    while (param<nParams) {
      infile >> func[side][param];
      cout << func[side][param] << " ";
      param++;
    }
    cout << endl;
  }
  infile.close();
  return func;
}

Double_t triggerProbability(vector <Double_t> params, Double_t En) {
  Double_t prob = params[0]+params[1]*TMath::Erf((En-params[2])/params[3])
    + params[4]*TMath::Gaus(En,params[5],params[6]);
  return prob;
}

UInt_t getSrcRunPeriod(Int_t runNumber) {
  UInt_t calibrationPeriod=0;
  if (runNumber <= 17297) {
    calibrationPeriod=1;
  }
  else if (runNumber <= 17439) {
    calibrationPeriod=2;
  }
  else if (runNumber <= 17734) {
    calibrationPeriod=3;
  }
  else if (runNumber <= 17955) {
    calibrationPeriod=4;
  }
  else if (runNumber <= 18386) {
    calibrationPeriod=5;
  }
  else if (runNumber <= 18683) {
    calibrationPeriod=6;
  }
  else if (runNumber <= 18994) {
    calibrationPeriod=7;
  }
  else if (runNumber <= 19239) {
    calibrationPeriod=8;
  }
  else if (runNumber <= 19544) {
    calibrationPeriod=9;
  }
  else if (runNumber <= 20000) {
    calibrationPeriod=11;
  }
  return calibrationPeriod;
}

UInt_t getXeRunPeriod(Int_t runNumber) {
  UInt_t calibrationPeriod=0;
  if (runNumber <= 18080) {
    calibrationPeriod=2;
  }
  else if (runNumber <= 18389) {
    calibrationPeriod=3;
  }
  else if (runNumber <= 18711) {
    calibrationPeriod=4;
  }
  else if (runNumber <= 19238) {
    calibrationPeriod=5;
  }
  //else if (runNumber <= 19872) {
  //calibrationPeriod=6;
  //}
  else if (runNumber <= 20000) {
    calibrationPeriod=7;
  }
  else {
    cout << "Bad run number\n";
    exit(0);}

  return calibrationPeriod;
}

vector <Int_t> getPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"../residuals/PMT_runQuality_master.dat"); 
  ifstream pmt;
  pmt.open(temp);
  Int_t run_hold;
  while (pmt >> run_hold >> pmtQuality[0] >> pmtQuality[1] >> pmtQuality[2]
	 >> pmtQuality[3] >> pmtQuality[4] >> pmtQuality[5]
	 >> pmtQuality[6] >> pmtQuality[7]) {
    if (run_hold==runNumber) break;
    if (pmt.fail()) break;
  }
  pmt.close();
  if (run_hold!=runNumber) {
    cout << "Run not found in PMT quality file!" << endl;
    exit(0);
  }
  return pmtQuality;
}

vector < Double_t > GetAlphaValues(Int_t runPeriod)
{
  Char_t temp[500];
  vector < Double_t > alphas (8,0.);
  sprintf(temp,"../smeared_peaks/nPE_Kev_%i.dat",runPeriod);
  ifstream infile;
  infile.open(temp);
  Int_t i = 0;

  while (infile >> alphas[i]) i++;
  return alphas;
}

void SetUpTree(TTree *tree) {
  tree->Branch("PID", &PID, "PID/I");
  tree->Branch("side", &side, "side/I");
  tree->Branch("type", &type, "type/I");
  tree->Branch("EvisTot", &EvisTot,"EvisTot/D");
  
  tree->Branch("Evis",&evis,"EvisE/D:EvisW");
  tree->Branch("Edep",&edep,"EdepE/D:EdepW");
  tree->Branch("EdepQ",&edepQ,"EdepQE/D:EdepQW");
  
  tree->Branch("time",&Time,"timeE/D:timeW");
  tree->Branch("MWPCEnergy",&mwpcE,"MWPCEnergyE/D:MWPCEnergyW");
  tree->Branch("MWPCPos",&mwpc_pos,"MWPCPosE/D:MWPCPosW");
  tree->Branch("ScintPos",&scint_pos,"ScintPosE/D:ScintPosW");
  tree->Branch("PMT_Evis",&pmt_Evis,"Evis0/D:Evis1:Evis2:Evis3:Evis4:Evis5:Evis6:Evis7:weight0:weight1:weight2:weight3:weight4:weight5:weight6:weight7");
  
}
  

void revCalSimulation (Int_t runNumber, string source) 
{
  //First get the number of total electron events from the data file
  Char_t temp[500];
  sprintf(temp,"%s/replay_pass4_%i.root",getenv("REPLAY_PASS4"),runNumber);
  TFile *dataFile = new TFile(temp,"READ");
  TTree *data = (TTree*)(dataFile->Get("pass4"));
  sprintf(temp,"type_pass4==0 || type_pass4==1 || type_pass4==2");
  UInt_t BetaEvents = data->GetEntries(temp);
  cout << "Electron Events in Data file: " << BetaEvents << endl;
  delete data;
  dataFile->Close();

  //Create simulation output file
  Char_t outputfile[500];
  sprintf(outputfile,"%s/revCalSim_%i.root",getenv("REVCALSIM"),runNumber);
  //sprintf(outputfile,"revCalSim_%i.root",runNumber);
  TFile *outfile = new TFile(outputfile, "RECREATE")
;
  vector <Int_t> pmtQuality = getPMTQuality(runNumber); // Get the quality of the PMTs for that run
  UInt_t calibrationPeriod = getSrcRunPeriod(runNumber); // retrieve the calibration period for this run
  UInt_t XePeriod = getXeRunPeriod(runNumber); // Get the proper Xe run period for the Trigger functions
  vector <Double_t> alphas = GetAlphaValues(calibrationPeriod); // fill vector with the alpha (nPE/keV) values for this run period
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7); // 2D vector with trigger function for East side and West side in that order

  //Setup the output tree
  TTree *tree = new TTree("revCalSim", "revCalSim");
  SetUpTree(tree); //Setup the output tree and branches

  //Histograms of event types for quick checks
  vector <TH1D*> finalEn (6,NULL);
  finalEn[0] = new TH1D("finalE0", "Simulated Weighted Sum East Type 0", 400, 0., 1200.);
  finalEn[1] = new TH1D("finalW0", "Simulated Weighted Sum West Type 0", 400, 0., 1200.);
  finalEn[2] = new TH1D("finalE1", "Simulated Weighted Sum East Type 1", 400, 0., 1200.);
  finalEn[3] = new TH1D("finalW1", "Simulated Weighted Sum West Type 1", 400, 0., 1200.);
  finalEn[4] = new TH1D("finalE23", "Simulated Weighted Sum East Type 2/3", 400, 0., 1200.);
  finalEn[5] = new TH1D("finalW23", "Simulated Weighted Sum West Type 2/3", 400, 0., 1200.);


  //Read in simulated data and put in a TChain
  TChain *chain = new TChain("anaTree");
  for (int i=0; i<250; i++) {
    sprintf(temp,"/extern/mabrow05/ucna/geant4work/output/10mil_2011-2012/%s/analyzed_%i.root",source.c_str(),i);
    //sprintf(temp,"../../../data/analyzed_%i.root",i);
    chain->AddFile(temp);
  }
  
  // Set the addresses of the information read in from the simulation file
  chain->SetBranchAddress("MWPCEnergy",&mwpcE);
  chain->SetBranchAddress("time",&Time);
  chain->SetBranchAddress("Edep",&edep);
  chain->SetBranchAddress("EdepQ",&edepQ);
  chain->SetBranchAddress("MWPCPos",&mwpc_pos);
  chain->SetBranchAddress("ScintPos",&scint_pos);
  //chain->GetBranch("EdepQ")->GetLeaf("EdepQE")->SetAddress(&EdepQE);
  //chain->GetBranch("EdepQ")->GetLeaf("EdepQW")->SetAddress(&EdepQW);
  //chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&MWPCEnergyE);
  //chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&MWPCEnergyW);

  //Trigger booleans
  bool EastScintTrigger, WestScintTrigger, EMWPCTrigger, WMWPCTrigger;
  Double_t MWPCThreshold=0.0;

  //Set random number generator
  TRandom3 *rand = new TRandom3(0);
  
  //Get total number of events in TChain
  UInt_t nevents = chain->GetEntries();
  cout << "events = " << nevents << endl;
 
  //Start from random position in evt sequence
  UInt_t evtStart = rand->Rndm()*nevents;
  UInt_t evtTally = 0; //To keep track ot the number of events 
  UInt_t evt = evtStart; //current event number

  //Read in events and determine evt type based on triggers
  while (evtTally<=BetaEvents) {
    if (evt>nevents) evt=0; //Wrapping the events back to the beginning
    EastScintTrigger = WestScintTrigger = EMWPCTrigger = WMWPCTrigger = false; //Resetting triggers each event

    chain->GetEvent(evt);
   
    //MWPC triggers
    if (mwpcE.MWPCEnergyE>MWPCThreshold) EMWPCTrigger=true;
    if (mwpcE.MWPCEnergyW>MWPCThreshold) WMWPCTrigger=true;

    //East Side smeared PMT energies
    for (UInt_t p=0; p<4; p++) {
      if (pmtQuality[p]) { //Check to make sure PMT was functioning
	pmt_Evis.Evis[p] = -1.; // set the Evis for the while loop to return positive first time
	while (pmt_Evis.Evis[p]<0.) { //Must have positive energies since we force the PMTs to intercept 0
	  pmt_Evis.Evis[p] = rand->Gaus(edepQ.EdepQE,sqrt(edepQ.EdepQE/alphas[p]));
	}
	if (pmt_Evis.Evis[p]>0.001) {
	  pmt_Evis.weight[p] = alphas[p]/pmt_Evis.Evis[p];
	}
	else {pmt_Evis.weight[p]=alphas[p]/0.001;} //This sets a hard cut on the weight of an event so that events with zero energy still carry weight.. Probably not the most effective way to do this, but it should address low energy behavior
      }
      // If PMT quality failed, set the evis and the weight to zero for this PMT
      else {
	pmt_Evis.weight[p]=0.;
	pmt_Evis.Evis[p]=0.;
      }
    }

    //Calculate the weighted energy on a side
    Double_t numer=0., denom=0.;
    for (UInt_t p=0;p<4;p++) {
      numer+=pmt_Evis.Evis[p]*pmt_Evis.weight[p];
      denom+=pmt_Evis.weight[p];
    }
    //Now we apply the trigger probability
    Double_t totalEnE = numer/denom;
    evis.EvisE = totalEnE;
    Double_t triggProb = triggerProbability(triggerFunc[0],totalEnE);
 
    //Set East Scint Trigger to true if event passes triggProb
    if (rand->Rndm(0)<triggProb && evis.EvisE>0.) {
      EastScintTrigger=true;
    }
    	
    
    //West Side
    for (UInt_t p=4; p<8; p++) {
      if (pmtQuality[p]) { //Check to make sure PMT was functioning
	pmt_Evis.Evis[p] = -1.; // set the Evis for the while loop to return positive first time
	while (pmt_Evis.Evis[p]<0.) { //Must have positive energies since we force the PMTs to intercept 0
	  pmt_Evis.Evis[p] = rand->Gaus(edepQ.EdepQW,sqrt(edepQ.EdepQW/alphas[p]));
	}
	if (pmt_Evis.Evis[p]>0.001) {
	  pmt_Evis.weight[p] = alphas[p]/pmt_Evis.Evis[p];
	}
	else {pmt_Evis.weight[p]=alphas[p]/0.001;} //This sets a hard cut on the weight of an event so that events with zero energy still carry weight.. Probably not the most effective way to do this, but it should address low energy behavior
      }
      // If PMT quality failed, set the evis and the weight to zero for this PMT
      else {
	pmt_Evis.weight[p]=0.;
	pmt_Evis.Evis[p]=0.;
      }
    }

    //Calculate the total weighted energy
    numer=denom=0.;
    for (UInt_t p=4;p<8;p++) {
      numer+=pmt_Evis.Evis[p]*pmt_Evis.weight[p];
      denom+=pmt_Evis.weight[p];
    }
    //Now we apply the trigger probability
    Double_t totalEnW = numer/denom;
    evis.EvisW = totalEnW;
    triggProb = triggerProbability(triggerFunc[1],totalEnW);
    //Fill histograms if event passes trigger function
    if (rand->Rndm(0)<triggProb && evis.EvisW>0.) {
      WestScintTrigger = true;      
    }
    
    //Fill total Energy loss
    EvisTot = evis.EvisW+evis.EvisE+mwpcE.MWPCEnergyE+mwpcE.MWPCEnergyW;

    //Fill proper total event histogram based on event type

    //Type 0 East
    if (EastScintTrigger && EMWPCTrigger && !WestScintTrigger && !WMWPCTrigger) {
      PID=1;
      type=0;
      side=0;
      finalEn[0]->Fill(evis.EvisE);
      //cout << "Type 0 East E = " << totalEnE << endl;
    }
    //Type 0 West
    else if (WestScintTrigger && WMWPCTrigger && !EastScintTrigger && !EMWPCTrigger) {
      PID=1;
      type=0;
      side=1;
      finalEn[1]->Fill(totalEnW);
    }
    //Type 1 
    else if (EastScintTrigger && EMWPCTrigger && WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=1;
      //East
      if (Time.timeE<Time.timeW) {
	finalEn[2]->Fill(totalEnE);
	side=0;
      }
      //West
      else if (Time.timeE>Time.timeW) {
	finalEn[3]->Fill(totalEnW);
	side=1;
      }
    }
    else if (EastScintTrigger && EMWPCTrigger && !WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=2;
      side=0;
      finalEn[4]->Fill(totalEnE);
      //cout << "Type 2/3 East E = " << totalEnE << endl;
    }
    else if (!EastScintTrigger && EMWPCTrigger && WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=2;
      side=1;
      finalEn[5]->Fill(totalEnW);
      //cout << "Type 2/3 East W = " << totalEnW << endl;
    }   
    //Gamma events and missed events (Type 4)
    else {
      if (!WMWPCTrigger && !EMWPCTrigger) {
	if (EastScintTrigger && !WestScintTrigger) {
	  PID=0;
	  type=0;
	  side=0;
	}
	else if (!EastScintTrigger && WestScintTrigger) {
	  PID=0;
	  type=0;
	  side=1;
	}
	else if (EastScintTrigger && WestScintTrigger) {
	  PID=0;
	  type=1;
	  if (Time.timeE<Time.timeW) {
	    side=0;
	  }
	  else if (Time.timeE>Time.timeW) {
	    side=1;
	  }
	}
      }
      else {
	PID=6;
	type=4;
	side=2; //Unknown Side
      }
    }
	
    // Increment the event tally if the event was PID = 1 (electron)
    if (PID==1) evtTally++;
    evt++;
    tree->Fill();
    
    //cout << "filled event " << evt << endl;
  }
  delete chain;
  outfile->Write();
  outfile->Close();
}

int main(int argc, char *argv[]) {
  string source = string(argv[2]);
  revCalSimulation(atoi(argv[1]),source);

  //tests
  /*UInt_t XePeriod = getXeRunPeriod(atoi(argv[1]));
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7);
  Double_t triggProbE = triggerProbability(triggerFunc[0],25.);
  Double_t triggProbW = triggerProbability(triggerFunc[1],25.);
  cout << triggProbE << " " << triggProbW << endl;*/


}
  
