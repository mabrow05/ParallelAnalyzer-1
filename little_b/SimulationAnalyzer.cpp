/* Code to take a run number, retrieve it's runperiod, and construct the 
weighted spectra which would be seen as a reconstructed energy on one side
of the detector. Also applies the trigger functions */

#include "SimulationAnalyzer.hh"
#include "posMapReader.h"
//#include "sourcePeaks.h"
//#include "runInfo.h"

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
//#include <TSQLServer.h>
//#include <TSQLResult.h>
//#include <TSQLRow.h>



//Function to return the trigger function for each side in a vector in the form vec[side][param]
// where side==0 is East and side==1 is West
std::vector < std::vector < Double_t > > getTriggerFunctionParams(Int_t XeRunPeriod, Int_t nParams = 7) {
  Char_t file[500];
  sprintf(file,"trigger_functions/trigger_functions_XePeriod_%i.dat",XeRunPeriod);
  ifstream infile(file);
  std::vector < std::vector <Double_t> > func;
  func.resize(2,std::vector <Double_t> (nParams,0.));
  //cout << "made it here\n";
  for (Int_t side = 0; side<2; side++) {
    Int_t param = 0;
    while (param<nParams) {
      infile >> func[side][param];
      //std::cout << func[side][param] << " ";
      param++;
    }
    //cout << endl;
  }
  infile.close();
  std::cout << "Loaded trigger functions.\n\n";
  return func;
}

//Function which returns the trigger probability for an erf with 7 parameters.
// These parameters are passed into the function with the quenched energy of the event
Double_t triggerProbability(std::vector <Double_t> params, Double_t En) {
  Double_t prob = params[0]+params[1]*TMath::Erf((En-params[2])/params[3])
    + params[4]*TMath::Gaus(En,params[5],params[6]);
  return prob;
}

// Returns a vector which holds the alpha value for each PMT
std::vector < Double_t > getAlphaValues(Int_t fileNum)
{
  Char_t temp[500];
  std::vector < Double_t > alphas (8,0.);
  sprintf(temp,"nPE_keV/nPE_keV_%i.dat",fileNum);
  ifstream infile;
  infile.open(temp);
  Int_t i = 0;

  while (infile >> alphas[i]) i++;
  return alphas;
}

//Get the conversion from EQ2Etrue
std::vector < std::vector < std::vector <double> > > getEQ2EtrueParams(int runNumber=0) {
  ifstream infile;
  if (runNumber<20000) infile.open("../simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat");
  else infile.open("../simulation_comparison/EQ2EtrueConversion/2012-2013_EQ2EtrueFitParams.dat");
  
  std::vector < std::vector < std::vector < double > > > params;
  params.resize(2,std::vector < std::vector < double > > (3, std::vector < double > (6,0.)));

  char holdType[10];
  int side=0, type=0;
  while (infile >> holdType >> params[side][type][0] >> params[side][type][1] >> params[side][type][2] >> params[side][type][3] >> params[side][type][4] >> params[side][type][5]) { 
    cout << holdType << " " << params[side][type][0] << " " << params[side][type][1] << " " << params[side][type][2] << " " << params[side][type][3] << " " << params[side][type][4] << " " << params[side][type][5] << endl;
    type+=1;
    if (type==3) {type=0; side=1;}
  }
  return params;
};
  
std::vector < std::vector <Double_t> > getLinearityCurve (UInt_t linearityFileIndex) {

}

// Sets up the output tree for the simulated data
void SetUpOutputTree(TTree *tree) {
  tree->Branch("PID", &PID, "PID/I");
  tree->Branch("side", &side, "side/I");
  tree->Branch("type", &type, "type/I");
  tree->Branch("Erecon", &Erecon,"Erecon/D");
  
  tree->Branch("Evis",&evis,"EvisE/D:EvisW");
  tree->Branch("Edep",&edep,"EdepE/D:EdepW");
  tree->Branch("EdepQ",&edepQ,"EdepQE/D:EdepQW");
  tree->Branch("Eprim",&Eprim,"Eprim/D");
  tree->Branch("AsymWeight",&AsymWeight,"AsymWeight/D");
  
  tree->Branch("time",&Time,"timeE/D:timeW");
  tree->Branch("MWPCEnergy",&mwpcE,"MWPCEnergyE/D:MWPCEnergyW");
  tree->Branch("MWPCPos",&mwpc_pos,"MWPCPosE[3]/D:MWPCPosW[3]");
  tree->Branch("ScintPos",&scint_pos,"ScintPosE[3]/D:ScintPosW[3]");
  tree->Branch("ScintPosAdjusted",&scint_pos_adj,"ScintPosAdjE[3]/D:ScintPosAdjW[3]");
  tree->Branch("PMT_Evis",&pmt_Evis,"Evis0/D:Evis1:Evis2:Evis3:Evis4:Evis5:Evis6:Evis7:weight0:weight1:weight2:weight3:weight4:weight5:weight6:weight7");
  
}
  

void revCalSimulation (std::string source, UInt_t numEvents, UInt_t linParamsFile) 
{
  std::cout << "Running reverse calibration for " << source << std::endl;

  UInt_t alphaFileIndex = 0; // the index on the nPE_keV file
  UInt_t XeMapPeriod = 2; // The xenon position map to be used
  UInt_t BetaEvents = numEvents;
  // Set the location of the source and it's spread in mm. srcPos[0=east,1=west][0=x,1=y,2=sigma]
  std::vector < std::vector <Double_t> > srcPos(2, std::vector <Double_t> (3,0.));
  srcPos[0][0] = srcPos[1][0] = 0.;
  srcPos[0][1] = srcPos[1][1] = 0.;
  srcPos[0][2] = srcPos[1][2] = 5.;

  //Decide which simulation to use...
  std::string simLocation;
  if (runNumber<16258) simLocation = string(getenv("SIM_2011_2012"));
  if (runNumber<20000) simLocation = string(getenv("SIM_2011_2012"));
  else if (runNumber>21087 && runNumber<21679) simLocation = string(getenv("SIM_2012_2013_ISOBUTANE"));
  else simLocation = string(getenv("SIM_2012_2013"));

  std::cout << "Using simulation from " << simLocation << "...\n";

  //Create simulation output file
  Char_t outputfile[500];
  sprintf(outputfile,"analyzed_files/SimAnalyzed_%s_%u.root",source.c_str(),linParamsFile);
  TFile *outfile = new TFile(outputfile, "RECREATE");


  // Load information for applying detector effects
  std::vector <Double_t> alphas = getAlphaValues(alphaFileIndex); // fill vector with the alpha (nPE/keV) values
  std::vector < std::vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7); // 2D vector with trigger function for East side and West side in that order
  GetPositionMap(XeMapPeriod); //Loads position map via posMapReader.h methods

  //Setup the output tree
  TTree *tree = new TTree("SimAnalyzed", "SimAnalyzed");
  SetUpOuputTree(tree); //Setup the output tree and branches

  //Histograms of event types for quick checks
  std::vector <TH1D*> finalEn (6,NULL);
  finalEn[0] = new TH1D("finalE0", "Simulated Weighted Sum East Type 0", 400, 0., 1200.);
  finalEn[1] = new TH1D("finalW0", "Simulated Weighted Sum West Type 0", 400, 0., 1200.);
  finalEn[2] = new TH1D("finalE1", "Simulated Weighted Sum East Type 1", 400, 0., 1200.);
  finalEn[3] = new TH1D("finalW1", "Simulated Weighted Sum West Type 1", 400, 0., 1200.);
  finalEn[4] = new TH1D("finalE23", "Simulated Weighted Sum East Type 2/3", 400, 0., 1200.);
  finalEn[5] = new TH1D("finalW23", "Simulated Weighted Sum West Type 2/3", 400, 0., 1200.);

  

  //Read in simulated data and put in a TChain
  int numFiles = source=="Beta"?400:250; //Sets the number of analyzed files
  TChain *chain = new TChain("anaTree");
  for (int i=0; i<numFiles; i++) {
    sprintf(temp,"%s/%s/analyzed_%i.root",simLocation.c_str(),source.c_str(),i);
    chain->AddFile(temp);
  }
  
  // Set the addresses of the information read in from the simulation file
  chain->SetBranchAddress("MWPCEnergy",&mwpcE);
  chain->SetBranchAddress("time",&Time);
  chain->SetBranchAddress("Edep",&edep);
  chain->SetBranchAddress("EdepQ",&edepQ);
  chain->SetBranchAddress("MWPCPos",&mwpc_pos);
  chain->SetBranchAddress("ScintPos",&scint_pos);
  chain->SetBranchAddress("primKE",&Eprim);
  
  //Trigger booleans
  bool EastScintTrigger, WestScintTrigger, EMWPCTrigger, WMWPCTrigger;
  Double_t MWPCThreshold=0.001; // generic MWPC Threshold

  //Set random number generators for sampling of gaussians in smearing
  TRandom3 *rand = new TRandom3(0);
  TRandom3 *rand2 = new TRandom3(0);
  
  //Get total number of events in TChain
  UInt_t nevents = chain->GetEntries();
  cout << "events in chain = " << nevents << endl;
 
  //Start from random position in evt sequence
  UInt_t evtStart = rand->Rndm()*nevents;
  UInt_t evtTally = 0; //To keep track of the number of events 
  UInt_t evt = evtStart; //current event number

  std::vector < std::vector <Int_t> > gridPoint; // Will hold the grid location for the position map

  //Read in events and determine evt type based on triggers
  while (evtTally<=BetaEvents) {
    if (evt>=nevents) evt=0; //Wrapping the events back to the beginning
    EastScintTrigger = WestScintTrigger = EMWPCTrigger = WMWPCTrigger = false; //Resetting triggers each event

    chain->GetEvent(evt);

    //Checking that the event occurs within the fiducial volume in the simulation to minimize
    // contamination from edge effects and interactions with detector walls
    Double_t fidCut = 45.;
    if (source=="Beta") fidCut=50.; //Nominal fiducial cut for betas
    else if (sqrt(scint_pos.ScintPosE[0]*scint_pos.ScintPosE[0]+scint_pos.ScintPosE[1]+scint_pos.ScintPosE[1])*sqrt(0.6)*10.>fidCut
	     || sqrt(scint_pos.ScintPosW[0]*scint_pos.ScintPosW[0]+scint_pos.ScintPosW[1]+scint_pos.ScintPosW[1])*sqrt(0.6)*10.>fidCut) {evt++; continue;} //For source events, 
    // We don't want edge contamination, and I use a cut of 45 mm when selecting what sources are present in calibration runs.

    //Setting adjusted positions when we are simulating a source
    scint_pos_adj.ScintPosAdjE[0] = source=="Beta"?scint_pos.ScintPosE[0]*sqrt(0.6)*10.:rand2->Gaus(srcPos[0][0], fabs(srcPos[0][2]));
    scint_pos_adj.ScintPosAdjE[1] = source=="Beta"?scint_pos.ScintPosE[1]*sqrt(0.6)*10.:rand2->Gaus(srcPos[0][1], fabs(srcPos[0][2]));
    scint_pos_adj.ScintPosAdjW[0] = source=="Beta"?scint_pos.ScintPosW[0]*sqrt(0.6)*10.:rand2->Gaus(srcPos[1][0], fabs(srcPos[1][2]));//sqrt(0.6)*10.*scint_pos.ScintPosW[0]+displacementX;
    scint_pos_adj.ScintPosAdjW[1] = source=="Beta"?scint_pos.ScintPosW[1]*sqrt(0.6)*10.:rand2->Gaus(srcPos[1][1], fabs(srcPos[1][2]));//sqrt(0.6)*10.*scint_pos.ScintPosW[1]+displacementY;
    scint_pos_adj.ScintPosAdjE[2] = scint_pos.ScintPosE[2]*10.;
    scint_pos_adj.ScintPosAdjW[2] = scint_pos.ScintPosW[2]*10.;
    
      
    //retrieve point on grid for each side of detector [E/W][x/y]
    gridPoint = getGridPoint(scint_pos_adj.ScintPosAdjE[0],scint_pos_adj.ScintPosAdjE[1],scint_pos_adj.ScintPosAdjW[0],scint_pos_adj.ScintPosAdjW[1]);

    Int_t intEastBinX = gridPoint[0][0];
    Int_t intEastBinY = gridPoint[0][1];
    Int_t intWestBinX = gridPoint[1][0];
    Int_t intWestBinY = gridPoint[1][1];
      
    //MWPC triggers
    if (mwpcE.MWPCEnergyE>MWPCThreshold) EMWPCTrigger=true;
    if (mwpcE.MWPCEnergyW>MWPCThreshold) WMWPCTrigger=true;

    Double_t pmtEnergyLowerLimit = 1.; //To put a hard cut on the weight
    

    //East Side smeared PMT energies
    for (UInt_t p=0; p<4; p++) {
      //cout << p << " " << positionMap[p][intEastBinX][intEastBinY] << endl;
      Double_t posCorrAlpha = alphas[p]*positionMap[p][intEastBinX][intEastBinY];
      if (posCorrAlpha==0.) posCorrAlpha=1.; //This occurs when the event occurs outside the position map... This is to prevent the Gaussian from having inf sigma
      //Sample a gaussian centered on the EQ of the event with a spread related to alpha
      pmt_Evis.Evis[p] = rand->Gaus(edepQ.EdepQE,sqrt(edepQ.EdepQE/posCorrAlpha));
     
      if (pmt_Evis.Evis[p]>pmtEnergyLowerLimit) {
	pmt_Evis.weight[p] = posCorrAlpha/pmt_Evis.Evis[p];
      }
      else {pmt_Evis.weight[p]=posCorrAlpha/pmtEnergyLowerLimit;} //This sets a hard cut on the weight of an event so that events with zero energy still carry weight.. Probably not the most effective way to do this, but it should address low energy behavior 
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
    if (rand->Rndm(0)<triggProb && evis.EvisE>0.) EastScintTrigger=true;
   
      
    //West Side
    for (UInt_t p=4; p<8; p++) {
      Double_t posCorrAlpha = alphas[p]*positionMap[p][intWestBinX][intWestBinY];
      if (posCorrAlpha==0.) posCorrAlpha=1.; //This occurs when the event occurs outside the position map... This is to prevent the Gaussian from having inf sigma
      pmt_Evis.Evis[p] = rand->Gaus(edepQ.EdepQW,sqrt(edepQ.EdepQW/posCorrAlpha));
      
      if (pmt_Evis.Evis[p]>pmtEnergyLowerLimit) {
	pmt_Evis.weight[p] = posCorrAlpha/pmt_Evis.Evis[p];
      }
      else {pmt_Evis.weight[p]=posCorrAlpha/pmtEnergyLowerLimit;} //This sets a hard cut on the weight of an event so that events with zero energy still carry weight.. Probably not the most effective way to do this, but it should address low energy behavior
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
    if (rand->Rndm(0)<triggProb && evis.EvisW>0.) WestScintTrigger = true;     

      
    //Fill proper total event histogram based on event type
    PID=-1;
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
    //Type 2/3 East
    else if (EastScintTrigger && EMWPCTrigger && !WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=2;
      side=0;
      finalEn[4]->Fill(totalEnE);
      //cout << "Type 2/3 East E = " << totalEnE << endl;
    }
    //Type 2/3 West
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
	  else {
	    side=1;
	  }
	}
	else {
	  PID=6;
	  type=4;
	  side=2;
	}
      }
      else {
	PID=1;
	type=4;
	side = (WMWPCTrigger && EMWPCTrigger) ? 2 : (WMWPCTrigger ? 1 : 0);
      }
    }

    //Calculate Erecon
    Double_t totalEvis = 0.; 
    if (side==0) {
	totalEvis = type==1 ? (evis.EvisE+evis.EvisW):evis.EvisE;
	if (totalEvis>0.) {
	  Erecon = EQ2Etrue[0][typeIndex][0]+EQ2Etrue[0][typeIndex][1]*totalEvis+EQ2Etrue[0][typeIndex][2]/(totalEvis+EQ2Etrue[0][typeIndex][3])+EQ2Etrue[0][typeIndex][4]/((totalEvis+EQ2Etrue[0][typeIndex][5])*(totalEvis+EQ2Etrue[0][typeIndex][5]));
	}
	else Erecon=-1.;
      }
      if (side==1) {
	totalEvis = type==1 ? (evis.EvisE+evis.EvisW):EvisW;
	if (totalEvis>0.) {
	  Erecon = EQ2Etrue[1][typeIndex][0]+EQ2Etrue[1][typeIndex][1]*totalEvis+EQ2Etrue[1][typeIndex][2]/(totalEvis+EQ2Etrue[1][typeIndex][3])+EQ2Etrue[1][typeIndex][4]/((totalEvis+EQ2Etrue[1][typeIndex][5])*(totalEvis+EQ2Etrue[1][typeIndex][5]));
	}
	else Erecon=-1.;
      }
    
      
    // Increment the event tally if the event was PID = 1 (electron) and the event was inside the fiducial radius used to determine num of events in data file
    if (PID==1 && Erecon>0. && sqrt(scint_pos.ScintPosE[0]*scint_pos.ScintPosE[0]+scint_pos.ScintPosE[1]+scint_pos.ScintPosE[1])*sqrt(0.6)*10.<fidCut && sqrt(scint_pos.ScintPosW[0]*scint_pos.ScintPosW[0]+scint_pos.ScintPosW[1]+scint_pos.ScintPosW[1])*sqrt(0.6)*10.<fidCut)
      evtTally++;

    evt++;
    if (PID>=0 && Erecon>0.) tree->Fill();
    //cout << evtTally << endl;
    if (evtTally%1000==0) {cout << "*";}//cout << "filled event " << evt << endl;
  }
  cout << endl;
  delete chain;
  delete rand;
  delete rand2;
  outfile->Write();
  outfile->Close();
  
}

int main(int argc, char *argv[]) {
  if (argc!=4) {
    std::cout << "Usage: ./SimulationAnalyzer [Source] [numEvents] [parameterFile]\n";
    exit(0);
  }
  std::string source = string(argv[1]);
  UInt_t numEvents = (unsigned int)atoi(argv[2]);
  UInt_t paramFile = (unsigned int)atoi(argv[3]);
  revCalSimulation(source, numEvents,paramFile);

  //tests
  /*UInt_t XePeriod = getXeRunPeriod(atoi(argv[1]));
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7);
  Double_t triggProbE = triggerProbability(triggerFunc[0],25.);
  Double_t triggProbW = triggerProbability(triggerFunc[1],25.);
  cout << triggProbE << " " << triggProbW << endl;*/


}
  
