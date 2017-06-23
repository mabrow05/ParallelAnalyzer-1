 /* Code to take a run number, retrieve it's runperiod, and construct the 
weighted spectra which would be seen as a reconstructed energy on one side
of the detector. Also applies the trigger functions */

#include "revCalSim.h"
#include "posMapReader.h"
#include "positionMapHandler.hh"
#include "sourcePeaks.h"
#include "runInfo.h"
#include "calibrationTools.hh"
#include "TriggerMap.hh"
//#include "../Asymmetry/SQLinterface.hh"
#include "MWPCPositionResponse.hh"

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
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>


using namespace std;

std::string getRunTypeFromOctetFile(int octet, int run) {
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  int numRuns = 0;
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {
    if ( runNumberHold == run ) { infile.close(); return runTypeHold; }
  }
  return "BAD";
  
};

int getPolarization(int run) {

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
  
  else {
    std::cout << "You chose a Depol Run.. no polarization\n";
    return 0;
  }
};



vector <vector <double> > returnSourcePosition (Int_t runNumber, string src) {
  Char_t temp[500];
  sprintf(temp,"%s/source_list_%i.dat",getenv("SOURCE_LIST"),runNumber);
  ifstream file(temp);
  cout << src << endl;
  int num = 0;
  file >> num;
  cout << num << endl;
  int srcNum = 0;
  string src_check;
  for (int i=0; i<num;srcNum++,i++) {
    file >> src_check;
    cout << src_check << endl;
    if (src_check==src) break;   
  }
  cout << "The source Number is: " << srcNum << endl;
  if (srcNum==num) {
    cout << "Didn't find source in that run\n"; exit(0);
  }
  file.close();
  
  sprintf(temp,"%s/source_positions_%i.dat",getenv("SOURCE_POSITIONS"),runNumber);
  file.open(temp);
  
  vector < vector < double > > srcPos;
  srcPos.resize(2,vector <double> (3,0.));
  
  for (int i=0; i<srcNum+1; i++) {
    for (int j=0; j<2; j++) {
      for (int jj=0; jj<3; jj++) {
	file >> srcPos[j][jj];
      }
    }
  }
  return srcPos;
}

vector <Int_t> getEreconPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"%s/residuals/PMT_EreconQuality_master.dat",getenv("ANALYSIS_CODE")); 
  ifstream pmt;
  std::cout << temp << std::endl;
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
  sprintf(temp,"%s/simulation_comparison/nPE_per_keV/nPE_per_keV_%i.dat",getenv("ANALYSIS_CODE"),runPeriod);
  ifstream infile;
  infile.open(temp);
  Int_t i = 0;

  while (infile >> alphas[i]) { std::cout << alphas[i] << std::endl; i++; }
  return alphas;
}

  

std::vector < std::vector <Double_t> > loadPMTpedestals(Int_t runNumber) {

  Char_t temp[500];
  std::vector < std::vector < Double_t > > peds (8,std::vector <Double_t> (2,0.));
  sprintf(temp,"%s/PMT_pedestals_%i.dat",getenv("PEDESTALS"),runNumber);
  ifstream infile;
  infile.open(temp);

  Int_t i = 0;
  Int_t run;

  while (infile >> run >> peds[i][0] >> peds[i][1]) { std::cout << "Pedestal " << i << ": " << peds[i][0] << " " << peds[i][1] << std::endl; i++; }
  return peds;

};

std::vector <Double_t> loadGainFactors(Int_t runNumber) {

  // Read gain corrections file                                                                                                                
  char tempFileGain[500];
  sprintf(tempFileGain, "%s/gain_bismuth_%i.dat",getenv("GAIN_BISMUTH"), runNumber);
  std::cout << "... Reading: " << tempFileGain << std::endl;

  std::vector <Double_t> gainCorrection(8,1.);
  std::vector <Double_t> fitMean(8,1.);
  ifstream fileGain(tempFileGain);
  for (int i=0; i<8; i++) {
    fileGain >> fitMean[i] >> gainCorrection[i];
  }
  std::cout << "...   PMT E1: " << gainCorrection[0] << std::endl;
  std::cout << "...   PMT E2: " << gainCorrection[1] << std::endl;
  std::cout << "...   PMT E3: " << gainCorrection[2] << std::endl;
  std::cout << "...   PMT E4: " << gainCorrection[3] << std::endl;
  std::cout << "...   PMT W1: " << gainCorrection[4] << std::endl;
  std::cout << "...   PMT W2: " << gainCorrection[5] << std::endl;
  std::cout << "...   PMT W3: " << gainCorrection[6] << std::endl;
  std::cout << "...   PMT W4: " << gainCorrection[7] << std::endl;

  return gainCorrection;
};


void SetUpTree(TTree *tree) {
  tree->Branch("PrimaryID", &primaryID, "primaryID/I");
  tree->Branch("PID", &PID, "PID/I");
  tree->Branch("side", &side, "side/I");
  tree->Branch("type", &type, "type/I");
  tree->Branch("Erecon", &Erecon,"Erecon/D");
  //  tree->Branch("old_Erecon", &old_Erecon,"old_Erecon/D");
  //tree->Branch("gaus_Erecon", &gaus_Erecon,"gaus_Erecon/D");

  tree->Branch("Evis",&evis,"EvisE/D:EvisW");
  tree->Branch("Edep",&edep,"EdepE/D:EdepW");
  tree->Branch("EdepQ",&edepQ,"EdepQE/D:EdepQW");
  tree->Branch("primaryKE",&primKE,"primKE/D");
  tree->Branch("primaryMomentum",&primMom,"primTheta/D");
  //tree->Branch("AsymWeight",&AsymWeight,"AsymWeight/D");
  //tree->Branch("hitCountSD",&hitCountSD,"hitCountSD[24]/I");
  
  tree->Branch("time",&Time,"timeE/D:timeW");
  tree->Branch("MWPCEnergy",&mwpcE,"MWPCEnergyE/D:MWPCEnergyW");
  tree->Branch("MWPCPos",&mwpc_pos,"MWPCPosE[3]/D:MWPCPosW[3]");
  //tree->Branch("Cath_EX",Cath_EX,"Cath_EX[16]/F");
  //tree->Branch("Cath_EY",Cath_EY,"Cath_EY[16]/F");
  //tree->Branch("Cath_WX",Cath_WX,"Cath_WX[16]/F");
  //tree->Branch("Cath_WY",Cath_WY,"Cath_WY[16]/F");
  //tree->Branch("nClipped_EX",&nClipped_EX,"nClipped_EX/I");
  //tree->Branch("nClipped_EY",&nClipped_EY,"nClipped_EY/I");
  //tree->Branch("nClipped_WX",&nClipped_WX,"nClipped_WX/I");
  //tree->Branch("nClipped_WY",&nClipped_WY,"nClipped_WY/I");
  tree->Branch("ScintPos",&scint_pos,"ScintPosE[3]/D:ScintPosW[3]");
  tree->Branch("ScintPosAdjusted",&scint_pos_adj,"ScintPosAdjE[3]/D:ScintPosAdjW[3]");
  //tree->Branch("cathRespPos",&cathResp_pos,"cathRespPosE[3]/D:cathRespPosW[3]");
  tree->Branch("PMT",&pmt,"Evis0/D:Evis1:Evis2:Evis3:Evis4:Evis5:Evis6:Evis7:etaEvis0/D:etaEvis1:etaEvis2:etaEvis3:etaEvis4:etaEvis5:etaEvis6:etaEvis7:nPE0/D:nPE1:nPE2:nPE3:nPE4:nPE5:nPE6:nPE7");

  //tree->Branch("xE",&xE,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  //tree->Branch("yE",&yE,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  //tree->Branch("xW",&xW,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  //tree->Branch("yW",&yW,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  
  //tree->Branch("old_xE",&old_xE,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  //tree->Branch("old_yE",&old_yE,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  //tree->Branch("old_xW",&old_xW,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  //tree->Branch("old_yW",&old_yW,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");

  //tree->Branch("gaus_xE",&gaus_xE,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  //tree->Branch("gaus_yE",&gaus_yE,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  //tree->Branch("gaus_xW",&gaus_xW,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  //tree->Branch("gaus_yW",&gaus_yW,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  
}
  

void revCalSimulation (Int_t runNumber, string b, string A) 
{
  bool allEvtsTrigg = false; //This just removes the use of the trigger function for an initial calibration. 
                            // Once a calibration is established (or if running on Betas), you can keep this false

  bool simProperStatistics = false; // If True, this uses the actual data run to determine the number of Type 0s to simulate

  
  cout << "Running reverse calibration for run " << runNumber << endl;  

  string outputBase;
  UInt_t BetaEvents = 0; //Holds the number of electron like Type 0 events to process

  Char_t temp[500],tempE[500],tempW[500];
  
  //First get the number of total electron events around the source position from the data file and also the length of the run
  
  //sprintf(temp,"%s/replay_pass4_%i.root",getenv("REPLAY_PASS4"),runNumber);
  sprintf(temp,"%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),runNumber);
  TFile *dataFile = new TFile(temp,"READ");
  TTree *data = (TTree*)(dataFile->Get("pass3"));
  
  sprintf(tempE,"Type<4 && Type>=0 && PID==1 && Side==0 && (xE.center*xE.center+yE.center*yE.center)<2500. && (xW.center*xW.center+yW.center*yW.center)<2500.");
  sprintf(tempW,"Type<4 && Type>=0 && PID==1 && Side==1 && (xW.center*xW.center+yW.center*yW.center)<2500. && (xE.center*xE.center+yE.center*yE.center)<2500.");
 
  BetaEvents = data->GetEntries(tempE) + data->GetEntries(tempW);
  cout << "Electron Events in Data file: " << BetaEvents << endl;
  cout << "East: " << data->GetEntries(tempE) << "\tWest: " << data->GetEntries(tempW) << endl;
  delete data;
  dataFile->Close(); 
  
  
  //If we have a Beta run and we don't want exact statistics, simulate 16 times the number of events
  if ( !simProperStatistics  ) {
    BetaEvents = 16*BetaEvents;
    outputBase = string(getenv("REVCALSIM")) + "beta_highStatistics/";
  } 

  std::cout << "Processing " << BetaEvents << " events...\n";



  ///////////////////////// SETTING GAIN OF FIRST and second DYNYODE
  Double_t g_d = 16.;
  Double_t g_rest = 12500.;

  /////// Loading other run dependent quantities
  vector <Int_t> pmtQuality = getEreconPMTQuality(runNumber); // Get the quality of the PMTs for that run
  UInt_t calibrationPeriod = getSrcRunPeriod(runNumber); // retrieve the calibration period for this run
  UInt_t XePeriod = getXeRunPeriod(runNumber); // Get the proper Xe run period for the Trigger functions
  //GetPositionMap(XePeriod);
  PositionMap posmap(5.0,50.); //Load position map with 5 mm bins
  posmap.readPositionMap(XePeriod,"endpoint");
  vector <Double_t> alpha = GetAlphaValues(calibrationPeriod); // fill vector with the alpha (nPE/keV) values for this run period


  /////////////// Load trigger functions //////////////////
  TriggerFunctions trigger(runNumber);

  std::vector < std::vector <Double_t> > pedestals = loadPMTpedestals(runNumber);
  std::vector < Double_t > PMTgain = loadGainFactors(runNumber);

  //Load the simulated relationship between EQ and Etrue
  EreconParameterization eRecon(runNumber);
  // Int_t pol = source=="Beta" ? getPolarization(runNumber) : 1;

  LinearityCurve linCurve(calibrationPeriod,false);
  std::cout << "Loaded Linearity Curves\n";

  //Decide which simulation to use...
  TChain *chain = new TChain("anaTree"); 
  
  TString fileLocation = TString::Format("/extern/mabrow05/ucna/xuan_stuff/A_%s_b_%s/",A.c_str(),b.c_str());
 
  std::cout << "Using simulation from " << fileLocation << "...\n";

  //Read in simulated data and put in a TChain
  TRandom3 *randFile = new TRandom3(runNumber*2);
  int numFiles = 10;
  int fileNum = (int)(randFile->Rndm()*numFiles);
  delete randFile;
  for (int i=0; i<numFiles; i++) {
    if (fileNum==numFiles) fileNum=0;
    chain->AddFile(TString::Format("%s/analyzed_%i.root",fileLocation.Data(),fileNum));
    fileNum++;
  }

  // Set the addresses of the information read in from the simulation file
  /*chain->SetBranchAddress("MWPCEnergy",&mwpcE);       x
  chain->SetBranchAddress("time",&Time);                x
  chain->SetBranchAddress("Edep",&edep);                x
  chain->SetBranchAddress("EdepQ",&edepQ);              x
  chain->SetBranchAddress("MWPCPos",&mwpc_pos);         x 
  chain->SetBranchAddress("ScintPos",&scint_pos);       x
  chain->SetBranchAddress("primKE",&primKE);            x
  chain->SetBranchAddress("primTheta",&primTheta);
  chain->SetBranchAddress("Cath_EX",Cath_EX);
  chain->SetBranchAddress("Cath_EY",Cath_EY);
  chain->SetBranchAddress("Cath_WX",Cath_WX);
  chain->SetBranchAddress("Cath_WY",Cath_WY);
  chain->SetBranchAddress("hitCountSD",hitCountSD);*/
  
  //These are for feeding in Xuan's simulations... this needs to be updated so that I can pass a flag and change these on the fly
  chain->SetBranchAddress("PrimaryParticleSpecies",&primaryID);
  chain->SetBranchAddress("mwpcEnergy",&mwpcE);
  chain->SetBranchAddress("scintTimeToHit",&Time);
  chain->SetBranchAddress("scintillatorEdep",&edep);
  chain->SetBranchAddress("scintillatorEdepQuenched",&edepQ);
  chain->SetBranchAddress("MWPCPos",&mwpc_pos);
  chain->SetBranchAddress("ScintPos",&scint_pos);
  chain->SetBranchAddress("primaryKE",&primKE);
  chain->SetBranchAddress("primaryMomentum",primMom);


  //Set random number generator

  TRandom3 *seed = new TRandom3( 3*runNumber ); // seed generator
  TRandom3 *rand0 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
  TRandom3 *rand1 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
  TRandom3 *rand2 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
  TRandom3 *rand3 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );

  //Initialize random numbers for 8 pmt trigger probabilities
  TRandom3 *randPMTE1 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
  TRandom3 *randPMTE2 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
  TRandom3 *randPMTE3 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
  TRandom3 *randPMTE4 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
  TRandom3 *randPMTW1 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
  TRandom3 *randPMTW2 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
  TRandom3 *randPMTW3 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );
  TRandom3 *randPMTW4 = new TRandom3( (unsigned int) (seed->Rndm()*10000.) );

  std::vector <Double_t> triggRandVec(4,0.);

  // Wirechamber information
  bool EastScintTrigger, WestScintTrigger, EMWPCTrigger, WMWPCTrigger; //Trigger booleans
  Double_t MWPCAnodeThreshold=0.; // keV dep in the wirechamber.. 


  //Get total number of events in TChain
  UInt_t nevents = chain->GetEntries();
  cout << "events = " << nevents << endl;
 
  //Start from random position in evt sequence if we aren't simulating veryHighStatistics
  UInt_t evtStart = seed->Rndm()*nevents;
    
  
  UInt_t evtTally = 0; //To keep track of the number of events 
  UInt_t evt = evtStart; //current event number


  vector < vector <Int_t> > gridPoint;

  
  //Create simulation output file
  
  TFile *outfile = new TFile(TString::Format("%s/revCalSim_%i.root",fileLocation.Data(),runNumber), "RECREATE");

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


  //Read in events and determine evt type based on triggers
  while (evtTally<=BetaEvents) {
    if (evt>=nevents) evt=0; //Wrapping the events back to the beginning
    EastScintTrigger = WestScintTrigger = EMWPCTrigger = WMWPCTrigger = false; //Resetting triggers each event

    chain->GetEvent(evt);

    //Checking that the event occurs within the fiducial volume in the simulation to minimize
    // contamination from edge effects and interactions with detector walls
    // This is also the fiducial cut used in extracting asymmetries and doing calibrations
    Double_t fidCut = 50.;


    /////////////// Do position and Wirechamber stuff ///////////////////
    
    // Following negative sign is to turn x-coordinate into global coordinate on East
    scint_pos_adj.ScintPosAdjE[0] = scint_pos.ScintPosE[0]*sqrt(0.6)*1000.;
    scint_pos_adj.ScintPosAdjE[1] = scint_pos.ScintPosE[1]*sqrt(0.6)*1000.;
    scint_pos_adj.ScintPosAdjW[0] = scint_pos.ScintPosW[0]*sqrt(0.6)*1000.;//sqrt(0.6)*10.*scint_pos.ScintPosW[0]+displacementX;
    scint_pos_adj.ScintPosAdjW[1] = scint_pos.ScintPosW[1]*sqrt(0.6)*1000.;//sqrt(0.6)*10.*scint_pos.ScintPosW[1]+displacementY;
    scint_pos_adj.ScintPosAdjE[2] = scint_pos.ScintPosE[2]*1000.;
    scint_pos_adj.ScintPosAdjW[2] = scint_pos.ScintPosW[2]*1000.;

    
    Double_t mwpcAdjE[3] = {0.,0.,0.};
    Double_t mwpcAdjW[3] = {0.,0.,0.};
   
    mwpcAdjE[0] = mwpc_pos.MWPCPosE[0] * sqrt(0.6) * 1000.; 
    mwpcAdjE[1] = mwpc_pos.MWPCPosE[1] * sqrt(0.6) * 1000.; 
    mwpcAdjW[0] = mwpc_pos.MWPCPosW[0] * sqrt(0.6) * 1000. ;
    mwpcAdjW[1] = mwpc_pos.MWPCPosW[1] * sqrt(0.6) * 1000. ;
    mwpcAdjE[2] = mwpc_pos.MWPCPosE[2] * 1000. ;
    mwpcAdjW[2] = mwpc_pos.MWPCPosW[2] * 1000. ;
    

    
    //std::cout << mwpcAdjE[0] << "\t" << mwpcAdjW[0] << std::endl;
    //if (evtTally==5) exit(0);


    /////////////////////////////////////////////////////////////////////
    // Creating all extra structs to hold the alternate position reconstruction
    // crap before writing it all out
    /////////////////////////////////////////////////////////////////////
    

    std::vector <Double_t> eta; 
    std::vector <Double_t> trueEta; // This is the position map value of the actual event position, not the position as determined
                                    // by cathode reconstruction
    

    trueEta = posmap.getInterpolatedEta(scint_pos_adj.ScintPosAdjE[0],
					scint_pos_adj.ScintPosAdjE[1],
					scint_pos_adj.ScintPosAdjW[0],
					scint_pos_adj.ScintPosAdjW[1]);
   
    eta = posmap.getInterpolatedEta(scint_pos_adj.ScintPosAdjE[0],
				    scint_pos_adj.ScintPosAdjE[1],
				    scint_pos_adj.ScintPosAdjW[0],
				    scint_pos_adj.ScintPosAdjW[1]);

      
      
    //MWPC triggers
    if (mwpcE.MWPCEnergyE>MWPCAnodeThreshold) EMWPCTrigger=true;
    if (mwpcE.MWPCEnergyW>MWPCAnodeThreshold) WMWPCTrigger=true;
    
    std::vector<Double_t> ADCvecE(4,0.);
    std::vector<Double_t> ADCvecW(4,0.);

    //East Side smeared PMT energies
    
    
    for (UInt_t p=0; p<4; p++) {
      if ( edep.EdepE>0. ) { //Check to make sure that there is light to see in the scintillator
	
	if (eta[p]>0.) {

	  pmt.etaEvis[p] = (1./(alpha[p]*g_d*g_rest)) * (rand3->Poisson(g_rest*rand2->Poisson(g_d*rand1->Poisson(alpha[p]*trueEta[p]*edepQ.EdepQE))));
	  Double_t ADC = linCurve.applyInverseLinCurve(p, pmt.etaEvis[p]) + rand0->Gaus(0.,pedestals[p][1]); //Take into account non-zero width of the pedestal
	  ADCvecE[p] = ADC;
	  pmt.etaEvis[p] = linCurve.applyLinCurve(p, ADC);
	  pmt.Evis[p] = pmt.etaEvis[p]/eta[p];

	}

	else { //To avoid dividing by zero.. these events won't be used in analysis since they are outside the fiducial cut
	  
	  pmt.etaEvis[p] = (1./(alpha[p]*g_d*g_rest)) * (rand3->Poisson(g_rest*rand2->Poisson(g_d*rand1->Poisson(alpha[p]*edepQ.EdepQE))));
	  Double_t ADC = linCurve.applyInverseLinCurve(p, pmt.etaEvis[p]);// + rand0->Gaus(0.,pedestals[p][1]); //Take into account non-zero width of the pedestal
	  ADCvecE[p] = ADC;
	  pmt.etaEvis[p] = linCurve.applyLinCurve(p, ADC);
	  pmt.Evis[p] = pmt.etaEvis[p];

	}


	
	pmt.nPE[p] = alpha[p]*pmt.etaEvis[p] ;

      }

    // If eQ is 0...
      else {
	pmt.etaEvis[p] = 0.;
	pmt.Evis[p] = 0.;
	pmt.nPE[p] = 0.;
	
      }
    }
  
    Double_t numer=0., denom=0.;
    
    for (UInt_t p=0;p<4;p++) {
      numer += ( pmtQuality[p] ) ? pmt.nPE[p] : 0.;
      denom += ( pmtQuality[p] ) ? eta[p]*alpha[p] : 0.;
      //numer += ( pmtQuality[p] && pmt.etaEvis[p]>0. ) ? pmt.nPE[p] : 0.;
      //denom += ( pmtQuality[p] && pmt.etaEvis[p]>0. ) ? eta[p]*alpha[p] : 0.;
    }

    
    Double_t totalEnE = denom>0. ? numer/denom : 0.;
    evis.EvisE = totalEnE;

    //Now we apply the trigger probability
    triggRandVec[0] = randPMTE1->Rndm();
    triggRandVec[1] = randPMTE2->Rndm();
    triggRandVec[2] = randPMTE3->Rndm();
    triggRandVec[3] = randPMTE4->Rndm();
    
    if ( (allEvtsTrigg && totalEnE>0.) || ( trigger.decideEastTrigger(ADCvecE,triggRandVec) ) ) EastScintTrigger = true; 
      
    //West Side
    for (UInt_t p=4; p<8; p++) {
      if ( !(p==5 && runNumber>16983 && runNumber<17249) &&  edep.EdepW>0. ) { //Check to make sure that there is light to see in the scintillator and that run isn't one where PMTW2 was dead
	
	if (eta[p]>0.) {

	  pmt.etaEvis[p] = (1./(alpha[p]*g_d*g_rest)) * (rand3->Poisson(g_rest*rand2->Poisson(g_d*rand1->Poisson(alpha[p]*trueEta[p]*edepQ.EdepQW))));
	  Double_t ADC = linCurve.applyInverseLinCurve(p, pmt.etaEvis[p]) + rand0->Gaus(0.,pedestals[p][1]); //Take into account non-zero width of the pedestal
	  ADCvecW[p-4] = ADC;
	  //cout << ADCvecW[p-4] << endl;
	  pmt.etaEvis[p] = linCurve.applyLinCurve(p, ADC);	  
	  pmt.Evis[p] = pmt.etaEvis[p]/eta[p];

	}

	else { //To avoid dividing by zero.. these events won't be used in analysis since they are outside the fiducial cut

	  pmt.etaEvis[p] = (1./(alpha[p]*g_d*g_rest)) * (rand3->Poisson(g_rest*rand2->Poisson(g_d*rand1->Poisson(alpha[p]*edepQ.EdepQW))));
	  Double_t ADC = linCurve.applyInverseLinCurve(p, pmt.etaEvis[p]);// + rand0->Gaus(0.,pedestals[p][1]); //Take into account non-zero width of the pedestal
	  ADCvecW[p-4] = ADC;
	  pmt.etaEvis[p] = linCurve.applyLinCurve(p, ADC);
	  pmt.Evis[p] = pmt.etaEvis[p];
	  
	}

	pmt.nPE[p] = alpha[p]*pmt.etaEvis[p];
	
      }
      // If PMT is dead and EQ=0...
      else {
	pmt.etaEvis[p] = 0.;
	pmt.Evis[p] = 0.;
	pmt.nPE[p] = 0.;
      }
    }

     //Calculate the total weighted energy
    numer=denom=0.;
    for (UInt_t p=4;p<8;p++) {
      numer += ( pmtQuality[p] ) ? pmt.nPE[p] : 0.;
      denom += ( pmtQuality[p] ) ? eta[p]*alpha[p] : 0.;
      //numer += ( pmtQuality[p] && pmt.etaEvis[p]>0. ) ? pmt.nPE[p] : 0.;
      //denom += ( pmtQuality[p] && pmt.etaEvis[p]>0. ) ? eta[p]*alpha[p] : 0.;
    }
    //Now we apply the trigger probability
    Double_t totalEnW = denom>0. ? numer/denom : 0.;
    evis.EvisW = totalEnW;

    triggRandVec[0] = randPMTW1->Rndm();
    triggRandVec[1] = randPMTW2->Rndm();
    triggRandVec[2] = randPMTW3->Rndm();
    triggRandVec[3] = randPMTW4->Rndm();
    
    if ( (allEvtsTrigg && totalEnW>0.) || ( trigger.decideWestTrigger(ADCvecW,triggRandVec) ) ) WestScintTrigger = true;

    //if (totalEnW<25.) std::cout << totalEnW << " " << WestScintTrigger << "\n";

            
    //Fill proper total event histogram based on event type
    PID=6; //This is an unidentified particle
    type=4; //this is a lost event
    side=2; //This means there are no scintillator triggers

    //Type 0 East
    if (EastScintTrigger && EMWPCTrigger && !WestScintTrigger && !WMWPCTrigger) {
      PID=1;
      type=0;
      side=0;
      //finalEn[0]->Fill(evis.EvisE);
      //cout << "Type 0 East E = " << totalEnE << endl;
    }
    //Type 0 West
    else if (WestScintTrigger && WMWPCTrigger && !EastScintTrigger && !EMWPCTrigger) {
      PID=1;
      type=0;
      side=1;
      //finalEn[1]->Fill(totalEnW);
    }
    //Type 1 
    else if (EastScintTrigger && EMWPCTrigger && WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=1;
      //East
      if (Time.timeE<Time.timeW) {
	//finalEn[2]->Fill(totalEnE);
	side=0;
      }
      //West
      else if (Time.timeE>Time.timeW) {
	//finalEn[3]->Fill(totalEnW);
	side=1;
      }
    }
    //Type 2/3 East
    else if (EastScintTrigger && EMWPCTrigger && !WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=2;
      side=0;
      //finalEn[4]->Fill(totalEnE);
      //cout << "Type 2/3 East E = " << totalEnE << endl;
    }
    //Type 2/3 West
    else if (!EastScintTrigger && EMWPCTrigger && WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=2;
      side=1;
      //finalEn[5]->Fill(totalEnW);
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
	side = (WMWPCTrigger && EMWPCTrigger) ? 2 : (WMWPCTrigger ? 1 : 0); //Side==2 means the event triggered both wirechambers, but neither scint
      }
    }
    
    //Calculate Erecon
    Erecon = -1.;
    Int_t typeIndex = (type==0 || type==4) ? 0:(type==1 ? 1:2); //for retrieving the parameters from EQ2Etrue
    if (side==0) {
      Double_t totalEvis = type==1 ? (evis.EvisE+evis.EvisW):evis.EvisE;
      if (evis.EvisE>0. && totalEvis>0.) {
	Erecon = eRecon.getErecon(0,typeIndex,totalEvis);
	if (type==0) finalEn[0]->Fill(Erecon); 
	else if (type==1) finalEn[2]->Fill(Erecon); 
	else if(type==2 ||type==3) finalEn[4]->Fill(Erecon);
      }
      else Erecon=-1.;
      
    }
    if (side==1) {
      Double_t totalEvis = type==1 ? (evis.EvisE+evis.EvisW):evis.EvisW;
      if (evis.EvisW>0. && totalEvis>0.) {
	Erecon = eRecon.getErecon(1,typeIndex,totalEvis);	
	if (type==0) finalEn[1]->Fill(Erecon); 
	else if (type==1) finalEn[3]->Fill(Erecon); 
	else if(type==2 ||type==3) finalEn[5]->Fill(Erecon);
      }
      else Erecon=-1.;
    }    
  
    if ( PID==1 && Erecon>0. &&  sqrt( scint_pos_adj.ScintPosAdjE[0]*scint_pos_adj.ScintPosAdjE[0]
				       + scint_pos_adj.ScintPosAdjE[1]*scint_pos_adj.ScintPosAdjE[1] )<fidCut &&
	 sqrt( scint_pos_adj.ScintPosAdjW[0]*scint_pos_adj.ScintPosAdjW[0]
				       + scint_pos_adj.ScintPosAdjW[1]*scint_pos_adj.ScintPosAdjW[1] )<fidCut ) evtTally++;
    
    evt++;
    
    tree->Fill();
    if (evtTally%10000==0) {std::cout << evtTally  << "    PID=" << PID << "    Erecon=" << Erecon << std::endl;}//cout << "filled event " << evt << endl;
  }
  cout << endl;

  delete chain; delete seed; delete rand0; delete rand1; delete rand2; delete rand3;
  delete randPMTE1; delete randPMTE2; delete randPMTE3; delete randPMTE4; delete randPMTW1; delete randPMTW2; delete randPMTW3; delete randPMTW4;

  outfile->Write();
  outfile->Close();
  

}

int main(int argc, char *argv[]) {

  if (argc!=4) std::cout << "Usage: ./revCalSim.exe [runNumber] [b_val = 0,inf] [A_val=0,A0]\n";
  
  string b = string(argv[2]);
  string A = string(argv[3]);


  revCalSimulation(atoi(argv[1]),b,A);

  //tests
  /*UInt_t XePeriod = getXeRunPeriod(atoi(argv[1]));
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7);
  Double_t triggProbE = triggerProbability(triggerFunc[0],25.);
  Double_t triggProbW = triggerProbability(triggerFunc[1],25.);
  cout << triggProbE << " " << triggProbW << endl;*/


}
  
