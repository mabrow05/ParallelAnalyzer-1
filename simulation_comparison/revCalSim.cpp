/* Code to take a run number, retrieve it's runperiod, and construct the 
weighted spectra which would be seen as a reconstructed energy on one side
of the detector. Also applies the trigger functions */

#include "revCalSim.h"
#include "posMapReader.h"
#include "sourcePeaks.h"
#include "runInfo.h"

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

//Returns the average eta value from the position map for the Sn source used to calculate alpha
vector < Double_t > getMeanEtaForAlpha(Int_t run) 
{
  Char_t temp[500];
  vector < Double_t > eta (8,0.);
  sprintf(temp,"%s/nPE_meanEtaVal_%i.dat",getenv("NPE_WEIGHTS"),run);
  ifstream infile;
  infile.open(temp);
  Int_t i = 0;

  while (infile >> eta[i]) i++;
  return eta;
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
  tree->Branch("MWPCPos",&mwpc_pos,"MWPCPosE[3]/D:MWPCPosW[3]");
  tree->Branch("ScintPos",&scint_pos,"ScintPosE[3]/D:ScintPosW[3]");
  tree->Branch("ScintPosAdjusted",&scint_pos_adj,"ScintPosAdjE[3]/D:ScintPosAdjW[3]");
  tree->Branch("PMT_Evis",&pmt_Evis,"Evis0/D:Evis1:Evis2:Evis3:Evis4:Evis5:Evis6:Evis7:weight0:weight1:weight2:weight3:weight4:weight5:weight6:weight7");
  
}
  

void revCalSimulation (Int_t runNumber, string source) 
{
  cout << "Running reverse calibration for run " << runNumber << " and source " << source << endl;
  //Start by getting the source position if not a Beta decay run
  vector < vector <double> > srcPos;
  if (source!="Beta") {
    string srcShort = source;
    srcShort.erase(2);
    srcPos = returnSourcePosition(runNumber, srcShort);
  }
  //If the source is In114, checking which side the Indium was facing
  if (source=="In114") {
    string side = getIndiumSide(runNumber);
    if (side=="East") source+="E";
    else if (side=="West") source+="W";
    else {cout << "can't find Indium in run requested\n"; exit(0);}
  }

  //First get the number of total electron events around the source position from the data file
  Char_t temp[500],tempE[500],tempW[500];
  //sprintf(temp,"%s/replay_pass4_%i.root",getenv("REPLAY_PASS4"),runNumber);
  sprintf(temp,"%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),runNumber);
  TFile *dataFile = new TFile(temp,"READ");
  TTree *data = (TTree*)(dataFile->Get("pass3"));
  string outputBase;
  if (source!="Beta") {
    sprintf(tempE,"Type<3 && PID==1 && Side==0 && (xE.center>(%f-2.*%f) && xE.center<(%f+2.*%f) && yE.center>(%f-2.*%f) && yE.center<(%f+2.*%f))",srcPos[0][0],fabs(srcPos[0][2]),srcPos[0][0],fabs(srcPos[0][2]),srcPos[0][1],fabs(srcPos[0][2]),srcPos[0][1],fabs(srcPos[0][2]));
    sprintf(tempW,"Type<3 && PID==1 && Side==1 && (xW.center>(%f-2.*%f) && xW.center<(%f+2.*%f) && yW.center>(%f-2.*%f) && yW.center<(%f+2.*%f))",srcPos[1][0],fabs(srcPos[1][2]),srcPos[1][0],fabs(srcPos[1][2]),srcPos[1][1],fabs(srcPos[1][2]),srcPos[1][1],fabs(srcPos[1][2]));
    outputBase = string(getenv("REVCALSIM")) + "sources/";
  }
  else if (source=="Beta") {
    sprintf(tempE,"Type<3 && PID==1 && Side==0 && (xE.center*xE.center+yE.center*yE.center)<2500.");
    sprintf(tempW,"Type<3 && PID==1 && Side==1 && (xW.center*xW.center+yW.center*yW.center)<2500.");
    outputBase = string(getenv("REVCALSIM")) + "beta/";
  }
  UInt_t BetaEvents = data->GetEntries(tempE) + data->GetEntries(tempW);
  cout << "Electron Events in Data file: " << BetaEvents << endl;
  cout << "East: " << data->GetEntries(tempE) << "\tWest: " << data->GetEntries(tempW) << endl;
  delete data;
  dataFile->Close();

  //Create simulation output file
  Char_t outputfile[500];
  sprintf(outputfile,"%s/revCalSim_%i_%s.root",outputBase.c_str(),runNumber,source.c_str());
  //sprintf(outputfile,"revCalSim_%i.root",runNumber);
  TFile *outfile = new TFile(outputfile, "RECREATE")
;
  vector <Int_t> pmtQuality = getPMTQuality(runNumber); // Get the quality of the PMTs for that run
  UInt_t calibrationPeriod = getSrcRunPeriod(runNumber); // retrieve the calibration period for this run
  UInt_t XePeriod = getXeRunPeriod(runNumber); // Get the proper Xe run period for the Trigger functions
  GetPositionMap(XePeriod);
  vector <Double_t> alphas = GetAlphaValues(calibrationPeriod); // fill vector with the alpha (nPE/keV) values for this run period
  vector <Double_t> meanEta = getMeanEtaForAlpha(runNumber); // fill vector with the calculated mean position correction for the source used to calculate alpha
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

  //Decide which simulation to use...
  std::string simLocation;
  TChain *chain = new TChain("anaTree");
  
  if (runNumber<20000) simLocation = string(getenv("SIM_2011_2012"));
  else if (runNumber>21087 && runNumber<21679) simLocation = string(getenv("SIM_2012_2013_ISOBUTANE"));
  else simLocation = string(getenv("SIM_2012_2013"));

  std::cout << "Using simulation from " << simLocation << "...\n";

  //Read in simulated data and put in a TChain
 
  int numFiles = source=="Beta"?400:250;
  for (int i=0; i<numFiles; i++) {
    sprintf(temp,"%s/%s/analyzed_%i.root",simLocation.c_str(),source.c_str(),i);
    //sprintf(temp,"../../../data/analyzed_%i.root",i);
    chain->AddFile(temp);
  }

  // Determine the center of the source in simulation in order to construct the displacement vector between
  // the source in data and the source in simulation. This way we can move the sim events to data positions
  /*vector < vector <double> > simSourcePos;
  simSourcePos.resize(2,vector <double> (2,0.));
  TF1 *f1 = new TF1("f1","gaus", -50., 50.);
  TH1F *simPos = new TH1F("simPos","simPos",800, -25., 25.);

  chain->Draw("ScintPosE[0]>>simPos","EdepE>0. && MWPCEnergyE>0. && EdepW==0. && MWPCEnergyW==0."); //East x
  Double_t mBin = simPos->GetMaximumBin();
  Double_t mBinCenter = simPos->GetXaxis()->GetBinCenter(mBin);
  f1->SetParameter(1,mBinCenter);
  simPos->Fit(f1,"L","",mBinCenter-10., mBinCenter+10.);
  simSourcePos[0][0] = f1->GetParameter(1);
  
  chain->Draw("ScintPosE[1]>>simPos","EdepE>0. && MWPCEnergyE>0. && EdepW==0. && MWPCEnergyW==0."); //East y
  mBin = simPos->GetMaximumBin();
  mBinCenter = simPos->GetXaxis()->GetBinCenter(mBin);
  f1->SetParameter(1,mBinCenter);
  simPos->Fit(f1,"L","",mBinCenter-10., mBinCenter+10.);
  simSourcePos[0][1] = f1->GetParameter(1);
  
  chain->Draw("ScintPosW[0]>>simPos","EdepW>0. && MWPCEnergyW>0. && EdepE==0. && MWPCEnergyE==0."); //West x
  mBin = simPos->GetMaximumBin();
  mBinCenter = simPos->GetXaxis()->GetBinCenter(mBin);
  f1->SetParameter(1,mBinCenter);
  simPos->Fit(f1,"L","",mBinCenter-10., mBinCenter+10.);
  simSourcePos[1][0] = f1->GetParameter(1);

  chain->Draw("ScintPosW[1]>>simPos","EdepW>0. && MWPCEnergyW>0. && EdepE==0. && MWPCEnergyE==0."); //West y
  mBin = simPos->GetMaximumBin();
  mBinCenter = simPos->GetXaxis()->GetBinCenter(mBin);
  f1->SetParameter(1,mBinCenter);
  simPos->Fit(f1,"L","",mBinCenter-10., mBinCenter+10.);
  simSourcePos[1][1] = f1->GetParameter(1);

  cout << "Simulated Source Positions:\n"
  << "xE = " << simSourcePos[0][0] <<  "\nyE = " << simSourcePos[0][1] <<  "\nxW = " << simSourcePos[1][0] <<  "\nyW = " << simSourcePos[1][1] << endl; */

  
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
  Double_t MWPCThreshold=0.001;

  //Set random number generator
  TRandom3 *rand = new TRandom3(0);
  TRandom3 *rand2 = new TRandom3(0);
  
  //Get total number of events in TChain
  UInt_t nevents = chain->GetEntries();
  cout << "events = " << nevents << endl;
 
  //Start from random position in evt sequence
  UInt_t evtStart = rand->Rndm()*nevents;
  UInt_t evtTally = 0; //To keep track of the number of events 
  UInt_t evt = evtStart; //current event number
  vector < vector <Int_t> > gridPoint;
  

  //Read in events and determine evt type based on triggers
  while (evtTally<=BetaEvents) {
    if (evt>=nevents) evt=0; //Wrapping the events back to the beginning
    EastScintTrigger = WestScintTrigger = EMWPCTrigger = WMWPCTrigger = false; //Resetting triggers each event

    chain->GetEvent(evt);

    //Checking that the event occurs within the fiducial volume in the simulation to minimize
    // contamination from edge effects and interactions with detector walls
    Double_t fidCut = 45.;
    if (source=="Beta") fidCut=50.; //Since I cut at 50 when determining the number of beta events
    if (sqrt(scint_pos.ScintPosE[0]*scint_pos.ScintPosE[0]+scint_pos.ScintPosE[1]+scint_pos.ScintPosE[1])*sqrt(0.6)*10.>fidCut
	|| sqrt(scint_pos.ScintPosW[0]*scint_pos.ScintPosW[0]+scint_pos.ScintPosW[1]+scint_pos.ScintPosW[1])*sqrt(0.6)*10.>fidCut) {evt++; continue;}
    //cout << "evt passed fiducial cut\n";
    //calculate adjusted event position by sampling a gaussian centered on source position. Do it for primary event side.
    /*Int_t primSide=0., primType=0;
    if (edep.EdepE>0. && edep.EdepW>0.) {
      if (Time.timeE>Time.timeW) {primSide=1;primType=1;}
      else {primSide=0;primType=1;}
    }
    else if (edep.EdepE>0.) primSide=0;
    else if (edep.EdepW>0.) primSide=1;
    else cout << "What else is there?\n";*/
    
    //Double_t primX, primY;
    //primX = rand2->Gaus(srcPos[primSide][0],srcPos[primSide][2]);
    //primY = rand2->Gaus(srcPos[primSide][1],srcPos[primSide][2]);
    
    scint_pos_adj.ScintPosAdjE[0] = source=="Beta"?scint_pos.ScintPosE[0]*sqrt(0.6)*10.:rand2->Gaus(srcPos[0][0], fabs(srcPos[0][2]));
    scint_pos_adj.ScintPosAdjE[1] = source=="Beta"?scint_pos.ScintPosE[1]*sqrt(0.6)*10.:rand2->Gaus(srcPos[0][1], fabs(srcPos[0][2]));
    scint_pos_adj.ScintPosAdjW[0] = source=="Beta"?scint_pos.ScintPosW[0]*sqrt(0.6)*10.:rand2->Gaus(srcPos[1][0], fabs(srcPos[1][2]));//sqrt(0.6)*10.*scint_pos.ScintPosW[0]+displacementX;
    scint_pos_adj.ScintPosAdjW[1] = source=="Beta"?scint_pos.ScintPosW[1]*sqrt(0.6)*10.:rand2->Gaus(srcPos[1][1], fabs(srcPos[1][2]));//sqrt(0.6)*10.*scint_pos.ScintPosW[1]+displacementY;
    scint_pos_adj.ScintPosAdjE[2] = scint_pos.ScintPosE[2]*10.;
    scint_pos_adj.ScintPosAdjW[2] = scint_pos.ScintPosW[2]*10.;
 
    /*Double_t displacementX=0., displacementY=0.;
    if (primType==1) {
      if (primSide==0) {
	displacementX = primX - sqrt(0.6)*10.*scint_pos.ScintPosE[0];
	displacementY = primY - sqrt(0.6)*10.*scint_pos.ScintPosE[1];    
      }
      else {
	displacementX = primX - sqrt(0.6)*10.*scint_pos.ScintPosW[0];
	displacementY = primY - sqrt(0.6)*10.*scint_pos.ScintPosW[1];    
      }
    }
    //for (UInt_t i=0;i<2;i++) {
    //scint_pos_adj.ScintPosAdjE[0] = sqrt(0.6)*10.*scint_pos.ScintPosE[0];//+displacement[0][i]; 
    // scint_pos_adj.ScintPosAdjW[0] = sqrt(0.6)*10.*scint_pos.ScintPosW[0];//+displacement[1][i];
    //}
    if (primSide==0) {
      scint_pos_adj.ScintPosAdjE[0] = primX;
      scint_pos_adj.ScintPosAdjE[1] = primY;
      scint_pos_adj.ScintPosAdjW[0] = sqrt(0.6)*10.*scint_pos.ScintPosW[0]+displacementX;
      scint_pos_adj.ScintPosAdjW[1] = sqrt(0.6)*10.*scint_pos.ScintPosW[1]+displacementY;
    }
    else {
      scint_pos_adj.ScintPosAdjW[0] = primX;
      scint_pos_adj.ScintPosAdjW[1] = primY;
      scint_pos_adj.ScintPosAdjE[0] = sqrt(0.6)*10.*scint_pos.ScintPosE[0]+displacementX;
      scint_pos_adj.ScintPosAdjE[1] = sqrt(0.6)*10.*scint_pos.ScintPosE[1]+displacementY;
      }*/

    
      
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
      if (pmtQuality[p]) { //Check to make sure PMT was functioning
	//cout << p << " " << positionMap[p][intEastBinX][intEastBinY] << endl;
	Double_t posCorrAlpha = alphas[p]*positionMap[p][intEastBinX][intEastBinY]/meanEta[p];
	pmt_Evis.Evis[p] = -1.; // set the Evis for the while loop to return positive first time
	while (pmt_Evis.Evis[p]<0.) { //Must have positive energies since we force the PMTs to intercept 0
	  pmt_Evis.Evis[p] = rand->Gaus(edepQ.EdepQE,sqrt(edepQ.EdepQE/posCorrAlpha));
	}
	if (pmt_Evis.Evis[p]>pmtEnergyLowerLimit) {
	  pmt_Evis.weight[p] = posCorrAlpha/pmt_Evis.Evis[p];
	}
	else {pmt_Evis.weight[p]=posCorrAlpha/pmtEnergyLowerLimit;} //This sets a hard cut on the weight of an event so that events with zero energy still carry weight.. Probably not the most effective way to do this, but it should address low energy behavior
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
	Double_t posCorrAlpha = alphas[p]*positionMap[p][intWestBinX][intWestBinY]/meanEta[p];
	pmt_Evis.Evis[p] = -1.; // set the Evis for the while loop to return positive first time
	while (pmt_Evis.Evis[p]<0.) { //Must have positive energies since we force the PMTs to intercept 0
	  pmt_Evis.Evis[p] = rand->Gaus(edepQ.EdepQW,sqrt(edepQ.EdepQW/posCorrAlpha));
	}
	if (pmt_Evis.Evis[p]>pmtEnergyLowerLimit) {
	  pmt_Evis.weight[p] = posCorrAlpha/pmt_Evis.Evis[p];
	}
	else {pmt_Evis.weight[p]=posCorrAlpha/pmtEnergyLowerLimit;} //This sets a hard cut on the weight of an event so that events with zero energy still carry weight.. Probably not the most effective way to do this, but it should address low energy behavior
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
      
    // Increment the event tally if the event was PID = 1 (electron)
    if (PID==1) evtTally++;
    evt++;
    if (PID>=0) tree->Fill();
    //cout << evtTally << endl;
    if (evtTally%1000==0) {cout << "*";}//cout << "filled event " << evt << endl;
  }
  cout << endl;
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
  
