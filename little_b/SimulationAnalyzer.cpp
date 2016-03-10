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
#include <utility>

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


//Setting global peak values for simulated sources
/*std::map < std::string, std::pair < std::string, std::vector < Double_t > > > peaks;
std::vector < Double_t > peak_hold(2,0.);

peak_hold[0] = 0.;
peak_hold[1] = 0.;
peaks["2010"] = std::make_pair("Ce139",peak_hold);

peak_hold[0] = 0.;
peak_hold[1] = 0.;
peaks["2010"] = std::make_pair("Sn113",peak_hold);

peak_hold[0] = 0.;
peak_hold[1] = 0.;
peaks["2010"] = std::make_pair("Bi207",peak_hold);
*/

//Function to return the x value of the max bin in a histogram 
Double_t GetXatMax(TH1D* hist)
{
  Double_t xBinAtMax = -1;
  Int_t binValue = -1;
  for(int i = 1; i < hist -> GetNbinsX(); i++)
  {
    if((hist -> GetBinContent(i)) > binValue)
    {
      binValue = hist -> GetBinContent(i);
      xBinAtMax = i;
    }
  }
  Double_t xAtMax = hist -> GetBinCenter(xBinAtMax);
  return xAtMax;
}

//Function to return a 2 element vector holding the fit mean and sigma of a peak
std::vector <Double_t> FitGaus(TH1D* histToFit, Double_t gausMean, Double_t min, Double_t max)
{
  TF1 *f1 = new TF1("Gaussian Fit", "gaus", min, max);
  histToFit -> Fit("Gaussian Fit", "R");

  std::vector <Double_t> functionFit;

  functionFit.push_back(f1->GetParameter(1));	// mean
  functionFit.push_back(f1->GetParameter(2));	// sigma

  return functionFit;

  // if you wanted to draw the lines or something.
  histToFit -> GetFunction("Gaussian Fit") -> SetLineColor(4);
  histToFit -> GetFunction("Gaussian Fit") -> SetLineStyle(1);
  histToFit -> GetFunction("Gaussian Fit") -> SetLineWidth(4);
  histToFit -> Draw();
}

bool CheckPeakValues2011(vector <int> parameters, std::string sourceName, std::string side)
{
  Double_t mean, sigma, diffMean, diffSigma;

  if(sourceName == "113Sn")
  {
    if(side == "east")
    {
      mean = 365.629;	// keV for all
      sigma = 33.7571;
    }
    else if(side == "west")
    {
      mean = 365.394;
      sigma = 34.2878;
    }
    else { 
    std::cout << "Bad input for side into CheckPeakValue\n";
    exit(0);
    }
  }
  else if(sourceName == "139Ce")
  {
    if(side == "east")
    {
      mean = 0;
      sigma = 0;
    }
    else if(side == "west")
    {
      mean = 0;
      sigma = 0;
    }
    else { 
    std::cout << "Bad input for side into CheckPeakValue\n";
    exit(0);
    }
  }
  else { 
    std::cout << "Bad input for Src into CheckPeakValue\n";
    exit(0);
  }

  // we enforce absolute value as a check
  diffMean = abs(parameters[0]) - mean;
  diffSigma = abs(parameters[1]) - sigma;

  if(diffMean < 0 && diffSigma < 0)
  {
    cout << "Set of parameters are GOOD: difference in mean is " << abs(diffMean)
    << ", difference in sigma is " << abs(diffSigma) << endl;


    return true;
  }
  if(diffMean > 0)
  {
    cout << "Parameters BAD: mean is outside range." << endl;
  }
  if(diffSigma > 0)
  {
    cout << "Parameters BAD: sigma is outside range." << endl;
  }

  return false;
}


bool CheckPeakValues2010(std::vector <int> parameters, std::string sourceName, std::string side)
{
  Double_t mean, sigma, diffMean, diffSigma;

  if(sourceName == "113Sn")
  {
    if(side == "east")
    {
      mean = 0;
      sigma = 0;
    }
    else if(side == "west")
    {
      mean = 0;
      sigma = 0;
    }
    else { 
    std::cout << "Bad input for side into CheckPeakValue\n";
    exit(0);
    }
  }
  else if(sourceName == "139Ce")
  {
    if(side == "east")
    {
      mean = 0;
      sigma = 0;
    }
    else if(side == "west")
    {
      mean = 0;
      sigma = 0;
    }
    else { 
    std::cout << "Bad input for side into CheckPeakValue\n";
    exit(0);
    }
  }
  else { 
    std::cout << "Bad input for Src into CheckPeakValue\n";
    exit(0);
  }
  
  // we enforce absolute value just as a check
  diffMean = abs(parameters[0]) - mean;
  diffSigma = abs(parameters[1]) - sigma;

  if(diffMean < 0 && diffSigma < 0)
  {
    cout << "Set of parameters are GOOD: difference in mean is " << abs(diffMean)
    << ", difference in sigma is " << abs(diffSigma) << endl;


    return true;
  }
  if(diffMean > 0)
  {
    cout << "Parameters BAD: mean is outside range." << endl;
  }
  if(diffSigma > 0)
  {
    cout << "Parameters BAD: sigma is outside range." << endl;
  }

  return false;
}


//Function to return the trigger function for each side in a std::vector in the form vec[side][param]
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
  return params[0]+params[1]*TMath::Erf((En-params[2])/params[3])
    + params[4]*TMath::Gaus(En,params[5],params[6]);
}

// Returns a vector which holds the alpha value (nPE/keV) for each PMT
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
std::vector < std::vector < std::vector <double> > > getEQ2EtrueParams(std::string geometry) {
  ifstream infile;
  if (geometry=="2010") infile.open("../simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat");
  else if (geometry=="2011/2012") infile.open("../simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat");
  else if (geometry=="2012/2013") infile.open("../simulation_comparison/EQ2EtrueConversion/2012-2013_EQ2EtrueFitParams.dat");
  else {
    std::cout << "Bad geometry passed to getEQ2EtrueParams\n";
    exit(0);
  }
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
  
//Function to return the correction to linearity offset
Double_t applyLinearityTwiddle (std::vector <Double_t> &params, Double_t EQ) {
  return params[0]+(1+params[1])*EQ+params[2]*EQ*EQ+params[3]*EQ*EQ*EQ;
}

// Sets up the output tree for the simulated data
void SetUpOutputTree(TTree& tree) {
  tree.Branch("PID", &PID, "PID/I");
  tree.Branch("side", &side, "side/I");
  tree.Branch("type", &type, "type/I");
  tree.Branch("Erecon", &Erecon,"Erecon/D");
  
  tree.Branch("Evis",&evis,"EvisE/D:EvisW");
  tree.Branch("Edep",&edep,"EdepE/D:EdepW");
  tree.Branch("EdepQ",&edepQ,"EdepQE/D:EdepQW");
  tree.Branch("Eprim",&Eprim,"Eprim/D");
  tree.Branch("AsymWeight",&AsymWeight,"AsymWeight/D");
  
  tree.Branch("time",&Time,"timeE/D:timeW");
  tree.Branch("MWPCEnergy",&mwpcE,"MWPCEnergyE/D:MWPCEnergyW");
  tree.Branch("MWPCPos",&mwpc_pos,"MWPCPosE[3]/D:MWPCPosW[3]");
  tree.Branch("ScintPos",&scint_pos,"ScintPosE[3]/D:ScintPosW[3]");
  tree.Branch("ScintPosAdjusted",&scint_pos_adj,"ScintPosAdjE[3]/D:ScintPosAdjW[3]");
  tree.Branch("PMT_Evis",&pmt_Evis,"Evis0/D:Evis1:Evis2:Evis3:Evis4:Evis5:Evis6:Evis7:weight0:weight1:weight2:weight3:weight4:weight5:weight6:weight7");
  
};
  

void revCalSimulation (std::string source, UInt_t numEvents, bool linCorr, std::vector <Double_t> params) 
{
  std::cout << "Running reverse calibration for " << source << std::endl;

  UInt_t alphaFileIndex = 0; // the index on the nPE_keV file
  UInt_t XeMapPeriod = 2; // The xenon position map to be used
  bool printTree = true; //Boolean to force printing of of TTree. The tree will be printed regardless if the source is Beta
  UInt_t BetaEvents = numEvents;
  std::string geometry = "2010"; // "2010","2011/2012","2012/2013"
  bool doLinCorr = linCorr;
  std::vector <Double_t> linCorrParams = params;

  // Set the location of the source and it's spread in mm. srcPos[0=east,1=west][0=x,1=y,2=sigma]
  bool moveSource = false; // If this is true, set the position of the source to the following location
                           // and sample a gaussian of position with sigma at srcPos[][2].
  std::vector < std::vector <Double_t> > srcPos(2, std::vector <Double_t> (3,0.));
  srcPos[0][0] = srcPos[1][0] = 0.;
  srcPos[0][1] = srcPos[1][1] = 0.;
  srcPos[0][2] = srcPos[1][2] = 5.;

  //Decide which simulation to use...
  std::string simLocation;
  if (geometry=="2010") simLocation = "../../data/XuanSim/"; //This needs to be set to wherever you are storing your sims..
  else if (geometry=="2011/2012") simLocation = string(getenv("SIM_2011_2012"));
  else if (geometry=="2012/2013") simLocation = string(getenv("SIM_2012_2013"));
  else if (geometry=="2012/2013_isobutane") simLocation = string(getenv("SIM_2012_2013_ISOBUTANE"));
  else {
    std::cout << "The geometry is set to an invalid option\n\n";
    exit(0);
  }

  std::cout << "Using simulation from " << simLocation << "...\n";

  //Create simulation output file
  //Char_t outputfile[500];
  //sprintf(outputfile,"analyzed_files/SimAnalyzed_%s.root",source.c_str());
  //TFile outfile;
  //if (source=="Beta" || printTree) outfile.Open(outputfile, "RECREATE");


  // Load information for applying detector effects
  std::vector <Double_t> alphas = getAlphaValues(alphaFileIndex); // fill vector with the alpha (nPE/keV) values
  std::vector < std::vector <Double_t> > triggerFunc = getTriggerFunctionParams(XeMapPeriod,7); // 2D vector with trigger function for East side and West side in that order
  GetPositionMap(XeMapPeriod); //Loads position map via posMapReader.h methods
  std::vector < std::vector < std::vector <double> > > EQ2Etrue = getEQ2EtrueParams(geometry);

  //Setup the output tree if the flag is set above, or we have beta decay source
  TTree tree("SimAnalyzed","SimAnalyzed");
  SetUpOutputTree(tree); //Setup the output tree and branches
  
  //Histograms of event types for quick checks
  std::vector <TH1D> finalEn;// (6,NULL);
  finalEn.push_back(TH1D("finalE0", "Simulated Weighted Sum East Type 0", 400, 0., 1200.));
  finalEn.push_back(TH1D("finalW0", "Simulated Weighted Sum West Type 0", 400, 0., 1200.));
  finalEn.push_back(TH1D("finalE1", "Simulated Weighted Sum East Type 1", 400, 0., 1200.));
  finalEn.push_back(TH1D("finalW1", "Simulated Weighted Sum West Type 1", 400, 0., 1200.));
  finalEn.push_back(TH1D("finalE23", "Simulated Weighted Sum East Type 2/3", 400, 0., 1200.));
  finalEn.push_back(TH1D("finalW23", "Simulated Weighted Sum West Type 2/3", 400, 0., 1200.));
  
  
  //Read in simulated data and put in a TChain
  int numFiles = 1; // This is for testing! can be more later depending on how Xuan formats his output of simulated data
  Char_t temp[500];
  TChain chain("anaTree");
  for (int i=0; i<numFiles; i++) {
    sprintf(temp,"%s/%s/xuan_analyzed.root",simLocation.c_str(),source.c_str());
    chain.AddFile(temp);
  }
  
  // Set the addresses of the information read in from the simulation file
  Int_t primaryID;
  Double_t initialMomentum[3]; //For reconstructing which detector type 1s initially struck
  // The rest of the variables are defined in SimulationAnalyzer.hh for no good reason (legacy)
  chain.SetBranchAddress("primaryParticleSpecies",&primaryID);
  chain.SetBranchAddress("primaryMomentum",&initialMomentum);
  chain.SetBranchAddress("mwpcEnergy",&mwpcE);
  //chain.SetBranchAddress("time",&Time);
  chain.SetBranchAddress("scintillatorEdep",&edep);
  chain.SetBranchAddress("scintillatorEdepQuenched",&edepQ);
  chain.SetBranchAddress("MWPCPos",&mwpc_pos);
  chain.SetBranchAddress("ScintPos",&scint_pos);
  chain.SetBranchAddress("primaryKE",&Eprim);
  
  //Trigger booleans
  bool EastScintTrigger, WestScintTrigger, EMWPCTrigger, WMWPCTrigger;
  Double_t MWPCThreshold=0.001; // generic MWPC Threshold

  //Set random number generators for sampling of gaussians in smearing
  TRandom3 rand(0);
  TRandom3 rand2(0);
  
  //Get total number of events in TChain
  UInt_t nevents = chain.GetEntries();
  cout << "events in chain = " << nevents << endl;
 
  //Start from random position in evt sequence
  UInt_t evtStart = rand.Rndm()*nevents;
  UInt_t evtTally = 0; //To keep track of the number of events 
  UInt_t evt = evtStart; //current event number

  std::vector < std::vector <Int_t> > gridPoint; // Will hold the grid location for the position map

  //Read in events and determine evt type based on triggers
  while (evtTally<BetaEvents) {
    if (evt>=nevents) evt=0; //Wrapping the events back to the beginning
    EastScintTrigger = WestScintTrigger = EMWPCTrigger = WMWPCTrigger = false; //Resetting triggers each event

    chain.GetEvent(evt);
    
    if (primaryID!=11) {evt++;continue;}
    
    //Dividing out the weighted energy dependence which is currently in the position branch of input tree
    scint_pos.ScintPosE[0] /= edepQ.EdepQE>0. ? edepQ.EdepQE : 1.;
    scint_pos.ScintPosW[0] /= edepQ.EdepQW>0. ? edepQ.EdepQW : 1.;
    scint_pos.ScintPosE[1] /= edepQ.EdepQE>0. ? edepQ.EdepQE : 1.;
    scint_pos.ScintPosW[1] /= edepQ.EdepQW>0. ? edepQ.EdepQW : 1.;
    scint_pos.ScintPosE[2] /= edepQ.EdepQE>0. ? edepQ.EdepQE : 1.;
    scint_pos.ScintPosW[2] /= edepQ.EdepQW>0. ? edepQ.EdepQW : 1.;
   
    //Setting adjusted positions when we are simulating a source and move source flag is set to true
    if (moveSource && source!="Beta") {
      scint_pos_adj.ScintPosAdjE[0] = rand2.Gaus(srcPos[0][0], fabs(srcPos[0][2]));
      scint_pos_adj.ScintPosAdjE[1] = rand2.Gaus(srcPos[0][1], fabs(srcPos[0][2]));
      scint_pos_adj.ScintPosAdjW[0] = rand2.Gaus(srcPos[1][0], fabs(srcPos[1][2]));
      scint_pos_adj.ScintPosAdjW[1] = rand2.Gaus(srcPos[1][1], fabs(srcPos[1][2]));
    }
    else {
      scint_pos_adj.ScintPosAdjE[0] = scint_pos.ScintPosE[0]*sqrt(0.6)*1000.;
      scint_pos_adj.ScintPosAdjE[1] = scint_pos.ScintPosE[1]*sqrt(0.6)*1000.;
      scint_pos_adj.ScintPosAdjW[0] = scint_pos.ScintPosW[0]*sqrt(0.6)*1000.;
      scint_pos_adj.ScintPosAdjW[1] = scint_pos.ScintPosW[1]*sqrt(0.6)*1000.;
    }
    scint_pos_adj.ScintPosAdjE[2] = scint_pos.ScintPosE[2]*1000.;
    scint_pos_adj.ScintPosAdjW[2] = scint_pos.ScintPosW[2]*1000.;
    
    
    //retrieve point on grid for each side of detector [E/W][x/y] to implement position map
    gridPoint = getGridPoint(scint_pos_adj.ScintPosAdjE[0],scint_pos_adj.ScintPosAdjE[1],scint_pos_adj.ScintPosAdjW[0],scint_pos_adj.ScintPosAdjW[1]);

    Int_t intEastBinX = gridPoint[0][0];
    Int_t intEastBinY = gridPoint[0][1];
    Int_t intWestBinX = gridPoint[1][0];
    Int_t intWestBinY = gridPoint[1][1];
      
    //MWPC triggers checked against threshold set above
    if (mwpcE.MWPCEnergyE>MWPCThreshold) EMWPCTrigger=true;
    if (mwpcE.MWPCEnergyW>MWPCThreshold) WMWPCTrigger=true;

    Double_t pmtEnergyLowerLimit = 1.; //To put a hard cut on the weight to avoid dividing by zero
    

    //East Side smeared PMT energies
    for (UInt_t p=0; p<4; p++) {
      Double_t posCorrAlpha = alphas[p]*positionMap[p][intEastBinX][intEastBinY];
      if (posCorrAlpha==0.) posCorrAlpha=1.; //This occurs when the event occurs outside the position map... This is to prevent the Gaussian from having inf sigma
      //Sample a gaussian centered on the EQ of the event with a spread related to alpha
      pmt_Evis.Evis[p] = rand.Gaus(edepQ.EdepQE,sqrt(edepQ.EdepQE/posCorrAlpha));
      //Set the weight of the event based on the number of photoelectrons as seen in the detector
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
    if (rand.Rndm(0)<triggProb && evis.EvisE>0.) EastScintTrigger=true;
   
      
    //West Side
    for (UInt_t p=4; p<8; p++) {
      Double_t posCorrAlpha = alphas[p]*positionMap[p][intWestBinX][intWestBinY];
      if (posCorrAlpha==0.) posCorrAlpha=1.; //This occurs when the event occurs outside the position map... This is to prevent the Gaussian from having inf sigma
      pmt_Evis.Evis[p] = rand.Gaus(edepQ.EdepQW,sqrt(edepQ.EdepQW/posCorrAlpha));
      
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
    if (rand.Rndm(0)<triggProb && evis.EvisW>0.) WestScintTrigger = true;     

    //////////////////////////////////////////////////////////////////////////////////////////////  
    //Applying the twiddle from the linearity curves not being perfect
    if (doLinCorr) {
      evis.EvisE = applyLinearityTwiddle(linCorrParams,evis.EvisE); //Applies the twiddle
      evis.EvisW = applyLinearityTwiddle(linCorrParams,evis.EvisW); //Applies the twiddle
    }
    // std::cout << evis.EvisE << "\t" << evis.EvisW << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////

   
    // Do event ID selection on data (to mimic hardware and software event selection)
    PID=-1;

    //Type 0 East
    if (EastScintTrigger && EMWPCTrigger && !WestScintTrigger && !WMWPCTrigger) {
      PID=1;
      type=0;
      side=0;
    }
    //Type 0 West
    else if (WestScintTrigger && WMWPCTrigger && !EastScintTrigger && !EMWPCTrigger) {
      PID=1;
      type=0;
      side=1;
    }
    //Type 1 
    else if (EastScintTrigger && EMWPCTrigger && WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=1;
      //East
      if (initialMomentum[2]<0.) {//(Time.timeE<Time.timeW) {
	side=0;
      }
      //West
      else if (initialMomentum[2]>=0.) {//(Time.timeE>Time.timeW) {
	side=1;
      }
    }
    //Type 2/3 East
    else if (EastScintTrigger && EMWPCTrigger && !WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=2;
      side=0;
    }
    //Type 2/3 West
    else if (!EastScintTrigger && EMWPCTrigger && WestScintTrigger && WMWPCTrigger) {
      PID=1;
      type=2;
      side=1;
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
    Int_t typeIndex = (type==0 || type==4) ? 0:(type==1 ? 1:2); //for retrieving the parameters from EQ2Etrue
    if (side==0) {
      Double_t totalEvis = type==1 ? (evis.EvisE+evis.EvisW):evis.EvisE;
      if (totalEvis>0.) {
	Erecon = EQ2Etrue[0][typeIndex][0]+EQ2Etrue[0][typeIndex][1]*totalEvis+EQ2Etrue[0][typeIndex][2]/(totalEvis+EQ2Etrue[0][typeIndex][3])+EQ2Etrue[0][typeIndex][4]/((totalEvis+EQ2Etrue[0][typeIndex][5])*(totalEvis+EQ2Etrue[0][typeIndex][5]));
	if (type==0) finalEn[0].Fill(Erecon); 
	else if (type==1) finalEn[2].Fill(Erecon); 
	else if(type==2 ||type==3) finalEn[4].Fill(Erecon);
      }
      else Erecon=-1.;
    }
    if (side==1) {
      Double_t totalEvis = type==1 ? (evis.EvisE+evis.EvisW):evis.EvisW;
      if (totalEvis>0.) {
	Erecon = EQ2Etrue[1][typeIndex][0]+EQ2Etrue[1][typeIndex][1]*totalEvis+EQ2Etrue[1][typeIndex][2]/(totalEvis+EQ2Etrue[1][typeIndex][3])+EQ2Etrue[1][typeIndex][4]/((totalEvis+EQ2Etrue[1][typeIndex][5])*(totalEvis+EQ2Etrue[1][typeIndex][5]));
	
	if (type==0) finalEn[1].Fill(Erecon); 
	else if (type==1) finalEn[3].Fill(Erecon); 
	else if(type==2 ||type==3) finalEn[5].Fill(Erecon);
      }
      else Erecon=-1.;
    }
      
    // Increment the event tally if the event was PID = 1 (electron) and the Erecon was valid
    if (PID==1 && Erecon!=-1.)
      {evtTally++;}

    evt++;
    if (PID>=0 && Erecon!=-1.) tree.Fill();
   
    if (evtTally%1000==0) {std::cout << evtTally << std::endl;}//cout << "filled event " << evt << endl;
  }
  std::cout << endl;

  if (source!="Beta") {
    //Fit the histograms
    Double_t width = source=="Bi207" ? 60. : (source=="Sn113" ? 45. : 30);
    Double_t mean = GetXatMax(&finalEn[0]);
    std::vector <Double_t> EastMeanAndSig = FitGaus(&finalEn[0],mean, mean-width, mean+width);
    std::cout << "East mean = " << EastMeanAndSig[0] << "    East sigma = " << EastMeanAndSig[1] << endl;
    
    mean = GetXatMax(&finalEn[1]);
    std::vector <Double_t> WestMeanAndSig = FitGaus(&finalEn[1],mean, mean-width, mean+width);
    std::cout << "West mean = " << WestMeanAndSig[0] << "    West sigma = " << WestMeanAndSig[1] << endl;

    //bool printToFile= false;
    //if (geometry=="2010") printToFile = CheckPeakValues2010(
    sprintf(temp,"passingParams_%s.dat",source.c_str());
    ofstream ofile(temp,ios::app);
    ofile << linCorrParams[0] << "\t" << linCorrParams[1] << "\t" << linCorrParams[2] << "\t" << linCorrParams[3]
	  << "\t" << EastMeanAndSig[0] << "\t" << EastMeanAndSig[1] << "\t" 
	  << WestMeanAndSig[0] << "\t" << WestMeanAndSig[1] << std::endl;
    ofile.close();
  }
  TFile *outfile;

  if (source=="Beta" || printTree)  {
    //Create simulation output file
    Char_t outputfile[500];
    sprintf(outputfile,"analyzed_files/SimAnalyzed_%s.root",source.c_str());
    outfile = new TFile(outputfile, "RECREATE");
    //outfile.Open(outputfile, "RECREATE");
    finalEn[0].Write(); 
    finalEn[1].Write();
    tree.Write();
    outfile->Close();
  }
}

int main(int argc, char *argv[]) {
  if (argc!=4 && argc!=8) {
    std::cout << "Usage: ./SimulationAnalyzer [Source] [numEvents] [bool linCorrections] if linCorrections, [p0] [p1] [p2] [p3]\n";
    exit(0);
  }
  std::string source = std::string(argv[1]); //Input source (Beta, Sn113, Bi207, Ce139)
  UInt_t numEvents = (unsigned int)atoi(argv[2]); //Number of electrons to accumulate
  std::string linCorrString = std::string(argv[3]); //False if you want no linearity corrections applied
  bool linCorrBool = false;
  std::vector <Double_t> params(4,0.);
  if (linCorrString=="true" || linCorrString=="1") {
    linCorrBool = true;
    params[0]=atof(argv[4]);
    params[1]=atof(argv[5]);
    params[2]=atof(argv[6]);
    params[3]=atof(argv[7]);
  }
  cout << "made it here\n";
  revCalSimulation(source, numEvents, linCorrBool, params);

  //tests
  /*UInt_t XePeriod = getXeRunPeriod(atoi(argv[1]));
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7);
  Double_t triggProbE = triggerProbability(triggerFunc[0],25.);
  Double_t triggProbW = triggerProbability(triggerFunc[1],25.);
  cout << triggProbE << " " << triggProbW << endl;*/


}
  
