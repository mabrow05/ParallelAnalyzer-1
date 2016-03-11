/* Written by: Michael A. Brown (UKy)
 Contributors: Xuan Sun (Caltech)

This code takes source or beta simulation data, applies detector effects (nPE smearing, trigger functions),
has the ability to twiddle non-linear terms in the calibration parameters, and fits source peaks.

If the simulation is beta decay, it outputs a tree with spectra 
*/

#include "SimAnalyzerUtils.hh"
#include "posMapReader.h"
#include <map>
  

// This function holds all the meat
void revCalSimulation (std::string source, UInt_t numEvents, bool linCorr);


int main(int argc, char *argv[]) {

  if (argc!=4) {
    std::cout << "Usage: ./SimulationAnalyzer [Source] [numEvents] [bool runParameterSets]\n";
    exit(0);
  }

  std::string source = std::string(argv[1]); //Input source (Beta, Sn113, Bi207, Ce139)
  UInt_t numEvents = (unsigned int)atoi(argv[2]); //Number of electrons to accumulate
  std::string linCorrString = std::string(argv[3]); //False if you want no linearity corrections applied
  bool linCorrBool = false;
  if (linCorrString=="true" || linCorrString=="1") {
    linCorrBool = true;
  }
  //std::vector <Double_t> checker = {5., 0.2, 0.0005,-3.e-6};
  //std::cout << checkParamSetStatus(checker,source.substr(0,2),"2010") << std::endl;
  
  revCalSimulation(source, numEvents, linCorrBool);

}


void revCalSimulation (std::string source, UInt_t numEvents, bool doParamSets) 
{
  std::cout << "Running reverse calibration for " << source << std::endl;

  UInt_t alphaFileIndex = 0; // the index on the nPE_keV file
  UInt_t XeMapPeriod = 2; // The xenon position map to be used
  bool printTree = true; //Boolean to force printing of of TTree. The tree will be printed regardless if the source is Beta
  UInt_t BetaEvents = numEvents;
  std::string geometry = "2010"; // "2010","2011/2012","2012/2013"
  bool runParamSets = doParamSets;
  //std::vector <Double_t> linCorrParams = params;

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

  
  //Histograms of event types for quick checks
  /*std::vector <TH1D> finalEn;// (6,NULL);
  finalEn.push_back(TH1D("finalE0", "Simulated Weighted Sum East Type 0", 400, 0., 1200.));
  finalEn.push_back(TH1D("finalW0", "Simulated Weighted Sum West Type 0", 400, 0., 1200.));
  finalEn.push_back(TH1D("finalE1", "Simulated Weighted Sum East Type 1", 400, 0., 1200.));
  finalEn.push_back(TH1D("finalW1", "Simulated Weighted Sum West Type 1", 400, 0., 1200.));
  finalEn.push_back(TH1D("finalE23", "Simulated Weighted Sum East Type 2/3", 400, 0., 1200.));
  finalEn.push_back(TH1D("finalW23", "Simulated Weighted Sum West Type 2/3", 400, 0., 1200.));
  */


  // Setting the parameter sets. If sources, generate parameters. If beta, read in good parameters
  //std::map <std::string,Double_t> paramDeltas = {{"p0",0.}, {"p1",0.}, {"p2",0.}, {"p3",0.}};
  std::map <std::string, std::pair <Double_t,Double_t> > paramDeltaRanges = {{"p0",std::make_pair(0.,0.)}, 
								      {"p1",std::make_pair(0.,0.)}, 
								      {"p2",std::make_pair(-1.e-4,1.e-4)}, 
								      {"p3",std::make_pair(-1.e-6,1.e-6)}};

  std::map <std::string,Double_t> paramInc = {{"p0",1.}, {"p1",0.1}, {"p2",1.e-7}, {"p3",1.e-9}};
  std::map <std::string,Int_t> paramSteps = {{"p0",0}, {"p1",0}, {"p2",21}, {"p3",21}};

  /*Int_t nParamSets = (int)((paramDeltaRanges.at("p0").second-paramDeltaRanges.at("p0").first)/paramInc.at("p0")+1.)*
    ((paramDeltaRanges.at("p1").second-paramDeltaRanges.at("p1").first)/paramInc.at("p1")+1.)*
    ((paramDeltaRanges.at("p2").second-paramDeltaRanges.at("p2").first)/paramInc.at("p2")+1.)*
    ((paramDeltaRanges.at("p3").second-paramDeltaRanges.at("p3").first)/paramInc.at("p3")+1.);
  */

  std::vector <Double_t> param0;
  std::vector <Double_t> param1;
  std::vector <Double_t> param2;
  std::vector <Double_t> param3;
  
  /*for (Double_t p0 = paramDeltaRanges.at("p0").first; p0<=paramDeltaRanges.at("p0").second; p0+=paramInc.at("p0")) {
      for (Double_t p1 = paramDeltaRanges.at("p1").first; p1<=paramDeltaRanges.at("p1").second; p1+=paramInc.at("p1")) {
	for (Double_t p2 = paramDeltaRanges.at("p2").first; p2<=paramDeltaRanges.at("p2").second; p2+=paramInc.at("p2")) {
	for (Double_t p3 = paramDeltaRanges.at("p3").first; p3<=paramDeltaRanges.at("p3").second; p3+=paramInc.at("p3")) {*/
  Double_t p0,p1,p2,p3;
  if (source!="Beta" && runParamSets==true) {

    for (Int_t i0 = 0; i0<paramSteps.at("p0")+1; i0++) {
      p0 = paramInc.at("p0")*(double)i0-paramInc.at("p0")*(double)paramSteps.at("p0")/2.;

      for (Int_t i1 = 0; i1<paramSteps.at("p1")+1; i1++) {
	p1 = paramInc.at("p1")*(double)i1-paramInc.at("p1")*(double)paramSteps.at("p1")/2.; 

	for (Int_t i2 = 0; i2<paramSteps.at("p2")+1; i2++) {
	  p2 = paramInc.at("p2")*(double)i2-paramInc.at("p2")*(double)paramSteps.at("p2")/2.;

	  for (Int_t i3 = 0; i3<paramSteps.at("p3")+1; i3++) {
	    p3 = paramInc.at("p3")*(double)i3-paramInc.at("p3")*(double)paramSteps.at("p3")/2.;
	    
	    std::cout << p0 << " " << p1 << " " <<  p2 << " " << p3 << std::endl; 
	    std::vector < Double_t > p = {p0,p1,p2,p3};
	    std::string src_hold = source.substr(0,2);

	    if (source=="Bi207") {
	      if (checkParamSetStatus(p,src_hold+std::string("1"),geometry) && checkParamSetStatus(p,src_hold+std::string("2"),geometry)) {
		param0.push_back(p0);
		param1.push_back(p1);
		param2.push_back(p2);
		param3.push_back(p3);
	      }
	    }
	    else {
	      if (checkParamSetStatus(p,src_hold,geometry)) {
		param0.push_back(p0);
		param1.push_back(p1);
		param2.push_back(p2);
		param3.push_back(p3);
	      }
	    }
	  }
	}
      }
    }
  }
  else if (runParamSets==false) {
    param0.push_back(0.);
    param1.push_back(0.);
    param2.push_back(0.);
    param3.push_back(0.);
  }
  else {
    //code later
  }


  //Vectors of histograms and trees for holding source plots for all parameter combinations  
  std::vector <TH1D> Etype0 (param0.size(), TH1D());
  std::vector <TH1D> Wtype0 (param0.size(), TH1D());
  std::vector <TTree*> tree (param0.size(), NULL);

  
  

  Char_t temp[500];
  for (UInt_t i = 0; i<param0.size(); i++) {

    sprintf(temp,"Tree%i",i);
    //tree[i].SetObject(temp,"SimAnalyzed");
    tree[i] = new TTree(temp,"SimAnalyzed");
    SetUpOutputTree(*tree[i]); //Setup the output tree and branches

    sprintf(temp,"East p%i",i);
    Etype0[i] = TH1D(temp,temp,200,0.,1200.);
    sprintf(temp,"West p%i",i);
    Wtype0[i] = TH1D(temp,temp,200,0.,1200.);
    std::cout << param0[i] << " " << param1[i] << " " <<  param2[i] << " " << param3[i] << std::endl; 
  }
 
  
  //Read in simulated data and put in a TChain
  int numFiles = 1; // This is for testing! can be more later depending on how Xuan formats his output of simulated data
  TChain chain("anaTree");
  for (int i=0; i<numFiles; i++) {
    sprintf(temp,"%s/%s/xuan_analyzed.root",simLocation.c_str(),source.c_str());
    chain.AddFile(temp);
  }
  
  // Set the addresses of the information read in from the simulation file
 
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



    for (UInt_t i=0; i<param0.size(); i++) {

      //////////////////////////////////////////////////////////////////////////////////////////////  
      //Applying the twiddle from the linearity curves not being perfect
      std::vector <Double_t> paramHold = {param0[i],param1[i],param2[i],param3[i]};
      evis.EvisE = applyLinearityTwiddle(paramHold,evis.EvisE); //Applies the twiddle
      evis.EvisW = applyLinearityTwiddle(paramHold,evis.EvisW); //Applies the twiddle
      
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
	  if (type==0) Etype0[i].Fill(Erecon); 
	  //else if (type==1) finalEn[2].Fill(Erecon); 
	  //else if(type==2 ||type==3) finalEn[4].Fill(Erecon);
	}
	else Erecon=-1.;
      }
      if (side==1) {
	Double_t totalEvis = type==1 ? (evis.EvisE+evis.EvisW):evis.EvisW;
	if (totalEvis>0.) {
	  Erecon = EQ2Etrue[1][typeIndex][0]+EQ2Etrue[1][typeIndex][1]*totalEvis+EQ2Etrue[1][typeIndex][2]/(totalEvis+EQ2Etrue[1][typeIndex][3])+EQ2Etrue[1][typeIndex][4]/((totalEvis+EQ2Etrue[1][typeIndex][5])*(totalEvis+EQ2Etrue[1][typeIndex][5]));
	  
	  if (type==0) Wtype0[i].Fill(Erecon); 
	  //else if (type==1) finalEn[3].Fill(Erecon); 
	  //else if(type==2 ||type==3) finalEn[5].Fill(Erecon);
	}
	else Erecon=-1.;
      }
      
      if (source=="Beta") tree[i]->Fill();
      
    }

    // Increment the event tally if the event was PID = 1 (electron) and the Erecon was valid
    if (PID==1) evtTally++;
    evt++;

    if (evtTally%1000==0) {std::cout << evtTally << std::endl;}//cout << "filled event " << evt << endl;

  }

  std::cout << endl;
  
  if (source!="Beta") {

    sprintf(temp,"passingParams_%s.dat",source.c_str());
    ofstream ofile(temp);

    for (UInt_t i=0;i<param0.size();i++) {
      //Fit the histograms
      Double_t width = source=="Bi207" ? 60. : (source=="Sn113" ? 45. : 30);
      if (source!="Bi207") {
	Double_t mean = GetXatMax(&Etype0[i]);
	std::vector < std::vector < Double_t > > MeanAndSig;// (2, std::vector<Double_t>);
	MeanAndSig.push_back(FitGaus(&Etype0[i],mean, mean-width, mean+width));
	std::cout << "East mean = " << MeanAndSig[0][0] << "    East sigma = " << MeanAndSig[0][1] << endl;
	
	mean = GetXatMax(&Wtype0[i]);
	MeanAndSig.push_back(FitGaus(&Wtype0[i],mean, mean-width, mean+width));
	std::cout << "West mean = " << MeanAndSig[1][0] << "    West sigma = " << MeanAndSig[1][1] << endl;
	
	//bool printToFile= false;
	//if (geometry=="2010") printToFile = CheckPeakValues2010(

	if (checkPeakStatus(MeanAndSig,source.substr(0,2),geometry)) {
      
	  ofile << i << "\t" << param0[i] << "\t" << param1[i] << "\t" << param2[i] << "\t" << param3[i]
		<< "\t" << MeanAndSig[0][0] << "\t" << MeanAndSig[0][1] << "\t" 
		<< MeanAndSig[1][0] << "\t" << MeanAndSig[1][1] << std::endl;
	}
      }
      else  {
	// CODE BISMUTH PEAK FITTERS
      }
    }

    ofile.close();
    
  }

  TFile *outfile;
  Char_t outputfile[500];

  if (printTree && source!="Beta") {
    //Create simulation output file    
    sprintf(outputfile,"analyzed_files/SimAnalyzed_%s.root",source.c_str());
    outfile = new TFile(outputfile, "RECREATE");
    for (UInt_t i=0; i<param0.size();i++) {
      Etype0[i].Write(); 
      Wtype0[i].Write();
    }
    outfile->Close();
  }
  else if (source=="Beta")  {
    //Create simulation output file
    for (UInt_t i=0; i<param0.size();i++) {
      sprintf(outputfile,"analyzed_files/SimAnalyzed_%s_%i.root",source.c_str(),i);
      outfile = new TFile(outputfile, "RECREATE");
      tree[i]->Write();
      outfile->Close();
    }
  }
 
}

  
