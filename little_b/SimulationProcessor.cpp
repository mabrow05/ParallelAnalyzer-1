/* Written by: Michael A. Brown (UKy)
 Contributors: Xuan Sun (Caltech)

This code takes source or beta simulation data, applies detector effects (nPE smearing, trigger functions),
has the ability to twiddle non-linear terms in the calibration parameters, and fits source peaks.

If the simulation is beta decay, it outputs a tree with spectra 
*/

#include "SimulationProcessor.hh"
#include "posMapReader.h"
#include <map>
  



int main(int argc, char *argv[]) {

  if (argc!=5 && argc!=13 && argc!=14) {
    std::cout << "Usage: ./SimulationAnalyzer [Source] [geometry] [numEvents] [bool doLinCorr] if doLinCorr [Ep0] [Ep1] [Ep2] [Ep3] [Wp0] [Wp1] [Wp2] [Wp3]\n";
    exit(0);
  }
  Int_t paramSetIndex = -1;
  if (argc==10) paramSetIndex = atoi(argv[9]);

  std::string source = std::string(argv[1]); //Input source (Beta, Beta_fierz, Sn113, Bi207, Ce139)
  UInt_t numEvents = (unsigned int)atoi(argv[3] ); //Number of electrons to accumulate
  std::string linCorrString = std::string(argv[4]); //False if you want no linearity corrections applied
  bool linCorrBool = false;
  std::vector < std::vector <Double_t> > params(2, std::vector <Double_t>(4,0.));	// a 2-vector of 4-vector where
  if (linCorrString=="True" || linCorrString=="true" || linCorrString=="1") {		// the 4-vector is a vec of polynomial coeff
    linCorrBool = true;
    params[0][0] = atof(argv[5]);
    params[0][1] = atof(argv[6]);
    params[0][2] = atof(argv[7]);
    params[0][3] = atof(argv[8]);
    params[1][0] = atof(argv[9]);
    params[1][1] = atof(argv[10]);
    params[1][2] = atof(argv[11]);
    params[1][3] = atof(argv[12]);
  }
  //std::vector <Double_t> checker = {5., 0.2, 0.0005,-3.e-6};
  //std::cout << checkParamSetStatus(checker,source.substr(0,2),"2010") << std::endl;
  
  revCalSimulation(source,std::string(argv[2]), numEvents, linCorrBool, params, paramSetIndex);

}


void revCalSimulation (std::string source, std::string geom, UInt_t numEvents, bool doParamSets, std::vector < std::vector <Double_t> > params, Int_t paramSetIndex) 
{
  std::cout << "Running reverse calibration for " << source << std::endl;

  std::string geometry = geom; // "2010","2011/2012","2012/2013"
  UInt_t alphaFileIndex = 0; // the index on the nPE_keV file
  UInt_t XeMapPeriod = 2; // The xenon position map to be used
  UInt_t triggerMap = geometry=="2010"?-1:(geometry=="2011/2012"?2:10);
  bool printTree = true; //Boolean to force printing of of TTree. The tree will be printed regardless if the source is Beta
  UInt_t BetaEvents = numEvents;
  
  std::vector < std::vector <Double_t> > linCorrParams = params;

  /*
  //Checking that the parameters have a chance
  if (source.substr(0,4)!="Beta" && paramSetIndex!=-1) { 
    if (source!="Bi207") {
      if (!checkParamSetStatus(linCorrParams,source.substr(0,2),geometry)) {
	std::cout << "Sn or Ce source parameters didn't pass envelope check!\n";
	return;
      }
    }
    else if (!checkParamSetStatus(linCorrParams,source.substr(0,2)+"2",geometry) || !checkParamSetStatus(linCorrParams,source.substr(0,2)+"1",geometry)) {
      std::cout << "Bi source parameters didn't pass envelope check!\n";
      return;
    }
  }
  */

  // Set the location of the source and it's spread in mm. srcPos[0=east,1=west][0=x,1=y,2=sigma]
  bool moveSource = false; // If this is true, set the position of the source to the following location
                           // and sample a gaussian of position with sigma at srcPos[][2].
  std::vector < std::vector <Double_t> > srcPos(2, std::vector <Double_t> (3,0.));
  srcPos[0][0] = srcPos[1][0] = 0.;
  srcPos[0][1] = srcPos[1][1] = 0.;
  srcPos[0][2] = srcPos[1][2] = 5.;

  //Decide which simulation to use...
  std::string simLocation;
  if (geometry=="2010") simLocation = std::string(getenv("CALTECH_SIMS"))+"2010/"; //This needs to be set to wherever you are storing your sims..
  else if (geometry=="2011/2012") simLocation = std::string(getenv("CALTECH_SIMS"))+"2011-2012/";
  else if (geometry=="2012/2013") simLocation = std::string(getenv("CALTECH_SIMS"))+"2012-2013/";
  else if (geometry=="2012/2013_isobutane") simLocation = std::string(getenv("CALTECH_SIMS"))+"2012-2013_ISOBUTANE";
  else {
    std::cout << "The geometry is set to an invalid option\n\n";
    exit(0);
  }

  std::cout << "Using simulation from " << simLocation << "...\n";


  // Load information for applying detector effects
  std::vector <Double_t> alphas = getAlphaValues(alphaFileIndex); // fill vector with the alpha (nPE/keV) values
  std::vector < std::vector <Double_t> > triggerFunc = getTriggerFunctionParams(triggerMap,8); // 2D vector with trigger function for East side and West side in that order
  GetPositionMap(XeMapPeriod); //Loads position map via posMapReader.h methods
  std::vector < std::vector < std::vector <double> > > EQ2Etrue = getEQ2EtrueParams(geometry);

  // Histograms for plotting the peaks to fit
  TH1D Etype0("Etype0","East type 0",200,0.,1200.);
  TH1D Wtype0("Wtype0","West type 0",200,0.,1200.);

  TTree *tree = new TTree("SimAnalyzed","SimAnalyzed");
  SetUpOutputTree(*tree); //Setup the output tree and branches
  
  
  //Read in simulated data and put in a TChain
  Char_t temp[500];
  int numFiles = 1; // This is for testing! can be more later depending on how Xuan formats his output of simulated data
  TChain chain("anaTree");
  for (int i=0; i<numFiles; i++) {
    if (source.substr(0,4)!="Beta") sprintf(temp,"%s/%s/xuan_analyzed.root",simLocation.c_str(),source.c_str());
    else if (source.length()>4) sprintf(temp,"%s/%s/%s/xuan_analyzed.root",simLocation.c_str(),source.substr(0,4).c_str(),"fierz");
    else sprintf(temp,"%s/%s/%s/xuan_analyzed.root",simLocation.c_str(),source.substr(0,4).c_str(),"base");
    chain.AddFile(temp);
  }
  
  // Set the addresses of the information read in from the simulation file
 
  chain.SetBranchAddress("primaryParticleSpecies",&primaryID);
  chain.SetBranchAddress("primaryMomentum",&initialMomentum);
  chain.SetBranchAddress("mwpcEnergy",&mwpcE);
  chain.SetBranchAddress("scintTimeToHit",&Time);	// we can save this variable, Xuan sim has it
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
/*    // this is no longer needed since Xuan explicitly divided out the energy dependence in simulation
    scint_pos.ScintPosE[0] /= edepQ.EdepQE>0. ? edepQ.EdepQE : 1.;
    scint_pos.ScintPosW[0] /= edepQ.EdepQW>0. ? edepQ.EdepQW : 1.;
    scint_pos.ScintPosE[1] /= edepQ.EdepQE>0. ? edepQ.EdepQE : 1.;
    scint_pos.ScintPosW[1] /= edepQ.EdepQW>0. ? edepQ.EdepQW : 1.;
    scint_pos.ScintPosE[2] /= edepQ.EdepQE>0. ? edepQ.EdepQE : 1.;
    scint_pos.ScintPosW[2] /= edepQ.EdepQW>0. ? edepQ.EdepQW : 1.;
*/
   
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
    Double_t totalEnE = numer/denom;
    evis.EvisE = totalEnE;
    
      
    //West Side
    for (UInt_t p=4; p<8; p++) {
      Double_t posCorrAlpha = alphas[p]*positionMap[p][intWestBinX][intWestBinY];
      //std::cout << posCorrAlpha << std::endl;
      if (posCorrAlpha==0.) posCorrAlpha=1.; //This occurs when the event occurs outside the position map... This is to prevent the Gaussian from having inf sigma
      pmt_Evis.Evis[p] = rand.Gaus(edepQ.EdepQW,sqrt(edepQ.EdepQW/posCorrAlpha));
      
      //std::cout << edepQ.EdepQW << " " << pmt_Evis.Evis[p] << std::endl;

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
    //std::cout << numer << std::endl;
    Double_t totalEnW = numer/denom;
    evis.EvisW = totalEnW;
    

    //////////////////////////////////////////////////////////////////////////////////////////////  
    //Applying the twiddle from the linearity curves not being perfect
    std::vector <Double_t> paramHold = {linCorrParams[0][0],linCorrParams[0][1],linCorrParams[0][2],linCorrParams[0][3]}; //East parameters
    evis.EvisE = applyLinearityTwiddle(paramHold,evis.EvisE); //Applies the twiddle
    paramHold = {linCorrParams[1][0],linCorrParams[1][1],linCorrParams[1][2],linCorrParams[1][3]}; //West parameters
    evis.EvisW = applyLinearityTwiddle(paramHold,evis.EvisW); //Applies the twiddle
    
    // std::cout << evis.EvisE << "\t" << evis.EvisW << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////
    

    // Now we can apply the trigger probability

    //EAST Trigger
    Double_t triggProb = triggerProbability(triggerFunc[0],evis.EvisE);
    if (rand.Rndm(0)<triggProb && evis.EvisE>0.)  {/*std::cout << "East Trigger\n";*/EastScintTrigger=true;}
    //WEST Trigger
    triggProb = triggerProbability(triggerFunc[1],evis.EvisW);
    if (rand.Rndm(0)<triggProb && evis.EvisW>0.) {WestScintTrigger = true;}     

    
      
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
      if (Time.timeE<Time.timeW) {
	side=0;
      }
      //West
      else {//if (Time.timeE>Time.timeW) { 
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
	if (type==0) Etype0.Fill(Erecon); 
	//else if (type==1) finalEn[2].Fill(Erecon); 
	//else if(type==2 ||type==3) finalEn[4].Fill(Erecon);
      }
      else Erecon=-1.;
    }
    if (side==1) {
      Double_t totalEvis = type==1 ? (evis.EvisE+evis.EvisW):evis.EvisW;
      if (totalEvis>0.) {
	Erecon = EQ2Etrue[1][typeIndex][0]+EQ2Etrue[1][typeIndex][1]*totalEvis+EQ2Etrue[1][typeIndex][2]/(totalEvis+EQ2Etrue[1][typeIndex][3])+EQ2Etrue[1][typeIndex][4]/((totalEvis+EQ2Etrue[1][typeIndex][5])*(totalEvis+EQ2Etrue[1][typeIndex][5]));
	
	if (type==0) Wtype0.Fill(Erecon); 
	//else if (type==1) finalEn[3].Fill(Erecon); 
	//else if(type==2 ||type==3) finalEn[5].Fill(Erecon);
      }
      else Erecon=-1.;
    }
    
    if (source.substr(0,4)=="Beta") tree->Fill();
    
    // Increment the event tally if the event was PID = 1 (electron) and the Erecon was valid
    if (PID==1) evtTally++;
    evt++;
    
    if (evtTally%1000==0) {std::cout << evtTally << std::endl;}//cout << "filled event " << evt << endl;
  }
  std::cout << endl;
  
  if (source.substr(0,4)!="Beta") {
    
    sprintf(temp,"linCurves/passingParams_%s_%s.dat",geometry.substr(0,4).c_str(),source.c_str());
    ofstream ofile(temp,ios::app);
    
    //Fit the histograms
    Double_t width = source=="Bi207" ? 60. : (source=="Sn113" ? 45. : 30);
    if (source!="Bi207") {
      Double_t mean = GetXatMax(&Etype0);
      std::vector < std::vector < Double_t > > MeanAndSig;// (2, std::vector<Double_t>);
      MeanAndSig.push_back(FitGaus(&Etype0,mean, mean-width, mean+width));
      std::cout << "East mean = " << MeanAndSig[0][0] << "    East sigma = " << MeanAndSig[0][1] << endl;
      
      mean = GetXatMax(&Wtype0);
      MeanAndSig.push_back(FitGaus(&Wtype0,mean, mean-width, mean+width));
      std::cout << "West mean = " << MeanAndSig[1][0] << "    West sigma = " << MeanAndSig[1][1] << endl;
      
      //bool printToFile= false;
      
      if (paramSetIndex!=-1 && checkPeakStatus(MeanAndSig,source.substr(0,2),geometry)) {
	
	ofile << paramSetIndex << "\t" << linCorrParams[0][0] << "\t" << linCorrParams[0][1] << "\t" << linCorrParams[0][2] << "\t" << linCorrParams[0][3]
	      << "\t" << linCorrParams[1][0] << "\t" << linCorrParams[1][1] << "\t" << linCorrParams[1][2] << "\t" << linCorrParams[1][3]
	      << "\t" << MeanAndSig[0][0] << "\t" << MeanAndSig[0][1] << "\t" 
	      << MeanAndSig[1][0] << "\t" << MeanAndSig[1][1] << std::endl;
      }
    }
    else  {
      // CODE BISMUTH PEAK FITTERS
      Double_t mean = GetXatMax(&Etype0);
      std::vector < std::vector < Double_t > > MeanAndSig1;// (2, std::vector<Double_t>);
      MeanAndSig1.push_back(FitGaus(&Etype0,mean, mean-width, mean+width));
      std::cout << "East mean = " << MeanAndSig1[0][0] << "    East sigma = " << MeanAndSig1[0][1] << endl;
      
      mean = GetXatMax(&Wtype0);
      MeanAndSig1.push_back(FitGaus(&Wtype0,mean, mean-width, mean+width));
      std::cout << "West mean = " << MeanAndSig1[1][0] << "    West sigma = " << MeanAndSig1[1][1] << endl;
      
      
      mean = GetXatMax(&Etype0, 400., 650.);
      std::vector < std::vector < Double_t > > MeanAndSig2;// (2, std::vector<Double_t>);
      MeanAndSig2.push_back(FitGaus(&Etype0,mean, mean-width, mean+width));
      std::cout << "East mean = " << MeanAndSig2[0][0] << "    East sigma = " << MeanAndSig2[0][1] << endl;
      
      mean = GetXatMax(&Wtype0, 400., 650.);
      MeanAndSig2.push_back(FitGaus(&Wtype0,mean, mean-width, mean+width));
      std::cout << "West mean = " << MeanAndSig2[1][0] << "    West sigma = " << MeanAndSig2[1][1] << endl;

      if (paramSetIndex!=-1 && checkPeakStatus(MeanAndSig1,source.substr(0,2)+"1",geometry) && checkPeakStatus(MeanAndSig2,source.substr(0,2)+"2",geometry)) {
	
	ofile << paramSetIndex << "\t" << linCorrParams[0][0] << "\t" << linCorrParams[0][1] << "\t" << linCorrParams[0][2] << "\t" << linCorrParams[0][3]
	      << "\t" << linCorrParams[1][0] << "\t" << linCorrParams[1][1] << "\t" << linCorrParams[1][2] << "\t" << linCorrParams[1][3]
	      << "\t" << MeanAndSig1[0][0] << "\t" << MeanAndSig1[0][1] << "\t" 
	      << MeanAndSig1[1][0] << "\t" << MeanAndSig1[1][1] 
	      << "\t" << MeanAndSig2[0][0] << "\t" << MeanAndSig2[0][1] << "\t" 
	      << MeanAndSig2[1][0] << "\t" << MeanAndSig2[1][1] << std::endl;
      }
    }
    ofile.close();
  }  
  TFile *outfile;
  Char_t outputfile[500];
  
  if (printTree && source.substr(0,4)!="Beta") {
    //Create simulation output file    
    sprintf(outputfile,"analyzed_files/SimAnalyzed_%s_%s.root",source.c_str(),geometry.c_str());
    outfile = new TFile(outputfile, "RECREATE");
    Etype0.Write(); 
    Wtype0.Write();
    
    outfile->Close();
  }
  else if (source.substr(0,4)=="Beta")  {
    //Create simulation output file
    sprintf(outputfile,"analyzed_files/SimAnalyzed_%s.root",source.c_str());
    outfile = new TFile(outputfile, "RECREATE");
    tree->Write();
    outfile->Close();
  } 
}

  


//////////////////////////////////////////////////////////////////////////////////
// Function Definitions
//////////////////////////////////////////////////////////////////////////////////


//Function to return the x value of the max bin in a histogram 
Double_t GetXatMax(TH1D* hist, Double_t xmin, Double_t xmax)
{
  Double_t xBinAtMax = -1;
  //Int_t binValue = -1;

  if (xmin>0. && xmax>0.)  hist->GetXaxis()->SetRangeUser(xmin,xmax);
  xBinAtMax = hist->GetMaximumBin();
  
  Double_t xAtMax = hist -> GetBinCenter(xBinAtMax);
  hist->GetXaxis()->SetRange(0,hist->GetNbinsX());
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

bool checkParamSetStatus(std::vector <Double_t> params, std::string source, std::string geometry) 
{
  std::map <std::string,std::pair<Double_t,Double_t> > peaks;
  if (geometry=="2010") peaks = peaks2010;
  else if (geometry=="2011/2012") peaks = peaks2011;
  else if (geometry=="2012/2013") peaks = peaks2012;
  else {
    std::cout << "Bad geometry given to checkParamStatus!!!!\n";
    exit(0);
  }

  Double_t absSlope=0., absYIntercept=0.; //These are for the error envelope at this point... 
  for (auto it = envelopes.at(geometry).begin(); it!=envelopes.at(geometry).end(); it++) {
    if (it->first > peaks.at(source).first) {
      absSlope = (it->second-std::prev(it)->second)/(it->first-std::prev(it)->first);
      absYIntercept = it->second-absSlope*(it->first);
      //std::cout << "slope = " << absSlope << " y-intercept = " << absYIntercept << std::endl;
      break;
    }
  }
  
  Double_t resid = applyLinearityTwiddle(params,peaks.at(source).first) - peaks.at(source).first;
  if (resid<(absSlope*peaks.at(source).first+absYIntercept) && resid>(-absSlope*peaks.at(source).first-absYIntercept)) return true;
  else return false;

}

bool checkPeakStatus(std::vector < std::vector < Double_t > > meanAndSig, std::string source, std::string geometry) 
{
  if (meanAndSig[0][1]<(source=="Ce"?10.:20.) || meanAndSig[1][1]<(source=="Ce"?10.:20.)) return false; //Check that peak has real sigma

  std::map <std::string,std::pair<Double_t,Double_t> > peaks;
  if (geometry=="2010") peaks = peaks2010;
  else if (geometry=="2011/2012") peaks = peaks2011;
  else if (geometry=="2012/2013") peaks = peaks2012;
  else {
    std::cout << "Bad geometry given to checkPeakStatus!!!!\n";
    exit(0);
  }

  Double_t absSlope=0., absYIntercept=0.; //These are for the error envelope at this point... 
  for (auto it = envelopes.at(geometry).begin(); it!=envelopes.at(geometry).end(); it++) {
    if (it->first > peaks.at(source).first) {
      absSlope = (it->second-std::prev(it)->second)/(it->first-std::prev(it)->first);
      absYIntercept = it->second-absSlope*(it->first);
      //std::cout << "slope = " << absSlope << " y-intercept = " << absYIntercept << std::endl;
      break;
    }
  }
  
  Double_t residE = meanAndSig[0][0] - peaks.at(source).first;
  Double_t residW = meanAndSig[1][0] - peaks.at(source).second;
  if ((residE<(absSlope*peaks.at(source).first+absYIntercept) && residE>(-absSlope*peaks.at(source).first-absYIntercept)) && (residW<(absSlope*peaks.at(source).second+absYIntercept) && residW>(-absSlope*peaks.at(source).second-absYIntercept))) return true;
  
  else return false;
  
}

//Function to return the trigger function for each side in a std::vector in the form vec[side][param]
// where side==0 is East and side==1 is West
std::vector < std::vector < Double_t > > getTriggerFunctionParams(Int_t XeRunPeriod, Int_t nParams=8) {
  Char_t file[500];
  sprintf(file,"trigger_functions/trigger_functions_XePeriod_%i.dat",XeRunPeriod);
  ifstream infile(file);
  std::vector < std::vector <Double_t> > func;
  func.resize(2,std::vector <Double_t> (nParams,0.));
  //std::cout << "made it here\n";
  for (Int_t side = 0; side<2; side++) {
    Int_t param = 0;
    while (param<nParams) {
      infile >> func[side][param];
      //std::cout << func[side][param] << " ";
      param++;
    }
    //std::cout << std::endl;
  }
  infile.close();
  std::cout << "Loaded trigger functions.\n\n";
  return func;
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

Double_t triggerProbability(std::vector <Double_t> params, Double_t En) {
  //return params[0]+params[1]*TMath::Erf((En-params[2])/params[3])*TMath::TanH(params[7]*En) + params[4]*TMath::Gaus(En,params[5],params[6]);
  return (params[0]+params[1]*TMath::Erf((En-params[2])/params[3]))*(0.5-0.5*TMath::TanH((En-params[2])/params[4])) + 
    (0.5+0.5*TMath::TanH((En-params[2])/params[4]))*(params[5]+params[6]*TMath::TanH((En-params[2])/params[7]));
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
    std::cout << holdType << " " << params[side][type][0] << " " << params[side][type][1] << " " << params[side][type][2] << " " << params[side][type][3] << " " << params[side][type][4] << " " << params[side][type][5] << std::endl;
    type+=1;
    if (type==3) {type=0; side=1;}
  }
  return params;
}

//Function to return the correction to linearity offset
Double_t applyLinearityTwiddle (std::vector <Double_t> &params, Double_t EQ) {
  return params[0]+(1+params[1])*EQ+params[2]*EQ*EQ+params[3]*EQ*EQ*EQ;
}

// Sets up the output tree for the simulated data
void SetUpOutputTree(TTree& tree) {
  tree.Branch("PID", &PID, "PID/I");
  tree.Branch("primMomentum", &initialMomentum, "primMo[3]/D");
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
  
}
  

