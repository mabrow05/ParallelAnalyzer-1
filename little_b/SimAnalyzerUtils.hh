// Utility code holding all pertinent functions, variables, and structs for SimulationAnalyzer


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


// This is the default error envelope as given by Michael M. (2008,2010) and Michael B. (2011,2012)

std::vector < std::pair <Double_t,Double_t> > envelope2010 = {std::make_pair(0,2.5),std::make_pair(200,200*0.0125),std::make_pair(500,500*0.0125),std::make_pair(1000,500*0.0125)};
std::vector < std::pair <Double_t,Double_t> > envelope2011 = {std::make_pair(0,0.025*130.3),std::make_pair(130.3,130.3*0.025),std::make_pair(368.4938,0.018*368.4938),std::make_pair(993.789,993.789*0.013),std::make_pair(1000,1000.*0.013)};
std::vector < std::pair <Double_t,Double_t> > envelope2012 = {std::make_pair(0,0.025*130.3),std::make_pair(130.3,130.3*0.025),std::make_pair(368.4938,0.018*368.4938),std::make_pair(993.789,993.789*0.013),std::make_pair(1000,1000.*0.013)};

std::map <std::string, std::vector < std::pair<Double_t, Double_t> > > envelopes = {{"2010",envelope2010},{"2011/2012",envelope2011},{"2012/2013",envelope2012}};


// These are the nominal peak values to compare against for each source
std::map <std::string,std::pair<Double_t,Double_t> > peaks2010 = {{"Sn",std::make_pair(365.629,365.394)},{"Ce",std::make_pair(0.,0.)},{"Bi1",std::make_pair(0.,0.)},{"Bi2",std::make_pair(0.,0.)}};
std::map <std::string,std::pair<Double_t,Double_t> > peaks2011 = {{"Sn",std::make_pair(365.629,365.394)},{"Ce",std::make_pair(0.,0.)},{"Bi1",std::make_pair(0.,0.)},{"Bi2",std::make_pair(0.,0.)}};
std::map <std::string,std::pair<Double_t,Double_t> > peaks2012 = {{"Sn",std::make_pair(365.629,365.394)},{"Ce",std::make_pair(0.,0.)},{"Bi1",std::make_pair(0.,0.)},{"Bi2",std::make_pair(0.,0.)}};

//std::map <std::string, std::vector < std::pair<Double_t, Double_t> > > nominalPeaks;
//envelope["2010"] = peaks2010;
//envelope["2011"] = peaks2011;
//envelope["2012"] = peaks2012;

Int_t PID, type, side, primaryID; // basic analysis tags

Double_t initialMomentum[3]; //For reconstructing which detector type 1 events initially struck

Double_t Eprim, AsymWeight; // initial energy from simulation, and weight of event based on A*beta*cos(theta)

Double_t Erecon; // smeared, weighted, and trigger func corrected energy with the EQ2Etrue conversion applied


//Nominal peak values for simulated data which has detector effects applied
// TODO


///////// Function declarations /////////////////

//Function to return the x value of the max bin in a histogram 
Double_t GetXatMax(TH1D* hist);

//Function to make initial check on state of parameter set.. i.e. directly comparing against error envelope
bool checkParamSetStatus(std::vector <Double_t> params); // TODO

//Function to return a 2 element vector holding the fit mean and sigma of a peak
std::vector <Double_t> FitGaus(TH1D* histToFit, Double_t gausMean, Double_t min, Double_t max);

//Function to check if the peak value for a source falls within the error envelope
bool CheckPeakValues2010(std::vector <int> parameters, std::string sourceName, std::string side);
bool CheckPeakValues2011(std::vector <int> parameters, std::string sourceName, std::string side);
bool CheckPeakValues2012(std::vector <int> parameters, std::string sourceName, std::string side);

//Function to return the trigger function for each side in a std::vector in the form vec[side][param]
// where side==0 is East and side==1 is West
std::vector < std::vector < Double_t > > getTriggerFunctionParams(Int_t XeRunPeriod, Int_t nParams = 7);

//Function to return the probability of a trigger at a given energy given the parameters
Double_t triggerProbability(std::vector <Double_t> params, Double_t En) {
  return params[0]+params[1]*TMath::Erf((En-params[2])/params[3])
    + params[4]*TMath::Gaus(En,params[5],params[6]);
}

//Get the conversion from EQ2Etrue
std::vector < std::vector < std::vector <double> > > getEQ2EtrueParams(std::string geometry);
  
//Function to return the correction to linearity offset
Double_t applyLinearityTwiddle (std::vector <Double_t> &params, Double_t EQ) {
  return params[0]+(1+params[1])*EQ+params[2]*EQ*EQ+params[3]*EQ*EQ*EQ;
}

// Sets up the output tree for the simulated data
void SetUpOutputTree(TTree& tree);
  

//////////////////////// Data Structures /////////////////////////

struct Evis {
  double EvisE;
  double EvisW;
} evis;

struct Edep {
  double EdepE;
  double EdepW;
} edep;

struct EdepQ {
  double EdepQE;
  double EdepQW;
} edepQ;

struct time {
  double timeE;
  double timeW;
} Time;

struct MWPCEnergy {
  double MWPCEnergyE; 
  double MWPCEnergyW;
} mwpcE;

struct MWPCPos {
  double MWPCPosE[3];
  double MWPCPosW[3];
} mwpc_pos;

struct ScintPos {
  double ScintPosE[3];
  double ScintPosW[3];
} scint_pos;

struct ScintPosAdjusted {
  double ScintPosAdjE[3];
  double ScintPosAdjW[3];
} scint_pos_adj;

struct PMT_Evis {
  double Evis[8];
  double weight[8]; 
} pmt_Evis;




//////////////////////////////////////////////////////////////////////////////////
// Function Definitions
//////////////////////////////////////////////////////////////////////////////////

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

bool CheckPeakValues2011(std::vector <int> parameters, std::string sourceName, std::string side)
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
    std::cout << "Set of parameters are GOOD: difference in mean is " << abs(diffMean)
    << ", difference in sigma is " << abs(diffSigma) << std::endl;


    return true;
  }
  if(diffMean > 0)
  {
    std::cout << "Parameters BAD: mean is outside range." << std::endl;
  }
  if(diffSigma > 0)
  {
    std::cout << "Parameters BAD: sigma is outside range." << std::endl;
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
    std::cout << "Set of parameters are GOOD: difference in mean is " << abs(diffMean)
    << ", difference in sigma is " << abs(diffSigma) << std::endl;


    return true;
  }
  if(diffMean > 0)
  {
    std::cout << "Parameters BAD: mean is outside range." << std::endl;
  }
  if(diffSigma > 0)
  {
    std::cout << "Parameters BAD: sigma is outside range." << std::endl;
  }

  return false;
}


//Function to return the trigger function for each side in a std::vector in the form vec[side][param]
// where side==0 is East and side==1 is West
std::vector < std::vector < Double_t > > getTriggerFunctionParams(Int_t XeRunPeriod, Int_t nParams) {
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
  
}
  
