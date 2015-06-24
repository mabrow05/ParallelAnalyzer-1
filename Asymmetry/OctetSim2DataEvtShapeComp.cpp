//This root script loads the replayed data files for a given octet and produces
// histograms of the different event types. Then it also takes the reverse
// calibrated simulation files for each run, weighs them by the physics 
// asymmetry, and normalizes the Type 0 events before producing histograms of
// the simulated event types

#include "OctetSim2DataEvtShapeComp.hh"

Sim2Data::Sim2Data(int r, double fid)
{
  run = r;
  fiducialCut = fid;
  char temp[200];
  sprintf(temp," ");
  outfile = new TFile(temp);
};

Sim2Data::~Sim2Data(){outfile->Write(); outfile->Close();};

void Sim2Data::histoCreator()
{
  eastData.resize(3,0);
  westData.resize(3,0);
  eastSim.resize(3,0);
  westSim.resize(3,0);
  
  eastData[0] = new TH1D("eastData0", "Type 0 Data East", 100, 0., 800.);
  westData[0] = new TH1D("westData0", "Type 0 Data West", 100, 0., 800.);
  eastData[1] = new TH1D("eastData1", "Type 1 Data East", 100, 0., 800.);
  westData[1] = new TH1D("westData1", "Type 1 Data West", 100, 0., 800.);
  eastData[2] = new TH1D("eastData23", "Type 2/3 Data East", 100, 0., 800.);
  westData[2] = new TH1D("westData23", "Type 2/3 Data West", 100, 0., 800.);
  eastSim[0] = new TH1D("eastSim0", "Type 0 Sim East", 100, 0., 800.);
  westSim[0] = new TH1D("westSim0", "Type 0 Sim West", 100, 0., 800.);
  eastSim[1] = new TH1D("eastSim1", "Type 1 Sim East", 100, 0., 800.);
  westSim[1] = new TH1D("westSim1", "Type 1 Sim West", 100, 0., 800.);
  eastSim[2] = new TH1D("eastSim23", "Type 2/3 Sim East", 100, 0., 800.);
  westSim[2] = new TH1D("westSim23", "Type 2/3 Sim West", 100, 0., 800.);
};

void Sim2Data::dataReader()
{
  
  Char_t temp[200];
  Char_t *dataDir = getenv("UCNAOUTPUTDIR");
 
  cout << "Reading in data file" << endl;
  //Read in file  
  sprintf(temp, "%s/hists/spec_%i.root",dataDir, run);
  TFile *f1 = new TFile(temp,"READ");
  TTree *t1 = (TTree*)(f1->Get("phys"));
  TH1D *hold = new TH1D("hold", "hold", 100, 0., 800.); //histogram to hold the data prior to adding to
  
  //Begin by extracting the number of type 0 events from the Type Branch. Then draw the energy 
  // distribution and extract the integral for possible use in normalization of the simulated files
  sprintf(temp,"Type==0 && PID==1 && Side==0 && (xEmpm.center*xEmpm.center+yEmpm.center*yEmpm.center)<(%f*%f)",fiducialCut,fiducialCut);
  t1->Draw("Erecon>>hold",temp,"goff");
  sprintf(temp,"Type==0 && PID==1 && Side==1 && (xWmpm.center*xWmpm.center+yWmpm.center*yWmpm.center)<(%f*%f)",fiducialCut,fiducialCut);
  t1->Draw("Erecon>>+hold",temp,"goff");
  type0evts = hold->GetEntries();
  cout << type0evts << endl;
  type0integral = hold->Integral(0,100,"width");
  cout << type0integral << endl;

  //Draw from the TTree the type 0 events to histogram called hold and add them to eastData[0] and westData[0], etc
  sprintf(temp,"Type==0 && PID==1 && Side==0 && (xEmpm.center*xEmpm.center+yEmpm.center*yEmpm.center)<(%f*%f)",fiducialCut,fiducialCut);
  t1->Draw("Erecon>>hold",temp,"goff");
  Int_t eastEvts=hold->GetEntries();
  eastData[0]->Add(hold);
  sprintf(temp,"Type==0 && PID==1 && Side==1 && (xWmpm.center*xWmpm.center+yWmpm.center*yWmpm.center)<(%f*%f)",fiducialCut,fiducialCut);
  t1->Draw("Erecon>>hold",temp,"goff");
  Int_t westEvts=hold->GetEntries();
  westData[0]->Add(hold);
  sprintf(temp,"Type==1 && PID==1 && Side==0 && (xEmpm.center*xEmpm.center+yEmpm.center*yEmpm.center)<(%f*%f)",fiducialCut,fiducialCut);
  t1->Draw("Erecon>>hold",temp,"goff");
  eastData[1]->Add(hold);
  sprintf(temp,"Type==1 && PID==1 && Side==1 && (xWmpm.center*xWmpm.center+yWmpm.center*yWmpm.center)<(%f*%f)",fiducialCut,fiducialCut);
  t1->Draw("Erecon>>hold",temp,"goff");
  westData[1]->Add(hold);
  sprintf(temp,"(Type==2 || Type==3) && PID==1 && Side==0 && (xEmpm.center*xEmpm.center+yEmpm.center*yEmpm.center)<(%f*%f)",fiducialCut,fiducialCut);
  t1->Draw("Erecon>>hold",temp,"goff");
  eastData[2]->Add(hold);
  sprintf(temp,"(Type==2 || Type==3) && PID==1 && Side==1 && (xWmpm.center*xWmpm.center+yWmpm.center*yWmpm.center)<(%f*%f)",fiducialCut,fiducialCut);
  t1->Draw("Erecon>>hold",temp,"goff");
  westData[2]->Add(hold);

  if (hold) delete hold;
  f1->Close();
  if (f1) delete f1;
  if (t1) delete t1;

  //determining polarization of data...
  pol = eastEvts>westEvts?1.:-1.;
      
};

void Sim2Data::reverseCalSims()
{
  Char_t temp[200];
  cout << "Reading in Simulated files" << endl;
      //Read in file
  sprintf(temp,"/extern/mabrow05/ucna/SimAnalysisOutput/hists/Octet_24_18433-18453/%i_sim_fid%0.f.root",run,fiducialCut);
  TFile *f1 = new TFile(temp,"READ");
  TTree *t1 = (TTree*)(f1->Get("SimOut"));
  cout << "Read in file for run " << run << endl;
  TH1D *hold = new TH1D("hold", "hold", 100, 0., 800.); //histogram to hold the data prior to adding to
  
  //Extract the number of type 0 events and calculate the normalization factor, eta, based on the ratio of 
  // data/sim type 0 events
  t1->Draw("Erecon>>hold","physicsWeight*(type==0)","goff");
  //sprintf(temp,"(1-%f*0.12*sqrt(1-(1/((ePrim/500.+1.)*(ePrim/500.+1.))))*costheta)*(type==0)",pol);
  //t1->Draw("Erecon>>hold",temp);
  Double_t etaE = (double)type0evts/((double)(hold->GetEntries()));
  Double_t etaIntegralsE = type0integral/(hold->Integral(0,100,"width"));
  Double_t etaW = (double)type0evts/((double)(hold->GetEntries()));     
  Double_t etaIntegralsW = type0integral/(hold->Integral(0,100,"width"));
  //Draw the appropriate type of events, weigh them, and then normalize them, before adding them to the 
  // correct histogram from above
  
  t1->Draw("Erecon>>hold","physicsWeight*(type==0 && side==0)");
  eastSim[0]->Add(hold,etaIntegralsE);
  
  t1->Draw("Erecon>>hold","physicsWeight*(type==0 && side==1)");
  westSim[0]->Add(hold, etaIntegralsW);
  
  t1->Draw("Erecon>>hold","physicsWeight*(type==1 && side==0)");
  eastSim[1]->Add(hold,etaIntegralsE);
  t1->Draw("Erecon>>hold","physicsWeight*(type==1 && side==1)");
  westSim[1]->Add(hold,etaIntegralsW);
  t1->Draw("Erecon>>hold","physicsWeight*((type==2 || type==3) && side==0)");
  eastSim[2]->Add(hold,etaIntegralsE);
  t1->Draw("Erecon>>hold","physicsWeight*((type==2 || type==3) && side==1)");
  westSim[2]->Add(hold,etaIntegralsW);
  
  if (hold) delete hold;
  f1->Close();
  if (f1) delete f1;
  if (t1) delete t1;
};

int main(int argc, char *argv[])
{

  if (argc < 3)
    {
      cout << "Usage: ./OctetSim2DataShapeComp [run] [fiducialCut]\n"; 
      return -1;
    }
  Char_t t[100];
  for (int i = 0; i<argc; i++)
    {
      sprintf(t, "Argument %i: ", i);
      cout << t << argv[i] << endl;
    }

  Sim2Data *S2D = new Sim2Data(atoi(argv[1]),atof(argv[2]));
  S2D->histoCreator();
  S2D->dataReader();
  S2D->reverseCalSims();

  if (S2D) delete S2D;
}
  

