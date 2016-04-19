#include "Etrue_Edep.hh"
#include "cstdlib"
#include "fstream"
#include "iomanip"
#include "math.h"
#include "iostream"
#include "string"
#include "vector"

using std::ifstream;
using std::ofstream;
using namespace std;

DataTree::DataTree() : inputFile(NULL), inputTree(NULL) { };

DataTree::~DataTree() {


  if (inputTree) delete inputTree;
  if (inputFile) { inputFile->Close(); delete inputFile;}
};

void DataTree::setupInputTree(string inputFileName, string inputTreeName) {

  inputFile = new TFile(inputFileName.c_str(),"READ");
  inputTree = (TTree*)(inputFile->Get(inputTreeName.c_str()));
  
  inputTree->SetBranchAddress("primKE", &primKE);
  inputTree->SetBranchAddress("primTheta", &primTheta);
  inputTree->SetBranchAddress("Edep", &Edep);
  inputTree->SetBranchAddress("EdepQ", &EdepQ); 
  inputTree->SetBranchAddress("MWPCPos", &MWPCPos); 
  inputTree->SetBranchAddress("MWPCEnergy", &MWPCEnergy); 
  inputTree->SetBranchAddress("time", &time); 

  cout << "Prepared input tree " << inputTreeName << " in " << inputFileName << "...\n";
};


int main(int argc, char *argv[]) {

  TFile *outputHists = new TFile("Hists.root","RECREATE");

  ofstream outdata;
  outdata.open ("HistMeans.dat");
  outdata << setw(5) << "Nhist" << setw(7) << "Emin" << setw(7) << "Emax" 
                     << setw(11) << "EQType0E" << setw(7) << "Neve" 
                     << setw(11) << "EQType0W" << setw(7) << "Neve"  
                     << setw(11) << "EQType1E" << setw(7) << "Neve" 
                     << setw(11) << "EQType1W" << setw(7) << "Neve"  
                     << setw(11) << "EQType2E" << setw(7) << "Neve" 
                     << setw(11) << "EQType2W" << setw(7) << "Neve" 
                     << setw(11) << "EQType3E" << setw(7) << "Neve" 
                     << setw(11) << "EQType3W" << setw(7) << "Neve"
                     << setw(12) << "EQType23E" << setw(7) << "Neve" 
                     << setw(12) << "EQType23W" << setw(7) << "Neve" << endl; 

  char temp[100];
  int nFile = 1000;
  double MaxE = 800.0;
  double DeltaE = 10.0;
  #define pi M_PI

  vector<double> Emin;
  vector<double> Emax;

  vector<TH1D*> EDepQType0E;
  vector<TH1D*> EDepQType1E;
  vector<TH1D*> EDepQType2E;
  vector<TH1D*> EDepQType3E;
  vector<TH1D*> EDepQType23E;
  vector<TH1D*> EDepQType0W;
  vector<TH1D*> EDepQType1W;
  vector<TH1D*> EDepQType2W;
  vector<TH1D*> EDepQType3W;
  vector<TH1D*> EDepQType23W;

  int Nhist = 0;
  for (double E=DeltaE/2.; E<=MaxE; E+=DeltaE) Nhist++;
  int MaxNhist = Nhist;

  EDepQType0E.resize(MaxNhist);
  EDepQType1E.resize(MaxNhist);
  EDepQType2E.resize(MaxNhist);
  EDepQType3E.resize(MaxNhist);
  EDepQType23E.resize(MaxNhist);
  EDepQType0W.resize(MaxNhist);
  EDepQType1W.resize(MaxNhist);
  EDepQType2W.resize(MaxNhist);
  EDepQType3W.resize(MaxNhist);
  EDepQType23W.resize(MaxNhist);

  Nhist = 0;
  for (double E=DeltaE/2.; E<=MaxE; E+=DeltaE) 
      {
       Emin.push_back(E-DeltaE/2.);
       Emax.push_back(E+DeltaE/2.);

       sprintf(temp, "EDepQType0E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType0E[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 0 East", 200, 0. , 800.);

       sprintf(temp, "EDepQType0W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType0W[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 0 West", 200, 0. , 800.);

       sprintf(temp, "EDepQType1E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType1E[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 1 East", 200, 0. , 800.);

       sprintf(temp, "EDepQType1W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType1W[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 1 West", 200, 0. , 800.);

       sprintf(temp, "EDepQType2E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType2E[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 2 East", 200, 0. , 800.);

       sprintf(temp, "EDepQType2W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType2W[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 2 West", 200, 0. , 800.);

       sprintf(temp, "EDepQType3E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType3E[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 3 East", 200, 0. , 800.);

       sprintf(temp, "EDepQType3W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType3W[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 3 West", 200, 0. , 800.);

       sprintf(temp, "EDepQType23E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType23E[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 23 East", 200, 0. , 800.);

       sprintf(temp, "EDepQType23W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType23W[Nhist] = new TH1D(temp, "Energy Deposition Quenched, TYPE 23 West", 200, 0. , 800.);

       Nhist++;
       }
 
  double FidRad = 50.0;
  char path[] = "/extern/mabrow05/ucna/geant4work/output/bigSim/flatField_2011-2012/F0/";//"/extern/mabrow05/ucna/geant4work/output/10mil_2011-2012/Beta/"; 
  for (int nF = 0; nF<nFile; nF++) 
      {
       DataTree Tin;
       sprintf(temp, "%s/analyzed_%d.root", path, nF);
       Tin.setupInputTree (temp, "anaTree");
       int Nevt = Tin.getEntries();
       for (Int_t evt = 0; evt<Nevt; evt++) 
           {
            Tin.getEvent(evt);
            if (Tin.MWPCPos.MWPCPosE[0]*Tin.MWPCPos.MWPCPosE[0]+Tin.MWPCPos.MWPCPosE[1]*Tin.MWPCPos.MWPCPosE[1] > FidRad*FidRad) continue; 

            for (Nhist=0; Nhist<MaxNhist; Nhist++) //Sets Nhist to the proper histogram based on the true energy of the event
                 if (Emin[Nhist]<=Tin.primKE && Tin.primKE<Emax[Nhist]) break; 

            //TYPE 0E EVENTS
            if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.Edep.EdepE > 0 &&
                Tin.MWPCEnergy.MWPCEnergyW == 0 && Tin.Edep.EdepW == 0) 
                EDepQType0E[Nhist]->Fill(Tin.EdepQ.EdepQE);

            //TYPE 0W EVENTS
            else if (Tin.MWPCEnergy.MWPCEnergyE == 0 && Tin.Edep.EdepE == 0 &&
                     Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.Edep.EdepW > 0)
                     EDepQType0W[Nhist]->Fill(Tin.EdepQ.EdepQW);

            //TYPE 1E & 1W EVENTS
            else if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.Edep.EdepE > 0 &&
                     Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.Edep.EdepW > 0)
                     {
                      if (Tin.time.timeE < Tin.time.timeW) EDepQType1E[Nhist]->Fill(Tin.EdepQ.EdepQE+Tin.EdepQ.EdepQW);
                      else EDepQType1W[Nhist]->Fill(Tin.EdepQ.EdepQW+Tin.EdepQ.EdepQE);
                      }

            //TYPE 3E & 2E EVENTS
            else if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.Edep.EdepE > 0 &&
                     Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.Edep.EdepW == 0)
                     {
                      EDepQType23E[Nhist]->Fill(Tin.EdepQ.EdepQE);
                      if (Tin.primTheta>pi/2.) EDepQType3E[Nhist]->Fill(Tin.EdepQ.EdepQE);
                      else EDepQType2E[Nhist]->Fill(Tin.EdepQ.EdepQE);
                      }

            //TYPE 2W & 3W EVENTS
            else if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.Edep.EdepE == 0 &&
                     Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.Edep.EdepW > 0)
                     {
                      EDepQType23W[Nhist]->Fill(Tin.EdepQ.EdepQW);
                      if (Tin.primTheta<pi/2) EDepQType3W[Nhist]->Fill(Tin.EdepQ.EdepQW);
                      else EDepQType2W[Nhist]->Fill(Tin.EdepQ.EdepQW);
                      }

            }
       cout << "******* analyzed_" << nF << ".root done *******" << endl;
       cout << endl;
       }

  for (Nhist=0; Nhist<MaxNhist; Nhist++)
      {
       outdata << setw(5)  << Nhist << setw(7) << Emin[Nhist] << setw(7) << Emax[Nhist]
               << setw(11) << EDepQType0E[Nhist]->GetMean() << setw(7) << EDepQType0E[Nhist]->GetEntries() 
               << setw(11) << EDepQType0W[Nhist]->GetMean() << setw(7) << EDepQType0W[Nhist]->GetEntries()
               << setw(11) << EDepQType1E[Nhist]->GetMean() << setw(7) << EDepQType1E[Nhist]->GetEntries() 
               << setw(11) << EDepQType1W[Nhist]->GetMean() << setw(7) << EDepQType1W[Nhist]->GetEntries()
               << setw(11) << EDepQType2E[Nhist]->GetMean() << setw(7) << EDepQType2E[Nhist]->GetEntries() 
               << setw(11) << EDepQType2W[Nhist]->GetMean() << setw(7) << EDepQType2W[Nhist]->GetEntries()
               << setw(11) << EDepQType3E[Nhist]->GetMean() << setw(7) << EDepQType3E[Nhist]->GetEntries() 
               << setw(11) << EDepQType3W[Nhist]->GetMean() << setw(7) << EDepQType3W[Nhist]->GetEntries()  
               << setw(12) << EDepQType23E[Nhist]->GetMean() << setw(7) << EDepQType23E[Nhist]->GetEntries()
               << setw(12) << EDepQType23W[Nhist]->GetMean() << setw(7) << EDepQType23W[Nhist]->GetEntries() << endl;  
       }
  
  outputHists->Write();
  outputHists->Close();
}
  
