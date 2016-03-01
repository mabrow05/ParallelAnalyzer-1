#include "Etrue_Edep.hh"
#include "cstdlib"
#include "fstream"
#include "iomanip"
#include "math.h"
#include "iostream"
#include "string"

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
  int nFile = 3;
  double MaxE = 800.0;
  double DeltaE = 100.0;
  #define pi M_PI

  TObjArray HistList(0);          // create an array of Histograms
  TH1D* EDepQType0E;              // create a pointer to a histogram
  TH1D* EDepQType0W;              // create a pointer to a histogram
  TH1D* EDepQType1E;              // create a pointer to a histogram
  TH1D* EDepQType1W;              // create a pointer to a histogram
  TH1D* EDepQType2E;              // create a pointer to a histogram
  TH1D* EDepQType2W;              // create a pointer to a histogram
  TH1D* EDepQType3E;              // create a pointer to a histogram
  TH1D* EDepQType3W;              // create a pointer to a histogram
  TH1D* EDepQType23E;             // create a pointer to a histogram
  TH1D* EDepQType23W;             // create a pointer to a histogram

  int Nhist = 0;
  for (double E=DeltaE/2.; E<=MaxE; E+=DeltaE) 
      {
       sprintf(temp, "EDepQType0E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType0E = new TH1D(temp, "Energy Deposition Quenched, TYPE 0 East", 200, 0. , 800.);
       HistList.Add(EDepQType0E);

       sprintf(temp, "EDepQType0W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType0W = new TH1D(temp, "Energy Deposition Quenched, TYPE 0 West", 200, 0. , 800.);
       HistList.Add(EDepQType0W);

       sprintf(temp, "EDepQType1E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType1E = new TH1D(temp, "Energy Deposition Quenched, TYPE 1 East", 200, 0. , 800.);
       HistList.Add(EDepQType1E);

       sprintf(temp, "EDepQType1W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType1W = new TH1D(temp, "Energy Deposition Quenched, TYPE 1 West", 200, 0. , 800.);
       HistList.Add(EDepQType1W);

       sprintf(temp, "EDepQType2E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType2E = new TH1D(temp, "Energy Deposition Quenched, TYPE 2 East", 200, 0. , 800.);
       HistList.Add(EDepQType2E);

       sprintf(temp, "EDepQType2W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType2W = new TH1D(temp, "Energy Deposition Quenched, TYPE 2 West", 200, 0. , 800.);
       HistList.Add(EDepQType2W);

       sprintf(temp, "EDepQType3E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType3E = new TH1D(temp, "Energy Deposition Quenched, TYPE 3 East", 200, 0. , 800.);
       HistList.Add(EDepQType3E);

       sprintf(temp, "EDepQType3W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType3W = new TH1D(temp, "Energy Deposition Quenched, TYPE 3 West", 200, 0. , 800.);
       HistList.Add(EDepQType3W);

       sprintf(temp, "EDepQType23E_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType23E = new TH1D(temp, "Energy Deposition Quenched, TYPE 23 East", 200, 0. , 800.);
       HistList.Add(EDepQType23E);

       sprintf(temp, "EDepQType23W_%.2f_%.2f", E-DeltaE/2., E+DeltaE/2.);
       EDepQType23W = new TH1D(temp, "Energy Deposition Quenched, TYPE 23 West", 200, 0. , 800.);
       HistList.Add(EDepQType23W);


       double FidRad = 50.0;
       char path[] = "/home/nima/nima/UCNA/Etrue_Edep/optimized_Etrue_Edep_code/"; 
       for (int nF = 0; nF<nFile; nF++) 
           {
            DataTree Tin;
            sprintf(temp, "%s/analyzed_%d.root", path, nF);
            Tin.setupInputTree (temp, "anaTree");
            int Nevt = Tin.getEntries();
            for (Int_t evt = 0; evt<Nevt; evt++) 
                {
                 Tin.getEvent(evt);

                 if (Tin.primKE < E-DeltaE/2. || Tin.primKE > E+DeltaE/2.) continue;
                 if (Tin.MWPCPos.MWPCPosE[0]*Tin.MWPCPos.MWPCPosE[0]+Tin.MWPCPos.MWPCPosE[1]*Tin.MWPCPos.MWPCPosE[1] > FidRad*FidRad) continue; 


                 //TYPE 0E EVENTS
                 if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.Edep.EdepE > 0 &&
                     Tin.MWPCEnergy.MWPCEnergyW == 0 && Tin.Edep.EdepW == 0) 
                     EDepQType0E->Fill(Tin.EdepQ.EdepQE);

                 //TYPE 0W EVENTS
                 else if (Tin.MWPCEnergy.MWPCEnergyE == 0 && Tin.Edep.EdepE == 0 &&
                          Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.Edep.EdepW > 0)
                          EDepQType0W->Fill(Tin.EdepQ.EdepQW);

                 //TYPE 1E & 1W EVENTS
                 else if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.Edep.EdepE > 0 &&
                          Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.Edep.EdepW > 0)
                          {
                           if (Tin.time.timeE < Tin.time.timeW) EDepQType1E->Fill(Tin.EdepQ.EdepQE);
                           else EDepQType1W->Fill(Tin.EdepQ.EdepQW);
                           }

                 //TYPE 3E & 2E EVENTS
                 else if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.Edep.EdepE > 0 &&
                          Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.Edep.EdepW == 0)
                          {
                           EDepQType23E->Fill(Tin.EdepQ.EdepQE);
                           if (Tin.primTheta>pi/2.) EDepQType3E->Fill(Tin.EdepQ.EdepQE);
                           else EDepQType2E->Fill(Tin.EdepQ.EdepQE);
                           }

                 //TYPE 2W & 3W EVENTS
                 else if (Tin.MWPCEnergy.MWPCEnergyE > 0 && Tin.Edep.EdepE == 0 &&
                          Tin.MWPCEnergy.MWPCEnergyW > 0 && Tin.Edep.EdepW > 0)
                          {
                           EDepQType23W->Fill(Tin.EdepQ.EdepQW);
                           if (Tin.primTheta<pi/2) EDepQType3W->Fill(Tin.EdepQ.EdepQW);
                           else EDepQType2W->Fill(Tin.EdepQ.EdepQW);
                          }
                 }
            cout << "******* analyzed_" << nF << ".root done *******" << endl;
            cout << endl;
            }
       outdata << setw(5)  << Nhist << setw(7) << E-DeltaE/2. << setw(7) << E+DeltaE/2.
               << setw(11) << EDepQType0E->GetMean() << setw(7) << EDepQType0E->GetEntries() 
               << setw(11) << EDepQType0W->GetMean() << setw(7) << EDepQType0W->GetEntries()
               << setw(11) << EDepQType1E->GetMean() << setw(7) << EDepQType1E->GetEntries() 
               << setw(11) << EDepQType1W->GetMean() << setw(7) << EDepQType1W->GetEntries()
               << setw(11) << EDepQType2E->GetMean() << setw(7) << EDepQType2E->GetEntries() 
               << setw(11) << EDepQType2W->GetMean() << setw(7) << EDepQType2W->GetEntries()
               << setw(11) << EDepQType3E->GetMean() << setw(7) << EDepQType3E->GetEntries() 
               << setw(11) << EDepQType3W->GetMean() << setw(7) << EDepQType3W->GetEntries()  
               << setw(12) << EDepQType23E->GetMean() << setw(7) << EDepQType23E->GetEntries()
               << setw(12) << EDepQType23W->GetMean() << setw(7) << EDepQType23W->GetEntries() << endl;  

       cout << "Nhist= " << Nhist << ", Range " << E-DeltaE/2. << " - " << E+DeltaE/2. << " done" <<endl<<endl;

       Nhist++;
       }

  TFile *outputHists = new TFile("Hists.root","RECREATE");
  HistList.Write(); 
  //outputHists->Write();
  outputHists->Close();
}
  
