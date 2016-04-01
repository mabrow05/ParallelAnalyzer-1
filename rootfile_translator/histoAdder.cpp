/*////////////////////////////////////////////////////////////////////
This code appends certain histograms to the final replay file. Michael
used these histograms for parts of his analysis and they weren't initially 
included in the UK analysis. 

Through the translator, these will be added to the spec files as well
*/////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

// ROOT libraries
#include "TRandom3.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH1D.h>

#include "fullTreeVariables.h"
//#include "MWPCGeometry.h"
//#include "pedestals.h"
//#include "cuts.h"
//#include "basic_reconstruction.h"
//#include "WireChamberResponse.h"
//#include "DataTree.hh"



int main(int argc, char *argv[]) {

  if (argc!=2) {
    std::cout << "Usage: ./histoAdder runNumber\n";
    exit(0);
  }
  
  int runNumber = atoi(argv[1]);
  std::cout << "Running histoAdder on run " << runNumber << std::endl;

  //For adding the UCN monitor rates, we need to open the raw data files.
  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "/extern/mabrow05/ucna/rawdata/full%s.root", runNumber);

  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("h1"));

  // Variables
  Tin->SetBranchAddress("Pdc30",  &Pdc30);
  Tin->SetBranchAddress("Pdc34",  &Pdc34);

  Tin->SetBranchAddress("Tdc016", &Tdc016);
  Tin->SetBranchAddress("Tdc017", &Tdc017);
  Tin->SetBranchAddress("Tdc00", &Tdc00);
  Tin->SetBranchAddress("Tdc01", &Tdc01);
  Tin->SetBranchAddress("Tdc02", &Tdc02);
  Tin->SetBranchAddress("Tdc03", &Tdc03);
  Tin->SetBranchAddress("Tdc08", &Tdc08);
  Tin->SetBranchAddress("Tdc09", &Tdc09);
  Tin->SetBranchAddress("Tdc010", &Tdc010);
  Tin->SetBranchAddress("Tdc011", &Tdc011);

  Tin->SetBranchAddress("Sis00",  &Sis00);

  Tin->SetBranchAddress("Qadc0",  &Qadc[0]);
  Tin->SetBranchAddress("Qadc1",  &Qadc[1]);
  Tin->SetBranchAddress("Qadc2",  &Qadc[2]);
  Tin->SetBranchAddress("Qadc3",  &Qadc[3]);
  Tin->SetBranchAddress("Qadc4",  &Qadc[4]);
  Tin->SetBranchAddress("Qadc5",  &Qadc[5]);
  Tin->SetBranchAddress("Qadc6",  &Qadc[6]);
  Tin->SetBranchAddress("Qadc7",  &Qadc[7]);

  Tin->SetBranchAddress("Pdc20",  &Pdc2[0]);
  Tin->SetBranchAddress("Pdc21",  &Pdc2[1]);
  Tin->SetBranchAddress("Pdc22",  &Pdc2[2]);
  Tin->SetBranchAddress("Pdc23",  &Pdc2[3]);
  Tin->SetBranchAddress("Pdc24",  &Pdc2[4]);
  Tin->SetBranchAddress("Pdc25",  &Pdc2[5]);
  Tin->SetBranchAddress("Pdc26",  &Pdc2[6]);
  Tin->SetBranchAddress("Pdc27",  &Pdc2[7]);
  Tin->SetBranchAddress("Pdc28",  &Pdc2[8]);
  Tin->SetBranchAddress("Pdc29",  &Pdc2[9]);
  Tin->SetBranchAddress("Pdc210", &Pdc2[10]);
  Tin->SetBranchAddress("Pdc211", &Pdc2[11]);
  Tin->SetBranchAddress("Pdc212", &Pdc2[12]);
  Tin->SetBranchAddress("Pdc213", &Pdc2[13]);
  Tin->SetBranchAddress("Pdc214", &Pdc2[14]);
  Tin->SetBranchAddress("Pdc215", &Pdc2[15]);
  Tin->SetBranchAddress("Pdc216", &Pdc2[16]);
  Tin->SetBranchAddress("Pdc217", &Pdc2[17]);
  Tin->SetBranchAddress("Pdc218", &Pdc2[18]);
  Tin->SetBranchAddress("Pdc219", &Pdc2[19]);
  Tin->SetBranchAddress("Pdc220", &Pdc2[20]);
  Tin->SetBranchAddress("Pdc221", &Pdc2[21]);
  Tin->SetBranchAddress("Pdc222", &Pdc2[22]);
  Tin->SetBranchAddress("Pdc223", &Pdc2[23]);
  Tin->SetBranchAddress("Pdc224", &Pdc2[24]);
  Tin->SetBranchAddress("Pdc225", &Pdc2[25]);
  Tin->SetBranchAddress("Pdc226", &Pdc2[26]);
  Tin->SetBranchAddress("Pdc227", &Pdc2[27]);
  Tin->SetBranchAddress("Pdc228", &Pdc2[28]);
  Tin->SetBranchAddress("Pdc229", &Pdc2[29]);
  Tin->SetBranchAddress("Pdc230", &Pdc2[30]);
  Tin->SetBranchAddress("Pdc231", &Pdc2[31]);

  Tin->SetBranchAddress("Padc0",  &Padc[0]);
  Tin->SetBranchAddress("Padc1",  &Padc[1]);
  Tin->SetBranchAddress("Padc2",  &Padc[2]);
  Tin->SetBranchAddress("Padc3",  &Padc[3]);
  Tin->SetBranchAddress("Padc4",  &Padc[4]);
  Tin->SetBranchAddress("Padc5",  &Padc[5]);
  Tin->SetBranchAddress("Padc6",  &Padc[6]);
  Tin->SetBranchAddress("Padc7",  &Padc[7]);
  Tin->SetBranchAddress("Padc8",  &Padc[8]);
  Tin->SetBranchAddress("Padc9",  &Padc[9]);
  Tin->SetBranchAddress("Padc10", &Padc[10]);
  Tin->SetBranchAddress("Padc11", &Padc[11]);
  Tin->SetBranchAddress("Padc12", &Padc[12]);
  Tin->SetBranchAddress("Padc13", &Padc[13]);
  Tin->SetBranchAddress("Padc14", &Padc[14]);
  Tin->SetBranchAddress("Padc15", &Padc[15]);
  Tin->SetBranchAddress("Padc16", &Padc[16]);
  Tin->SetBranchAddress("Padc17", &Padc[17]);
  Tin->SetBranchAddress("Padc18", &Padc[18]);
  Tin->SetBranchAddress("Padc19", &Padc[19]);
  Tin->SetBranchAddress("Padc20", &Padc[20]);
  Tin->SetBranchAddress("Padc21", &Padc[21]);
  Tin->SetBranchAddress("Padc22", &Padc[22]);
  Tin->SetBranchAddress("Padc23", &Padc[23]);
  Tin->SetBranchAddress("Padc24", &Padc[24]);
  Tin->SetBranchAddress("Padc25", &Padc[25]);
  Tin->SetBranchAddress("Padc26", &Padc[26]);
  Tin->SetBranchAddress("Padc27", &Padc[27]);
  Tin->SetBranchAddress("Padc28", &Padc[28]);
  Tin->SetBranchAddress("Padc29", &Padc[29]);
  Tin->SetBranchAddress("Padc30", &Padc[30]);
  Tin->SetBranchAddress("Padc31", &Padc[31]);

  Tin->SetBranchAddress("S83028", &S83028);
  Tin->SetBranchAddress("S8200", &S8200);
  Tin->SetBranchAddress("Clk0", &Clk0);
  Tin->SetBranchAddress("Clk1", &Clk1);
  Tin->SetBranchAddress("Clk2", &Clk2);
  Tin->SetBranchAddress("Clk3", &Clk3);

  Tin->SetBranchAddress("Pdc38",  &Pdc38);
  Tin->SetBranchAddress("Pdc39",  &Pdc39);
  Tin->SetBranchAddress("Pdc310", &Pdc310);
  Tin->SetBranchAddress("Pdc311", &Pdc311);

  Tin->SetBranchAddress("Qadc9",  &Qadc9);
  Tin->SetBranchAddress("Tdc019", &Tdc019);
  Tin->SetBranchAddress("Pdc313", &Pdc313);
  Tin->SetBranchAddress("Pdc315", &Pdc315);
  Tin->SetBranchAddress("Qadc8",  &Qadc8);
  Tin->SetBranchAddress("Tdc018", &Tdc018);
  Tin->SetBranchAddress("Qadc10", &Qadc10);
  Tin->SetBranchAddress("Tdc020", &Tdc020);

  Tin->SetBranchAddress("Number", &Number);
  Tin->SetBranchAddress("Delt0", &Delt0);

  Tin->SetBranchAddress("Evnb0", &Evnb[0]);
  Tin->SetBranchAddress("Evnb1", &Evnb[1]);
  Tin->SetBranchAddress("Evnb2", &Evnb[2]);
  Tin->SetBranchAddress("Evnb3", &Evnb[3]);
  Tin->SetBranchAddress("Evnb4", &Evnb[4]);
  Tin->SetBranchAddress("Bkhf0", &Bkhf[0]);
  Tin->SetBranchAddress("Bkhf1", &Bkhf[1]);
  Tin->SetBranchAddress("Bkhf2", &Bkhf[2]);
  Tin->SetBranchAddress("Bkhf3", &Bkhf[3]);
  Tin->SetBranchAddress("Bkhf4", &Bkhf[4]);



  


  int nEvents = Tin->GetEntries();
  std::cout << "... Processing nEvents = " << nEvents << std::endl;

    //Get length of run from UNBLINDED TIME
  Tin->GetEvent(nEvents-1);
  float time = S83028*scalerCountsToTime;
  //Histograms to hold UCNMonRate
  int nbins = (int)((time+10.)/10.);
  TH1F *UCN_Mon_1_Rate = new TH1F("UCN_Mon_1_Rate","UCN Mon 1 Rate",nbins, 0., nbins*10.);
  TH1F *UCN_Mon_2_Rate = new TH1F("UCN_Mon_2_Rate","UCN Mon 2 Rate",nbins, 0., nbins*10.);
  TH1F *UCN_Mon_3_Rate = new TH1F("UCN_Mon_3_Rate","UCN Mon 3 Rate",nbins, 0., nbins*10.);
  TH1F *UCN_Mon_4_Rate = new TH1F("UCN_Mon_4_Rate","UCN Mon 4 Rate",nbins, 0., nbins*10.);

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    Tin->GetEvent(i);
    Int_t iSis00 = (int) Sis00;
    time = S83028*scalerCountsToTime;
    if (iSis00==260) UCN_Mon_1_Rate->Fill(time);
    else if (iSis00==516) UCN_Mon_2_Rate->Fill(time);
    else if (iSis00==1028) UCN_Mon_3_Rate->Fill(time);
    else if (iSis00==2052) UCN_Mon_4_Rate->Fill(time);
  }

  fileIn->Close();
  

  //Append to appropriate files
  char temp[500];
  sprintf(temp,"%s/spec_%i.root",getenv("UK_SPEC_REPLAY"), runNumber);
  TFile *spec = new TFile(temp,"UPDATE");
  spec->cd();
  UCN_Mon_1_Rate->Write();
  UCN_Mon_2_Rate->Write();
  UCN_Mon_3_Rate->Write();
  UCN_Mon_4_Rate->Write();

  spec->Close();

  sprintf(temp,"%s/replay_pass4_%i.root",getenv("REPLAY_PASS4"), runNumber);
  TFile *rep4 = new TFile(temp,"UPDATE");
  rep4->cd();
  UCN_Mon_1_Rate->Write();
  UCN_Mon_2_Rate->Write();
  UCN_Mon_3_Rate->Write();
  UCN_Mon_4_Rate->Write();

  rep4->Close();

  delete UCN_Mon_1_Rate;
  delete UCN_Mon_2_Rate;
  delete UCN_Mon_3_Rate;
  delete UCN_Mon_4_Rate;

  return 0;
}
