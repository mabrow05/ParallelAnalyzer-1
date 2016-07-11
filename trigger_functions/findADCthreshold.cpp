#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <TString.h>
#include <TH2D.h>
#include <TH1.h>
#include <TF1.h>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TMath.h>


struct MWPC {
  Double_t center;
  Double_t width;
  Double_t cathSum;
  Double_t maxValue;
  Int_t maxWire;
  Int_t mult;
  Int_t nClipped;
  Int_t err;
  Double_t rawCenter;
  Double_t height;
};

struct cuts {
  double cutBeamBurstTime; // Beam Burst T0
  int nCutsTimeWindows; // Number of Time Window Cuts
  std::vector <Double_t> cutTimeWindowLower;
  std::vector <Double_t> cutTimeWindowUpper;
  Double_t cutEastAnode;
  Double_t cutWestAnode;
  Double_t cutEastTwoFold;
  Double_t cutWestTwoFold;
  Double_t cutEastTopVetoQADC;
  Double_t cutEastTopVetoTDC;
  Double_t cutEastDriftTubeTAC;
  Double_t cutWestDriftTubeTAC;
  Double_t cutEastBackingVetoQADC;
  Double_t cutEastBackingVetoTDC;
  Double_t cutWestBackingVetoQADC;
  Double_t cutWestBackingVetoTDC;
};

std::vector < std::vector <Double_t> > loadPMTpedestals(Int_t runNumber);

void loadCuts(Int_t runNumber, cuts* Cuts);

void findDiscriminatorThresh(TString octetORxenon, Int_t runListIndex); 

std::vector<Double_t> fitThreshold(TH1F *trigg, std::vector <Double_t> paramGuesses = std::vector <Double_t>(0));

int main(int argc, char* argv[]) {
  
  std::vector <Int_t> evtTypes;

  if (argc!=3) { std::cout << "USAGE: ./findADCthreshold.exe [octet or xenon] [runListIndex]\n"; exit(0); }

  findDiscriminatorThresh(TString(argv[1]),atoi(argv[2]));

}


std::vector < std::vector <Double_t> > loadPMTpedestals(Int_t runNumber) {

  Char_t temp[500];
  std::vector < std::vector < Double_t > > peds (8,std::vector <Double_t> (2,0.));
  sprintf(temp,"%s/pedestal_widths_%i.dat",getenv("PEDESTALS"),runNumber);
  ifstream infile;
  infile.open(temp);

  Int_t i = 0;
  Int_t run;

  while (infile >> run >> peds[i][0] >> peds[i][1]) { std::cout << "Pedestal " << i << ": " << peds[i][0] << " " << peds[i][1] << std::endl; i++; }
  return peds;

};


void findDiscriminatorThresh(TString octetORxenon, Int_t runListIndex) {

  // Read in the Beta/octet runs in this Octet
  std::vector <Int_t> runs;
  //runs.push_back(18090);
  

  Int_t rn;
  UInt_t numRuns=0;
  Char_t temp[200];
  ifstream runList;

  if (octetORxenon==TString("octet")) {

    Char_t type[5];
    sprintf(temp,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),runListIndex);
    runList.open(temp);
    while (runList >> type >> rn) {
      runs.push_back(rn);
      std::cout << runs[numRuns] << std::endl;
      numRuns++;
      //if (numRuns>30) break;
    }
    runList.close();
  }
  else { 
    
    sprintf(temp,"../run_lists/Xenon_Calibration_Run_Period_%i.dat",runListIndex);
    runList.open(temp);
    while (runList >> rn) {
      runs.push_back(rn);
      std::cout << runs[numRuns] << std::endl;
      numRuns++;
      //if (numRuns>30) break;
    }
    runList.close();
  }
  
  if (numRuns!=runs.size()) {
    std::cout << "Runs!=runs.size()\n";
    exit(0);
  }

  std::cout << "There are " << numRuns << " in this Octet/XePeriod\n";
  
  
  TFile *f = NULL;
  TTree *t = NULL;

  Double_t lower_limit = -30.;
  Int_t hisBinWidth = 2;
  Double_t upper_limit = 150.;
  Int_t nBins = (int)(upper_limit/hisBinWidth);

  TH1F *totalE1 = new TH1F("totalE1", "total events E1", nBins, lower_limit, upper_limit);
  TH1F *triggerE1 = new TH1F("triggerE1", "trigger events E1", nBins, lower_limit, upper_limit);
  TH1F *probE1 = new TH1F("probE1", "probability E1", nBins, lower_limit, upper_limit);
  
  TH1F *totalE2 = new TH1F("totalE2", "total events E2", nBins, lower_limit, upper_limit);
  TH1F *triggerE2 = new TH1F("triggerE2", "trigger events E2", nBins, lower_limit, upper_limit);
  TH1F *probE2 = new TH1F("probE2", "probability E2", nBins, lower_limit, upper_limit);

  TH1F *totalE3 = new TH1F("totalE3", "total events E3", nBins, lower_limit, upper_limit);
  TH1F *triggerE3 = new TH1F("triggerE3", "trigger events E3", nBins, lower_limit, upper_limit);
  TH1F *probE3 = new TH1F("probE3", "probability E3", nBins, lower_limit, upper_limit);

  TH1F *totalE4 = new TH1F("totalE4", "total events E4", nBins, lower_limit, upper_limit);
  TH1F *triggerE4 = new TH1F("triggerE4", "trigger events E4", nBins, lower_limit, upper_limit);
  TH1F *probE4 = new TH1F("probE4", "probability E4", nBins, lower_limit, upper_limit);

  TH1F *totalW1 = new TH1F("totalW1", "total events W1", nBins, lower_limit, upper_limit);
  TH1F *triggerW1 = new TH1F("triggerW1", "trigger events W1", nBins, lower_limit, upper_limit);
  TH1F *probW1 = new TH1F("probW1", "probability W1", nBins, lower_limit, upper_limit);
  
  TH1F *totalW2 = new TH1F("totalW2", "total events W2", nBins, lower_limit, upper_limit);
  TH1F *triggerW2 = new TH1F("triggerW2", "trigger events W2", nBins, lower_limit, upper_limit);
  TH1F *probW2 = new TH1F("probW2", "probability W2", nBins, lower_limit, upper_limit);

  TH1F *totalW3 = new TH1F("totalW3", "total events W3", nBins, lower_limit, upper_limit);
  TH1F *triggerW3 = new TH1F("triggerW3", "trigger events W3", nBins, lower_limit, upper_limit);
  TH1F *probW3 = new TH1F("probW3", "probability W3", nBins, lower_limit, upper_limit);

  TH1F *totalW4 = new TH1F("totalW4", "total events W4", nBins, lower_limit, upper_limit);
  TH1F *triggerW4 = new TH1F("triggerW4", "trigger events W4", nBins, lower_limit, upper_limit);
  TH1F *probW4 = new TH1F("probW4", "probability W4", nBins, lower_limit, upper_limit);
  

  // All pertinent variables to be read in from data files
  Float_t QadcE1, QadcE2, QadcE3, QadcE4, QadcW1, QadcW2, QadcW3, QadcW4;
  Float_t PdcWest, PdcEast; //Wirechamber Anode values
  Float_t TdcE1, TdcE2, TdcE3, TdcE4, TdcW1, TdcW2, TdcW3, TdcW4; 
  Float_t TdcEast2fold, TdcWest2fold;
  Float_t Sis00;
  

  // pedestal variables
  Double_t pedE1, pedE2, pedE3, pedE4, pedW1, pedW2, pedW3, pedW4;
  Double_t pedE1_width, pedE2_width, pedE3_width, pedE4_width, pedW1_width, pedW2_width, pedW3_width, pedW4_width;

  //Cut variables  
  cuts Cuts;

  std::vector < std::vector < Double_t > > peds(8, std::vector <Double_t> (2,0.));

  // Cycling through all the input files and filling the histograms
  for (UInt_t i=0; i<runs.size(); i++) {
    Int_t run = runs[i];
    std::cout << "MADE IT\n";

    loadCuts(run, &Cuts);
    peds = loadPMTpedestals(run);    

    pedE1 = peds[0][0];
    pedE2 = peds[1][0];
    pedE3 = peds[2][0];
    pedE4 = peds[3][0];
    pedW1 = peds[4][0];
    pedW2 = peds[5][0];
    pedW3 = peds[6][0];
    pedW4 = peds[7][0];

    pedE1_width = peds[0][1];
    pedE2_width = peds[1][1];
    pedE3_width = peds[2][1];
    pedE4_width = peds[3][1];
    pedW1_width = peds[4][1];
    pedW2_width = peds[5][1];
    pedW3_width = peds[6][1];
    pedW4_width = peds[7][1];

    std::cout << "Ped width: " << pedE1 << std::endl;
    
	 
    f = new TFile(TString::Format("%s/full%i.root",getenv("UCNA_RAW_DATA"),run),"READ");
    t = (TTree*)(f->Get("h1"));

    // Variables
    t->SetBranchAddress("Pdc30",  &PdcEast);
    t->SetBranchAddress("Pdc34",  &PdcWest);
    
    t->SetBranchAddress("Tdc016", &TdcEast2fold);
    t->SetBranchAddress("Tdc017", &TdcWest2fold);
    t->SetBranchAddress("Tdc00", &TdcE1);
    t->SetBranchAddress("Tdc01", &TdcE2);
    t->SetBranchAddress("Tdc02", &TdcE3);
    t->SetBranchAddress("Tdc03", &TdcE4);
    t->SetBranchAddress("Tdc08", &TdcW1);
    t->SetBranchAddress("Tdc09", &TdcW2);
    t->SetBranchAddress("Tdc014", &TdcW3);
    t->SetBranchAddress("Tdc011", &TdcW4);
    
    t->SetBranchAddress("Sis00",  &Sis00);
    
    t->SetBranchAddress("Qadc0",  &QadcE1);
    t->SetBranchAddress("Qadc1",  &QadcE2);
    t->SetBranchAddress("Qadc2",  &QadcE3);
    t->SetBranchAddress("Qadc3",  &QadcE4);
    t->SetBranchAddress("Qadc4",  &QadcW1);
    t->SetBranchAddress("Qadc5",  &QadcW2);
    t->SetBranchAddress("Qadc6",  &QadcW3);
    t->SetBranchAddress("Qadc7",  &QadcW4);
    

    Int_t nevents = t->GetEntries();
    
    std::cout << "nEvents = " << nevents << std::endl;

    for (int evt=0; evt<nevents; evt++) {
      
      t->GetEvent(evt);
      

      if (evt%100000==0) std::cout << evt << " in run " << run << std::endl;

      
      if ((int)Sis00>3) continue; //Not a beta trigger

      //Int_t side = (int)Sis00 - 1;

      if ( ((int)Sis00==1 || (int)Sis00==3) && PdcEast<Cuts.cutEastAnode) continue; //not an electron event
      if ( ((int)Sis00==2  || (int)Sis00==3) && PdcWest<Cuts.cutWestAnode) continue; //not an electron event
      

      // East Side Triggers (including backscatters)
      if ( TdcEast2fold > 0. ) {

	//Fill West side total histograms
	totalW1->Fill(QadcW1-pedW1);
	totalW2->Fill(QadcW2-pedW2);
	totalW3->Fill(QadcW3-pedW3);
	totalW4->Fill(QadcW4-pedW4);

	if ( TdcW1>0. ) triggerW1->Fill(QadcW1-pedW1);
	if ( TdcW2>0. ) triggerW2->Fill(QadcW2-pedW2);
	if ( TdcW3>0. ) triggerW3->Fill(QadcW3-pedW3);
	if ( TdcW4>0. ) triggerW4->Fill(QadcW4-pedW4);

	//3 PMTs trigger
	if ( (TdcE1>0. && TdcE2>0. && TdcE3>0.) || (TdcE2>0. && TdcE3>0. && TdcE4>0.) || (TdcE1>0. && TdcE3>0. && TdcE4>0.) || (TdcE1>0. && TdcE2>0. && TdcE4>0.) ) { 
	  
	  //East side 3-fold trigger events
	  totalE1->Fill(QadcE1-pedE1);
	  totalE2->Fill(QadcE2-pedE2);
	  totalE3->Fill(QadcE3-pedE3);
	  totalE4->Fill(QadcE4-pedE4);

	  if ( TdcE1>0. ) triggerE1->Fill(QadcE1-pedE1);
	  if ( TdcE2>0. ) triggerE2->Fill(QadcE2-pedE2);
	  if ( TdcE3>0. ) triggerE3->Fill(QadcE3-pedE3);
	  if ( TdcE4>0. ) triggerE4->Fill(QadcE4-pedE4);
	}
	
      }

      // West Side Triggers (including backscatters)
      if ( TdcWest2fold > 0. ) {

	//Fill East side total histograms
	totalE1->Fill(QadcE1-pedE1);
	totalE2->Fill(QadcE2-pedE2);
	totalE3->Fill(QadcE3-pedE3);
	totalE4->Fill(QadcE4-pedE4);

	if ( TdcE1>0. ) triggerE1->Fill(QadcE1-pedE1);
	if ( TdcE2>0. ) triggerE2->Fill(QadcE2-pedE2);
	if ( TdcE3>0. ) triggerE3->Fill(QadcE3-pedE3);
	if ( TdcE4>0. ) triggerE4->Fill(QadcE4-pedE4);

	
	if ( (TdcW1>0. && TdcW2>0. && TdcW3>0.) || (TdcW2>0. && TdcW3>0. && TdcW4>0.) || (TdcW1>0. && TdcW3>0. && TdcW4>0.) || (TdcW1>0. && TdcW2>0. && TdcW4>0.) ) { 
	  
	  //West side events
	  totalW1->Fill(QadcW1-pedW1);
	  totalW2->Fill(QadcW2-pedW2);
	  totalW3->Fill(QadcW3-pedW3);
	  totalW4->Fill(QadcW4-pedW4);
	  
	  if ( TdcW1>0. ) triggerW1->Fill(QadcW1-pedW1);
	  if ( TdcW2>0. ) triggerW2->Fill(QadcW2-pedW2);
	  if ( TdcW3>0. ) triggerW3->Fill(QadcW3-pedW3);
	  if ( TdcW4>0. ) triggerW4->Fill(QadcW4-pedW4);
	}
      }
      
    }
    f->Close();
  }

  probE1->SetMarkerStyle(20);
  probE2->SetMarkerStyle(20);
  probE3->SetMarkerStyle(20);
  probE4->SetMarkerStyle(20);
  probW1->SetMarkerStyle(20);
  probW2->SetMarkerStyle(20);
  probW3->SetMarkerStyle(20);
  probW4->SetMarkerStyle(20);

  triggerE1->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerE2->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerE3->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerE4->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerW1->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerW2->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerW3->GetXaxis()->SetTitle("ADC channels above pedestal");
  triggerW4->GetXaxis()->SetTitle("ADC channels above pedestal");

  totalE1->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalE2->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalE3->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalE4->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalW1->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalW2->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalW3->GetXaxis()->SetTitle("ADC channels above pedestal");
  totalW4->GetXaxis()->SetTitle("ADC channels above pedestal");

  probE1->GetXaxis()->SetTitle("ADC channels above pedestal");
  probE2->GetXaxis()->SetTitle("ADC channels above pedestal");
  probE3->GetXaxis()->SetTitle("ADC channels above pedestal");
  probE4->GetXaxis()->SetTitle("ADC channels above pedestal");
  probW1->GetXaxis()->SetTitle("ADC channels above pedestal");
  probW2->GetXaxis()->SetTitle("ADC channels above pedestal");
  probW3->GetXaxis()->SetTitle("ADC channels above pedestal");
  probW4->GetXaxis()->SetTitle("ADC channels above pedestal");

  probE1->SetMinimum(0.);
  probE2->SetMinimum(0.);
  probE3->SetMinimum(0.);
  probE4->SetMinimum(0.);
  probW1->SetMinimum(0.);
  probW2->SetMinimum(0.);
  probW3->SetMinimum(0.);
  probW4->SetMinimum(0.);

  
  triggerE1->SetLineColor(kRed);
  triggerE1->SetLineStyle(7);
  triggerE2->SetLineColor(kRed);
  triggerE2->SetLineStyle(7);
  triggerE3->SetLineColor(kRed);
  triggerE3->SetLineStyle(7);
  triggerE4->SetLineColor(kRed);
  triggerE4->SetLineStyle(7);
  triggerW1->SetLineColor(kRed);
  triggerW1->SetLineStyle(7);
  triggerW2->SetLineColor(kRed);
  triggerW2->SetLineStyle(7);
  triggerW3->SetLineColor(kRed);
  triggerW3->SetLineStyle(7);
  triggerW4->SetLineColor(kRed);
  triggerW4->SetLineStyle(7);

  probE1->Divide(triggerE1, totalE1);
  probE2->Divide(triggerE2, totalE2);
  probE3->Divide(triggerE3, totalE3);
  probE4->Divide(triggerE4, totalE4);

  probW1->Divide(triggerW1, totalW1);
  probW2->Divide(triggerW2, totalW2);
  probW3->Divide(triggerW3, totalW3);
  probW4->Divide(triggerW4, totalW4);

  //Set all the error bars in prob histos
  for (Int_t bin=1; bin<=nBins; bin++) {
    probE1->SetBinError(bin, (totalE1->GetBinContent(bin) > 0. ? probE1->GetBinContent(bin)
			      *TMath::Sqrt(1./(triggerE1->GetBinContent(bin)>0.?triggerE1->GetBinContent(bin):1.) + 1./totalE1->GetBinContent(bin)) : 1.));
    probE2->SetBinError(bin, (totalE2->GetBinContent(bin) > 0. ? probE2->GetBinContent(bin)
			      *TMath::Sqrt(1./(triggerE2->GetBinContent(bin)>0.?triggerE2->GetBinContent(bin):1.) + 1./totalE2->GetBinContent(bin)) : 1.));
    probE3->SetBinError(bin, (totalE3->GetBinContent(bin) > 0. ? probE3->GetBinContent(bin)
			      *TMath::Sqrt(1./(triggerE3->GetBinContent(bin)>0.?triggerE3->GetBinContent(bin):1.) + 1./totalE3->GetBinContent(bin)) : 1.));
    probE4->SetBinError(bin, (totalE4->GetBinContent(bin) > 0. ? probE4->GetBinContent(bin)
			      *TMath::Sqrt(1./(triggerE4->GetBinContent(bin)>0.?triggerE4->GetBinContent(bin):1.) + 1./totalE4->GetBinContent(bin)) : 1.));
    probW1->SetBinError(bin, (totalW1->GetBinContent(bin) > 0. ? probW1->GetBinContent(bin)
			      *TMath::Sqrt(1./(triggerW1->GetBinContent(bin)>0.?triggerW1->GetBinContent(bin):1.) + 1./totalW1->GetBinContent(bin)) : 1.));
    probW2->SetBinError(bin, (totalW2->GetBinContent(bin) > 0. ? probW2->GetBinContent(bin)
			      *TMath::Sqrt(1./(triggerW2->GetBinContent(bin)>0.?triggerW2->GetBinContent(bin):1.) + 1./totalW2->GetBinContent(bin)) : 1.));
    probW3->SetBinError(bin, (totalW3->GetBinContent(bin) > 0. ? probW3->GetBinContent(bin)
			      *TMath::Sqrt(1./(triggerW3->GetBinContent(bin)>0.?triggerW3->GetBinContent(bin):1.) + 1./totalW3->GetBinContent(bin)) : 1.));
    probW4->SetBinError(bin, (totalW4->GetBinContent(bin) > 0. ? probW4->GetBinContent(bin)
			      *TMath::Sqrt(1./(triggerW4->GetBinContent(bin)>0.?triggerW4->GetBinContent(bin):1.) + 1./totalW4->GetBinContent(bin)) : 1.));
  }
	

  //Fit the probabaility histograms and write to file
  TString filename = TString::Format("%s/%s/%s_%i_thresholds.dat",getenv("TRIGGER_FUNC"),octetORxenon.Data(),octetORxenon.Data(), runListIndex);
  ofstream funcfile(filename.Data());
  
  std::vector <Double_t> params(8,0);
  
  params = fitThreshold(probE1);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probE2);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probE3);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params =   fitThreshold(probE4);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probW1);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probW2);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probW3);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;
  params = fitThreshold(probW4);
  funcfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] 
	   << " " << params[4] << " " << params[5] << " " << params[6] << " " << params[7] << std::endl;

  funcfile.close();

  TString pdf_file = TString::Format("%s/%s/%s_%i_Thresholds.pdf",getenv("TRIGGER_FUNC"),octetORxenon.Data(),octetORxenon.Data(), runListIndex);

  TCanvas *cEtrigg = new TCanvas("cEtrigg"," ",2000,800);
  cEtrigg->Divide(4,1);

  cEtrigg->cd(1);
  gPad->SetLogy();
  totalE1->Draw();
  triggerE1->Draw("SAME");
  cEtrigg->cd(2);
  gPad->SetLogy();
  totalE2->Draw();
  triggerE2->Draw("SAME");
  cEtrigg->cd(3);
  gPad->SetLogy();
  totalE3->Draw();
  triggerE3->Draw("SAME");
  cEtrigg->cd(4);
  gPad->SetLogy();
  totalE4->Draw();
  triggerE4->Draw("SAME");

  cEtrigg->Print(pdf_file+"(");

  TCanvas *cWtrigg = new TCanvas("cWtrigg"," ",2000,800);
  cWtrigg->Divide(4,1);

  cWtrigg->cd(1);
  gPad->SetLogy();
  totalW1->Draw();
  triggerW1->Draw("SAME");
  cWtrigg->cd(2);
  gPad->SetLogy();
  totalW2->Draw();
  triggerW2->Draw("SAME");
  cWtrigg->cd(3);
  gPad->SetLogy();
  totalW3->Draw();
  triggerW3->Draw("SAME");
  cWtrigg->cd(4);
  gPad->SetLogy();
  totalW4->Draw();
  triggerW4->Draw("SAME");

  cWtrigg->Print(pdf_file);

  TCanvas *cEprob = new TCanvas("cEprob"," ",2000,800);
  cEprob->Divide(4,1);

  cEprob->cd(1);
  probE1->Draw("EP");
  cEprob->cd(2);
  probE2->Draw("EP");
  cEprob->cd(3);
  probE3->Draw("EP");
  cEprob->cd(4);
  probE4->Draw("EP");

  cEprob->Print(pdf_file);

  TCanvas *cWprob = new TCanvas("cWprob"," ",2000,800);
  cWprob->Divide(4,1);

  cWprob->cd(1);
  probW1->Draw("EP");
  cWprob->cd(2);
  probW2->Draw("EP");
  cWprob->cd(3);
  probW3->Draw("EP");
  cWprob->cd(4);
  probW4->Draw("EP");

  cWprob->Print(pdf_file+")");

  //Clean-up
  delete cEtrigg; delete cWtrigg; delete cEprob; delete cWprob;

  delete triggerE1; delete triggerE2; delete triggerE3; delete triggerE4;
  delete triggerW1; delete triggerW2; delete triggerW3; delete triggerW4;

  delete totalE1; delete totalE2; delete totalE3; delete totalE4;
  delete totalW1; delete totalW2; delete totalW3; delete totalW4;

  delete probE1; delete probE2; delete probE3; delete probE4;
  delete probW1; delete probW2; delete probW3; delete probW4;


}


void loadCuts(Int_t runNumber, cuts* Cuts) {

  // Read cuts file
  char tempFileCuts[500];
  sprintf(tempFileCuts, "%s/cuts_%i.dat", getenv("CUTS"),runNumber);
  std::cout << "... Reading: " << tempFileCuts << std::endl;

  ifstream fileCuts(tempFileCuts);
  Char_t comment[500];
  fileCuts >> Cuts->cutBeamBurstTime >> comment;
  fileCuts >> Cuts->nCutsTimeWindows >> comment;
  
  Double_t cutLowerHold, cutUpperHold;

  if (Cuts->nCutsTimeWindows > 0) {
    for (int i=0; i<Cuts->nCutsTimeWindows; i++) {
      fileCuts >> cutLowerHold >> cutUpperHold;
      Cuts->cutTimeWindowLower.push_back(cutLowerHold);
      Cuts->cutTimeWindowUpper.push_back(cutUpperHold);
    }
  }
  fileCuts >> Cuts->cutEastAnode >> comment;
  fileCuts >> Cuts->cutWestAnode >> comment;
  fileCuts >> Cuts->cutEastTwoFold >> comment;
  fileCuts >> Cuts->cutWestTwoFold >> comment;
  fileCuts >> Cuts->cutEastTopVetoQADC >> comment;
  fileCuts >> Cuts->cutEastTopVetoTDC >> comment;
  fileCuts >> Cuts->cutEastDriftTubeTAC >> comment;
  fileCuts >> Cuts->cutWestDriftTubeTAC >> comment;
  fileCuts >> Cuts->cutEastBackingVetoQADC >> comment;
  fileCuts >> Cuts->cutEastBackingVetoTDC >> comment;
  fileCuts >> Cuts->cutWestBackingVetoQADC >> comment;
  fileCuts >> Cuts->cutWestBackingVetoTDC >> comment;

  std::cout << "... Beam Burst T0 Cut: " << Cuts->cutBeamBurstTime << std::endl;
  std::cout << "... Number of Time Windows Cuts: " << Cuts->nCutsTimeWindows << std::endl;
  if (Cuts->nCutsTimeWindows > 0) {
    for (int i=0; i<Cuts->nCutsTimeWindows; i++) {
      std::cout << "        [" << Cuts->cutTimeWindowLower[i] << ", " << Cuts->cutTimeWindowUpper[i] << "]" << std::endl;
    }
  }
  std::cout << "... East MWPC Anode Cut: " << Cuts->cutEastAnode << std::endl;
  std::cout << "... West MWPC Anode Cut: " << Cuts->cutWestAnode << std::endl;
  std::cout << "... East Scintillator Two-Fold Trigger Cut: " << Cuts->cutEastTwoFold << std::endl;
  std::cout << "... West Scintillator Two-Fold Trigger Cut: " << Cuts->cutWestTwoFold << std::endl;
  std::cout << "... East Top Veto QADC Cut: " << Cuts->cutEastTopVetoQADC << std::endl;
  std::cout << "... East Top Veto TDC Cut: " << Cuts->cutEastTopVetoTDC << std::endl;
  std::cout << "... East Drift Tube TAC Cut: " << Cuts->cutEastDriftTubeTAC << std::endl;
  std::cout << "... West Drift Tube TAC Cut: " << Cuts->cutWestDriftTubeTAC << std::endl;
  std::cout << "... East Backing Veto QADC Cut: " << Cuts->cutEastBackingVetoQADC << std::endl;
  std::cout << "... East Backing Veto TDC Cut: " << Cuts->cutEastBackingVetoTDC << std::endl;
  std::cout << "... West Backing Veto QADC Cut: " << Cuts->cutWestBackingVetoQADC << std::endl;
  std::cout << "... West Backing Veto TDC Cut: " << Cuts->cutWestBackingVetoTDC << std::endl;

}

std::vector <Double_t> fitThreshold(TH1F *trigg, std::vector <Double_t> params) {

  bool defaultParams = true;
  if (params.size()==8) defaultParams = false;

  TF1 *erf = new TF1("erf","([0]+[1]*TMath::Erf((x-[2])/[3]))*(0.5-.5*TMath::TanH((x-[2])/[4]))+(0.5+.5*TMath::TanH((x-[2])/[4]))*([5]+[6]*TMath::TanH((x-[2])/[7]))",-20.,150.);

  if (defaultParams) {
    params.resize(8);
    erf->SetParameter(0,0.5); //Constant offset of erf                                                                                       
    erf->SetParameter(1,0.5); //Scaling of erf                                                                                               
    erf->SetParameter(2,25.); //Mean of gaussian integrated for erf                                                                          
    erf->SetParameter(3,7.); //std. dev. of gaussian integrated for erf                                                                      
    erf->SetParameter(4,9.); //severity of transition function "turn on"                                                                     
    erf->SetParameter(5,0.5); //constant offset of second tanh                                                                               
    erf->SetParameter(6,0.5); //Scaling of second tanh                                                                                        
    erf->SetParameter(7,12.); // stretching factor of second tanh                                                                            
    
  }
  else {
    erf->SetParameter(0,params[0]); //Constant offset of erf                                                                            
    erf->SetParameter(1,params[1]); //Scaling of erf                                                                                   
    erf->SetParameter(2,params[2]); //Mean of gaussian integrated for erf                                                               
    erf->SetParameter(3,params[3]); //std. dev. of gaussian integrated for erf                                                         
    erf->SetParameter(4,params[4]); //severity of transition function "turn on"                                                        
    erf->SetParameter(5,params[5]); //constant offset of second tanh                                                                   
    erf->SetParameter(6,params[6]); //Scaling of second tanh                                                                           
    erf->SetParameter(7,params[7]); // stretching factor of second tanh
  }
  
  erf->FixParameter(0,0.5);
  erf->FixParameter(1,0.5);
  // erf->SetParLimits(4,1.,25.);
  //erf->SetParLimits(4,0.,2.);
  erf->FixParameter(5,0.5);                                                                                                              
  erf->FixParameter(6,0.5);
  //erf->SetParLimits(7,0.,50.); //Constant offset of erf     
  
    
  trigg->Fit("erf","R");

  for (UInt_t i=0; i<params.size();i++) params[i] = erf->GetParameter(i);

  delete erf;

  return params;
}
