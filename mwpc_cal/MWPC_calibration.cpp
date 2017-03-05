#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include <sstream>

// ROOT libraries
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TCanvas.h>

#include "runInfo.h"
#include "DataTree.hh"
#include "MWPCPositionResponse.hh"
#include "peaks.hh"
#include "positionMapHandler.hh"

std::vector <int>  readOctetFile(int octet) {

  std::vector <int> runs;
  
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  int numRuns = 0;
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    if (runTypeHold=="A2" || runTypeHold=="A5" || runTypeHold=="A7" || runTypeHold=="A10" || 
	runTypeHold=="B2" || runTypeHold=="B5" || runTypeHold=="B7" || runTypeHold=="B10" )  {
     
      runs.push_back(runNumberHold);

    }
    numRuns++;
  }

  infile.close();
 
  std::cout << "Read in octet file for octet " << octet << "\n";
  return runs;

};


int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  int octetNumber = atoi(argv[1]);
 
  cout << "Octet Number = " << octetNumber << endl;

  std::vector <int> betaRuns = readOctetFile( octetNumber );

  // Output root file
  TString path = TString::Format("%s/MWPC_Cal_octet_%i.root",getenv("MWPC_CALIBRATION"),octetNumber);
  TFile *fileOut = new TFile(path,"RECREATE");
  
  // Output histograms
  TH1D *hisE_data[8];
  TH1D *hisE_sim[8];
  TH1D *hisW_data[8];
  TH1D *hisW_sim[8];

  double EnergyBinWidth = 100.;
  int nEnergyBins = 7;
  double enBinStart = 100.;
  int nBinHist = 200;

  for ( int i = 0; i<nEnergyBins; ++i ) {

    double elow = enBinStart + i*EnergyBinWidth;
    double ehigh = elow + EnergyBinWidth;
    
    hisE_data[i] = new TH1D(TString::Format("hisE_data%i",i),TString::Format("East Anode Data %0.0f-%0.0f keV",elow,ehigh), nBinHist, 0.,2000.);
    hisE_sim[i] = new TH1D(TString::Format("hisE_sim%i",i),TString::Format("East Anode Sim %0.0f-%0.0f keV",elow,ehigh), nBinHist, 0.,20.);
    hisW_data[i] = new TH1D(TString::Format("hisW_data%i",i),TString::Format("West Anode Data %0.0f-%0.0f keV",elow,ehigh), nBinHist, 0.,2000.);
    hisW_sim[i] = new TH1D(TString::Format("hisW_sim%i",i),TString::Format("West Anode Sim %0.0f-%0.0f keV",elow,ehigh), nBinHist, 0.,20.);

  }	
  
  //Loop over all runs in octet
  
  for ( auto rn : betaRuns ) {
    
    // DataTree structure
    DataTree *t = new DataTree();

    // Input ntuple
    char tempIn[500];
    TString infile = TString::Format("%s/replay_pass3_%i.root", getenv("REPLAY_PASS3"),rn);
    t->setupInputTree(std::string(infile.Data()),"pass3");
 
    int nEvents = t->getEntries();
    std::cout << "... Processing nEvents = " << nEvents << std::endl;

    // Load Anode position map
    MWPCPositionMap posmap(5., 50.);
    
    int XePeriod = getXeRunPeriod(rn);
    posmap.readMWPCPositionMap(XePeriod,250.,300.); // Using 250-300 keV because that's the most probable range
	  
    // Loop over all electron events 
    for ( int evt=0; evt<nEvents; ++evt ) {

      t->getEvent(evt);

      if ( evt%10000==0 ) std::cout << "Data Run " << rn << " event " << evt << std::endl;
      
      if ( t->PID != 1 || t->Type != 0 ) continue; // only use type 0 electrons

      // Check that the event is well within the fiducial volume
      double r2E = t->xE.center*t->xE.center + t->yE.center*t->yE.center;
      double r2W = t->xW.center*t->xW.center + t->yW.center*t->yW.center;

      if ( r2E > 40. || r2W > 40. ) continue; 
      
      // Get the position response
      std::vector <Double_t> eta = posmap.getInterpolatedEta(t->xE.center,t->yE.center,
							     t->xW.center,t->yW.center);
   
      // East side first
      if ( t->Side == 0 && t->xE.nClipped==0 && t->yE.nClipped==0 && t->Erecon>enBinStart ) { 

	int hist = 0;
	while ( t->Erecon > ( enBinStart + EnergyBinWidth*(hist+1) ) ) hist++;

	hisE_data[hist]->Fill( t->AnodeE / eta[0] );
	
      }

      // West Side Next
      if ( t->Side == 1 && t->xW.nClipped==0 && t->yW.nClipped==0 && t->Erecon>enBinStart ) { 

	int hist = 0;
	while ( t->Erecon > ( enBinStart + EnergyBinWidth*(hist+1) ) ) hist++;

	hisW_data[hist]->Fill( t->AnodeE / eta[1] );
	
      }

    }

    delete t;
  }


  ////////// Now we do the same thing for the simulated beta runs //////////////////

  int PID, Side, Type, nClippedEX, nClippedEY, nClippedWX, nClippedWY;
  double Erecon, mwpcEX, mwpcEY, mwpcWX, mwpcWY;
  double MWPCEnergy[2];
  
  for ( auto rn : betaRuns ) {

    
    TString infile = TString::Format("%s/beta/revCalSim_%i_Beta.root",getenv("REVCALSIM"),rn);
    TFile *input = new TFile(infile, "READ");
    TTree *Tin = (TTree*)input->Get("revCalSim");
    //std::cout << "made it here in dataReader\n";
    //Set branch addresses
    
    Tin->SetBranchAddress("PID", &PID);
    Tin->SetBranchAddress("type", &Type);
    Tin->SetBranchAddress("side", &Side); 
    Tin->SetBranchAddress("Erecon",&Erecon);
    Tin->SetBranchAddress("MWPCEnergy",MWPCEnergy);
    Tin->GetBranch("xE")->GetLeaf("center")->SetAddress(&mwpcEX);
    Tin->GetBranch("yE")->GetLeaf("center")->SetAddress(&mwpcEY);
    Tin->GetBranch("xW")->GetLeaf("center")->SetAddress(&mwpcWX);
    Tin->GetBranch("yW")->GetLeaf("center")->SetAddress(&mwpcWY);
    Tin->GetBranch("xE")->GetLeaf("nClipped")->SetAddress(&nClippedEX);
    Tin->GetBranch("yE")->GetLeaf("nClipped")->SetAddress(&nClippedEY);
    Tin->GetBranch("xW")->GetLeaf("nClipped")->SetAddress(&nClippedWX);
    Tin->GetBranch("yW")->GetLeaf("nClipped")->SetAddress(&nClippedWY);

    unsigned int nevents = Tin->GetEntriesFast();
    std::cout << "... Processing nEvents = " << nevents << std::endl;
    
    // Loop over all electron events
    for ( unsigned int evt=0; evt<nevents; ++evt ) {

      Tin->GetEvent(evt);

      if ( evt%10000==0 ) std::cout << "Sim Run " << rn << " event " << evt << std::endl;
      
      if ( PID != 1 || Type != 0 ) continue; // only use type 0 electrons

      // Check that the event is well within the fiducial volume
      double r2E = mwpcEX*mwpcEX + mwpcEY*mwpcEY;
      double r2W = mwpcWX*mwpcWX + mwpcWY*mwpcWY;

      if ( r2E > 40. || r2W > 40. ) continue; 
   
      // East side first
      if ( Side == 0 && nClippedEX==0 && nClippedEY==0 && Erecon>enBinStart ) { 

	int hist = 0;
	while ( Erecon > ( enBinStart + EnergyBinWidth*(hist+1) ) ) hist++;

	hisE_sim[hist]->Fill( MWPCEnergy[0] );
	
      }

      // West Side Next
      if ( Side == 1 && nClippedWX==0 && nClippedWY==0 && Erecon>enBinStart ) { 

	int hist = 0;
	while ( Erecon > ( enBinStart + EnergyBinWidth*(hist+1) ) ) hist++;

	hisW_sim[hist]->Fill( MWPCEnergy[1] );
	
      }
    }
    delete input;
  }

  
  /////////////////////////////////////////////////
  // Now fit the histograms with landaus and record the MPV

  double EdataMPV[nEnergyBins];
  double EdataMPVerr[nEnergyBins];
  double EsimMPV[nEnergyBins];
  double EsimMPVerr[nEnergyBins];

  double WdataMPV[nEnergyBins];
  double WdataMPVerr[nEnergyBins];
  double WsimMPV[nEnergyBins];
  double WsimMPVerr[nEnergyBins];

  SinglePeakHist *peak;
  
  for ( int i=0; i<nEnergyBins; ++i ) {

    peak = new SinglePeakHist(hisE_data[i],0., 2000., true, 5, 2., 7., true);
    EdataMPV[i] = peak->ReturnMean();
    EdataMPVerr[i] = peak->ReturnMeanError();
    if ( !peak->isGoodFit() ) { EdataMPV[i] = hisE_data[i]->GetMean(); EdataMPVerr[i] = hisE_data[i]->GetMeanError(); }
    delete peak;
  
    peak = new SinglePeakHist(hisW_data[i],0., 2000., true, 5, 2., 7., true);
    WdataMPV[i] = peak->ReturnMean();
    WdataMPVerr[i] = peak->ReturnMeanError();
    if ( !peak->isGoodFit() ) { WdataMPV[i] = hisW_data[i]->GetMean(); WdataMPVerr[i] = hisW_data[i]->GetMeanError(); }
    delete peak;

    peak = new SinglePeakHist(hisE_sim[i],0., 20., true, 5, 2., 7., true);
    EsimMPV[i] = peak->ReturnMean();
    EsimMPVerr[i] = peak->ReturnMeanError();
    if ( !peak->isGoodFit() ) { EsimMPV[i] = hisE_sim[i]->GetMean(); EsimMPVerr[i] = hisE_sim[i]->GetMeanError(); }
    delete peak;
  
    peak = new SinglePeakHist(hisW_sim[i],0., 20., true, 5, 2., 7., true);
    WsimMPV[i] = peak->ReturnMean();
    WsimMPVerr[i] = peak->ReturnMeanError();
    if ( !peak->isGoodFit() ) { WsimMPV[i] = hisW_sim[i]->GetMean(); WsimMPVerr[i] = hisW_sim[i]->GetMeanError(); }
    delete peak;

  }

  fileOut->Write();
  delete fileOut;

  // Now fit the data to extract the conversion factor

  double eastFactor, eastFactorErr, westFactor, westFactorErr;
  
  TF1 *f = new TF1("f","[0]*x",0.,1000.);
  f->SetParameter(0,0.01);
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->cd(1);
  
  TGraphErrors *gEast = new TGraphErrors(nEnergyBins, EdataMPV, EsimMPV, EdataMPVerr, EsimMPVerr);
  gEast->SetTitle(TString::Format("Octet %i East Wirechamber Calibration",octetNumber));
  gEast->SetMarkerStyle(20);
  gEast->SetMarkerColor(kBlue);
  gEast->SetLineColor(kRed);
  gEast->GetXaxis()->SetTitle("Position Correction Anode (ADC)");
  gEast->GetXaxis()->CenterTitle();
  gEast->GetYaxis()->SetTitle("Simulated MWPC Energy Deposition");
  gEast->GetYaxis()->CenterTitle();
  gEast->Draw("AC");
  
  gEast->Fit(f,"R");

  eastFactor = f->GetParameter(0);
  eastFactorErr = f->GetParError(0);

  c1->Update();

  c1->cd(2);
  
  TGraphErrors *gWest = new TGraphErrors(nEnergyBins, WdataMPV, WsimMPV, WdataMPVerr, WsimMPVerr);
  gWest->SetTitle(TString::Format("Octet %i West Wirechamber Calibration",octetNumber));
  gWest->SetMarkerStyle(20);
  gWest->SetMarkerColor(kBlue);
  gWest->SetLineColor(kRed);
  gWest->GetXaxis()->SetTitle("Position Correction Anode (ADC)");
  gWest->GetXaxis()->CenterTitle();
  gWest->GetYaxis()->SetTitle("Simulated MWPC Energy Deposition");
  gWest->GetYaxis()->CenterTitle();
  gWest->Draw("AC");
  
  gWest->Fit(f,"R");

  westFactor = f->GetParameter(0);
  westFactorErr = f->GetParError(0);

  c1->Update();

  c1->Print(TString::Format("%s/MWPC_cal_octet_%i.pdf()",getenv("MWPC_CAL"),octetNumber));

  // Write out calibration information
    
  std::ofstream mwpcFile(TString::Format("%s/MWPC_cal_octet_%i.dat",getenv("MWPC_CAL"),octetNumber).Data());

  mwpcFile << "CalFactorE:\t" << eastFactor << "\t" << eastFactorErr << std::endl
	   << "CalFactorW:\t" << westFactor << "\t" << westFactorErr << std::endl << std::endl;
  
  mwpcFile << "#EastData" << std::setw(12)
	   << "EastDataErr"  << std::setw(12)
	   << "WestData"  << std::setw(12)
	   << "WestDataErr" << std::setw(12)
	   << "EastSim" << std::setw(12)
	   << "EastSimErr"  << std::setw(12)
	   << "WestSim"  << std::setw(12)
	   << "WestSimErr" << std::endl;
  
  mwpcFile << std::setprecision(7);
  
  for ( int h = 0; h<nEnergyBins; ++h ) {
    
    mwpcFile << EdataMPV[h] << std::setw(12)
	     << EdataMPVerr[h] << std::setw(12)
	     << WdataMPV[h] << std::setw(12)
	     << WdataMPVerr[h] << std::setw(12)
	     << EsimMPV[h] << std::setw(12)
	     << EsimMPVerr[h] << std::setw(12)
	     << WsimMPV[h] << std::setw(12)
	     << WsimMPVerr[h];
    
    if ( h<nEnergyBins ) mwpcFile << std::endl;
    
  }

  mwpcFile.close();
  
  return 0;
}
