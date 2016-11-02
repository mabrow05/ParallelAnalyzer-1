/* 
    Produce the reference background spectra, as well as files which hold
    reference errors for low count bins
*/

#include "posMapReader.h"
#include "positionMapHandler.hh"
#include "sourcePeaks.h"
#include "runInfo.h"
#include "calibrationTools.hh"
#include "TriggerMap.hh"
#include "../Asymmetry/SQLinterface.hh"

#include "MBUtils.hh"

#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>

#include <TRandom3.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH1D.h>
#include <TChain.h>
#include <TMath.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TString.h>
#include <TStyle.h>


using namespace std;

std::map<Int_t,std::string> runType; //List of all runs in octet and their types
std::vector <Int_t> bgRuns_SFon; //list of bg runs with spin flipper on 
std::vector <Int_t> bgRuns_SFoff; //list of bg runs with spin flipper off




int separate23(int side, double mwpcEn) {
  int type = 2;
  if (side==0) {
    type = ( mwpcEn>4.14 ) ? 3 : 2;
  }
  
  if (side==1) {
    type = ( mwpcEn>4.14 ) ? 3 : 2;
  }
  return type;
};

void writeCorrectionFactorsToFile(Int_t octet) {
  TString fn_base = TString::Format("DeltaExp_OctetByOctetCorrections/ThOverProc_Octet-%i_Analysis-",octet); 
  ofstream oct; 
  ofstream quartA;
  ofstream quartB; 
  ofstream pairA1; 
  ofstream pairA2; 
  ofstream pairB1; 
  ofstream pairB2; 

  for (int i=1; i<11; i++) {
    oct.open(TString::Format("%s%i.txt",fn_base.Data(),i));
    //quartA.open(TString::Format("%s%i_quartA.txt",fn_base.Data(),i));
    //quartB.open(TString::Format("%s%i_quartB.txt",fn_base.Data(),i));
    //pairA1.open(TString::Format("%s%i_pairA1.txt",fn_base.Data(),i));
    //pairA2.open(TString::Format("%s%i_pairA2.txt",fn_base.Data(),i));
    //pairB1.open(TString::Format("%s%i_pairB1.txt",fn_base.Data(),i));
    //pairB2.open(TString::Format("%s%i_pairB2.txt",fn_base.Data(),i));

    for (int bin=0; bin<120; bin++) {
      Double_t binMidPoint = (double)bin*10.+5.;
      //Double_t A_theory = A0_PDG*asymmetryCorrectionFactor(binMidPoint)*beta(binMidPoint)/2.;
      //oct << binMidPoint << "\t" << ( fabs(oct_A_SR[i-1][bin] )>0.000001 ? A_theory/oct_A_SR[i-1][bin] : 1.) << "\t" 
      //	  << ( fabs(oct_A_SR[i-1][bin] )>0.000001 ? fabs(A_theory/power(oct_A_SR[i-1][bin],2)*oct_A_SR_err[i-1][bin]) : 1.) << "\n";
      /*quartA << binMidPoint << "\t" << A_theory/quartA_A_SR[i-1][bin] << "\t" << A_theory/power(quartA_A_SR[i-1][bin],2)*quartA_A_SR_err[i-1][bin] << "\n";
      quartB << binMidPoint << "\t" << A_theory/quartB_A_SR[i-1][bin] << "\t" << A_theory/power(quartB_A_SR[i-1][bin],2)*quartB_A_SR_err[i-1][bin] << "\n";
      pairA1 << binMidPoint << "\t" << A_theory/pairA1_A_SR[i-1][bin] << "\t" << A_theory/power(pairA1_A_SR[i-1][bin],2)*pairA1_A_SR_err[i-1][bin] << "\n";
      pairA2 << binMidPoint << "\t" << A_theory/pairA2_A_SR[i-1][bin] << "\t" << A_theory/power(pairA2_A_SR[i-1][bin],2)*pairA2_A_SR_err[i-1][bin] << "\n";
      pairB1 << binMidPoint << "\t" << A_theory/pairB1_A_SR[i-1][bin] << "\t" << A_theory/power(pairB1_A_SR[i-1][bin],2)*pairB1_A_SR_err[i-1][bin] << "\n";
      pairB2 << binMidPoint << "\t" << A_theory/pairB2_A_SR[i-1][bin] << "\t" << A_theory/power(pairB2_A_SR[i-1][bin],2)*pairB2_A_SR_err[i-1][bin] << "\n";
      */
    }
    oct.close();
    quartA.close();
    quartB.close();
    pairA1.close();
    pairA2.close();
    pairB1.close();
    pairB2.close();
    
  }
  std::cout << "Wrote All Corrections to file!\n";

};


void readOctetFileForBGruns(int octet) {
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  int runNumberHold;
  
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    runType[runNumberHold] = runTypeHold;

    if ( runTypeHold=="A1" || runTypeHold=="A12" || runTypeHold=="B4" || runTypeHold=="B9" ) {
      bgRuns_SFoff.push_back(runNumberHold);
    }
    if ( runTypeHold=="B1" || runTypeHold=="B12" || runTypeHold=="A4" || runTypeHold=="A9" )  {
      bgRuns_SFon.push_back(runNumberHold);
    }

  }

  infile.close();
 
  std::cout << "Read in octet file for octet " << octet << "\n";
};

std::string getRunTypeFromOctetFile(int run) {
  for (auto const& rt : runType) {
    if (rt.first == run) return rt.second;
  }
  
  return "BAD";
  
};


int getPolarization(int run) {
 
  std::string dbAddress = std::string(getenv("UCNADBADDRESS"));
  std::string dbname = std::string(getenv("UCNADB"));
  std::string dbUser = std::string(getenv("UCNADBUSER"));
  std::string dbPass = std::string(getenv("UCNADBPASS"));
  
  char cmd[200];
  sprintf(cmd,"SELECT flipper FROM run WHERE run_number=%i;",run);

  SQLdatabase *db = new SQLdatabase(dbname, dbAddress, dbUser, dbPass);
  db->fetchQuery(cmd);
  std::string flipperStatus = db->returnQueryEntry();
  delete db;

  std::cout << flipperStatus << std::endl;
  if (flipperStatus=="On") return 1;
  else if (flipperStatus=="Off") return -1;
  else {
    std::cout <<  "Polarization isn't applicaple or you chose a Depol Run";
    return 0;
  }
};




vector <Int_t> getPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"%s/residuals/PMT_runQuality_master.dat",getenv("ANALYSIS_CODE")); 
  ifstream pmt;
  std::cout << temp << std::endl;
  pmt.open(temp);
  Int_t run_hold;
  while (pmt >> run_hold >> pmtQuality[0] >> pmtQuality[1] >> pmtQuality[2]
	 >> pmtQuality[3] >> pmtQuality[4] >> pmtQuality[5]
	 >> pmtQuality[6] >> pmtQuality[7]) {
    if (run_hold==runNumber) break;
    if (pmt.fail()) break;
  }
  pmt.close();
  if (run_hold!=runNumber) {
    cout << "Run not found in PMT quality file!" << endl;
    exit(0);
  }
  return pmtQuality;
}



void doBackgroundSpectra (int octetMin, int octetMax) 
{

  gStyle->SetOptStat(0);
  double fiducialCut = 50.; //mm
 
  cout << "Calculating BG spectra..." << endl;

  //Load all the bg run numbers
  
  for ( int i=octetMin ; i<=octetMax ; i++ ) {

    readOctetFileForBGruns(i);

  }
   
  char temp[200];

  //Root file to output to
  TFile *outfile = new TFile(TString::Format("BackgroundSpectra_%i-%i.root",octetMin,octetMax),"RECREATE");

  //Make all the pertinent histograms

  TH1D* histOFF[3][2]; //We only go to type 2, then we separate them later and fill other histograms

  histOFF[0][0] = new TH1D("Type0E_sfOFF","Type0E",120,0.,1200.);
  histOFF[0][1] = new TH1D("Type0W_sfOFF","Type0W",120,0.,1200.);
  histOFF[1][0] = new TH1D("Type1E_sfOFF","Type1E",120,0.,1200.);
  histOFF[1][1] = new TH1D("Type1W_sfOFF","Type1W",120,0.,1200.);
  histOFF[2][0] = new TH1D("Type23E_sfOFF","Type23E",120,0.,1200.);
  histOFF[2][1] = new TH1D("Type23W_sfOFF","Type23W",120,0.,1200.);
  
  // Separated type 2
  TH1D* histOFF2[2];
  histOFF2[0] = new TH1D("Type2E_sfOFF","Type2E",120,0.,1200.);
  histOFF2[1] = new TH1D("Type2W_sfOFF","Type2W",120,0.,1200.);
  // Separated type 3
  TH1D* histOFF3[2];
  histOFF3[0] = new TH1D("Type3E_sfOFF","Type3E",120,0.,1200.);
  histOFF3[1] = new TH1D("Type3W_sfOFF","Type3W",120,0.,1200.);

  TH1D* histON[3][2]; //We only go to type 2, then we separate them later and fill other histograms

  histON[0][0] = new TH1D("Type0E_sfON","Type0E",120,0.,1200.);
  histON[0][1] = new TH1D("Type0W_sfON","Type0W",120,0.,1200.);
  histON[1][0] = new TH1D("Type1E_sfON","Type1E",120,0.,1200.);
  histON[1][1] = new TH1D("Type1W_sfON","Type1W",120,0.,1200.);
  histON[2][0] = new TH1D("Type23E_sfON","Type23E",120,0.,1200.);
  histON[2][1] = new TH1D("Type23W_sfON","Type23W",120,0.,1200.);
  
  // Separated type 2
  TH1D* histON2[2];
  histON2[0] = new TH1D("Type2E_sfON","Type2E",120,0.,1200.);
  histON2[1] = new TH1D("Type2W_sfON","Type2W",120,0.,1200.);
  // Separated type 3
  TH1D* histON3[2];
  histON3[0] = new TH1D("Type3E_sfON","Type3E",120,0.,1200.);
  histON3[1] = new TH1D("Type3W_sfON","Type3W",120,0.,1200.);
  

  double totalTimeOFF = 0.;
  double totalTimeON = 0.; //TODO

  //Process SF off runs first
  for ( auto rn : bgRuns_SFoff ) {
    
    sprintf(temp,"replay_pass3_%i.root",rn);
    std::string infile = getenv("REPLAY_PASS3")+std::string("/")+std::string(temp);
    TFile *input = new TFile(infile.c_str(), "READ");
    TTree *Tin = (TTree*)input->Get("pass3");


    double EmwpcX=0., EmwpcY=0., WmwpcX=0., WmwpcY=0., TimeE=0., TimeW=0., Time=0., Erecon=0.; //Branch Variables being read in
    int PID, type, side; // basic analysis tags

    Tin->SetBranchAddress("PID", &PID);
    Tin->SetBranchAddress("Type", &type);
    Tin->SetBranchAddress("Side", &side); 
    Tin->SetBranchAddress("Erecon",&Erecon);
    Tin->SetBranchAddress("TimeE",&TimeE);
    Tin->SetBranchAddress("TimeW",&TimeW);
    Tin->SetBranchAddress("Time",&Time);
    Tin->GetBranch("xE")->GetLeaf("center")->SetAddress(&EmwpcX);
    Tin->GetBranch("yE")->GetLeaf("center")->SetAddress(&EmwpcY);
    Tin->GetBranch("xW")->GetLeaf("center")->SetAddress(&WmwpcX);
    Tin->GetBranch("yW")->GetLeaf("center")->SetAddress(&WmwpcY);


    unsigned int nevents = Tin->GetEntriesFast();

    Tin->GetEvent(nevents-1);
    totalTimeOFF += Time;

    for (unsigned int n=0 ; n<nevents ; n++ ) {
      
      
      Tin->GetEvent(n);


      //Type 2/3 separation... ADD IN AFTER DOING MWPC CAL
      /*if (Erecon>0. && type==2) {
  
	if (side==0) {
	  type = separate23(side,mwpcE.MWPCEnergyE);
	  side = type==2 ? 1 : 0;
	}
	else if (side==1) {
	  type = separate23(side,mwpcE.MWPCEnergyW);
	  side = type==2 ? 0 : 1;
	}
	
	}*/

      
      //***********************************************************************************************************************
      // Filling rate histograms with "good" events to calculate the corrections

      
      if ( PID==1 && side<2 && type<4 && Erecon>0.) {

	double r2=EmwpcX*EmwpcX+EmwpcY*EmwpcY;

	if ( r2<(fiducialCut*fiducialCut) )	  {
		
	  //Type 0
	  if (type==0) histOFF[0][side]->Fill(Erecon);
	
	  //Type 1
	  if (type==1) histOFF[1][side]->Fill(Erecon);
	
	  //Type 23
	  if (type==2 || type==3) {
	    if (side==0) { 
	      histOFF[2][0]->Fill(Erecon);
	      //if (type==3) histOFF[2][0]->Fill(Erecon);
	      //else histOFF[2][1]->Fill(Erecon);
	    }
	    else if (side==1) {
	      histOFF[2][1]->Fill(Erecon);
	      //if (type==3) histOFF[2][1]->Fill(Erecon);
	      //else histOFF[2][0]->Fill(Erecon);
	    }
	  }
	
	  //Type 2
	  //if (type==2) histOFF2[side]->Fill(Erecon);
	
	  //Type 3
	  //if (type==3) histOFF3[side]->Fill(Erecon);

	}
      }
      
    }

    input->Close();
    if (input) delete input;
    cout << "Finished Run " << rn << endl;
  }

  //Process SF ON runs now
  for ( auto rn : bgRuns_SFon ) {
    
    sprintf(temp,"replay_pass3_%i.root",rn);
    std::string infile = getenv("REPLAY_PASS3")+std::string("/")+std::string(temp);
    TFile *input = new TFile(infile.c_str(), "READ");
    TTree *Tin = (TTree*)input->Get("pass3");


    double EmwpcX=0., EmwpcY=0., WmwpcX=0., WmwpcY=0., TimeE=0., TimeW=0., Time=0., Erecon=0.; //Branch Variables being read in
    int PID, type, side; // basic analysis tags

    Tin->SetBranchAddress("PID", &PID);
    Tin->SetBranchAddress("Type", &type);
    Tin->SetBranchAddress("Side", &side); 
    Tin->SetBranchAddress("Erecon",&Erecon);
    Tin->SetBranchAddress("TimeE",&TimeE);
    Tin->SetBranchAddress("TimeW",&TimeW);
    Tin->SetBranchAddress("Time",&Time);
    Tin->GetBranch("xE")->GetLeaf("center")->SetAddress(&EmwpcX);
    Tin->GetBranch("yE")->GetLeaf("center")->SetAddress(&EmwpcY);
    Tin->GetBranch("xW")->GetLeaf("center")->SetAddress(&WmwpcX);
    Tin->GetBranch("yW")->GetLeaf("center")->SetAddress(&WmwpcY);


    unsigned int nevents = Tin->GetEntriesFast();

    Tin->GetEvent(nevents-1);
    totalTimeON += Time;

    for (unsigned int n=0 ; n<nevents ; n++ ) {
      
      
      Tin->GetEvent(n);


      //Type 2/3 separation
      /*if (Erecon>0. && type==2) {
  
	if (side==0) {
	  type = separate23(side,mwpcE.MWPCEnergyE);
	  side = type==2 ? 1 : 0;
	}
	else if (side==1) {
	  type = separate23(side,mwpcE.MWPCEnergyW);
	  side = type==2 ? 0 : 1;
	}
	
	}*/

      
      //***********************************************************************************************************************
      // Filling rate histograms with "good" events to calculate the corrections

      
      if ( PID==1 && side<2 && type<4 && Erecon>0.) {

	double r2=EmwpcX*EmwpcX+EmwpcY*EmwpcY;

	if ( r2<(fiducialCut*fiducialCut) )	  {
		
	  //Type 0
	  if (type==0) histON[0][side]->Fill(Erecon);
	
	  //Type 1
	  if (type==1) histON[1][side]->Fill(Erecon);
	
	  //Type 23
	  if (type==2 || type==3) {
	    histON[2][side]->Fill(Erecon);
	    //if (side==0) { 
	    //histON[2][0]->Fill(Erecon);
	      //if (type==3) histON[2][0]->Fill(Erecon);
	      //else histON[2][1]->Fill(Erecon);
	    //}
	  //else if (side==1) {
	  //  histON[2][1]->Fill(Erecon);
	      //if (type==3) histON[2][1]->Fill(Erecon);
	      //else histON[2][0]->Fill(Erecon);
	    //}
	  }
	
	  //Type 2
	  //if (type==2) histON2[side]->Fill(Erecon);
	
	  //Type 3
	  //if (type==3) histON3[side]->Fill(Erecon);
    
	}
      }
      
    }

    input->Close();
    if (input) delete input;
    cout << "Finished Run " << rn << endl;
  }
  
  // Scale by 10^3mHz / (10 keV per bin) / total Time
  for (int s = 0; s<2; s++) {
    for (int t=0; t<3; t++) {

      histOFF[t][s]->SetYTitle("mHz/keV");
      histON[t][s]->SetYTitle("mHz/keV");
      
      histOFF[t][s]->Scale(1.e2/totalTimeOFF);
      histON[t][s]->Scale(1.e2/totalTimeON);
    }

    histON2[s]->SetYTitle("mHz/keV");
    histON3[s]->SetYTitle("mHz/keV");
    histOFF2[s]->SetYTitle("mHz/keV");
    histOFF3[s]->SetYTitle("mHz/keV");

    histON2[s]->Scale(1.e2/totalTimeON);
    histON3[s]->Scale(1.e2/totalTimeON);
    histOFF2[s]->Scale(1.e2/totalTimeOFF);
    histOFF3[s]->Scale(1.e2/totalTimeOFF);
  }
  
  
  

  outfile->Write();
  outfile->Close();
  cout << endl;

   
  /* for (int anaCh = 1; anaCh<11 ; anaCh++) {

      int type_low=-1, type_high=-1;
      bool sep23 = false;
      bool sep2 = false;
      bool sep3 = false;

      if (anaCh==1) { type_low=0; type_high=2; }                             //All types, 2/3 not separated
      else if (anaCh==3 || anaCh==5) { type_low=0; type_high=1; sep23=true;} // All event types, 2/3 separated
      else if (anaCh==2) { type_low=0; type_high=1;}                         // Type 0 and 1
      else if (anaCh==4) { type_low=0; type_high=0;}                         // Type 0
      else if (anaCh==6) { type_low=1; type_high=1;}                         // Type 1
      else if (anaCh==7) { type_low=type_high=2; }                           // Type 2/3 not separated
      else if (anaCh==8) { sep23 = sep2 = sep3 = true; }                 // Type 2/3, separated
      else if (anaCh==9) { sep23 = sep2 = true;}                   // Type 2
      else if (anaCh==10) { sep23 = sep3 = true;}                  // Type 3
      
      for (unsigned int side=0; side<2; side++) {
	for (unsigned int bin=1; bin<=120; bin++) {
	  for (int type=type_low; type<=type_high; type++) {	 

	    if ( runType == "A2" ) {
	      A2[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      A2_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A5" ) {
	      A5[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      A5_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A7" ) {
	      A7[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      A7_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A10" ) {
	      A10[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      A10_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B2" ) {
	      B2[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      B2_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B5" ) {
	      B5[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      B5_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B7" ) {
	      B7[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      B7_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B10" ) {
	      B10[anaCh-1][side][bin-1] += type>-1 ? hist[type][side]->GetBinContent(bin) : 0.;
	      B10_err[anaCh-1][side][bin-1] += type>-1 ? power(hist[type][side]->GetBinError(bin),2) : 0.;
	    }

	  }
	  
	  if ( sep23 ) {
	    if ( runType == "A2" ) {
	      A2[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      A2_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      A2[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      A2_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A5" ) {
	      A5[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      A5_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      A5[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      A5_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A7" ) {
	      A7[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      A7_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      A7[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      A7_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "A10" ) {
	      A10[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      A10_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      A10[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      A10_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B2" ) {
	      B2[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      B2_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      B2[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      B2_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B5" ) {
	      B5[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      B5_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      B5[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      B5_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B7" ) {
	      B7[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      B7_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      B7[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      B7_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	    if ( runType == "B10" ) {
	      B10[anaCh-1][side][bin-1] += sep2 ? hist2[side]->GetBinContent(bin) : 0.;
	      B10_err[anaCh-1][side][bin-1] += sep2 ? power(hist2[side]->GetBinError(bin),2) : 0.;
	      B10[anaCh-1][side][bin-1] += sep3 ? hist3[side]->GetBinContent(bin) : 0.;
	      B10_err[anaCh-1][side][bin-1] += sep3 ? power(hist3[side]->GetBinError(bin),2) : 0.;
	    }
	  }
	    
	  
	  A2_err[anaCh-1][side][bin-1] = sqrt(A2_err[anaCh-1][side][bin-1]);
	  A5_err[anaCh-1][side][bin-1] = sqrt(A5_err[anaCh-1][side][bin-1]);
	  A7_err[anaCh-1][side][bin-1] = sqrt(A7_err[anaCh-1][side][bin-1]);
	  A10_err[anaCh-1][side][bin-1] = sqrt(A10_err[anaCh-1][side][bin-1]);
	  B2_err[anaCh-1][side][bin-1] = sqrt(B2_err[anaCh-1][side][bin-1]);
	  B5_err[anaCh-1][side][bin-1] = sqrt(B5_err[anaCh-1][side][bin-1]);
	  B7_err[anaCh-1][side][bin-1] = sqrt(B7_err[anaCh-1][side][bin-1]);
	  B10_err[anaCh-1][side][bin-1]=sqrt(B10_err[anaCh-1][side][bin-1]);
	  //std::cout << anaChoice_A2[side][bin] << " " << anaChoice_A2_err[side][bin] << std::endl;
	}
      }

    }
  */
  /*for (int side = 0; side<2; side++) {
    for (int type=0; type<3; type++) {
      delete histOFF[type][side];
      delete histON[type][side];
    }
    delete histOFF2[side]; delete histOFF3[side]; delete histON2[side]; delete histON3[side];
    }*/
  

    
};
  

int main(int argc, char *argv[]) {

  int octetMin = atoi(argv[1]);
  int octetMax = atoi(argv[2]);

  doBackgroundSpectra(octetMin, octetMax);
  //writeCorrectionFactorsToFile(octet);

  //tests
  /*UInt_t XePeriod = getXeRunPeriod(atoi(argv[1]));
  vector < vector <Double_t> > triggerFunc = getTriggerFunctionParams(XePeriod,7);
  Double_t triggProbE = triggerProbability(triggerFunc[0],25.);
  Double_t triggProbW = triggerProbability(triggerFunc[1],25.);
  cout << triggProbE << " " << triggProbW << endl;*/


}
  
