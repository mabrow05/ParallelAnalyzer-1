/* Code to take a run number, retrieve it's runperiod, and contruct the 
weighted spectra which would be seen as a reconstructed energy on one side
of the detector */

#include <vector>

UInt_t getRunPeriod(Int_t runNumber) {
  UInt_t calibrationPeriod=0;
  if (runNumber <= 17297) {
    calibrationPeriod=1;
  }
  else if (runNumber <= 17439) {
    calibrationPeriod=2;
  }
  else if (runNumber <= 17734) {
    calibrationPeriod=3;
  }
  else if (runNumber <= 17955) {
    calibrationPeriod=4;
  }
  else if (runNumber <= 18386) {
    calibrationPeriod=5;
  }
  else if (runNumber <= 18683) {
    calibrationPeriod=6;
  }
  else if (runNumber <= 18994) {
    calibrationPeriod=7;
  }
  else if (runNumber <= 19239) {
    calibrationPeriod=8;
  }
  else if (runNumber <= 19544) {
    calibrationPeriod=9;
  }
  else if (runNumber <= 20000) {
    calibrationPeriod=11;
  }
  return calibrationPeriod;
}

vector <Int_t> getPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"../../residuals/PMT_runQuality_master.dat"); 
  ifstream pmt;
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

vector < Double_t > GetAlphaValues(Int_t runPeriod)
{
  Char_t temp[500];
  vector < Double_t > alphas (8,0.);
  sprintf(temp,"../nPE_Kev_%i.dat",runPeriod);
  ifstream infile;
  infile.open(temp);
  Int_t i = 0;

  while (infile >> alphas[i]) i++;
  return alphas;
}

void weightPeaks (Int_t runNumber, string source) 
{

  if (source!="Ce139" && source!="Sn113" && source!="Bi207") {
    cout << "Source Options are: \"Ce139\", \"Sn113\", \"Bi207\"\n";
    exit(0);}

  Char_t outputfile[500];
  sprintf(outputfile,"%i_%s_weightedSimPeaks.root",runNumber,source.c_str());
  TFile *outfile = new TFile(outputfile, "RECREATE");
  vector <Int_t> pmtQuality = getPMTQuality(runNumber); // Get the quality of the PMTs for that run
  UInt_t calibrationPeriod = getRunPeriod(runNumber); // retrieve the calibration period for this run
  vector <Double_t> alphas = GetAlphaValues(calibrationPeriod); // fill vector with the alpha (nPE/keV) values for this run period


  //Fill histograms with smeared and weighted energies, add them after
  vector <TH1D*> pmt (8,0);
  //Final histograms which are the addition of the previous histrograms
  TH1D *finalE = new TH1D("finalE", "Simulated Weighted Sum East", 400, 0., 1200.);
  TH1D *finalW = new TH1D("finalW", "Simulated Weighted Sum West", 400, 0., 1200.);
  Char_t name[200];
  Char_t file[500];

  Char_t temp[500];
  TChain *chain = new TChain("anaTree");
  for (int i=0; i<1; i++) {
    //  sprintf(temp,"/extern/mabrow05/ucna/geant4work/output/10mil_2011-2012/%s/analyzed_%i.root",source.c_str(),i);
    sprintf(temp,"../../../data/analyzed_%i.root",i);
    chain->AddFile(temp);
  }
  Int_t PID;
  //Double_t EdepQ[2], MWPCEnergy[2];
  Double_t EdepQE, EdepQW, MWPCEnergyE, MWPCEnergyW;
  Double_t E_sm[8]={0.}; //Hold the smeared energy
  Double_t weight[8]={0.}; // Hold the weight calculated based upon alpha
  //chain->SetBranchAddress("EdepQ",EdepQ);
  //chain->SetBranchAddress("MWPCEnergy",MWPCEnergy);
  chain->GetBranch("EdepQ")->GetLeaf("EdepQE")->SetAddress(&EdepQE);
  chain->GetBranch("EdepQ")->GetLeaf("EdepQW")->SetAddress(&EdepQW);
  chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyE")->SetAddress(&MWPCEnergyE);
  chain->GetBranch("MWPCEnergy")->GetLeaf("MWPCEnergyW")->SetAddress(&MWPCEnergyW);

  TRandom3 *rand = new TRandom3(0);
  Int_t nevents = chain->GetEntries();
  cout << "events = " << nevents << endl;
  
  //Create PMT histograms
  for (UInt_t p=0; p<8; p++) {
    sprintf(name,"PMT %i",p);
    pmt[p] = new TH1D(name,name, 400, 0., 1200.);
  }
  
  //Read in type 0 events and fill respective PMTs with smeared and weighted values
  for (Int_t evt=0; evt<nevents; evt++) {
    chain->GetEvent(evt);
    //Check for east type 0
    if (EdepQE>0. && EdepQW==0. && MWPCEnergyE>0. && MWPCEnergyW==0.) {
      for (UInt_t p=0; p<4; p++) {
	E_sm[p] = rand->Gaus(EdepQE,sqrt(EdepQE/alphas[p]));
	if (pmtQuality[p] && E_sm[p]>0.) {
	  weight[p] = alphas[p]/E_sm[p];
	  //cout << p << " " << E_sm << " " << weight << endl;
	  pmt[p]->Fill(E_sm[p]);
	}
	else {weight[p]=0.; E_sm[p]=0.;}
      }
      Double_t numer=0., denom=0.;
      for (UInt_t p=0;p<4;p++) {
	numer+=E_sm[p]*weight[p];
	denom+=weight[p];
      }
      finalE->Fill(numer/denom);
    }
    //Check for West type 0
    else if (EdepQW>0. && EdepQE==0. && MWPCEnergyW>0. && MWPCEnergyE==0.) {
      for (UInt_t p=4; p<8; p++) {
	E_sm[p] = rand->Gaus(EdepQW,sqrt(EdepQW/alphas[p]));
	if (pmtQuality[p] && E_sm[p]>0.) {
	  weight[p] = alphas[p]/E_sm[p];
	  //cout << p << " " << E_sm << " " << weight << endl;
	  pmt[p]->Fill(E_sm[p]);
	}
	else {weight[p]=0.; E_sm[p]=0.;}
      }
      Double_t numer=0., denom=0.;
      for (UInt_t p=4;p<8;p++) {
	numer+=E_sm[p]*weight[p];
	denom+=weight[p];
      }
      //cout << numer/denom << endl;
      finalW->Fill(numer/denom);
    }
    //cout << "filled event " << evt << endl;
  }
  
  outfile->Write();
  outfile->Close();
  
  //TCanvas *c1 = new TCanvas();
  //finalE->Draw();

  //TCanvas *c2 = new TCanvas();
  //finalW->Draw();

  /*    h->Draw();

    Int_t maxBin = h->GetMaximumBin();
    Double_t peak = h->GetXaxis()->GetBinCenter(maxBin);
    Double_t maxBinContent = h->GetBinContent(maxBin);
    Double_t high = 0., low=0.;
    
    for (int i=maxBin; i<400; i++) {
      if (h->GetBinContent(i+1) < 0.5*maxBinContent) {
	high= h->GetXaxis()->GetBinCenter(i+1);
	break;
      }
    }
    for (int i=maxBin; i<400; i--) {
      if (h->GetBinContent(i-1) < 0.5*maxBinContent) {
	low= h->GetXaxis()->GetBinCenter(i-1);
	break;
      }
    }
    
    
    TF1 *func = new TF1("func", "gaus", low, high);
    //h->SetParameter(1,peak);
    h->Fit("func", "LRQ");
    cout << alpha << " \t" << func->GetParameter(1) << endl;
    }*/
}
