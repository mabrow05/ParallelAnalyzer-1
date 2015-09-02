//Code to study the effect of changing alpha on the peaks in the 
// Bi spectrum
{
  Double_t alphaLow = 0.001;
  Double_t alphaHigh = 0.1;
  Double_t inc = 0.0005;
  Double_t alpha_comp = 5.; //used to determine peak to be compared against

  ofstream outfile;
  outfile.open("alphaTrack_0.001-0.1.dat");

  Char_t temp[500];
  TChain *chain = new TChain("Evts");
  for (int i=0; i<1000; i++) {
    sprintf(temp,"/extern/mabrow05/ucna/geant4work/output/10mil_src_events/Bi207/Evts_%i.root",i);
    chain->AddFile(temp);
  }
  Int_t PID; 
  Double_t KE;
  chain->SetBranchAddress("KE",&KE);
  chain->SetBranchAddress("PID",&PID);
  TRandom3 *rand = new TRandom3(0);
  Int_t nevents = chain->GetEntries();
  cout << "events = " << nevents << endl;
  
  //////////////////////////////////////////////////////////////////////
  //Determine a peak to compare against. Set alpha to high number
  TH1D *h = new TH1D("h", "h", 400, 0., 1600.);  
  for (Int_t evt=0; evt<nevents; evt++) {
    chain->GetEvent(evt);
    if (PID==11 && KE>200.) {
      h->Fill(rand->Gaus(KE,sqrt(KE/alpha_comp)));
      //cout << "filled event " << evt << endl;
    }
  }
  h->Draw();
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
  cout << alpha_comp << " \t" << func->GetParameter(1) << endl;
  Double_t comp_peak = 993.789;//func->GetParameter(1);
  delete func;
  delete h;
  /////////////////////////////////////////////////////////////////////

  for (double alpha = alphaLow; alpha<alphaHigh+inc; alpha+=inc) {
    TH1D *h = new TH1D("h", "h", 400, 0., 1600.);  
    for (Int_t evt=0; evt<nevents; evt++) {
      chain->GetEvent(evt);
      if (PID==11 && KE>200.) {
	h->Fill(rand->Gaus(KE,sqrt(KE/alpha)));
	//cout << "filled event " << evt << endl;
      }
    }
    h->Draw();
    maxBin = h->GetMaximumBin();
    peak = h->GetXaxis()->GetBinCenter(maxBin);
    maxBinContent = h->GetBinContent(maxBin);
    high = 0.; 
    low=0.;
    
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
    outfile << alpha << " \t" << comp_peak - func->GetParameter(1) << endl;
    delete func;
    delete h;
  }
  outfile.close();
}
