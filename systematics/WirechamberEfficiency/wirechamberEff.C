


void wirechamberEff() {


  int runmin=19899;
  int runmax=19960;

  TChain *chain = new TChain("pass3");

  for (int run=runmin; run<=runmax; ++run) {
    
    chain->Add(TString::Format("%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),run));

  }
  
  TH1D *triggE = new TH1D("triggE","East wirechamber trigger",100,280.,440.);
  TH1D *notriggE = new TH1D("notriggE","East no wirechamber trigger",100,280.,440.);

  TH1D *triggW = new TH1D("triggW","West wirechamber trigger",100,280.,440.);
  TH1D *notriggW = new TH1D("notriggW","West no wirechamber trigger",100,280.,440.);

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(2,2);
  
  c1->cd(1);
  chain->Draw("Erecon>>triggE","PID==1 && Side==0 && Erecon>320. && Erecon<400.");
  c1->cd(2);
  chain->Draw("Erecon>>notriggE","PID==0 && Side==0 && Erecon>320. && Erecon<400.");
  c1->cd(3);
  chain->Draw("Erecon>>triggW","PID==1 && Side==1 && Erecon>320. && Erecon<400.");
  c1->cd(4);
  chain->Draw("Erecon>>notriggW","PID==0 && Side==1 && Erecon>320. && Erecon<400.");

  Double_t nEastTrigg = triggE->Integral();
  Double_t nEastNoTrigg = notriggE->Integral();
  Double_t nEastTotal = nEastTrigg+nEastNoTrigg;
    
  cout << "East inefficiency = " << nEastNoTrigg/nEastTotal << endl;

  Double_t nWestTrigg = triggW->Integral();
  Double_t nWestNoTrigg = notriggW->Integral();
  Double_t nWestTotal = nWestTrigg+nWestNoTrigg;
    
  cout << "West inefficiency = " << nWestNoTrigg/nWestTotal << endl;



}
