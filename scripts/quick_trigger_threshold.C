

void quick_trigger_threshold(Int_t runNumber) {

  Char_t temp[500];
  sprintf(temp,"/extern/UCNA/replay_pass4_MB/replay_pass4_%i.root",runNumber);
  TFile *file = new TFile(temp,"READ");

  TTree *t = (TTree*)file->Get("pass4");

  Int_t nbins = 50;
  TH1F *Etype1 = new TH1F("Etype1","West Type 1: EvisE",nbins,0.,100.);
  TH1F *Etype23 = new TH1F("Etype23","West Type 2/3: EvisE",nbins,0.,100.);
  TH1F *Etotal = new TH1F("Etotal","West Type 1,2/3: EvisE",nbins,0.,100.);
  TH1F *Etrigg = new TH1F("Etrigg","East Trigger Probability",nbins,0.,100.);
  Etrigg->SetMarkerStyle(20);

  TH1F *Wtype1 = new TH1F("Wtype1","East Type 1: EvisW",nbins,0.,100.);
  TH1F *Wtype23 = new TH1F("Wtype23","East Type 2/3: EvisW",nbins,0.,100.);
  TH1F *Wtotal = new TH1F("Wtotal","East Type 1,2/3: EvisW",nbins,0.,100.);
  TH1F *Wtrigg = new TH1F("Wtrigg","West Trigger Probability",nbins,0.,100.);
  Wtrigg->SetMarkerStyle(20);


  TCanvas *c1 = new TCanvas("c1"," ",1200.,1600.);
  c1->Divide(2,2);
  c1->cd(1);
  t->Draw("EvisE>>Etype1","type_pass4==1 && side_pass4==1 && PID_pass4==1 && EvisE>0.");
  c1->cd(2);  
  t->Draw("EvisE>>Etype23","type_pass4==2 && side_pass4==1 && PID_pass4==1 && EvisE>0.");
  c1->cd(3);
  Etotal->Add(Etype1,Etype23);
  Etotal->Draw();
  c1->cd(4);
  Etrigg->Divide(Etype1,Etotal);
  Etrigg->SetStats(0);
  Etrigg->Draw("P");
  


  TCanvas *c2 = new TCanvas("c2"," ",1200.,1600.);
  c2->Divide(2,2);
  c2->cd(1);
  t->Draw("EvisW>>Wtype1","type_pass4==1 && side_pass4==0 && PID_pass4==1 && EvisW>0.");
  c2->cd(2);  
  t->Draw("EvisW>>Wtype23","type_pass4==2 && side_pass4==0 && PID_pass4==1 && EvisW>0.");
  c2->cd(3);
  Wtotal->Add(Wtype1,Wtype23);
  Wtotal->Draw();
  c2->cd(4);
  Wtrigg->Divide(Wtype1,Wtotal);
  Wtrigg->SetStats(0);
  Wtrigg->Draw("P");

}
