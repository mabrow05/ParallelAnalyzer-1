
Double_t getVud(bool bottle, Double_t x) {
  Double_t tau = bottle?878.:915.;
  Double_t chi = 1.0390*1.6887; // code rest of constants... make sense of units on G_f
  return 0;
}


void Vud_vs_lambda() {

  Double_t Vud = 0.974;
  Double_t Vud_err = 0.0001;

  Double_t Lambda1 = 1.2755;
  Double_t Lambda1_err = 0.0005;
  Double_t Lambda2 = 1.263;
  Double_t Lambda2_err = 0.002;
  

  Int_t fillStyle_vud = 3001;
  Int_t fillStyle_lambda1 = 3001;
  Int_t fillStyle_lambda2 = 3001;
  Int_t fillColor_vud = 1;
  Int_t fillColor_lambda1 = 3;
  Int_t fillColor_lambda2 = 3;

  
  gStyle->SetPadLeftMargin(1.2);
  //gStyle->SetPadRightMargin(1.2);
  gStyle->SetPadBottomMargin(1.2);
  
  TCanvas *c1 = new TCanvas("c1","c1");

  TPad *pad1 = new TPad("pad1","",0,0,1,1);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetFillStyle(4000);
  pad2->SetFrameFillStyle(0);
  pad2->Draw();

  pad1->cd();
  
  std::vector <Double_t> vud_x {0.,2.};
  std::vector <Double_t> vud_y {Vud,Vud};
  std::vector <Double_t> vud_yerr {Vud_err,Vud_err};
  TGraphErrors *vud = new TGraphErrors(2,&vud_x[0],&vud_y[0],0,&vud_yerr[0]);
  vud->SetFillStyle(fillStyle_vud);
  vud->SetFillColor(fillColor_vud);

  std::vector <Double_t> lambda1_x {Lambda1-Lambda1_err,Lambda1+Lambda1_err};
  std::vector <Double_t> lambda1_y {0.974,0.974};
  std::vector <Double_t> lambda1_yerr {1.,1.};
  TGraphErrors *lambda1 = new TGraphErrors(2,&lambda1_x[0],&lambda1_y[0],0,&lambda1_yerr[0]);
  lambda1->SetFillStyle(fillStyle_lambda1);
  lambda1->SetFillColor(fillColor_lambda1);


  std::vector <Double_t> lambda2_x {Lambda2-Lambda2_err,Lambda2+Lambda2_err};
  std::vector <Double_t> lambda2_y {0.974,0.974};
  std::vector <Double_t> lambda2_yerr {1.,1.};
  TGraphErrors *lambda2 = new TGraphErrors(2,&lambda2_x[0],&lambda2_y[0],0,&lambda2_yerr[0]);
  lambda2->SetFillStyle(fillStyle_lambda2);
  lambda2->SetFillColor(fillColor_lambda2);


 
  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(vud,"3");
  mg->Add(lambda1,"3");
  mg->Add(lambda2,"3");


  mg->Draw("A");
  mg->GetXaxis()->SetTitle("|#lambda|");
  mg->GetXaxis()->SetTitleOffset(1.2);
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("V_{ud}");
  mg->GetYaxis()->SetTitleOffset(1.4);
  mg->GetYaxis()->CenterTitle();

  mg->GetXaxis()->SetLimits(1.255,1.28);
  mg->SetMinimum(0.968);
  mg->SetMaximum(0.982);

  c1->Update();


  pad2->cd();
  
  std::vector <Double_t> lambdaMeas_x {1.2783};
  std::vector <Double_t> lambdaMeas_y {2017};
  std::vector <Double_t> lambdaMeas_xerr {0.0022};
  
  TGraphErrors *lambdaMeas = new TGraphErrors(1,&lambdaMeas_x[0],&lambdaMeas_y[0],&lambdaMeas_xerr[0],0);
  lambdaMeas->SetMarkerColor(fillColor_lambda1);
  lambdaMeas->SetMarkerStyle(24);
  lambdaMeas->Draw("APY+");
  lambdaMeas->GetXaxis()->SetLimits(1.255,1.28);
  lambdaMeas->GetYaxis()->SetTitle("Year");
  lambdaMeas->GetYaxis()->CenterTitle();
  lambdaMeas->GetYaxis()->SetTitleOffset(1.2);
  lambdaMeas->SetMinimum(1960);
  lambdaMeas->SetMaximum(2020);
  c1->Update();

}
