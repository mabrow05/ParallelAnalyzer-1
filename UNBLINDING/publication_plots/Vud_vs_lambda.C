#include <vector>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TAxis.h>
#include <TMultiGraph.h>

Double_t weightedAve(std::vector<Double_t>& val, std::vector<Double_t>& err) { 
  Double_t numer = 0.;
  Double_t denom = 0.;
  for (UInt_t i=0;i<val.size();++i) {
    numer+=(val[i]/err[i]/err[i]);
    denom+=(1./err[i]/err[i]);
  }
  return numer/denom;
}

Double_t weightedAveErr(std::vector<Double_t>& err) { 
  Double_t denom = 0.;
  for (UInt_t i=0;i<err.size();++i) {
    denom+=(1./err[i]/err[i]);
  }
  return TMath::Sqrt(1./denom);
}

Double_t getVud(Double_t x, Double_t tau) {
  //Double_t tau = bottle?878.:915.;
  Double_t chi = 4908.7;
  return TMath::Sqrt(chi/(tau*(1.+3.*x*x)));
}

Double_t getVudUncert(Double_t x, Double_t tau, Double_t dtau) {
  //Double_t tau = bottle?878.:915.;
  Double_t chi = 4908.7;
  //Double_t dtau = bottle?1.:1.;
  Double_t dchi = 1.9; // code rest of constants... make sense of units on G_f
  return TMath::Sqrt(TMath::Power(dtau/tau,2.)+TMath::Power(dchi/chi,2.))*getVud(x,tau)/2.;
}		     
		     


void Vud_vs_lambda() {

  bool color = true;
  
  Int_t fillStyle_vud = color?3002:3002;//3345;//3004;//3018;
  Int_t fillStyle_lambda1 = color?3002:3002;//3354;//3005;//3017;
  Int_t fillStyle_lambda2 = color?3002:3002;//3354//3005;//3017;
  Int_t fillStyle_tau1 = color?3002:3003;//3002;
  Int_t fillStyle_tau2 = color?3002:3003;//3002;
  Int_t fillColor_vud = color?2:13; 
  Int_t fillColor_lambda1 = color?1:1;//12;
  Int_t fillColor_lambda2 = color?1:1;//12;
  Int_t fillColor_tau1 = color?4:1;//14;
  Int_t fillColor_tau2 = color?4:1;//14;

  Double_t xmin = 1.253;
  Double_t xmax = 1.282;
  
  gStyle->SetPadLeftMargin(0.125);
  gStyle->SetPadRightMargin(0.125);
  gStyle->SetPadBottomMargin(0.115);
  gStyle->SetTitleSize(0.05,"xy");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetFillStyle(0);

  Double_t Vud = 0.97417;
  Double_t Vud_err = 0.00021;

  std::vector<TString> Lambda1name {"Abele et al.","Mund et al.","Mendenhall et al.","Brown et al.","Schumann et al."};
  std::vector<Double_t> Lambda1year {2002,2012,2013,2017,2008};
  std::vector<Double_t> Lambda1vec {1.2739,1.2761,1.2756,1.2783,1.275};
  std::vector<Double_t> Lambda1vecErr {0.0019,0.0017,0.0030,0.0022,0.016};
  Double_t Lambda1 = weightedAve(Lambda1vec,Lambda1vecErr);//1.2755;
  Double_t Lambda1_err = weightedAveErr(Lambda1vecErr);//0.0005;

  std::vector<TString> Lambda2name{"Bopp et al.","Yerozolimsky et al.","Liaud et al.","Mostovoi et al."};
  std::vector<Double_t> Lambda2year {1986,1997,1997.5,2001};
  std::vector<Double_t> Lambda2vec {1.262,1.2594,1.266,1.2686};
  std::vector<Double_t> Lambda2vecErr {0.005,0.0038,0.004,0.0047};
  Double_t Lambda2 = weightedAve(Lambda2vec,Lambda2vecErr);
  Double_t Lambda2_err = weightedAveErr(Lambda2vecErr);


  std::vector<Double_t> Tau1vec {878.5};
  std::vector<Double_t> Tau1vecErr {0.8}; 
  Double_t Tau1 = weightedAve(Tau1vec,Tau1vecErr);
  Double_t dTau1 = weightedAveErr(Tau1vecErr);

  std::vector<Double_t> Tau2vec {887.7};
  std::vector<Double_t> Tau2vecErr {2.2}; 
  Double_t Tau2 = weightedAve(Tau2vec,Tau2vecErr);
  Double_t dTau2 = weightedAveErr(Tau2vecErr);
  
  
  
  
  
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


  std::vector <Double_t> tau1_y;
  std::vector <Double_t> tau1_yerr;
  std::vector <Double_t> tau1_x;
  std::vector <Double_t> tau2_y;
  std::vector <Double_t> tau2_yerr;
  std::vector <Double_t> tau2_x;
  
  for (Double_t i=1.2;i<1.3;i+=0.0001) {
    tau1_y.push_back(getVud(i,Tau1));
    tau1_yerr.push_back(getVudUncert(i,Tau1,dTau1));
    tau1_x.push_back(i);
    tau2_y.push_back(getVud(i,Tau2));
    tau2_yerr.push_back(getVudUncert(i,Tau2,dTau2));
    tau2_x.push_back(i);
    //std::cout << i << "\t" << getVud(i,Tau1) << "\t" << getVudUncert(i,Tau1,dTau1) << "\n";
  }

  TGraphErrors *tau1 = new TGraphErrors(tau1_x.size(),&tau1_x[0],&tau1_y[0],0,&tau1_yerr[0]);
  tau1->SetFillStyle(fillStyle_tau1);
  tau1->SetFillColor(fillColor_tau1);

  TGraphErrors *tau2 = new TGraphErrors(tau2_x.size(),&tau2_x[0],&tau2_y[0],0,&tau2_yerr[0]);
  tau2->SetFillStyle(fillStyle_tau2);
  tau2->SetFillColor(fillColor_tau2);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("");
  mg->Add(vud,"3");
  mg->Add(lambda1,"3");
  mg->Add(lambda2,"3");
  mg->Add(tau1,"3");
  mg->Add(tau2,"3");

  mg->Draw("A");
  mg->GetXaxis()->SetTitle("|#lambda|");
  //mg->GetXaxis()->SetNdivisions(505);
  mg->GetXaxis()->SetTitleOffset(1);
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("V_{ud}");
  mg->GetYaxis()->SetTitleOffset(1.1);
  mg->GetYaxis()->CenterTitle();

  mg->GetXaxis()->SetLimits(xmin,xmax);
  mg->SetMinimum(0.967);
  mg->SetMaximum(0.982);

  c1->Update();


  pad2->cd();
  
  TGraphErrors *lambdaMeas10 = new TGraphErrors(1,&Lambda1vec[0],&Lambda1year[0],&Lambda1vecErr[0],0);
  lambdaMeas10->SetTitle("");
  lambdaMeas10->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas10->SetMarkerStyle(23);
  lambdaMeas10->Draw("APY+");
  lambdaMeas10->GetXaxis()->SetLimits(xmin,xmax);
  //lambdaMeas10->GetXaxis()->SetNdivisions(505);
  lambdaMeas10->GetYaxis()->SetTitle("Year of Publication");
  lambdaMeas10->GetYaxis()->CenterTitle();
  lambdaMeas10->GetYaxis()->SetTitleOffset(1);
  lambdaMeas10->SetMinimum(1960);
  lambdaMeas10->SetMaximum(2020);

  TGraphErrors *lambdaMeas11 = new TGraphErrors(1,&Lambda1vec[1],&Lambda1year[1],&Lambda1vecErr[1],0);
  lambdaMeas11->SetTitle("");
  lambdaMeas11->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas11->SetMarkerStyle(22);
  lambdaMeas11->Draw("P SAME");

  TGraphErrors *lambdaMeas12 = new TGraphErrors(1,&Lambda1vec[2],&Lambda1year[2],&Lambda1vecErr[2],0);
  lambdaMeas12->SetTitle("");
  lambdaMeas12->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas12->SetMarkerStyle(21);
  lambdaMeas12->Draw("P SAME");

  TGraphErrors *lambdaMeas13 = new TGraphErrors(1,&Lambda1vec[3],&Lambda1year[3],&Lambda1vecErr[3],0);
  lambdaMeas13->SetTitle("");
  lambdaMeas13->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas13->SetMarkerStyle(20);
  lambdaMeas13->Draw("P SAME");

  TGraphErrors *lambdaMeas14 = new TGraphErrors(1,&Lambda1vec[4],&Lambda1year[4],&Lambda1vecErr[4],0);
  lambdaMeas14->SetTitle("");
  lambdaMeas14->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas14->SetMarkerStyle(34);
  lambdaMeas14->Draw("P SAME");

  TGraphErrors *lambdaMeas20 = new TGraphErrors(1,&Lambda2vec[0],&Lambda2year[0],&Lambda2vecErr[0],0);
  lambdaMeas20->SetTitle("");
  lambdaMeas20->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas20->SetMarkerStyle(24);//33);
  lambdaMeas20->Draw("P SAME");

  TGraphErrors *lambdaMeas21 = new TGraphErrors(1,&Lambda2vec[1],&Lambda2year[1],&Lambda2vecErr[1],0);
  lambdaMeas21->SetTitle("");
  lambdaMeas21->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas21->SetMarkerStyle(25);//);
  lambdaMeas21->Draw("P SAME");

  TGraphErrors *lambdaMeas22 = new TGraphErrors(1,&Lambda2vec[2],&Lambda2year[2],&Lambda2vecErr[2],0);
  lambdaMeas22->SetTitle("");
  lambdaMeas22->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas22->SetMarkerStyle(26);//34);
  lambdaMeas22->Draw("P SAME");

  TGraphErrors *lambdaMeas23 = new TGraphErrors(1,&Lambda2vec[3],&Lambda2year[3],&Lambda2vecErr[3],0);
  lambdaMeas23->SetTitle("");
  lambdaMeas23->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas23->SetMarkerStyle(28);//25);
  lambdaMeas23->Draw("P SAME");

  
  
  TLegend *legA = new TLegend(0.16,0.16,0.46,0.40);
  legA->SetTextSize(0.030);
  legA->SetHeader("A_{0} measurements");
  legA->AddEntry(lambdaMeas13,Lambda1name[3]+" (this work)","p");
  legA->AddEntry(lambdaMeas12,Lambda1name[2],"p");
  legA->AddEntry(lambdaMeas11,Lambda1name[1],"p");
  legA->AddEntry(lambdaMeas10,Lambda1name[0],"p");
  legA->AddEntry(lambdaMeas22,Lambda2name[2],"p"); 
  legA->AddEntry(lambdaMeas21,Lambda2name[1],"p"); 
  legA->AddEntry(lambdaMeas20,Lambda2name[0],"p");
  //legA->AddEntry(lambdaMeas14,Lambda1name[4],"p"); 
  //legA->AddEntry(lambdaMeas23,Lambda2name[3],"p");
  legA->Draw("SAME");

  TLegend *legB = new TLegend(0.44,0.16,0.74,0.26);
  legB->SetTextSize(0.030);
  legB->SetHeader("Other measurements");
  legB->AddEntry(lambdaMeas14,Lambda1name[4],"p");
  legB->AddEntry(lambdaMeas23,Lambda2name[3],"p"); 
  legB->Draw("SAME");

  TPaveText *tVud = new TPaveText(0.13,0.505,0.3,0.535,"nbNDC");
  tVud->SetBorderSize(0);
  tVud->SetTextColor(color?fillColor_vud:1);
  tVud->AddText("PDG 0^{+}#rightarrow0^{+}");
  tVud->GetLine(0)->SetTextSize(0.035);
  tVud->Draw();

  TPaveText *tTau1 = new TPaveText(0.452,0.79,0.592,0.83,"nbNDC");
  tTau1->SetBorderSize(0);
  tTau1->SetTextColor(color?fillColor_tau2:1);
  tTau1->AddText("Bottle #tau_{n}");
  tTau1->GetLine(0)->SetTextAngle(-42.2);
  tTau1->GetLine(0)->SetTextSize(0.035);
  tTau1->Draw();

  TPaveText *tTau2 = new TPaveText(0.13,0.79,0.27,0.83,"nbNDC");
  tTau2->SetBorderSize(0);
  tTau2->SetTextColor(color?fillColor_tau2:1);
  tTau2->AddText("Beam #tau_{n}");
  tTau2->GetLine(0)->SetTextAngle(-42.2);
  tTau2->GetLine(0)->SetTextSize(0.035);
  tTau2->Draw();
  
  c1->Update();
  c1->Print(TString::Format("vud_vs_lambda%s.pdf",color?"_color":""));
}
