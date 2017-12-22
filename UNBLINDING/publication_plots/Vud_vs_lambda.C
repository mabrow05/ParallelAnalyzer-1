#include <vector>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TAxis.h>
#include <TMultiGraph.h>


std::vector<Double_t> calcLambda(Double_t A0,Double_t A0err) {
  Double_t l = (-1.-TMath::Sqrt(1.-3.*A0*A0-2.*A0))/(3.*A0+2.);
  Double_t lerr = TMath::Abs( ( (3.*A0+5.)/((3.*A0+2.)*(3.*A0+2.)*TMath::Sqrt(1.-3.*A0*A0-2.*A0)) + 3./(3.*A0+2.)/(3.*A0+2.) ) * A0err );
  return std::vector<Double_t>{l,lerr};
}

Double_t chi2nu(std::vector<Double_t>& val, std::vector<Double_t>& err) {
  Double_t numer = 0.;
  Double_t denom = 0.;
  for (UInt_t i=0;i<val.size();++i) {
    numer+=(val[i]/err[i]/err[i]);
    denom+=(1./err[i]/err[i]);
  }
  Double_t mu = numer/denom;
  Double_t c2 = 0;
  for (int i=0;i<val.size();++i) {
    //if (i!=0) std::cout << " + ";
    c2+= (val[i]-mu)*(val[i]-mu)/(err[i]*err[i]);
    //std::cout << "("<<val[i]<<"-"<<mu<<")^2/"<<err[i]<<"^2";
  }
  //std::cout << " = " << c2 << std::endl;
    
  return c2/(val.size()-1);
}

Double_t chi2(std::vector<Double_t>& val, std::vector<Double_t>& err) {
  Double_t numer = 0.;
  Double_t denom = 0.;
  for (UInt_t i=0;i<val.size();++i) {
    numer+=(val[i]/err[i]/err[i]);
    denom+=(1./err[i]/err[i]);
  }
  Double_t mu = numer/denom;
  Double_t c2 = 0;
  for (int i=0;i<val.size();++i) {
    if (i!=0) std::cout << " + ";
    c2+= (val[i]-mu)*(val[i]-mu)/(err[i]*err[i]);
    std::cout << "("<<val[i]<<"-"<<mu<<")^2/"<<err[i]<<"^2";
  }
  std::cout << " = " << c2 << std::endl;
  return c2;
}

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
  bool data = true;

  Double_t tau1scale = 1.38; //bottle
  Double_t tau2scale = 1.;   //beam
  
  Int_t fillStyle_vud = color?1001:1001;//3002;
  Int_t fillStyle_lambda1 = color?1001:1001;//3002;
  Int_t fillStyle_lambda2 = color?1001:1001;//3002;
  Int_t fillStyle_tau1 = color?1001:1001;//3003;
  Int_t fillStyle_tau2 = color?1001:1001;//3003;
  Int_t fillColor_vud = color?2:15;//13; 
  Int_t fillColor_lambda1 = color?14:16;//1;//12;
  Int_t fillColor_lambda2 = color?14:16;//1;//12;
  Int_t fillColor_tau1 = color?4:17;//1;//14;
  Int_t fillColor_tau2 = color?4:17;//1;//14;

  Double_t xmin = 1.253;
  Double_t xmax = 1.2805;
  
  gStyle->SetPadLeftMargin(0.130);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.125);
  gStyle->SetPadBottomMargin(0.115);
  gStyle->SetTitleSize(0.06,"xy");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetFillStyle(0);
  gStyle->SetLabelSize(0.042,"xyz");

  Double_t Vud = 0.97417;
  Double_t Vud_err = 0.00021;

  std::vector<TString> Lambda1name {"Mund #font[50]{et al.}","Brown #font[50]{et al.}","Schumann #font[50]{et al.}"};
  std::vector<Double_t> Lambda1year {2013,2017,2008};
  std::vector<Double_t> Lambda1vec {1.2748,1.2772,1.275}; //Brown result 1.2783 +/- 0.0022 // raw abele (1.2739) and Mund (1.2761)
  std::vector<Double_t> Lambda1vecErr {0.0014,0.0020,0.016};
  Double_t Lambda1 = weightedAve(Lambda1vec,Lambda1vecErr);//1.2755;
  Double_t Lambda1_err = weightedAveErr(Lambda1vecErr);//0.0005;

  std::vector<TString> Lambda2name{"Bopp #font[50]{et al.}","Yerozolimsky #font[50]{et al.}","Liaud #font[50]{et al.}","Mostovoi #font[50]{et al.}"};
  std::vector<Double_t> Lambda2year {1986,1997,1997.5,2001};
  std::vector<Double_t> Lambda2vec {1.262,1.2594,1.266,1.2686};
  std::vector<Double_t> Lambda2vecErr {0.005,0.0038,0.004,0.0047};
  Double_t Lambda2 = weightedAve(Lambda2vec,Lambda2vecErr);
  Double_t Lambda2_err = weightedAveErr(Lambda2vecErr);


  //UCN
  std::vector<TString> Tau1name {"Arzumanov #font[50]{et al.}","Steyerl #font[50]{et al.}","Pichlmaier #font[50]{et al.}","Serebrov #font[50]{et al.}","Mampe #font[50]{et al.}","Pattie #font[50]{et al.}"};
  std::vector<Double_t> Tau1vec {880.2,882.5,880.7,878.5,882.6,877.7};
  std::vector<Double_t> Tau1vecErr {1.2,2.1,1.8,0.8,2.7,0.8}; 
  Double_t Tau1 = weightedAve(Tau1vec,Tau1vecErr);
  Double_t dTau1 = weightedAveErr(Tau1vecErr);

  //beam
  std::vector<TString> Tau2name {"Yue #font[50]{et al.}","Byrne #font[50]{et al.}"};
  std::vector<Double_t> Tau2vec {887.7,889.2};
  std::vector<Double_t> Tau2vecErr {2.2,4.8}; 
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
  lambda1->SetFillColorAlpha(fillColor_lambda1,0.8);


  std::vector <Double_t> lambda2_x {Lambda2-Lambda2_err,Lambda2+Lambda2_err};
  std::vector <Double_t> lambda2_y {0.974,0.974};
  std::vector <Double_t> lambda2_yerr {1.,1.};
  TGraphErrors *lambda2 = new TGraphErrors(2,&lambda2_x[0],&lambda2_y[0],0,&lambda2_yerr[0]);
  lambda2->SetFillStyle(fillStyle_lambda2);
  lambda2->SetFillColorAlpha(fillColor_lambda2,0.8);


  std::vector <Double_t> tau1_y;
  std::vector <Double_t> tau1_yerr;
  std::vector <Double_t> tau1_x;
  std::vector <Double_t> tau2_y;
  std::vector <Double_t> tau2_yerr;
  std::vector <Double_t> tau2_x;
  
  for (Double_t i=1.2;i<1.3;i+=0.0001) {
    tau1_y.push_back(getVud(i,Tau1));
    tau1_yerr.push_back(getVudUncert(i,Tau1,dTau1)*tau1scale);
    tau1_x.push_back(i);
    tau2_y.push_back(getVud(i,Tau2));
    tau2_yerr.push_back(getVudUncert(i,Tau2,dTau2)*tau2scale);
    tau2_x.push_back(i);
    //std::cout << i << "\t" << getVud(i,Tau1) << "\t" << getVudUncert(i,Tau1,dTau1) << "\n";
  }

  TGraphErrors *tau1 = new TGraphErrors(tau1_x.size(),&tau1_x[0],&tau1_y[0],0,&tau1_yerr[0]);
  tau1->SetFillStyle(fillStyle_tau1);
  tau1->SetFillColorAlpha(fillColor_tau1,0.8);

  TGraphErrors *tau2 = new TGraphErrors(tau2_x.size(),&tau2_x[0],&tau2_y[0],0,&tau2_yerr[0]);
  tau2->SetFillStyle(fillStyle_tau2);
  tau2->SetFillColorAlpha(fillColor_tau2,0.8);
  
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
  mg->GetXaxis()->SetTitleOffset(0.85);
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("V_{ud}");
  mg->GetYaxis()->SetTitleOffset(1.05);
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
  if (data) {
    lambdaMeas10->Draw("APY+");
    lambdaMeas10->GetXaxis()->SetLimits(xmin,xmax);
    //lambdaMeas10->GetXaxis()->SetNdivisions(505);
    lambdaMeas10->GetYaxis()->SetTitle("Year of Publication");
    lambdaMeas10->GetYaxis()->CenterTitle();
    lambdaMeas10->GetYaxis()->SetTitleOffset(1);
    lambdaMeas10->SetMinimum(1960);
    lambdaMeas10->SetMaximum(2020);
  }

  TGraphErrors *lambdaMeas11 = new TGraphErrors(1,&Lambda1vec[1],&Lambda1year[1],&Lambda1vecErr[1],0);
  lambdaMeas11->SetTitle("");
  lambdaMeas11->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas11->SetMarkerStyle(20);
  if (data) lambdaMeas11->Draw("P SAME");

  TGraphErrors *lambdaMeas12 = new TGraphErrors(1,&Lambda1vec[2],&Lambda1year[2],&Lambda1vecErr[2],0);
  lambdaMeas12->SetTitle("");
  lambdaMeas12->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas12->SetMarkerStyle(21);
  if (data) lambdaMeas12->Draw("P SAME");

  /*TGraphErrors *lambdaMeas13 = new TGraphErrors(1,&Lambda1vec[3],&Lambda1year[3],&Lambda1vecErr[3],0);
  lambdaMeas13->SetTitle("");
  lambdaMeas13->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas13->SetMarkerStyle(20);
  if (data) lambdaMeas13->Draw("P SAME");

  TGraphErrors *lambdaMeas14 = new TGraphErrors(1,&Lambda1vec[4],&Lambda1year[4],&Lambda1vecErr[4],0);
  lambdaMeas14->SetTitle("");
  lambdaMeas14->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas14->SetMarkerStyle(34);
  if (data) lambdaMeas14->Draw("P SAME");*/
  
  TGraphErrors *lambdaMeas20 = new TGraphErrors(1,&Lambda2vec[0],&Lambda2year[0],&Lambda2vecErr[0],0);
  lambdaMeas20->SetTitle("");
  lambdaMeas20->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas20->SetMarkerStyle(24);//33);
  if (data) lambdaMeas20->Draw("P SAME");

  TGraphErrors *lambdaMeas21 = new TGraphErrors(1,&Lambda2vec[1],&Lambda2year[1],&Lambda2vecErr[1],0);
  lambdaMeas21->SetTitle("");
  lambdaMeas21->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas21->SetMarkerStyle(25);//);
  if (data) lambdaMeas21->Draw("P SAME");

  TGraphErrors *lambdaMeas22 = new TGraphErrors(1,&Lambda2vec[2],&Lambda2year[2],&Lambda2vecErr[2],0);
  lambdaMeas22->SetTitle("");
  lambdaMeas22->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas22->SetMarkerStyle(26);//34);
  if (data) lambdaMeas22->Draw("P SAME");

  TGraphErrors *lambdaMeas23 = new TGraphErrors(1,&Lambda2vec[3],&Lambda2year[3],&Lambda2vecErr[3],0);
  lambdaMeas23->SetTitle("");
  lambdaMeas23->SetMarkerColor(1);//(fillColor_lambda1);
  lambdaMeas23->SetMarkerStyle(28);//25);
  if (data) lambdaMeas23->Draw("P SAME");

  
  
  TLegend *legA = new TLegend(0.16,0.16,0.46,0.42);
  legA->SetTextSize(0.040);
  legA->SetHeader("A_{0} measurements");
  legA->AddEntry(lambdaMeas11,Lambda1name[1]+" (this work)","p");
  //legA->AddEntry(lambdaMeas11,Lambda1name[1],"p");
  legA->AddEntry(lambdaMeas10,Lambda1name[0],"p");
  //legA->AddEntry(lambdaMeas10,Lambda1name[0],"p");
  legA->AddEntry(lambdaMeas22,Lambda2name[2],"p"); 
  legA->AddEntry(lambdaMeas21,Lambda2name[1],"p"); 
  legA->AddEntry(lambdaMeas20,Lambda2name[0],"p");
  //legA->AddEntry(lambdaMeas14,Lambda1name[4],"p"); 
  //legA->AddEntry(lambdaMeas23,Lambda2name[3],"p");
  if (data) legA->Draw("SAME");

  TLegend *legB = new TLegend(0.48,0.16,0.78,0.30);
  legB->SetTextSize(0.040);
  legB->SetHeader("Other measurements");
  legB->AddEntry(lambdaMeas12,Lambda1name[2],"p");
  legB->AddEntry(lambdaMeas23,Lambda2name[3],"p"); 
  if (data) legB->Draw("SAME");

  TPaveText *tVud = new TPaveText(0.13,0.530,0.3,0.570,"nbNDC");
  tVud->SetBorderSize(0);
  tVud->SetTextColor(color?fillColor_vud:1);
  tVud->AddText("PDG 0^{+}#rightarrow0^{+}");
  tVud->GetLine(0)->SetTextSize(0.039);
  tVud->Draw();

  TPaveText *tTau1 = new TPaveText(0.46,0.84,0.61,0.88,"nbNDC");
  tTau1->SetBorderSize(0);
  tTau1->SetTextColor(color?fillColor_tau2:1);
  tTau1->AddText("UCN #tau_{n}");
  tTau1->GetLine(0)->SetTextAngle(-42.2);
  tTau1->GetLine(0)->SetTextSize(0.039);
  tTau1->Draw();

  TPaveText *tTau2 = new TPaveText(0.15,0.83,0.27,0.87,"nbNDC");
  tTau2->SetBorderSize(0);
  tTau2->SetTextColor(color?fillColor_tau2:1);
  tTau2->AddText("Beam #tau_{n}");
  tTau2->GetLine(0)->SetTextAngle(-42.2);
  tTau2->GetLine(0)->SetTextSize(0.039);
  tTau2->Draw();

  TPaveText *tLambda1 = new TPaveText(0.7,0.58,0.90,0.69,"nbNDC");
  tLambda1->SetBorderSize(0);
  tLambda1->SetTextColor(color?fillColor_lambda2:1);
  tLambda1->AddText("Post-2002 #lambda");
  tLambda1->GetLine(0)->SetTextAngle(-90);
  tLambda1->GetLine(0)->SetTextSize(0.039);
  tLambda1->Draw();

  TPaveText *tLambda2 = new TPaveText(0.234,0.57,0.434,0.67,"nbNDC");
  tLambda2->SetBorderSize(0);
  tLambda2->SetTextColor(color?fillColor_lambda2:1);
  tLambda2->AddText("Pre-2002 #lambda");
  tLambda2->GetLine(0)->SetTextAngle(-90);
  tLambda2->GetLine(0)->SetTextSize(0.039);
  tLambda2->Draw();
  
  c1->Update();
  c1->Print(TString::Format("vud_vs_lambda%s%s.pdf",data?"":"_noData",color?"_color":""));

  std::vector <Double_t> scaleVec;
  std::vector <Double_t> scaleVecErr;
  
  //Lambda pre-2002
  std::cout << "********** Pre-2002 Lambda ****************\n";
  Double_t l_pre2002 = weightedAve(Lambda2vec,Lambda2vecErr);
  Double_t lerr_pre2002 = weightedAveErr(Lambda2vecErr);
  Double_t d0 = 3.*TMath::Sqrt(Lambda2vec.size())*lerr_pre2002;
  for (int i=0;i<Lambda2vec.size(); ++i) {
    if (Lambda2vecErr[i]<d0) {
      scaleVec.push_back(Lambda2vec[i]);
      scaleVecErr.push_back(Lambda2vecErr[i]);
    }
    else {
      std::cout << "Measurement " << i << " (" << Lambda2name[i] << ") not included in scale\n";
    }
  }
  std::cout << "Scale Factor lambda pre-2002: " << TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)) << std::endl;
  std::cout << "Lambda pre-2002 w/ scaled error: " << l_pre2002 << " +/- "
	    << lerr_pre2002*(TMath::Sqrt(chi2nu(scaleVec,scaleVecErr))>1.?TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)):1.) << "\n";
  std::cout << "Lambda pre-2002 w/o scaled error: " << l_pre2002 << " +/- "
	    << lerr_pre2002 << "\n\n";
  scaleVec.resize(0);
  scaleVecErr.resize(0);

  //Lambda post-2002
  std::cout << "********** Post-2002 Lambda ****************\n";
  Double_t l_post2002 = weightedAve(Lambda1vec,Lambda1vecErr);
  Double_t lerr_post2002 = weightedAveErr(Lambda1vecErr);
  d0 = 3.*TMath::Sqrt(Lambda1vec.size())*lerr_post2002;
  for (int i=0;i<Lambda1vec.size(); ++i) {
    if (Lambda1vecErr[i]<d0) {
      scaleVec.push_back(Lambda1vec[i]);
      scaleVecErr.push_back(Lambda1vecErr[i]);
    }
    else {
      std::cout << "Measurement " << i << " (" << Lambda1name[i] << ") not included in scale\n";
    }
  }
  std::cout << "Scale Factor lambda post-2002: " << TMath::Sqrt(chi2nu(scaleVec,scaleVecErr))  << std::endl;
  std::cout << "Lambda post-2002 w/ scaled error: " << l_post2002 << " +/- "
	    << lerr_post2002*(TMath::Sqrt(chi2nu(scaleVec,scaleVecErr))>1.?TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)):1.) << "\n";
  std::cout << "Lambda post-2002 w/o scaled error: " << l_post2002 << " +/- "
	    << lerr_post2002 << "\n\n";
  scaleVec.resize(0);
  scaleVecErr.resize(0);

  //beam lifetime
  std::cout << "********** Beam Neutron Lifetime ****************\n";
  Double_t t_beam = weightedAve(Tau2vec,Tau2vecErr);
  Double_t terr_beam = weightedAveErr(Tau2vecErr);
  d0 = 3.*TMath::Sqrt(Tau2vec.size())*terr_beam;
  for (int i=0;i<Tau2vec.size(); ++i) {
    if (Tau2vecErr[i]<d0) {
      scaleVec.push_back(Tau2vec[i]);
      scaleVecErr.push_back(Tau2vecErr[i]);
    }
    else {
      std::cout << "Measurement " << i << "  not included in scale\n";
    }
  }
  std::cout << "Scale Factor beam lifetime: " << TMath::Sqrt(chi2nu(scaleVec,scaleVecErr))  << std::endl;
  std::cout << "beam lifetime w/ scaled error: " << t_beam << " +/- "
	    << terr_beam*(TMath::Sqrt(chi2nu(scaleVec,scaleVecErr))>1.?TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)):1.) << "\n";
  std::cout << "beam lifetime w/o scaled error: " << t_beam << " +/- "
	    << terr_beam << "\n\n";
  scaleVec.resize(0);
  scaleVecErr.resize(0);


  //bottle lifetime
  std::cout << "********** Bottle Neutron Lifetime ****************\n";
  Double_t t_bottle = weightedAve(Tau1vec,Tau1vecErr);
  Double_t terr_bottle = weightedAveErr(Tau1vecErr);
  d0 = 3.*TMath::Sqrt(Tau1vec.size())*terr_bottle;
  for (int i=0;i<Tau1vec.size(); ++i) {
    if (Tau1vecErr[i]<d0) {
      scaleVec.push_back(Tau1vec[i]);
      scaleVecErr.push_back(Tau1vecErr[i]);
    }
    else {
      std::cout << "Measurement " << i << "  not included in scale\n";
    }
  }
  std::cout << "Chi2/(n-1) bottle lifetime: " << chi2nu(scaleVec,scaleVecErr)*(scaleVec.size()-1.) << "/("<<scaleVec.size()<<"-1)" << std::endl;
  std::cout << "Scale Factor bottle lifetime: " << TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)) << std::endl;
  std::cout << "bottle lifetime w/ scaled error: " << t_bottle << " +/- "
	    << terr_bottle*(TMath::Sqrt(chi2nu(scaleVec,scaleVecErr))>1.?TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)):1.) << "\n";
  std::cout << "bottle lifetime w/o scaled error: " << t_bottle << " +/- "
	    << terr_bottle << "\n\n";
  
  scaleVec.resize(0);
  scaleVecErr.resize(0);


  std::vector<Double_t> combo;
  std::vector<Double_t> comboerr;
  std::vector<TString> comboname;
  
  //Combined lifetime
  std::cout << "\n********** Combined Lifetime ****************\n";
  for (auto i:Tau2vec) combo.push_back(i);
  for (auto i:Tau2vecErr) comboerr.push_back(i);
  for (auto i:Tau2name) comboname.push_back(i);
  for (auto i:Tau1vecErr) comboerr.push_back(i);
  for (auto i:Tau1vec) combo.push_back(i);
  for (auto i:Tau1name) comboname.push_back(i);
  Double_t t_combo = weightedAve(combo,comboerr);
  Double_t terr_combo = weightedAveErr(comboerr);
  d0 = 3.*TMath::Sqrt(combo.size())*terr_combo;
  for (int i=0;i<combo.size(); ++i) {
    if (comboerr[i]<d0) {
      scaleVec.push_back(combo[i]);
      scaleVecErr.push_back(comboerr[i]);
    }
    else {
      std::cout << "Measurement " << i << " (" << comboname[i] << ") not included in scale factor\n";
    }
  }
  std::cout << "Scale Factor combined lifetime: " << TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)) << std::endl;
  std::cout << "combo lifetime w/ scaled error: " << t_combo << " +/- "
	    << terr_combo*(TMath::Sqrt(chi2nu(scaleVec,scaleVecErr))>1.?TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)):1.) << "\n";
  std::cout << "combo lifetime w/o scaled error: " << t_combo << " +/- "
	    << terr_combo << "\n\n";
  
  scaleVec.resize(0);
  scaleVecErr.resize(0);
  combo.resize(0);
  comboerr.resize(0);
  comboname.resize(0);


  //Combined Lambda
  std::cout << "\n\n********** Combined Lambda ****************\n";
  for (auto i:Lambda2vec) combo.push_back(i);
  for (auto i:Lambda2vecErr) comboerr.push_back(i);
  for (auto i:Lambda2name) comboname.push_back(i);
  for (auto i:Lambda1vecErr) comboerr.push_back(i);
  for (auto i:Lambda1vec) combo.push_back(i);
  for (auto i:Lambda1name) comboname.push_back(i);
  Double_t l_combo = weightedAve(combo,comboerr);
  Double_t lerr_combo = weightedAveErr(comboerr);
  d0 = 3.*TMath::Sqrt(combo.size())*lerr_combo;
  for (int i=0;i<combo.size(); ++i) {
    if (comboerr[i]<d0) {
      scaleVec.push_back(combo[i]);
      scaleVecErr.push_back(comboerr[i]);
    }
    else {
      std::cout << "Measurement " << i << " (" << comboname[i] << ") not included in scale factor\n";
    }
  }
  std::cout << "Scale Factor combined lambda: " << TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)) << std::endl ;
  std::cout << "lambda combo w/ scaled error: " << l_combo << " +/- "
	    << lerr_combo*(TMath::Sqrt(chi2nu(scaleVec,scaleVecErr))>1.?TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)):1.) << "\n";
  std::cout << "lambda combo w/o scaled error: " << l_combo << " +/- "
	    << lerr_combo << "\n\n";
  
  scaleVec.resize(0);
  scaleVecErr.resize(0);
  
  //Combined Lambda with an extra point included
  std::cout << "\n\n********** Combined Lambda w/ extra measurements ****************\n";

  std::vector <Double_t> newl = calcLambda(-0.1184,0.001*0.1184);

  Double_t newLambdaOffsetFactor = 0.;//2.*0.00025; // 1 sigma (0.00025 or 0.025%)

  //Lambda1 = Lambda2;
  
  for (int i=0; i<1; ++i) {
    combo.push_back((Lambda1+Lambda1*newLambdaOffsetFactor));
    comboerr.push_back( (Lambda1+Lambda1*newLambdaOffsetFactor)*newl[1]/newl[0]);
    comboname.push_back(TString::Format("extra meas %i",i));
    Lambda1vec.push_back( (Lambda1+Lambda1*newLambdaOffsetFactor) );
    Lambda1vecErr.push_back( (Lambda1+Lambda1*newLambdaOffsetFactor)*newl[1]/newl[0]);
  }

  Double_t l_comboNew = weightedAve(combo,comboerr);
  Double_t lerr_comboNew = weightedAveErr(comboerr);
  d0 = 3.*TMath::Sqrt(combo.size())*lerr_comboNew;
  for (int i=0;i<combo.size(); ++i) {
    if (comboerr[i]<d0) {
      scaleVec.push_back(combo[i]);
      scaleVecErr.push_back(comboerr[i]);
    }
    else {
      std::cout << "Measurement " << i << " (" << comboname[i] << ") not included in scale factor\n";
    }
  }
  std::cout << "Scale Factor combined lambda: " << TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)) << std::endl << std::endl;
  
  Double_t New_lambda_combo = weightedAve(combo,comboerr);//1.2755;
  Double_t New_lambda_combo_err = weightedAveErr(comboerr);//0.0005;
  Double_t New_lambda_post2002 = weightedAve(Lambda1vec,Lambda1vecErr);//1.2755;
  Double_t New_lambda_post2002_err = weightedAveErr(Lambda1vecErr);//0.0005;
 
  std::cout << "\n\n\n";

  std::cout << "Combo Lambda w/ scaled error: " << New_lambda_combo << " +/- "
	    << New_lambda_combo_err*(TMath::Sqrt(chi2nu(scaleVec,scaleVecErr))>1.?TMath::Sqrt(chi2nu(scaleVec,scaleVecErr)):1.) << "\n";
  std::cout << "Post 2002 Lambda w/ scaled error: " << New_lambda_post2002 << " +/- "
	    << New_lambda_post2002_err*(TMath::Sqrt(chi2nu(Lambda1vec,Lambda1vecErr))>1.?TMath::Sqrt(chi2nu(Lambda1vec,Lambda1vecErr)):1.) << "\n";
  
  exit(0);
  
 

}
