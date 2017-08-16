
#include <fstream>
#include <iostream>
#include <vector>
#include <TString.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TMath.h>

Int_t groupBin=1;
Double_t enStart=40.;//10.*groupBin;

Double_t fitRangeMin = 180.;
Double_t fitRangeMax = 700.;

Double_t BSlimitLow = -1., BSlimitHigh = 3.;
Double_t AnglelimitLow = -5., AnglelimitHigh = 3.;

TString drawOpt = "0AL3";//"0AC4" "0AL3"

std::vector <std::vector<Double_t> > readBSCorr(TString corr,TString year,TString anaCh) {

  std::vector <std::vector<Double_t> > corrs;
  //  std::ifstream infile(TString::Format("../OldMCCorrection/deltaBS%s_anaCh%s_%s.txt",corr.Data(),anaCh.Data(),year.Data()));
  std::ifstream infile(TString::Format("../OldMCCorrection/%s_delta_2%s_anaCh%s.txt",year.Data(),corr!=std::string("ALL")?corr.Data():"",anaCh.Data()));
  std::vector<Double_t> enbin;
  std::vector<Double_t> c;
  std::vector<Double_t> cerr;
  Double_t en;
  TString hold1,hold2;
  Int_t i=0;
  Double_t aveEn=0., aveCorr=0.,aveErr=0.;
  while (infile >> en >> hold1 >> hold2 && en<800) {
    if (en<enStart) continue;
    std::cout << en << " " << hold1 << " " << hold2 << "\n";
    if (i<groupBin) {
      aveEn+=en;
      aveCorr+=( hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold1.Data());
      aveErr+=(hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold2.Data());
      i++;
    } else {
      enbin.push_back(aveEn/groupBin);
      c.push_back(aveCorr/groupBin*100.);
      cerr.push_back(aveErr/groupBin*100.);
      std::cout << aveErr/groupBin*100. << std::endl;
      i=0,aveEn=0.,aveCorr=0.,aveErr=0.;
      aveEn+=en;
      aveCorr+=( hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold1.Data());
      aveErr+=(hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold2.Data());
      i++;
    }
  }
  enbin.push_back(enbin[enbin.size()-1]+10.*groupBin);
  c.push_back(c[c.size()-1]);
  cerr.push_back(cerr[cerr.size()-1]);
  corrs.push_back(enbin);
  corrs.push_back(c);
  corrs.push_back(cerr);
  return corrs;
};

std::vector <std::vector<Double_t> > readAngleCorr(TString year, TString type) {

  std::vector <std::vector<Double_t> > corrs;
  std::ifstream infile(TString::Format("../AngleCorrections/%s_delta_3%s.txt",year.Data(),(type==TString("ALL")?"":type.Data())));
  std::vector<Double_t> enbin;
  std::vector<Double_t> c;
  std::vector<Double_t> cerr;
  TString hold1,hold2;
  Double_t en;
  Int_t i=0;
  Double_t aveEn=0., aveCorr=0., aveErr=0.;
  while (infile >> en >> hold1 >> hold2 && en<800.) {
    if (en<enStart) continue;
    //std::cout << en << " " << hold2 << " " << errhold << "\n";
    if (i<groupBin) {
      aveEn+=en;
      aveCorr+=( hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold1.Data());
      aveErr+=(hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold2.Data());
      i++;
    } else {
      enbin.push_back(aveEn/groupBin);
      c.push_back(aveCorr/groupBin*100.);
      cerr.push_back(aveErr/groupBin*100.);
      i=0,aveEn=0.,aveCorr=0., aveErr=0.;
      aveEn+=en;
      aveCorr+=( hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold1.Data());
      aveErr+=(hold1==TString("nan")||hold1==TString("-nan"))?0.:atof(hold2.Data());
      i++;
    }
  }
  corrs.push_back(enbin);
  corrs.push_back(c);
  corrs.push_back(cerr);
  return corrs;
};


void fitCorrections(TString anaCh) {

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.07,"t");
  gStyle->SetPadLeftMargin(0.12);
  //gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(1.);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleSize(0.04,"Z");
  gStyle->SetTitleOffset(1.6,"Z");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetFillStyle(0);
  gStyle->SetGridStyle(2);
  //  gStyle->SetHatchesSpacing(0.6);
  //gStyle->SetHatchesLineWidth(2);
  // read in the corrections and plot their deltaA/A for each bin
  // This can be done for both theory and MC

  //Int_t octmin = year==2011?0:60;
  //Int_t octmax = year==2011?59:121;

  int fill2011 = 3001;
  int fill2012 = 3001;

  int col2011 = 4;
  int col2012 = 3;

  int startPoint=0;
  
  TString year = "2011-2012";
  std::vector<std::vector<Double_t> > bs0_2011 = readBSCorr("0",year,anaCh);
  std::vector<std::vector<Double_t> > bs1_2011 = readBSCorr("1",year,anaCh);
  std::vector<std::vector<Double_t> > bs2_2011 = readBSCorr("2",year,anaCh);
  std::vector<std::vector<Double_t> > bs3_2011 = readBSCorr("3",year,anaCh);
  std::vector<std::vector<Double_t> > bsALL_2011 = readBSCorr("ALL",year,anaCh);

  std::vector<std::vector<Double_t> > delta30_2011 = readAngleCorr(year,"0");
  std::vector<std::vector<Double_t> > delta31_2011 = readAngleCorr(year,"1");
  std::vector<std::vector<Double_t> > delta32_2011 = readAngleCorr(year,"2");
  std::vector<std::vector<Double_t> > delta33_2011 = readAngleCorr(year,"3");
  std::vector<std::vector<Double_t> > delta3_2011 = readAngleCorr(year,"ALL");

  year = "2012-2013";
  std::vector<std::vector<Double_t> > bs0_2012 = readBSCorr("0",year,anaCh);
  std::vector<std::vector<Double_t> > bs1_2012 = readBSCorr("1",year,anaCh);
  std::vector<std::vector<Double_t> > bs2_2012 = readBSCorr("2",year,anaCh);
  std::vector<std::vector<Double_t> > bs3_2012 = readBSCorr("3",year,anaCh);
  std::vector<std::vector<Double_t> > bsALL_2012 = readBSCorr("ALL",year,anaCh);

  std::vector<std::vector<Double_t> > delta30_2012 = readAngleCorr(year,"0");
  std::vector<std::vector<Double_t> > delta31_2012 = readAngleCorr(year,"1");
  std::vector<std::vector<Double_t> > delta32_2012 = readAngleCorr(year,"2");
  std::vector<std::vector<Double_t> > delta33_2012 = readAngleCorr(year,"3");
  std::vector<std::vector<Double_t> > delta3_2012 = readAngleCorr(year,"ALL");
 
  TCanvas *c0 = new TCanvas("c0","c0");
  gPad->SetGrid(0,1);

  TMultiGraph *mg0 = new TMultiGraph();
  mg0->SetTitle("#Delta_{2,0} vs. Energy");  
  
  TGraphErrors *g_bs0_2011 = new TGraphErrors(bs0_2011[0].size()-startPoint,&bs0_2011[0][startPoint],&bs0_2011[1][startPoint],0,0);//&bs0_2011[2][startPoint]);
  g_bs0_2011->SetMarkerStyle(0);
  g_bs0_2011->SetLineWidth(3);
  g_bs0_2011->SetLineColor(col2011);
  g_bs0_2011->SetFillColor(col2011);
  g_bs0_2011->SetFillStyle(fill2011);

  TGraphErrors *g_bs0_2012 = new TGraphErrors(bs0_2012[0].size()-startPoint,&bs0_2012[0][startPoint],&bs0_2012[1][startPoint],0,0);//&bs0_2012[2][startPoint]);
  g_bs0_2012->SetMarkerStyle(0);
  g_bs0_2012->SetLineWidth(3);
  g_bs0_2012->SetLineColor(col2012);
  g_bs0_2012->SetFillColor(col2012);
  g_bs0_2012->SetFillStyle(fill2012);

  TF1 *f0 = new TF1("f0","[0]*TMath::Exp(-[1]*x)+pol3(2)",fitRangeMin,fitRangeMax);//+[2]+[3]*x+[4]*x*x",100.,700.);
  f0->SetParameters(1.,0.,1.,0.,0.,0.);
  g_bs0_2011->Fit(f0,"RQ");
  g_bs0_2011->Draw("AL");

  //calculate the err from fit 
  Double_t sigma = 0.;
  Double_t num=0.;
  for ( UInt_t i=0;i<bs0_2011[1].size();++i) {
    if (bs0_2011[0][i]>=fitRangeMin && bs0_2011[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(bs0_2011[1][i]-f0->Eval(bs0_2011[0][i]),2);
      //std::cout << bs0_2011[1][i] << " - " << f0->Eval(bs0_2011[0][i]) << " = " << bs0_2011[1][i]-f0->Eval(bs0_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA20_2011 = " << sigma << std::endl;

  g_bs0_2012->Fit(f0,"RQ");
  g_bs0_2012->Draw("AL");

  //calculate the err from fit 
  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<bs0_2012[1].size();++i) {
    if (bs0_2012[0][i]>=fitRangeMin && bs0_2012[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(bs0_2012[1][i]-f0->Eval(bs0_2012[0][i]),2);
      //std::cout << bs0_2011[1][i] << " - " << f0->Eval(bs0_2011[0][i]) << " = " << bs0_2011[1][i]-f0->Eval(bs0_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA20_2012 = " << sigma << std::endl;
  
  c0->cd();
  
  
  mg0->Add(g_bs0_2011);
  mg0->Add(g_bs0_2012);
  mg0->SetMinimum(BSlimitLow);
  mg0->SetMaximum(BSlimitHigh);
  
  mg0->Draw(drawOpt);
  mg0->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg0->GetYaxis()->CenterTitle();
  mg0->GetXaxis()->SetTitle("Energy (keV)");
  mg0->GetXaxis()->CenterTitle();
  mg0->GetXaxis()->SetLimits(0.,780.);
  gPad->Modified();

  TLegend *leg0 = new TLegend(0.57,0.7,0.87,0.8);
  leg0->AddEntry(g_bs0_2011,"2011-2012","lf");
  leg0->AddEntry(g_bs0_2012,"2012-2013","lf");
  leg0->SetTextSize(0.05);
  leg0->Draw("SAME");


  //////////////////// BS1///////////////////
  TCanvas *c1 = new TCanvas("c1","c1");
  gPad->SetGrid(0,1);

  TMultiGraph *mg1 = new TMultiGraph();
  mg1->SetTitle("#Delta_{2,1} vs. Energy");  
  
  TGraphErrors *g_bs1_2011 = new TGraphErrors(bs1_2011[0].size()-startPoint,&bs1_2011[0][startPoint],&bs1_2011[1][startPoint],0,0);//&bs1_2011[2][startPoint]);
  g_bs1_2011->SetMarkerStyle(0);
  g_bs1_2011->SetLineWidth(3);
  g_bs1_2011->SetLineColor(col2011);
  g_bs1_2011->SetFillColor(col2011);
  g_bs1_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_bs1_2012 = new TGraphErrors(bs1_2012[0].size()-startPoint,&bs1_2012[0][startPoint],&bs1_2012[1][startPoint],0,0);//&bs1_2012[2][startPoint]);
  g_bs1_2012->SetMarkerStyle(0);
  g_bs1_2012->SetLineWidth(3);
  g_bs1_2012->SetLineColor(col2012);
  g_bs1_2012->SetFillColor(col2012);
  g_bs1_2012->SetFillStyle(fill2012);

  g_bs1_2011->Fit(f0,"RQ");
  g_bs1_2011->Draw("AL");

  //calculate the err from fit 
  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<bs1_2011[1].size();++i) {
    if (bs1_2011[0][i]>=fitRangeMin && bs1_2011[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(bs1_2011[1][i]-f0->Eval(bs1_2011[0][i]),2);
      //std::cout << bs1_2011[1][i] << " - " << f0->Eval(bs1_2011[0][i]) << " = " << bs1_2011[1][i]-f0->Eval(bs1_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA21_2011 = " << sigma << std::endl;

  g_bs1_2012->Fit(f0,"RQ");
  g_bs1_2012->Draw("AL");

  //calculate the err from fit 
  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<bs1_2012[1].size();++i) {
    if (bs1_2012[0][i]>=fitRangeMin && bs1_2012[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(bs1_2012[1][i]-f0->Eval(bs1_2012[0][i]),2);
      //std::cout << bs1_2011[1][i] << " - " << f0->Eval(bs1_2011[0][i]) << " = " << bs1_2011[1][i]-f0->Eval(bs1_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA21_2012 = " << sigma << std::endl;
  
  c1->cd();  
  mg1->Add(g_bs1_2011);
  mg1->Add(g_bs1_2012);
  mg1->SetMinimum(BSlimitLow);
  mg1->SetMaximum(BSlimitHigh);
  
  mg1->Draw(drawOpt);
  mg1->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg1->GetYaxis()->CenterTitle();
  mg1->GetXaxis()->SetTitle("Energy (keV)");
  mg1->GetXaxis()->CenterTitle();
  mg1->GetXaxis()->SetLimits(0.,780.);

  gPad->Modified();

  TLegend *leg1 = new TLegend(0.57,0.7,0.87,0.8);
  leg1->AddEntry(g_bs1_2011,"2011-2012","lf");
  leg1->AddEntry(g_bs1_2012,"2012-2013","lf");
  leg1->SetTextSize(0.05);
  leg1->Draw("SAME");

  //////////////////// BS2///////////////////
  TCanvas *c2 = new TCanvas("c2","c2");
  gPad->SetGrid(0,1);

  TMultiGraph *mg2 = new TMultiGraph();
  mg2->SetTitle("#Delta_{2,2} vs. Energy");  
  
  TGraphErrors *g_bs2_2011 = new TGraphErrors(bs2_2011[0].size()-startPoint,&bs2_2011[0][startPoint],&bs2_2011[1][startPoint],0,0);//&bs2_2011[2][startPoint]);
  g_bs2_2011->SetMarkerStyle(0);
  g_bs2_2011->SetLineWidth(3);
  g_bs2_2011->SetLineColor(col2011);
  g_bs2_2011->SetFillColor(col2011);
  g_bs2_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_bs2_2012 = new TGraphErrors(bs2_2012[0].size()-startPoint,&bs2_2012[0][startPoint],&bs2_2012[1][startPoint],0,0);//&bs2_2012[2][startPoint]);
  g_bs2_2012->SetMarkerStyle(0);
  g_bs2_2012->SetLineWidth(3);
  g_bs2_2012->SetLineColor(col2012);
  g_bs2_2012->SetFillColor(col2012);
  g_bs2_2012->SetFillStyle(fill2012);


  //calculate the err from fit 
  g_bs2_2011->Fit(f0,"RQ");
  g_bs2_2011->Draw("AL");

  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<bs2_2011[1].size();++i) {
    if (bs2_2011[0][i]>=fitRangeMin && bs2_2011[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(bs2_2011[1][i]-f0->Eval(bs2_2011[0][i]),2);
      //std::cout << bs2_2011[1][i] << " - " << f0->Eval(bs2_2011[0][i]) << " = " << bs2_2011[1][i]-f0->Eval(bs2_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA22_2011 = " << sigma << std::endl;

  g_bs2_2012->Fit(f0,"RQ");
  g_bs2_2012->Draw("AL");

  //calculate the err from fit 
  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<bs2_2012[1].size();++i) {
    if (bs2_2012[0][i]>=fitRangeMin && bs2_2012[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(bs2_2012[1][i]-f0->Eval(bs2_2012[0][i]),2);
      //std::cout << bs2_2011[1][i] << " - " << f0->Eval(bs2_2011[0][i]) << " = " << bs2_2011[1][i]-f0->Eval(bs2_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA22_2012 = " << sigma << std::endl;
  

  c2->cd();
  
  mg2->Add(g_bs2_2011);
  mg2->Add(g_bs2_2012);
  mg2->SetMinimum(BSlimitLow);
  mg2->SetMaximum(BSlimitHigh);
  
  mg2->Draw(drawOpt);
  mg2->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg2->GetYaxis()->CenterTitle();
  mg2->GetXaxis()->SetTitle("Energy (keV)");
  mg2->GetXaxis()->CenterTitle();
  mg2->GetXaxis()->SetLimits(0.,780.);
    
  gPad->Modified();

  TLegend *leg2 = new TLegend(0.57,0.7,0.87,0.8);
  leg2->AddEntry(g_bs2_2011,"2011-2012","lf");
  leg2->AddEntry(g_bs2_2012,"2012-2013","lf");
  leg2->SetTextSize(0.05);
  leg2->Draw("SAME");

  //////////////////// BS3///////////////////
  TCanvas *c3 = new TCanvas("c3","c3");
  gPad->SetGrid(0,1);

  TMultiGraph *mg3 = new TMultiGraph();
  mg3->SetTitle("#Delta_{2,3} vs. Energy");  
  
  TGraphErrors *g_bs3_2011 = new TGraphErrors(bs3_2011[0].size()-startPoint,&bs3_2011[0][startPoint],&bs3_2011[1][startPoint],0,0);//&bs3_2011[2][startPoint]);
  g_bs3_2011->SetMarkerStyle(0);
  g_bs3_2011->SetLineWidth(3);
  g_bs3_2011->SetLineColor(col2011);
  g_bs3_2011->SetFillColor(col2011);
  g_bs3_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_bs3_2012 = new TGraphErrors(bs3_2012[0].size()-startPoint,&bs3_2012[0][startPoint],&bs3_2012[1][startPoint],0,0);//&bs3_2012[2][startPoint]);
  g_bs3_2012->SetMarkerStyle(0);
  g_bs3_2012->SetLineWidth(3);
  g_bs3_2012->SetLineColor(col2012);
  g_bs3_2012->SetFillColor(col2012);
  g_bs3_2012->SetFillStyle(fill2012);

  //calculate the err from fit 
  g_bs3_2011->Fit(f0,"RQ");
  g_bs3_2011->Draw("AL");

  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<bs3_2011[1].size();++i) {
    if (bs3_2011[0][i]>=fitRangeMin && bs3_2011[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(bs3_2011[1][i]-f0->Eval(bs3_2011[0][i]),2);
      //std::cout << bs3_2011[1][i] << " - " << f0->Eval(bs3_2011[0][i]) << " = " << bs3_2011[1][i]-f0->Eval(bs3_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA23_2011 = " << sigma << std::endl;

  g_bs3_2012->Fit(f0,"RQ");
  g_bs3_2012->Draw("AL");

  //calculate the err from fit 
  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<bs3_2012[1].size();++i) {
    if (bs3_2012[0][i]>=fitRangeMin && bs3_2012[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(bs3_2012[1][i]-f0->Eval(bs3_2012[0][i]),2);
      //std::cout << bs3_2011[1][i] << " - " << f0->Eval(bs3_2011[0][i]) << " = " << bs3_2011[1][i]-f0->Eval(bs3_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA23_2012 = " << sigma << std::endl;
  

  c3->cd();
  
  
  mg3->Add(g_bs3_2011);
  mg3->Add(g_bs3_2012);
  mg3->SetMinimum(BSlimitLow);
  mg3->SetMaximum(BSlimitHigh);
  
  mg3->Draw(drawOpt);
  mg3->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg3->GetYaxis()->CenterTitle();
  mg3->GetXaxis()->SetTitle("Energy (keV)");
  mg3->GetXaxis()->CenterTitle();
  mg3->GetXaxis()->SetLimits(0.,780.);

  gPad->Modified();

  TLegend *leg3 = new TLegend(0.57,0.7,0.87,0.8);
  leg3->AddEntry(g_bs3_2011,"2011-2012","lf");
  leg3->AddEntry(g_bs3_2012,"2012-2013","lf");
  leg3->SetTextSize(0.05);
  leg3->Draw("SAME");

  /////////////////// Angle corrections ////////////////////


  TCanvas *c30 = new TCanvas("c30","c30");
  gPad->SetGrid(0,1);

  TMultiGraph *mg30 = new TMultiGraph();
  mg30->SetTitle("#Delta_{3,0} vs. Energy");  
  
  TGraphErrors *g_delta30_2011 = new TGraphErrors(delta30_2011[0].size()-startPoint,&delta30_2011[0][startPoint],&delta30_2011[1][startPoint],0,0);//&delta30_2011[2][startPoint]);
  g_delta30_2011->SetMarkerStyle(0);
  g_delta30_2011->SetLineWidth(3);
  g_delta30_2011->SetLineColor(col2011);
  g_delta30_2011->SetFillColor(col2011);
  g_delta30_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_delta30_2012 = new TGraphErrors(delta30_2012[0].size()-startPoint,&delta30_2012[0][startPoint],&delta30_2012[1][startPoint],0,0);//&delta30_2012[2][startPoint]);
  g_delta30_2012->SetMarkerStyle(0);
  g_delta30_2012->SetLineWidth(3);
  g_delta30_2012->SetLineColor(col2012);
  g_delta30_2012->SetFillColor(col2012);
  g_delta30_2012->SetFillStyle(fill2012);

    //calculate the err from fit 
  g_delta30_2011->Fit(f0,"RQ");
  g_delta30_2011->Draw("AL");
  
  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<delta30_2011[1].size();++i) {
    if (delta30_2011[0][i]>=fitRangeMin && delta30_2011[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(delta30_2011[1][i]-f0->Eval(delta30_2011[0][i]),2);
      //std::cout << delta30_2011[1][i] << " - " << f0->Eval(delta30_2011[0][i]) << " = " << delta30_2011[1][i]-f0->Eval(delta30_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA30_2011 = " << sigma << std::endl;

  g_delta30_2012->Fit(f0,"RQ");
  g_delta30_2012->Draw("AL");

  //calculate the err from fit 
  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<delta30_2012[1].size();++i) {
    if (delta30_2012[0][i]>=fitRangeMin && delta30_2012[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(delta30_2012[1][i]-f0->Eval(delta30_2012[0][i]),2);
      //std::cout << delta30_2011[1][i] << " - " << f0->Eval(delta30_2011[0][i]) << " = " << delta30_2011[1][i]-f0->Eval(delta30_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA30_2012 = " << sigma << std::endl;

  
  mg30->Add(g_delta30_2011);
  mg30->Add(g_delta30_2012);
  mg30->SetMinimum(AnglelimitLow);
  mg30->SetMaximum(AnglelimitHigh);
  
  mg30->Draw("ALP3");
  mg30->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg30->GetYaxis()->CenterTitle();
  mg30->GetXaxis()->SetTitle("Energy (keV)");
  mg30->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg30 = new TLegend(0.57,0.7,0.87,0.8);
  leg30->AddEntry(g_delta30_2011,"2011-2012","lf");
  leg30->AddEntry(g_delta30_2012,"2012-2013","lf");
  leg30->SetTextSize(0.05);
  leg30->Draw("SAME");


  //////////////////// DELTA31///////////////////
  TCanvas *c31 = new TCanvas("c31","c31");
  gPad->SetGrid(0,1);

  TMultiGraph *mg31 = new TMultiGraph();
  mg31->SetTitle("#Delta_{3,1} vs. Energy");  
  
  TGraphErrors *g_delta31_2011 = new TGraphErrors(delta31_2011[0].size()-startPoint,&delta31_2011[0][startPoint],&delta31_2011[1][startPoint],0,0);//&delta31_2011[2][startPoint]);
  g_delta31_2011->SetMarkerStyle(0);
  g_delta31_2011->SetLineWidth(3);
  g_delta31_2011->SetLineColor(col2011);
  g_delta31_2011->SetFillColor(col2011);
  g_delta31_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_delta31_2012 = new TGraphErrors(delta31_2012[0].size()-startPoint,&delta31_2012[0][startPoint],&delta31_2012[1][startPoint],0,0);//&delta31_2012[2][startPoint]);
  g_delta31_2012->SetMarkerStyle(0);
  g_delta31_2012->SetLineWidth(3);
  g_delta31_2012->SetLineColor(col2012);
  g_delta31_2012->SetFillColor(col2012);
  g_delta31_2012->SetFillStyle(fill2012);

      //calculate the err from fit 
  g_delta31_2011->Fit(f0,"RQ");
  g_delta31_2011->Draw("AL");

  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<delta31_2011[1].size();++i) {
    if (delta31_2011[0][i]>=fitRangeMin && delta31_2011[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(delta31_2011[1][i]-f0->Eval(delta31_2011[0][i]),2);
      //std::cout << delta31_2011[1][i] << " - " << f0->Eval(delta31_2011[0][i]) << " = " << delta31_2011[1][i]-f0->Eval(delta31_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA31_2011 = " << sigma << std::endl;

  g_delta31_2012->Fit(f0,"RQ");
  g_delta31_2012->Draw("AL");

  //calculate the err from fit 
  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<delta31_2012[1].size();++i) {
    if (delta31_2012[0][i]>=fitRangeMin && delta31_2012[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(delta31_2012[1][i]-f0->Eval(delta31_2012[0][i]),2);
      //std::cout << delta31_2011[1][i] << " - " << f0->Eval(delta31_2011[0][i]) << " = " << delta31_2011[1][i]-f0->Eval(delta31_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA31_2012 = " << sigma << std::endl;

  mg31->Add(g_delta31_2011);
  mg31->Add(g_delta31_2012);
  mg31->SetMinimum(AnglelimitLow);
  mg31->SetMaximum(AnglelimitHigh);
  
  mg31->Draw("ALP3");
  mg31->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg31->GetYaxis()->CenterTitle();
  mg31->GetXaxis()->SetTitle("Energy (keV)");
  mg31->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg31 = new TLegend(0.57,0.7,0.87,0.8);
  leg31->AddEntry(g_delta31_2011,"2011-2012","lf");
  leg31->AddEntry(g_delta31_2012,"2012-2013","lf");
  leg31->SetTextSize(0.05);
  leg31->Draw("SAME");

  //////////////////// DELTA32///////////////////
  TCanvas *c32 = new TCanvas("c32","c32");
  gPad->SetGrid(0,1);

  TMultiGraph *mg32 = new TMultiGraph();
  mg32->SetTitle("#Delta_{3,2} vs. Energy");  
  
  TGraphErrors *g_delta32_2011 = new TGraphErrors(delta32_2011[0].size()-startPoint,&delta32_2011[0][startPoint],&delta32_2011[1][startPoint],0,0);//&delta32_2011[2][startPoint]);
  g_delta32_2011->SetMarkerStyle(0);
  g_delta32_2011->SetLineWidth(3);
  g_delta32_2011->SetLineColor(col2011);
  g_delta32_2011->SetFillColor(col2011);
  g_delta32_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_delta32_2012 = new TGraphErrors(delta32_2012[0].size()-startPoint,&delta32_2012[0][startPoint],&delta32_2012[1][startPoint],0,0);//&delta32_2012[2][startPoint]);
  g_delta32_2012->SetMarkerStyle(0);
  g_delta32_2012->SetLineWidth(3);
  g_delta32_2012->SetLineColor(col2012);
  g_delta32_2012->SetFillColor(col2012);
  g_delta32_2012->SetFillStyle(fill2012);

      //calculate the err from fit 
  g_delta32_2011->Fit(f0,"RQ");
  g_delta32_2011->Draw("AL");

  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<delta32_2011[1].size();++i) {
    if (delta32_2011[0][i]>=fitRangeMin && delta32_2011[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(delta32_2011[1][i]-f0->Eval(delta32_2011[0][i]),2);
      //std::cout << delta32_2011[1][i] << " - " << f0->Eval(delta32_2011[0][i]) << " = " << delta32_2011[1][i]-f0->Eval(delta32_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA32_2011 = " << sigma << std::endl;

  g_delta32_2012->Fit(f0,"RQ");
  g_delta32_2012->Draw("AL");

  //calculate the err from fit 
  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<delta32_2012[1].size();++i) {
    if (delta32_2012[0][i]>=fitRangeMin && delta32_2012[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(delta32_2012[1][i]-f0->Eval(delta32_2012[0][i]),2);
      //std::cout << delta32_2011[1][i] << " - " << f0->Eval(delta32_2011[0][i]) << " = " << delta32_2011[1][i]-f0->Eval(delta32_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA32_2012 = " << sigma << std::endl;

  
  mg32->Add(g_delta32_2011);
  mg32->Add(g_delta32_2012);
  mg32->SetMinimum(AnglelimitLow);
  mg32->SetMaximum(AnglelimitHigh);
  
  mg32->Draw("ALP3");
  mg32->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg32->GetYaxis()->CenterTitle();
  mg32->GetXaxis()->SetTitle("Energy (keV)");
  mg32->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg32 = new TLegend(0.57,0.7,0.87,0.8);
  leg32->AddEntry(g_delta32_2011,"2011-2012","lf");
  leg32->AddEntry(g_delta32_2012,"2012-2013","lf");
  leg32->SetTextSize(0.05);
  leg32->Draw("SAME");

  //////////////////// DELTA33///////////////////
  TCanvas *c33 = new TCanvas("c33","c33");
  gPad->SetGrid(0,1);

  TMultiGraph *mg33 = new TMultiGraph();
  mg33->SetTitle("#Delta_{3,3} vs. Energy");  
  
  TGraphErrors *g_delta33_2011 = new TGraphErrors(delta33_2011[0].size()-startPoint,&delta33_2011[0][startPoint],&delta33_2011[1][startPoint],0,0);//&delta33_2011[2][startPoint]);
  g_delta33_2011->SetMarkerStyle(0);
  g_delta33_2011->SetLineWidth(3);
  g_delta33_2011->SetLineColor(col2011);
  g_delta33_2011->SetFillColor(col2011);
  g_delta33_2011->SetFillStyle(fill2011);
  
  TGraphErrors *g_delta33_2012 = new TGraphErrors(delta33_2012[0].size()-startPoint,&delta33_2012[0][startPoint],&delta33_2012[1][startPoint],0,0);//&delta33_2012[2][startPoint]);
  g_delta33_2012->SetMarkerStyle(0);
  g_delta33_2012->SetLineWidth(3);
  g_delta33_2012->SetLineColor(col2012);
  g_delta33_2012->SetFillColor(col2012);
  g_delta33_2012->SetFillStyle(fill2012);

      //calculate the err from fit 
  g_delta33_2011->Fit(f0,"RQ");
  g_delta33_2011->Draw("AL");

  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<delta33_2011[1].size();++i) {
    if (delta33_2011[0][i]>=fitRangeMin && delta33_2011[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(delta33_2011[1][i]-f0->Eval(delta33_2011[0][i]),2);
      //std::cout << delta33_2011[1][i] << " - " << f0->Eval(delta33_2011[0][i]) << " = " << delta33_2011[1][i]-f0->Eval(delta33_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA33_2011 = " << sigma << std::endl;

  g_delta33_2012->Fit(f0,"RQ");
  g_delta33_2012->Draw("AL");

  //calculate the err from fit 
  sigma = 0.;
  num=0.;
  for ( UInt_t i=0;i<delta33_2012[1].size();++i) {
    if (delta33_2012[0][i]>=fitRangeMin && delta33_2012[0][i]<=fitRangeMax) {
      sigma+=TMath::Power(delta33_2012[1][i]-f0->Eval(delta33_2012[0][i]),2);
      //std::cout << delta33_2011[1][i] << " - " << f0->Eval(delta33_2011[0][i]) << " = " << delta33_2011[1][i]-f0->Eval(delta33_2011[0][i]) << "\n";
      num+=1.;
    }
  }
  sigma/=(num);
  sigma = TMath::Sqrt(sigma);
  std::cout << "SIGMA33_2012 = " << sigma << std::endl;

  
  mg33->Add(g_delta33_2011);
  mg33->Add(g_delta33_2012);
  mg33->SetMinimum(AnglelimitLow);
  mg33->SetMaximum(AnglelimitHigh);
  
  mg33->Draw("ALP3");
  mg33->GetYaxis()->SetTitle("#DeltaA/A (%)");
  mg33->GetYaxis()->CenterTitle();
  mg33->GetXaxis()->SetTitle("Energy (keV)");
  mg33->GetXaxis()->CenterTitle();
  gPad->Modified();

  TLegend *leg33 = new TLegend(0.57,0.7,0.87,0.8);
  leg33->AddEntry(g_delta33_2011,"2011-2012","lf");
  leg33->AddEntry(g_delta33_2012,"2012-2013","lf");
  leg33->SetTextSize(0.05);
  leg33->Draw("SAME");
  

  TString pdffile = TString::Format("MC_Fits_anaCh%s.pdf",anaCh.Data());

  c0->Print(TString::Format("%s(",pdffile.Data()));
  c1->Print(pdffile);
  c2->Print(pdffile);
  c3->Print(pdffile);
  c30->Print(pdffile);
  c31->Print(pdffile);
  c32->Print(pdffile);
  c33->Print(TString::Format("%s)",pdffile.Data()));
  
}
