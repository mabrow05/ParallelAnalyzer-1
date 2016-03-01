/*
Compilable code which houses the heart of the analysis. 
Here are classes which handle calculating asymmetries of either data 
or simulation. There will be flags to be passed dictating what type of 
asymmetry, what type of data, and what octet/runs to use. 

The code should be able to do BG subtraction, construct asymmetries, 
or do both depending on the flags passed.

Maybe even add in writing the final answer to the database if the user wants to

*/

#include "Asymmetries.hh"
#include "MBUtils.hh"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>

#include <TCanvas.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TPDF.h>
#include <TLegend.h>
#include <TPaveText.h>

int main(int argc, char *argv[])
{
  if (argc!=5) {
    std::cout << "Usage: ./UK2MPMcomparison [octet] [energy window low] [energy window high] [energy bin width]\n";
  }

  gStyle->SetOptStat(0);

  int octetNum = atoi(argv[1]);
  double enWinLow = atof(argv[2]);
  double enWinHigh = atof(argv[3]);
  int enBinWidth = atoi(argv[4]);

  std::string pdfFileBase = "calibrationRunComp_Octet"+itos(octetNum)+"_evis_"+itos((int)enWinLow)+"-"+itos((int)enWinHigh)+".pdf";

  int numBins = 800/enBinWidth;
  //Making histograms for comparison  
  TH1F *EvisEALL_mpm = new TH1F("EvisEALL_mpm","EvisE Type ALL", numBins, 0., 800.);
  TH1F *EvisWALL_mpm = new TH1F("EvisWALL_mpm","EvisW Type ALL", numBins, 0., 800.);
  TH1F *EvisEALL_uk = new TH1F("EvisEALL_uk","EvisE Type ALL", numBins, 0., 800.);
  EvisEALL_uk->GetXaxis()->SetTitle("Energy (keV)");
  TH1F *EvisWALL_uk = new TH1F("EvisWALL_uk","EvisW Type ALL", numBins, 0., 800.);
  EvisWALL_uk->GetXaxis()->SetTitle("Energy (keV)");

  TH1F *EvisE0_mpm = new TH1F("EvisE0_mpm","EvisE Type 0", numBins, 0., 800.);
  TH1F *EvisW0_mpm = new TH1F("EvisW0_mpm","EvisW Type 0", numBins, 0., 800.);
  TH1F *EvisE0_uk = new TH1F("EvisE0_uk","EvisE Type 0", numBins, 0., 800.);
  EvisE0_uk->GetXaxis()->SetTitle("Energy (keV)");
  TH1F *EvisW0_uk = new TH1F("EvisW0_uk","EvisW Type 0", numBins, 0., 800.);
  EvisW0_uk->GetXaxis()->SetTitle("Energy (keV)");

  TH1F *EvisE1_mpm = new TH1F("EvisE1_mpm","EvisE Type 1", numBins, 0., 800.);
  TH1F *EvisW1_mpm = new TH1F("EvisW1_mpm","EvisW Type 1", numBins, 0., 800.);
  TH1F *EvisE1_uk = new TH1F("EvisE1_uk","EvisE Type 1", numBins, 0., 800.);
  EvisE1_uk->GetXaxis()->SetTitle("Energy (keV)");
  TH1F *EvisW1_uk = new TH1F("EvisW1_uk","EvisW Type 1", numBins, 0., 800.);
  EvisW1_uk->GetXaxis()->SetTitle("Energy (keV)");

  TH1F *EvisE23_mpm = new TH1F("EvisE23_mpm","EvisE Type 2/3", numBins, 0., 800.);
  TH1F *EvisW23_mpm = new TH1F("EvisW23_mpm","EvisW Type 2/3", numBins, 0., 800.);
  TH1F *EvisE23_uk = new TH1F("EvisE23_uk","EvisE Type 2/3", numBins, 0., 800.);
  EvisE23_uk->GetXaxis()->SetTitle("Energy (keV)");
  TH1F *EvisW23_uk = new TH1F("EvisW23_uk","EvisW Type 2/3", numBins, 0., 800.); 
  EvisW23_uk->GetXaxis()->SetTitle("Energy (keV)");

  TH1I *Type_mpm = new TH1I("Type_mpm","Event Types", 4, 0, 4);
  TH1I *Type_uk = new TH1I("Type_uk","Event Types", 4, 0, 4);
  

  EvisEALL_mpm->SetLineColor(kRed);
  EvisWALL_mpm->SetLineColor(kRed);
  EvisEALL_uk->SetLineColor(kBlue);
  EvisWALL_uk->SetLineColor(kBlue);
  
  EvisE0_mpm->SetLineColor(kRed);
  EvisW0_mpm->SetLineColor(kRed);
  EvisE0_uk->SetLineColor(kBlue);
  EvisW0_uk->SetLineColor(kBlue);

  EvisE1_mpm->SetLineColor(kRed);
  EvisW1_mpm->SetLineColor(kRed);
  EvisE1_uk->SetLineColor(kBlue);
  EvisW1_uk->SetLineColor(kBlue);

  EvisE23_mpm->SetLineColor(kRed);
  EvisW23_mpm->SetLineColor(kRed);
  EvisE23_uk->SetLineColor(kBlue);
  EvisW23_uk->SetLineColor(kBlue);

  Type_mpm->SetLineColor(kRed);
  Type_uk->SetLineColor(kBlue);

  try {
    
    AsymmetryBase *UK = new AsymmetryBase(octetNum, enBinWidth, 50., true);
    AsymmetryBase *mpm = new AsymmetryBase(octetNum, enBinWidth, 50., false);
    
    UK->calcBGsubtractedEvts();
    mpm->calcBGsubtractedEvts();
    
    std::vector <double> evts;
    
    for (int bin=0; bin<numBins; bin++) {
      double aveBinEn = ((double)bin+.5)*(double)enBinWidth;
      
      for (int i=0;i<3;i++) {
	evts = UK->getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,i);
	EvisEALL_uk->Fill(aveBinEn,evts[0]);
	EvisWALL_uk->Fill(aveBinEn,evts[1]);
	
	evts = mpm->getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,i);
	EvisEALL_mpm->Fill(aveBinEn,evts[0]);
	EvisWALL_mpm->Fill(aveBinEn,evts[1]);
      }

      evts = UK->getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,0);
      EvisE0_uk->Fill(aveBinEn,evts[0]);
      EvisW0_uk->Fill(aveBinEn,evts[1]);
      
      evts = mpm->getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,0);
      EvisE0_mpm->Fill(aveBinEn,evts[0]);
      EvisW0_mpm->Fill(aveBinEn,evts[1]);

      evts = UK->getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,1);
      EvisE1_uk->Fill(aveBinEn,evts[0]);
      EvisW1_uk->Fill(aveBinEn,evts[1]);
      
      evts = mpm->getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,1);
      EvisE1_mpm->Fill(aveBinEn,evts[0]);
      EvisW1_mpm->Fill(aveBinEn,evts[1]);

      evts = UK->getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,2);
      EvisE23_uk->Fill(aveBinEn,evts[0]);
      EvisW23_uk->Fill(aveBinEn,evts[1]);
      
      evts = mpm->getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,2);
      EvisE23_mpm->Fill(aveBinEn,evts[0]);
      EvisW23_mpm->Fill(aveBinEn,evts[1]);

      evts = UK->getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,3);
      EvisE23_uk->Fill(aveBinEn,evts[0]);
      EvisW23_uk->Fill(aveBinEn,evts[1]);
      
      evts = mpm->getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,3);
      EvisE23_mpm->Fill(aveBinEn,evts[0]);
      EvisW23_mpm->Fill(aveBinEn,evts[1]);

    }
    
    evts = UK->getNumBGsubtrEvts(enWinLow,enWinHigh,0);
    Type_uk->SetBinContent(0, (int)(evts[0]+evts[1]));
    evts = UK->getNumBGsubtrEvts(enWinLow,enWinHigh,1);
    Type_uk->SetBinContent(1, (int)(evts[0]+evts[1]));
    evts = UK->getNumBGsubtrEvts(enWinLow,enWinHigh,2);
    Type_uk->SetBinContent(2, (int)(evts[0]+evts[1]));
    evts = UK->getNumBGsubtrEvts(enWinLow,enWinHigh,3);
    Type_uk->SetBinContent(3, (int)(evts[0]+evts[1]));

    evts = mpm->getNumBGsubtrEvts(enWinLow,enWinHigh,0);
    Type_mpm->SetBinContent(0, (int)(evts[0]+evts[1]));
    evts = mpm->getNumBGsubtrEvts(enWinLow,enWinHigh,1);
    Type_mpm->SetBinContent(1, (int)(evts[0]+evts[1]));
    evts = mpm->getNumBGsubtrEvts(enWinLow,enWinHigh,2);
    Type_mpm->SetBinContent(2, (int)(evts[0]+evts[1]));
    evts = mpm->getNumBGsubtrEvts(enWinLow,enWinHigh,3);
    Type_mpm->SetBinContent(3, (int)(evts[0]+evts[1]));

  }

  catch(const char* ex){
    std::cerr << "Error: " << ex << std::endl;
  }
    
 

  TCanvas *cEvisALL = new TCanvas("cEvisALL"," ", 1400., 600.);
  cEvisALL->Divide(2,1);
  cEvisALL->cd(1);
  
  EvisEALL_uk->Draw();
  EvisEALL_mpm->Draw("SAME");
  
  TLegend *legALL = new TLegend(0.55,0.75,0.85,0.875);
  legALL->AddEntry(EvisEALL_mpm,"  MPM","l");
  legALL->AddEntry(EvisEALL_uk,"  UK","l");
  legALL->Draw();
  
  cEvisALL->cd(2);
  
  EvisWALL_uk->Draw();
  EvisWALL_mpm->Draw("SAME");
  
  std::string pdfstart = pdfFileBase + "[";
  cEvisALL->Print(pdfstart.c_str());
  cEvisALL->Print(pdfFileBase.c_str());
  
  
  TCanvas *cEvis0 = new TCanvas("cEvis0"," ", 1400., 600.);
  cEvis0->Divide(2,1);
  cEvis0->cd(1);
  
  EvisE0_uk->Draw();
  EvisE0_mpm->Draw("SAME");
  
  TLegend *leg0 = new TLegend(0.55,0.75,0.85,0.875);
  leg0->AddEntry(EvisE0_mpm,"  MPM","l");
  leg0->AddEntry(EvisE0_uk,"  UK","l");
  leg0->Draw();
  
  cEvis0->cd(2);
  
  EvisW0_uk->Draw();
  EvisW0_mpm->Draw("SAME");
  
  cEvis0->Print(pdfFileBase.c_str());
  
  
  TCanvas *cEvis1 = new TCanvas("cEvis0"," ", 1400., 600.);
  cEvis1->Divide(2,1);
  cEvis1->cd(1);
  
  EvisE1_uk->Draw();
  EvisE1_mpm->Draw("SAME");
  
  TLegend *leg1 = new TLegend(0.55,0.75,0.85,0.875);
  leg1->AddEntry(EvisE1_mpm,"  MPM","l");
  leg1->AddEntry(EvisE1_uk,"  UK","l");
  leg1->Draw();
  
  cEvis1->cd(2);
  
  EvisW1_uk->Draw();
  EvisW1_mpm->Draw("SAME");
  
  cEvis1->Print(pdfFileBase.c_str());
  
  
  TCanvas *cEvis23 = new TCanvas("cEvis23"," ", 1400., 600.);
  cEvis23->Divide(2,1);
  cEvis23->cd(1);
  
  EvisE23_uk->Draw();
  EvisE23_mpm->Draw("SAME");
  
  TLegend *leg23 = new TLegend(0.55,0.75,0.85,0.875);
  leg23->AddEntry(EvisE23_mpm,"  MPM","l");
  leg23->AddEntry(EvisE23_uk,"  UK","l");
  leg23->Draw();
  
  cEvis23->cd(2);
  
  EvisW23_uk->Draw();
  EvisW23_mpm->Draw("SAME");
  
  cEvis23->Print(pdfFileBase.c_str());
  
  
  TCanvas *cTypes = new TCanvas("cTypes"," ", 1400., 600.);
  cTypes->Divide(2,1);
  cTypes->cd(1);
  gPad->SetLogy();
  
  Type_uk->Draw();
  Type_mpm->Draw("SAME");
  
  Type_uk->SetMinimum(1);
  
  TLegend *legTypes = new TLegend(0.55,0.75,0.85,0.875);
  legTypes->AddEntry(Type_mpm,"  MPM","l");
  legTypes->AddEntry(Type_uk,"  UK","l");
  legTypes->Draw();
  
  cTypes->Print(pdfFileBase.c_str());
  std::string pdfend = pdfFileBase + "]";
  cTypes->Print(pdfend.c_str());

  return 0;

}
    
	
 
