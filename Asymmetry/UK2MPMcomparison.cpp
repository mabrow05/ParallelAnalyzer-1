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
  if (argc!=6) {
    std::cout << "Usage: ./UK2MPMcomparison [octet start] [octet end] [energy window low] [energy window high] [energy bin width]\n";
    std::cout << "The code will produce comparisons for every octet in the range given, on an octet-by-octet basis\n";
    exit(0);
  }

  gStyle->SetOptStat(0);

  int octetNumStart = atoi(argv[1]);
  int octetNumEnd = atoi(argv[2]);
  double enWinLow = atof(argv[3]);
  double enWinHigh = atof(argv[4]);
  int enBinWidth = atoi(argv[5]);

  //OctetAsymmetry UK;
  //OctetAsymmetry mpm;

  TH1F EvisEALL_mpm;
  TH1F EvisWALL_mpm;
  TH1F EvisEALL_uk;
  TH1F EvisWALL_uk;

  TH1F EvisE0_mpm;
  TH1F EvisW0_mpm;
  TH1F EvisE0_uk;
  TH1F EvisW0_uk;

  TH1F EvisE1_mpm;
  TH1F EvisW1_mpm;
  TH1F EvisE1_uk;
  TH1F EvisW1_uk;

  TH1F EvisE23_mpm;
  TH1F EvisW23_mpm;
  TH1F EvisE23_uk;
  TH1F EvisW23_uk;

  TH1F Type_mpm;
  TH1F Type_uk;

  int numBins = 800/enBinWidth;

  double mpmAsym2;

  for (int octetNum=octetNumStart; octetNum<octetNumEnd+1; octetNum++) {
 
    std::string pdfFileBase = "comparison_plots/calibrationRunComp_Octet"+itos(octetNum)+"_evis_"+itos((int)enWinLow)+"-"+itos((int)enWinHigh)+".pdf";
 
    //Making histograms for comparison  
    EvisEALL_mpm = TH1F("EvisEALL_mpm","EvisE Type ALL", numBins, 0., 800.);
    EvisWALL_mpm = TH1F("EvisWALL_mpm","EvisW Type ALL", numBins, 0., 800.);
    EvisEALL_uk = TH1F("EvisEALL_uk","EvisE Type ALL", numBins, 0., 800.);
    EvisEALL_uk.GetXaxis()->SetTitle("Energy (keV)");
    EvisWALL_uk = TH1F("EvisWALL_uk","EvisW Type ALL", numBins, 0., 800.);
    EvisWALL_uk.GetXaxis()->SetTitle("Energy (keV)");

    EvisE0_mpm = TH1F("EvisE0_mpm","EvisE Type 0", numBins, 0., 800.);
    EvisW0_mpm = TH1F("EvisW0_mpm","EvisW Type 0", numBins, 0., 800.);
    EvisE0_uk = TH1F("EvisE0_uk","EvisE Type 0", numBins, 0., 800.);
    EvisE0_uk.GetXaxis()->SetTitle("Energy (keV)");
    EvisW0_uk = TH1F("EvisW0_uk","EvisW Type 0", numBins, 0., 800.);
    EvisW0_uk.GetXaxis()->SetTitle("Energy (keV)");

    EvisE1_mpm = TH1F("EvisE1_mpm","EvisE Type 1", numBins, 0., 800.);
    EvisW1_mpm = TH1F("EvisW1_mpm","EvisW Type 1", numBins, 0., 800.);
    EvisE1_uk = TH1F("EvisE1_uk","EvisE Type 1", numBins, 0., 800.);
    EvisE1_uk.GetXaxis()->SetTitle("Energy (keV)");
    EvisW1_uk = TH1F("EvisW1_uk","EvisW Type 1", numBins, 0., 800.);
    EvisW1_uk.GetXaxis()->SetTitle("Energy (keV)");

    EvisE23_mpm = TH1F("EvisE23_mpm","EvisE Type 2/3", numBins, 0., 800.);
    EvisW23_mpm = TH1F("EvisW23_mpm","EvisW Type 2/3", numBins, 0., 800.);
    EvisE23_uk = TH1F("EvisE23_uk","EvisE Type 2/3", numBins, 0., 800.);
    EvisE23_uk.GetXaxis()->SetTitle("Energy (keV)");
    EvisW23_uk = TH1F("EvisW23_uk","EvisW Type 2/3", numBins, 0., 800.); 
    EvisW23_uk.GetXaxis()->SetTitle("Energy (keV)");

    Type_mpm = TH1F("Type_mpm","Event Types", 4, 0., 4.);
    Type_uk = TH1F("Type_uk","Event Types", 4, 0., 4.);
  

    EvisEALL_mpm.SetLineColor(kRed);
    EvisWALL_mpm.SetLineColor(kRed);
    EvisEALL_uk.SetLineColor(kBlue);
    EvisWALL_uk.SetLineColor(kBlue);
  
    EvisE0_mpm.SetLineColor(kRed);
    EvisW0_mpm.SetLineColor(kRed);
    EvisE0_uk.SetLineColor(kBlue);
    EvisW0_uk.SetLineColor(kBlue);

    EvisE1_mpm.SetLineColor(kRed);
    EvisW1_mpm.SetLineColor(kRed);
    EvisE1_uk.SetLineColor(kBlue);
    EvisW1_uk.SetLineColor(kBlue);

    EvisE23_mpm.SetLineColor(kRed);
    EvisW23_mpm.SetLineColor(kRed);
    EvisE23_uk.SetLineColor(kBlue);
    EvisW23_uk.SetLineColor(kBlue);

    Type_mpm.SetLineColor(kRed);
    Type_uk.SetLineColor(kBlue);

    Int_t uk0E=0, mpm0E=0, uk1E=0, mpm1E=0, uk23E=0, mpm23E=0, ukTotE=0, mpmTotE=0, uk0W=0, mpm0W=0, uk1W=0, mpm1W=0, uk23W=0, mpm23W=0, ukTotW=0, mpmTotW=0; //Event totals
    double ukAsym = 0., mpmAsym=0., ukAsymError = 0., mpmAsymError=0.;
    try {
    
      OctetAsymmetry UK(octetNum, enBinWidth, 50., true);
      OctetAsymmetry mpm(octetNum, enBinWidth, 50., false);

      
      UK.calcBGsubtractedEvts();
      std::cout << "UK Asym " << ukAsym << std::endl;
      //exit(0);
      mpm.calcBGsubtractedEvts();
      std::cout << "UK Asym " << ukAsym << std::endl;
      
    
      std::vector <double> evts;
    
       for (int bin=0; bin<numBins; bin++) {
	double aveBinEn = ((double)bin+.5)*(double)enBinWidth;
      
	for (int i=0;i<3;i++) {
	  evts = UK.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,i);
	  EvisEALL_uk.Fill(aveBinEn,evts[0]);
	  EvisWALL_uk.Fill(aveBinEn,evts[1]);
	
	  evts = mpm.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,i);
	  EvisEALL_mpm.Fill(aveBinEn,evts[0]);
	  EvisWALL_mpm.Fill(aveBinEn,evts[1]);
	}

	evts = UK.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,0);
	EvisE0_uk.Fill(aveBinEn,evts[0]);
	EvisW0_uk.Fill(aveBinEn,evts[1]);
      
	evts = mpm.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,0);
	EvisE0_mpm.Fill(aveBinEn,evts[0]);
	EvisW0_mpm.Fill(aveBinEn,evts[1]);

	evts = UK.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,1);
	EvisE1_uk.Fill(aveBinEn,evts[0]);
	EvisW1_uk.Fill(aveBinEn,evts[1]);
      
	evts = mpm.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,1);
	EvisE1_mpm.Fill(aveBinEn,evts[0]);
	EvisW1_mpm.Fill(aveBinEn,evts[1]);

	evts = UK.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,2);
	EvisE23_uk.Fill(aveBinEn,evts[0]);
	EvisW23_uk.Fill(aveBinEn,evts[1]);
      
	evts = mpm.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,2);
	EvisE23_mpm.Fill(aveBinEn,evts[0]);
	EvisW23_mpm.Fill(aveBinEn,evts[1]);

	evts = UK.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,3);
	EvisE23_uk.Fill(aveBinEn,evts[0]);
	EvisW23_uk.Fill(aveBinEn,evts[1]);
      
	evts = mpm.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,3);
	EvisE23_mpm.Fill(aveBinEn,evts[0]);
	EvisW23_mpm.Fill(aveBinEn,evts[1]);

      }
    
      evts = UK.getNumBGsubtrEvts(enWinLow,enWinHigh,0);
      uk0E = (int)evts[0]; uk0W = (int)evts[1];
      Type_uk.SetBinContent(1, (int)(evts[0]+evts[1]));
      evts = UK.getNumBGsubtrEvts(enWinLow,enWinHigh,1);
      uk1E = (int)evts[0]; uk1W = (int)evts[1];
      Type_uk.SetBinContent(2, (int)(evts[0]+evts[1]));
      evts = UK.getNumBGsubtrEvts(enWinLow,enWinHigh,2);
      uk23E = (int)evts[0]; uk23W = (int)evts[1];
      Type_uk.SetBinContent(3, (int)(evts[0]+evts[1]));
      evts = UK.getNumBGsubtrEvts(enWinLow,enWinHigh,3);
      uk23E += (int)evts[0]; uk23W += (int)evts[1];
      Type_uk.SetBinContent(4, (int)(evts[0]+evts[1]));

      ukTotE = uk0E+uk1E+uk23E;
      ukTotW = uk0W+uk1W+uk23W;

      evts = mpm.getNumBGsubtrEvts(enWinLow,enWinHigh,0);
      mpm0E = (int)evts[0]; mpm0W = (int)evts[1];
      Type_mpm.SetBinContent(1, (int)(evts[0]+evts[1]));
      evts = mpm.getNumBGsubtrEvts(enWinLow,enWinHigh,1);
      mpm1E = (int)evts[0]; mpm1W = (int)evts[1];
      Type_mpm.SetBinContent(2, (int)(evts[0]+evts[1]));
      evts = mpm.getNumBGsubtrEvts(enWinLow,enWinHigh,2);
      mpm23E = (int)evts[0]; mpm23W = (int)evts[1];
      Type_mpm.SetBinContent(3, (int)(evts[0]+evts[1]));
      evts = mpm.getNumBGsubtrEvts(enWinLow,enWinHigh,3);
      mpm23E += (int)evts[0]; mpm23W += (int)evts[1];
      Type_mpm.SetBinContent(4, (int)(evts[0]+evts[1]));

      mpmTotE = mpm0E+mpm1E+mpm23E;
      mpmTotW = mpm0W+mpm1W+mpm23W;
      
      
      UK.calcTotalAsymmetry(225.,675.,1);
      ukAsym = UK.returnTotalAsymmetry();
      ukAsymError = UK.returnTotalAsymmetryError();
      std::cout << "UK Asym " << ukAsym << std::endl;
      

      mpm.calcTotalAsymmetry(225.,675.,1);
      mpmAsym = mpm.returnTotalAsymmetry();
      mpmAsymError = mpm.returnTotalAsymmetryError();
      std::cout << "MPM Asym " << mpmAsym << std::endl;

      
      std::cout << "UK Asym " << ukAsym << std::endl;
      
      //delete UK; delete mpm;
      std::cout << "mpm Asym " << mpmAsym << std::endl;
      //exit(0);

    }

    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
    

    std::cout << "mpm Asym " << mpmAsym << std::endl;
    exit(0);
      
    TCanvas *cEvisALL = new TCanvas("cEvisALL"," ", 1400., 600.);
    cEvisALL->Divide(2,1);
    cEvisALL->cd(1);
  
    EvisEALL_uk.Draw();
    EvisEALL_mpm.Draw("SAME");
  
    TLegend *legALL = new TLegend(0.55,0.75,0.85,0.875);
    legALL->AddEntry(&EvisEALL_mpm,"  MPM","l");
    legALL->AddEntry(&EvisEALL_uk,"  UK","l");
    legALL->Draw();
  
    cEvisALL->cd(2);
  
    EvisWALL_uk.Draw();
    EvisWALL_mpm.Draw("SAME");
  
    std::string pdfstart = pdfFileBase + "[";
    cEvisALL->Print(pdfstart.c_str());
    cEvisALL->Print(pdfFileBase.c_str());
  
  
    TCanvas *cEvis0 = new TCanvas("cEvis0"," ", 1400., 600.);
    cEvis0->Divide(2,1);
    cEvis0->cd(1);
  
    EvisE0_uk.Draw();
    EvisE0_mpm.Draw("SAME");
  
    TLegend *leg0 = new TLegend(0.55,0.75,0.85,0.875);
    leg0->AddEntry(&EvisE0_mpm,"  MPM","l");
    leg0->AddEntry(&EvisE0_uk,"  UK","l");
    leg0->Draw();
  
    cEvis0->cd(2);
  
    EvisW0_uk.Draw();
    EvisW0_mpm.Draw("SAME");
  
    cEvis0->Print(pdfFileBase.c_str());
  
  
    TCanvas *cEvis1 = new TCanvas("cEvis1"," ", 1400., 600.);
    cEvis1->Divide(2,1);
    cEvis1->cd(1);
  
    EvisE1_uk.Draw();
    EvisE1_mpm.Draw("SAME");
  
    TLegend *leg1 = new TLegend(0.55,0.75,0.85,0.875);
    leg1->AddEntry(&EvisE1_mpm,"  MPM","l");
    leg1->AddEntry(&EvisE1_uk,"  UK","l");
    leg1->Draw();
  
    cEvis1->cd(2);
  
    EvisW1_uk.Draw();
    EvisW1_mpm.Draw("SAME");
  
    cEvis1->Print(pdfFileBase.c_str());
    
  
    TCanvas *cEvis23 = new TCanvas("cEvis23"," ", 1400., 600.);
    cEvis23->Divide(2,1);
    cEvis23->cd(1);
  
    EvisE23_uk.Draw();
    EvisE23_mpm.Draw("SAME");
  
    TLegend *leg23 = new TLegend(0.55,0.75,0.85,0.875);
    leg23->AddEntry(&EvisE23_mpm,"  MPM","l");
    leg23->AddEntry(&EvisE23_uk,"  UK","l");
    leg23->Draw();
  
    cEvis23->cd(2);
  
    EvisW23_uk.Draw();
    EvisW23_mpm.Draw("SAME");
  
    cEvis23->Print(pdfFileBase.c_str());
  
  
    TCanvas *cTypes = new TCanvas("cTypes"," ", 1400., 600.);
    cTypes->Divide(2,1);
    cTypes->cd(1);
    gPad->SetLogy();
  
    Type_uk.Draw();
    Type_mpm.Draw("SAME");
  
    Type_uk.SetMinimum(1);
  
    TLegend *legTypes = new TLegend(0.55,0.75,0.85,0.875);
    legTypes->AddEntry(&Type_mpm,"  MPM","l");
    legTypes->AddEntry(&Type_uk,"  UK","l");
    legTypes->Draw();
  
    std::cout << "UK Asym " << ukAsym << std::endl;

    cTypes->cd(2);

    TPaveText *ukTxt = new TPaveText(0.1, 0.25, 0.45, 0.75);
    Char_t temp[500];
    ukTxt->AddText("UK Analyzer");
    //ukTxt->AddText(geometry.c_str());
    sprintf(temp, "Evis = %.0f-%.0f keV", enWinLow, enWinHigh);
    ukTxt->AddText(temp);
    sprintf(temp, "Events: East %i   West %i", ukTotE, ukTotW);
    ukTxt->AddText(temp);
    sprintf(temp, "Total:  %i",ukTotE+ukTotW);
    ukTxt->AddText(temp);
    ukTxt->AddText(" ");
    sprintf(temp, "Raw Asymmetry:  %.5f",float(ukTotE-ukTotW)/float(ukTotE+ukTotW));
    ukTxt->AddText(temp);
    sprintf(temp, "Octet Asymmetry:  %.5f +/- %.5f",ukAsym,ukAsymError);
    ukTxt->AddText(temp);
    ukTxt->AddText(" ");
    ukTxt->AddText("Type:     East:        West:       "); 
    sprintf(temp, "0         %.5f   %.5f",  (float)uk0E/(float)ukTotE, (float)uk0W/(float)ukTotW);
    ukTxt->AddText(temp); 
    sprintf(temp, "I        %.5f   %.5f", (float)uk1E/(float)ukTotE, (float)uk1W/(float)ukTotW);
    ukTxt->AddText(temp); 
    sprintf(temp, "II/III       %.5f   %.5f",  (float)uk23E/(float)ukTotE, (float)uk23W/(float)ukTotW);
    ukTxt->AddText(temp); 
    ukTxt->Draw();

    TPaveText *mpmTxt = new TPaveText(0.55, 0.25, 0.9, 0.75);
    mpmTxt->AddText("MPM Analyzer");
    //mpmTxt->AddText(geometry.c_str());
    sprintf(temp, "Evis = %.0f-%.0f keV", enWinLow, enWinHigh);
    mpmTxt->AddText(temp);
    sprintf(temp, "Events: East %i   West %i", mpmTotE, mpmTotW);
    mpmTxt->AddText(temp);
    sprintf(temp, "Total: %i",mpmTotE+mpmTotW);
    mpmTxt->AddText(temp);
    mpmTxt->AddText(" ");
    sprintf(temp, "Raw Asymmetry:  %.5f",float(mpmTotE-mpmTotW)/float(mpmTotE+mpmTotW));
    mpmTxt->AddText(temp);
    sprintf(temp, "Octet Asymmetry:  %.5f +/- %.5f",mpmAsym,mpmAsymError);
    mpmTxt->AddText(temp);
    mpmTxt->AddText(" ");
    mpmTxt->AddText("Type:     East:        West:       "); 
    sprintf(temp, "0         %.5f   %.5f",  (float)mpm0E/(float)mpmTotE, (float)mpm0W/(float)mpmTotW);
    mpmTxt->AddText(temp); 
    sprintf(temp, "I        %.5f   %.5f", (float)mpm1E/(float)mpmTotE, (float)mpm1W/(float)mpmTotW);
    mpmTxt->AddText(temp); 
    sprintf(temp, "II/III       %.5f   %.5f",  (float)mpm23E/(float)mpmTotE, (float)mpm23W/(float)mpmTotW);
    mpmTxt->AddText(temp); 
    mpmTxt->Draw();
  
    cTypes->Print(pdfFileBase.c_str());
    std::string pdfend = pdfFileBase + "]";
    cTypes->Print(pdfend.c_str());


    /*    delete EvisEALL_mpm;
    delete EvisWALL_mpm;
    delete EvisEALL_uk;
    delete EvisWALL_uk;

    delete EvisE0_mpm;
    delete EvisW0_mpm;
    delete EvisE0_uk;
    delete EvisW0_uk;

    delete EvisE1_mpm;
    delete EvisW1_mpm;
    delete EvisE1_uk;
    delete EvisW1_uk;

    delete EvisE23_mpm;
    delete EvisW23_mpm;
    delete EvisE23_uk;
    delete EvisW23_uk;

    delete Type_mpm;
    delete Type_uk;

    delete cEvisALL; delete cEvis0; delete cEvis1; delete cEvis23; delete cTypes;
    delete legALL; delete leg0; delete leg1; delete leg23; delete legTypes;
    delete ukTxt; delete mpmTxt;*/
  }


  

  return 0;

}
    
	
 
