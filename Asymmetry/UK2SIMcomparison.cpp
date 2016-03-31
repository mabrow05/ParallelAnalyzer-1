/*
 Compares the BG subtracted data spectra to simulation data which has 
 detector effects already taken into account
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
    std::cout << "Usage: ./UK2SIMcomparison [octet start] [octet end] [energy window low] [energy window high] [energy bin width]\n";
    std::cout << "The code will produce comparisons for every octet in the range given, on an octet-by-octet basis\n";
    exit(0);
  }

  gStyle->SetOptStat(0);

  int octetNumStart = atoi(argv[1]);
  int octetNumEnd = atoi(argv[2]);
  double enWinLow = (double) atof(argv[3]);
  double enWinHigh = (double) atof(argv[4]);
  int enBinWidth = atoi(argv[5]);

  int numBins = 800/enBinWidth;

  double normFactor=0.;

 
  for (int octetNum=octetNumStart; octetNum<octetNumEnd+1; octetNum++) {
 
    try {
      std::string pdfFileBase = "comparison_plots/SimulationComp_Octet"+itos(octetNum)+"_evis_"+itos((int)enWinLow)+"-"+itos((int)enWinHigh)+".pdf";
 
      //Making histograms for comparison  
      TH1D EreconEALL_sim("EreconEALL_sim","Erecon East Type ALL", numBins, 0., 800.);
      TH1D EreconWALL_sim("EreconWALL_sim","Erecon West Type ALL", numBins, 0., 800.);
      TH1D EreconEALL_uk("EreconEALL_uk","Erecon East Type ALL", numBins, 0., 800.);
      EreconEALL_uk.GetXaxis()->SetTitle("Energy (keV)");
      TH1D EreconWALL_uk("EreconWALL_uk","Erecon West Type ALL", numBins, 0., 800.);
      EreconWALL_uk.GetXaxis()->SetTitle("Energy (keV)");

      TH1D EreconE0_sim("EreconE0_sim","Erecon East Type 0", numBins, 0., 800.);
      TH1D EreconW0_sim("EreconW0_sim","Erecon West Type 0", numBins, 0., 800.);
      TH1D EreconE0_uk("EreconE0_uk","Erecon East Type 0", numBins, 0., 800.);
      EreconE0_uk.GetXaxis()->SetTitle("Energy (keV)");
      TH1D EreconW0_uk("EreconW0_uk","Erecon West Type 0", numBins, 0., 800.);
      EreconW0_uk.GetXaxis()->SetTitle("Energy (keV)");

      TH1D EreconE1_sim("EreconE1_sim","Erecon East Type 1", numBins, 0., 800.);
      TH1D EreconW1_sim("EreconW1_sim","Erecon West Type 1", numBins, 0., 800.);
      TH1D EreconE1_uk("EreconE1_uk","Erecon East Type 1", numBins, 0., 800.);
      EreconE1_uk.GetXaxis()->SetTitle("Energy (keV)");
      TH1D EreconW1_uk("EreconW1_uk","Erecon West Type 1", numBins, 0., 800.);
      EreconW1_uk.GetXaxis()->SetTitle("Energy (keV)");

      TH1D EreconE23_sim("EreconE23_sim","Erecon East Type 2/3", numBins, 0., 800.);
      TH1D EreconW23_sim("EreconW23_sim","Erecon West Type 2/3", numBins, 0., 800.);
      TH1D EreconE23_uk("EreconE23_uk","Erecon East Type 2/3", numBins, 0., 800.);
      EreconE23_uk.GetXaxis()->SetTitle("Energy (keV)");
      TH1D EreconW23_uk("EreconW23_uk","Erecon West Type 2/3", numBins, 0., 800.); 
      EreconW23_uk.GetXaxis()->SetTitle("Energy (keV)");

      TH1D Type_sim("Type_sim","Event Types", 4, 0., 4.);
      TH1D Type_uk("Type_uk","Event Types", 4, 0., 4.);
  

      EreconEALL_sim.SetLineColor(kRed);
      EreconWALL_sim.SetLineColor(kRed);
      EreconEALL_uk.SetLineColor(kBlue);
      EreconWALL_uk.SetLineColor(kBlue);
  
      EreconE0_sim.SetLineColor(kRed);
      EreconW0_sim.SetLineColor(kRed);
      EreconE0_uk.SetLineColor(kBlue);
      EreconW0_uk.SetLineColor(kBlue);

      EreconE1_sim.SetLineColor(kRed);
      EreconW1_sim.SetLineColor(kRed);
      EreconE1_uk.SetLineColor(kBlue);
      EreconW1_uk.SetLineColor(kBlue);

      EreconE23_sim.SetLineColor(kRed);
      EreconW23_sim.SetLineColor(kRed);
      EreconE23_uk.SetLineColor(kBlue);
      EreconW23_uk.SetLineColor(kBlue);

      Type_sim.SetLineColor(kRed);
      Type_uk.SetLineColor(kBlue);

      int uk0E=0, sim0E=0, uk1E=0, sim1E=0, uk23E=0, sim23E=0, ukTotE=0, simTotE=0, uk0W=0, sim0W=0, uk1W=0, sim1W=0, uk23W=0, sim23W=0, ukTotW=0, simTotW=0; //Event totals
      double ukAsym = 0., simAsym=0., ukAsymError = 0., simAsymError=0.;
   

      { //Setting scope for first OctetAsymmetry so that it will be deleted to clear memory
	OctetAsymmetry UK(octetNum, enBinWidth, 50., true);
	
	UK.calcTotalAsymmetry(225.,675.,1);
	ukAsym = UK.returnTotalAsymmetry();
	ukAsymError = UK.returnTotalAsymmetryError();
	std::cout << "UK Asym " << ukAsym << std::endl;
      
	UK.calcBGsubtractedEvts();
	
	std::vector <double> evts;
    
	for (int bin=0; bin<numBins; bin++) {

	  double aveBinEn = ((double)bin+.5)*(double)enBinWidth;

	  evts = UK.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,0);
	  EreconE0_uk.Fill(aveBinEn,evts[0]);
	  EreconW0_uk.Fill(aveBinEn,evts[1]);
	  EreconEALL_uk.Fill(aveBinEn,evts[0]);
	  EreconWALL_uk.Fill(aveBinEn,evts[1]);
      
	  evts = UK.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,1);
	  EreconE1_uk.Fill(aveBinEn,evts[0]);
	  EreconW1_uk.Fill(aveBinEn,evts[1]);
	  EreconEALL_uk.Fill(aveBinEn,evts[0]);
	  EreconWALL_uk.Fill(aveBinEn,evts[1]);
      
	  evts = UK.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,2);
	  EreconE23_uk.Fill(aveBinEn,evts[0]);
	  EreconW23_uk.Fill(aveBinEn,evts[1]);
	  EreconEALL_uk.Fill(aveBinEn,evts[0]);
	  EreconWALL_uk.Fill(aveBinEn,evts[1]);
      
	  evts = UK.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,3);
	  EreconE23_uk.Fill(aveBinEn,evts[0]);
	  EreconW23_uk.Fill(aveBinEn,evts[1]);
	  EreconEALL_uk.Fill(aveBinEn,evts[0]);
	  EreconWALL_uk.Fill(aveBinEn,evts[1]);
   
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

      }

      { //Setting scope for second OctetAsymmetry so that it will be deleted to clear memory
	OctetAsymmetry SIM(octetNum, enBinWidth, 50., true, true);
	std::cout << "Made it here\n";
	SIM.calcTotalAsymmetry(225.,675.,1);
	simAsym = SIM.returnTotalAsymmetry();
	simAsymError = SIM.returnTotalAsymmetryError();
	std::cout << "sim Asym " << simAsym << std::endl;
      
	SIM.calcBGsubtractedEvts();
	
	std::vector <double> evts;
	evts = SIM.getNumBGsubtrEvts(enWinLow,enWinHigh,0);
	normFactor = (uk0E+uk0W)/(evts[0]+evts[1]);
    
	for (int bin=0; bin<numBins; bin++) {

	  double aveBinEn = ((double)bin+.5)*(double)enBinWidth;
      
	  evts = SIM.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,0);
	  EreconE0_sim.Fill(aveBinEn,evts[0]*normFactor);
	  EreconW0_sim.Fill(aveBinEn,evts[1]*normFactor);
	  EreconEALL_sim.Fill(aveBinEn,evts[0]*normFactor);
	  EreconWALL_sim.Fill(aveBinEn,evts[1]*normFactor);
      
	  evts = SIM.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,1);
	  EreconE1_sim.Fill(aveBinEn,evts[0]*normFactor);
	  EreconW1_sim.Fill(aveBinEn,evts[1]*normFactor);
	  EreconEALL_sim.Fill(aveBinEn,evts[0]*normFactor);
	  EreconWALL_sim.Fill(aveBinEn,evts[1]*normFactor);
      
	  evts = SIM.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,2);
	  EreconE23_sim.Fill(aveBinEn,evts[0]*normFactor);
	  EreconW23_sim.Fill(aveBinEn,evts[1]*normFactor);
	  EreconEALL_sim.Fill(aveBinEn,evts[0]*normFactor);
	  EreconWALL_sim.Fill(aveBinEn,evts[1]*normFactor);
      
	  evts = SIM.getNumBGsubtrEvts(bin*enBinWidth,(bin+1)*enBinWidth,3);
	  EreconE23_sim.Fill(aveBinEn,evts[0]*normFactor);
	  EreconW23_sim.Fill(aveBinEn,evts[1]*normFactor);
	  EreconEALL_sim.Fill(aveBinEn,evts[0]*normFactor);
	  EreconWALL_sim.Fill(aveBinEn,evts[1]*normFactor);
   
	}
    
	evts = SIM.getNumBGsubtrEvts(enWinLow,enWinHigh,0);
	sim0E = (int)(evts[0]*normFactor); sim0W = (int)(evts[1]*normFactor);
	Type_sim.SetBinContent(1, (int)((evts[0]*normFactor)+(evts[1]*normFactor)));
	evts = SIM.getNumBGsubtrEvts(enWinLow,enWinHigh,1);
	sim1E = (int)(evts[0]*normFactor); sim1W = (int)(evts[1]*normFactor);
	Type_sim.SetBinContent(2, (int)((evts[0]*normFactor)+(evts[1]*normFactor)));
	evts = SIM.getNumBGsubtrEvts(enWinLow,enWinHigh,2);
	sim23E = (int)(evts[0]*normFactor); sim23W = (int)(evts[1]*normFactor);
	Type_sim.SetBinContent(3, (int)((evts[0]*normFactor)+(evts[1]*normFactor)));
	evts = SIM.getNumBGsubtrEvts(enWinLow,enWinHigh,3);
	sim23E += (int)(evts[0]*normFactor); sim23W += (int)(evts[1]*normFactor);
	Type_sim.SetBinContent(4, (int)((evts[0]*normFactor)+(evts[1]*normFactor)));
      
	simTotE = sim0E+sim1E+sim23E;
	simTotW = sim0W+sim1W+sim23W;

      }

    

    
    
      std::cout << "uk Asym " << ukAsym << std::endl;
      std::cout << "sim Asym " << simAsym << std::endl;
      
      TCanvas cEvisALL("cEvisALL"," ", 1400., 600.);
      cEvisALL.Divide(2,1);
      cEvisALL.cd(1);
  
      EreconEALL_uk.Draw();
      EreconEALL_sim.Draw("SAME");
  
      TLegend legALL(0.55,0.75,0.85,0.875);
      legALL.AddEntry(&EreconEALL_sim,"  SIM","l");
      legALL.AddEntry(&EreconEALL_uk,"  UK","l");
      legALL.Draw();
  
      cEvisALL.cd(2);
  
      EreconWALL_uk.Draw();
      EreconWALL_sim.Draw("SAME");
  
      std::string pdfstart = pdfFileBase + "[";
      cEvisALL.Print(pdfstart.c_str());
      cEvisALL.Print(pdfFileBase.c_str());
  
  
      TCanvas cEvis0("cEvis0"," ", 1400., 600.);
      cEvis0.Divide(2,1);
      cEvis0.cd(1);
  
      EreconE0_uk.Draw();
      EreconE0_sim.Draw("SAME");
  
      TLegend leg0(0.55,0.75,0.85,0.875);
      leg0.AddEntry(&EreconE0_sim,"  SIM","l");
      leg0.AddEntry(&EreconE0_uk,"  UK","l");
      leg0.Draw();
  
      cEvis0.cd(2);
  
      EreconW0_uk.Draw();
      EreconW0_sim.Draw("SAME");
  
      cEvis0.Print(pdfFileBase.c_str());
  
  
      TCanvas cEvis1("cEvis1"," ", 1400., 600.);
      cEvis1.Divide(2,1);
      cEvis1.cd(1);
  
      EreconE1_uk.Draw();
      EreconE1_sim.Draw("SAME");
  
      TLegend leg1(0.55,0.75,0.85,0.875);
      leg1.AddEntry(&EreconE1_sim,"  SIM","l");
      leg1.AddEntry(&EreconE1_uk,"  UK","l");
      leg1.Draw();
  
      cEvis1.cd(2);
  
      EreconW1_uk.Draw();
      EreconW1_sim.Draw("SAME");
  
      cEvis1.Print(pdfFileBase.c_str());
    
  
      TCanvas cEvis23("cEvis23"," ", 1400., 600.);
      cEvis23.Divide(2,1);
      cEvis23.cd(1);
  
      EreconE23_uk.Draw();
      EreconE23_sim.Draw("SAME");
  
      TLegend leg23(0.55,0.75,0.85,0.875);
      leg23.AddEntry(&EreconE23_sim,"  SIM","l");
      leg23.AddEntry(&EreconE23_uk,"  UK","l");
      leg23.Draw();
  
      cEvis23.cd(2);
  
      EreconW23_uk.Draw();
      EreconW23_sim.Draw("SAME");
  
      cEvis23.Print(pdfFileBase.c_str());
  
  
      TCanvas cTypes("cTypes"," ", 1400., 600.);
      cTypes.Divide(2,1);
      cTypes.cd(1);
      gPad->SetLogy();
  
      Type_uk.Draw();
      Type_sim.Draw("SAME");
  
      Type_uk.SetMinimum(1);
  
      TLegend legTypes(0.55,0.75,0.85,0.875);
      legTypes.AddEntry(&Type_sim,"  SIM","l");
      legTypes.AddEntry(&Type_uk,"  UK","l");
      legTypes.Draw();
  
      std::cout << "UK Asym " << ukAsym << std::endl;

      cTypes.cd(2);

      TPaveText ukTxt(0.1, 0.25, 0.45, 0.75);
      Char_t temp[500];
      ukTxt.AddText("Data");
      //ukTxt.AddText(geometry.c_str());
      sprintf(temp, "Erecon = %.0f-%.0f keV", enWinLow, enWinHigh);
      ukTxt.AddText(temp);
      sprintf(temp, "Events: East %i   West %i", ukTotE, ukTotW);
      ukTxt.AddText(temp);
      sprintf(temp, "Total:  %i",ukTotE+ukTotW);
      ukTxt.AddText(temp);
      ukTxt.AddText(" ");
      sprintf(temp, "Raw Asymmetry:  %.5f",float(ukTotE-ukTotW)/float(ukTotE+ukTotW));
      ukTxt.AddText(temp);
      sprintf(temp, "Octet Asymmetry:  %.5f +/- %.5f",ukAsym,ukAsymError);
      ukTxt.AddText(temp);
      ukTxt.AddText(" ");
      ukTxt.AddText("Type:     East:        West:       "); 
      sprintf(temp, "0         %.5f   %.5f",  (float)uk0E/(float)ukTotE, (float)uk0W/(float)ukTotW);
      ukTxt.AddText(temp); 
      sprintf(temp, "I        %.5f   %.5f", (float)uk1E/(float)ukTotE, (float)uk1W/(float)ukTotW);
      ukTxt.AddText(temp); 
      sprintf(temp, "II/III       %.5f   %.5f",  (float)uk23E/(float)ukTotE, (float)uk23W/(float)ukTotW);
      ukTxt.AddText(temp); 
      ukTxt.Draw();

      TPaveText simTxt(0.55, 0.25, 0.9, 0.75);
      simTxt.AddText("Simulation");
      //simTxt.AddText(geometry.c_str());
      sprintf(temp, "Erecon = %.0f-%.0f keV", enWinLow, enWinHigh);
      simTxt.AddText(temp);
      sprintf(temp, "Events: East %i   West %i", simTotE, simTotW);
      simTxt.AddText(temp);
      sprintf(temp, "Total: %i",simTotE+simTotW);
      simTxt.AddText(temp);
      simTxt.AddText(" ");
      sprintf(temp, "Raw Asymmetry:  %.5f",float(simTotE-simTotW)/float(simTotE+simTotW));
      simTxt.AddText(temp);
      sprintf(temp, "Octet Asymmetry:  %.5f +/- %.5f",simAsym,simAsymError);
      simTxt.AddText(temp);
      simTxt.AddText(" ");
      simTxt.AddText("Type:     East:        West:       "); 
      sprintf(temp, "0         %.5f   %.5f",  (float)sim0E/(float)simTotE, (float)sim0W/(float)simTotW);
      simTxt.AddText(temp); 
      sprintf(temp, "I        %.5f   %.5f", (float)sim1E/(float)simTotE, (float)sim1W/(float)simTotW);
      simTxt.AddText(temp); 
      sprintf(temp, "II/III       %.5f   %.5f",  (float)sim23E/(float)simTotE, (float)sim23W/(float)simTotW);
      simTxt.AddText(temp); 
      simTxt.Draw();
  
      cTypes.Print(pdfFileBase.c_str());
      std::string pdfend = pdfFileBase + "]";
      cTypes.Print(pdfend.c_str());

    }

    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }

  return 0;

}
    
	
 
