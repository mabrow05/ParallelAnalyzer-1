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
#include <iomanip>

#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TGraphErrors.h>
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
      //std::string pdfFileBase = "comparison_plots/SimulationComp_Octet"+itos(octetNum)+"_evis_"+itos((int)enWinLow)+"-"+itos((int)enWinHigh)+".pdf";
 
      TFile *f = new TFile("test.root","RECREATE");
      
      /*//Making graphs for comparison  
      TGraphErrors EreconALL_sim;
      //TGraphErrors EreconALL_uk;
      //EreconALL_uk.GetXaxis()->SetTitle("Energy (keV)");
      //EreconALL_uk.SetTitle("All Types Super-Sum");

      TGraphErrors Erecon0_sim;
      TGraphErrors Erecon0_uk;
      Erecon0_uk.GetXaxis()->SetTitle("Energy (keV)");
      Erecon0_uk.SetTitle("Type 0 Super-Sum");

      TGraphErrors Erecon1_sim;
      TGraphErrors Erecon1_uk;
      Erecon1_uk.GetXaxis()->SetTitle("Energy (keV)");
      Erecon1_uk.SetTitle("Type 1 Super-Sum");

      TGraphErrors Erecon23_sim;
      TGraphErrors Erecon23_uk;
      Erecon23_uk.GetXaxis()->SetTitle("Energy (keV)");
      Erecon23_uk.SetTitle("Type 2/3 Super-Sum");
    
      */
      TH1D Type_sim("Type_sim","Event Types", 4, 0., 4.);
      TH1D Type_uk("Type_uk","Event Types", 4, 0., 4.);
      

      /*EreconALL_sim.SetMarkerColor(kRed);
      EreconALL_uk.SetMarkerColor(kBlue);
     
      Erecon0_sim.SetMarkerColor(kRed);    
      Erecon0_uk.SetMarkerColor(kBlue);

      Erecon1_sim.SetMarkerColor(kRed); 
      Erecon1_uk.SetMarkerColor(kBlue);

      Erecon23_sim.SetMarkerColor(kRed);    
      Erecon23_uk.SetMarkerColor(kBlue);    

      Type_sim.SetLineColor(kRed);
      Type_uk.SetLineColor(kBlue);*/

      //set up energy binning vector
      std::vector <double> enBins(numBins,0.);
      
      for (int bin=0; bin<numBins; bin++) {
	enBins[bin] = ((double)bin+.5)*(double)enBinWidth;
      }

      int uk0E=0, sim0E=0, uk1E=0, sim1E=0, uk23E=0, sim23E=0, ukTotE=0, simTotE=0, uk0W=0, sim0W=0, uk1W=0, sim1W=0, uk23W=0, sim23W=0, ukTotW=0, simTotW=0; //Event totals
      double ukAsym = 0., simAsym=0., ukAsymError = 0., simAsymError=0.;
      
      double normFactor = 0., ukIntegral=0., simIntegral=0.;


      { //Setting scope for first OctetAsymmetry so that it will be deleted to clear memory
	OctetAsymmetry UK(octetNum, enBinWidth, 50., true, false);
	//UK.calcAsymmetryBinByBin(1);
	
	/*UK.calcTotalAsymmetry(225.,675.,1);
	ukAsym = UK.returnTotalAsymmetry();
	ukAsymError = UK.returnTotalAsymmetryError();
	std::cout << "UK Asym " << ukAsym << std::endl;
	*/
	//UK.calcBGsubtractedEvts();


	std::vector <double> evts;

	std::vector <double> superSum;
	std::vector <double> superSumError;
	std::vector <double> xErr(numBins,0.);

	
	
	//All event types
	UK.calcSuperSum(4);
	superSum = UK.returnSuperSum();
	superSumError = UK.returnSuperSumError();

	for (int n=0; n<numBins;n++) {
	  std::cout << enBins[n] << " " << superSum[n] << " " << superSumError[n] << std::endl; 
	  ukIntegral+=superSum[n];
	}

	f->cd();
	TGraphErrors *EreconALL_uk = new TGraphErrors(numBins,&enBins[0],&superSum[0],&xErr[0],&superSumError[0]);
	EreconALL_uk->SetMarkerColor(kBlue);
	EreconALL_uk->SetMarkerStyle(21);
	EreconALL_uk->GetXaxis()->SetTitle("Energy (keV)");
	EreconALL_uk->SetTitle("All Types Super-Sum");
	EreconALL_uk->Draw("AP");
	//EreconALL_uk->SetDirectory(f);
	
	EreconALL_uk->Write("EreconALL_uk");
	//f->Write();
	//f->Close();
	delete EreconALL_uk;
	//std::cout << ukIntegral << std::endl;
	//exit(0);

	
       
    
	/*evts = UK.getNumBGsubtrEvts(enWinLow,enWinHigh,0);
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
	*/
      }

      { //Setting scope for second OctetAsymmetry so that it will be deleted to clear memory
	OctetAsymmetry SIM(octetNum, enBinWidth, 50., true, true, false);
	/*std::cout << "Made it here\n";
	SIM.calcTotalAsymmetry(225.,675.,1);
	simAsym = SIM.returnTotalAsymmetry();
	simAsymError = SIM.returnTotalAsymmetryError();
	std::cout << "sim Asym " << simAsym << std::endl;
      
	SIM.calcBGsubtractedEvts();
	
	std::vector <double> evts;
	evts = SIM.getNumBGsubtrEvts(enWinLow,enWinHigh,0);
	normFactor = (uk0E+uk0W)/(evts[0]+evts[1]);
	*/
	std::vector <double> superSum;
	std::vector <double> superSumError;
	std::vector <double> xErr(numBins,0.);

	

	//All event types
	SIM.calcSuperSum(4);
	superSum = SIM.returnSuperSum();
	superSumError = SIM.returnSuperSumError();
	for (int n=0; n<numBins;n++) {
	  simIntegral+=superSum[n];
	  //std::cout << enBins[n] << " " << superSum[n] << " " << superSumError[n] << std::endl; 	  
	}
	normFactor = ukIntegral/simIntegral;
	for (int n=0; n<numBins;n++) {
	  superSum[n] = normFactor*superSum[n];
	  superSumError[n] = normFactor*superSumError[n];
	  std::cout << enBins[n] << " " << superSum[n] << " " << superSumError[n] << std::endl; 	  
	}
	f->cd();
	TGraphErrors *EreconALL_sim = new TGraphErrors(numBins,&enBins[0],&superSum[0],&xErr[0],&superSumError[0]);
	EreconALL_sim->SetMarkerColor(kBlue);
	EreconALL_sim->SetMarkerStyle(21);
	EreconALL_sim->GetXaxis()->SetTitle("Energy (keV)");
	EreconALL_sim->SetTitle("All Types Super-Sum");
	EreconALL_sim->Draw("AP");
	//EreconALL_uk->SetDirectory(f);
	
	EreconALL_sim->Write("EreconALL_sim");
	//f->Write();
	f->Close();
	delete EreconALL_sim;
	exit(0);
	
    
	/*evts = SIM.getNumBGsubtrEvts(enWinLow,enWinHigh,0);
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
	simTotW = sim0W+sim1W+sim23W;*/

      }

    }

    catch(const char* ex){
      std::cerr << "Error: " << ex << std::endl;
    }
  }

  return 0;

}

    
    
      /* std::cout << "uk Asym " << ukAsym << std::endl;
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
      */    
	
 
