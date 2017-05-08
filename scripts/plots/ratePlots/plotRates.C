#include <vector>
#include <algorithm>

std::vector < Float_t> readOctetFileForFGruns(int octet) {
  
  std::vector <Float_t> vec;
  
  char fileName[200];
  sprintf(fileName,"%s/All_Octets/octet_list_%i.dat",getenv("OCTET_LIST"),octet);
  std::ifstream infile(fileName);
  std::string runTypeHold;
  Float_t runNumberHold;
  
  // Populate the map with the runType and runNumber from this octet
  while (infile >> runTypeHold >> runNumberHold) {

    if ( runTypeHold=="A2" || runTypeHold=="A10" || 
	 runTypeHold=="B5" || runTypeHold=="B7" || 
	 runTypeHold=="B2" || runTypeHold=="B10" 
	 || runTypeHold=="A5" || runTypeHold=="A7" )  {
      
      vec.push_back(runNumberHold);
    
    }

  }

  infile.close();
 
  return vec;
};

void plotRates(TString year, TString anaCh) {

  gStyle->SetTitleSize(0.06,"t");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleYOffset(0.8);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetLabelSize(0.04,"xyz");
  //gStyle->SetOptFit(1111);
  gStyle->SetStatY(0.85);
  gStyle->SetStatX(0.975);
  gStyle->SetStatW(.09);

  gStyle->SetFillStyle(0000); 
  gStyle->SetStatStyle(0); 
  //gStyle->SetTitleStyle(0); 
  //gStyle->SetCanvasBorderSize(0); 
  //gStyle->SetFrameBorderSize(0); 
  //gStyle->SetLegendBorderSize(0); 
  //gStyle->SetStatBorderSize(0); 
  //gStyle->SetTitleBorderSize(0);

  

  //This will plot the rates as a function of Beta Runs

  Float_t binMin = 22;
  Float_t binMax = 68; // Not including this bin, so the sum goes from binMin->binMax-1

  
  std::vector <Int_t> badOct;// {7,60,61,62,63,64,65,66};
  Int_t octs[] = {7,60,61,62,63,64,65,66};

  for (int j=0; j<8; ++j) badOct.push_back(octs[j]);

  Int_t octMin, octMax;

  if (year==TString("2011-2012") ) {
    octMin = 0;
    octMax = 59;
  }

  else {
    octMin = 60;
    octMax = 121;
  }

  std::vector <Float_t> runs;
  
  for ( Int_t oct=octMin; oct<octMax+1; ++oct ) {
    if (std::find(badOct.begin(),badOct.end(),oct)!=badOct.end()) continue;
    std::vector<Float_t> octruns = readOctetFileForFGruns(oct);
    runs.insert(runs.end(),octruns.begin(),octruns.end());

  }

  std::vector < Float_t > runNumber;
  std::vector < std::vector<Float_t> > fgrate(2,std::vector<Float_t>(0));
  std::vector < std::vector<Float_t> > bgrate(2,std::vector<Float_t>(0));
  std::vector < std::vector<Float_t> > fgrateErr(2,std::vector<Float_t>(0));
  std::vector < std::vector<Float_t> > bgrateErr(2,std::vector<Float_t>(0));

  double eastMax=0.;
  double westMax=0.;

  for (unsigned i=0; i<runs.size(); ++i) {

    cout << "Processing Run : " << runs[i] << endl;
    
    // These are all double runs of one runtype. They are already accounted for in the previous 
    // run's rate.
    if (runs[i]==17167 || runs[i]==18963 || runs[i]==22200 || runs[i]==22823 || runs[i]==23073) continue;
    runNumber.push_back(runs[i]);

    Float_t fgsum = 0;
    Float_t bgsum = 0;

    TString fname = TString::Format("%s/Asymmetry/BinByBinComparison/UK_Erun%0.0f_anaCh%s.dat",
				    getenv("ANALYSIS_CODE"),runs[i],anaCh.Data());

    std::ifstream infile(fname.Data());

    std::string hold;
    Float_t timeFG, timeBG;

    infile >> hold >> hold >> hold >> timeFG
	   >> hold >> hold >> hold >> timeBG;
      

    int bin=binMin;

    Float_t en, fg, bg;

    while ( bin < binMax ) {
      
      infile >> en >> fg >> bg;
      //cout << " " << en << " " << fg << " " << bg << endl;
      fgsum+=fg;
      bgsum+=bg;
      bin++;

    }
  
    eastMax = fgsum>eastMax?fgsum:eastMax;
    
    fgrate[0].push_back(fgsum);
    bgrate[0].push_back(bgsum);

    fgrateErr[0].push_back(TMath::Sqrt(fgsum/timeFG));
    bgrateErr[0].push_back(TMath::Sqrt(bgsum/timeBG));

    infile.close();


    //West
    fgsum = 0;
    bgsum = 0;

    fname = TString::Format("%s/Asymmetry/BinByBinComparison/UK_Wrun%0.0f_anaCh%s.dat",
				    getenv("ANALYSIS_CODE"),runs[i],anaCh.Data());

    infile.open(fname.Data());

    infile >> hold >> hold >> hold >> timeFG
	   >> hold >> hold >> hold >> timeBG;

    bin=binMin;

    while ( bin <= binMax ) {
      
      infile >> en >> fg >> bg;
      fgsum+=fg;
      bgsum+=bg;
      bin++;

    }

    westMax = fgsum>westMax?fgsum:westMax;
    
    fgrate[1].push_back(fgsum);
    bgrate[1].push_back(bgsum);

    fgrateErr[1].push_back(TMath::Sqrt(fgsum/timeFG));
    bgrateErr[1].push_back(TMath::Sqrt(bgsum/timeBG));
  }

  //Plotting

  //East
  TCanvas *cEast = new TCanvas("cEast","East",1200, 600);
  //cEast->Divide(1,2);
  //cEast->cd(1);
  
  TGraphErrors *eastFG = new TGraphErrors(runNumber.size(),&runNumber[0], &fgrate[0][0],0,&fgrateErr[0][0]);
  eastFG->SetTitle(TString::Format("%s East %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));
  eastFG->GetYaxis()->SetTitle("Event Rate (Hz)");
  eastFG->GetXaxis()->SetTitle("Run Number");
  eastFG->SetMarkerStyle(20);
  eastFG->SetMarkerColor(kBlack);
  eastFG->SetMinimum(0.05);
  eastFG->SetMaximum(eastMax+10.);
  eastFG->Draw("AP");

  TF1 *eastFGfit = new TF1("eastFGfit","[0]",runNumber[0],runNumber[runNumber.size()-1]);
  eastFGfit->SetParameter(0,10.);
  eastFGfit->SetLineColor(kBlack);
  eastFGfit->SetLineStyle(2);
  eastFGfit->SetLineWidth(3);
  eastFG->Fit(eastFGfit,"R");

  //East->cd(2);

  TGraphErrors *eastBG = new TGraphErrors(runNumber.size(),&runNumber[0], &bgrate[0][0], 0, &bgrateErr[0][0]);
  eastBG->SetTitle(TString::Format("%s East %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));
  eastBG->GetYaxis()->SetTitle("Event Rate (Hz)");
  eastBG->GetXaxis()->SetTitle("Run Number");
  eastBG->SetMarkerStyle(26);
  eastBG->SetMarkerColor(kBlue);
  eastBG->SetLineColor(kBlue);
  eastBG->Draw("PSAME");

  TF1 *eastBGfit = new TF1("eastBGfit","[0]",runNumber[0],runNumber[runNumber.size()-1]);
  eastBGfit->SetParameter(0,0.1);
  eastBGfit->SetLineColor(kBlue);
  eastBGfit->SetLineStyle(6);
  eastBGfit->SetLineWidth(3);
  eastBG->Fit(eastBGfit,"R");

  TLegend *eastLeg = new TLegend(0.55,0.5,0.8,0.66);
  eastLeg->AddEntry(eastFG,"Foreground","ep");
  eastLeg->AddEntry(eastFGfit,TString::Format("Ave Rate = %0.2f Hz",eastFGfit->GetParameter(0)),"l");
  eastLeg->AddEntry(eastBG,"Background","ep");
  eastLeg->AddEntry(eastBGfit,TString::Format("Ave Rate = %0.2f Hz",eastBGfit->GetParameter(0)),"l");
  eastLeg->Draw("SAME");

  gPad->SetLogy();
  cEast->Update();


  //West
  TCanvas *cWest = new TCanvas("cWest","West",1200, 600);
  //cWest->Divide(1,2);
  //cWest->cd(1);
  
  TGraphErrors *westFG = new TGraphErrors(runNumber.size(),&runNumber[0], &fgrate[1][0],0,&fgrateErr[1][0]);
  westFG->SetTitle(TString::Format("%s West %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));  westFG->GetYaxis()->SetTitle("Event Rate (Hz)");
  westFG->GetXaxis()->SetTitle("Run Number");
  westFG->SetMarkerStyle(20);
  westFG->SetMarkerColor(kBlack);
  westFG->SetMinimum(0.05);
  westFG->SetMaximum(westMax+10.);
  westFG->Draw("AP");

  TF1 *westFGfit = new TF1("westFGfit","[0]",runNumber[0],runNumber[runNumber.size()-1]);
  westFGfit->SetParameter(0,10.);
  westFGfit->SetLineColor(kBlack);
  westFGfit->SetLineStyle(2);
  westFGfit->SetLineWidth(3);
  westFG->Fit(westFGfit,"R");

  //West->cd(2);

  TGraphErrors *westBG = new TGraphErrors(runNumber.size(),&runNumber[0], &bgrate[1][0], 0, &bgrateErr[1][0]);
  westBG->SetTitle(TString::Format("%s West %0.0f-%0.0f keV Integrated Rates",year.Data(),binMin*10., binMax*10.));
  westBG->GetYaxis()->SetTitle("Event Rate (Hz)");
  westBG->GetXaxis()->SetTitle("Run Number");
  westBG->SetMarkerStyle(26);
  westBG->SetMarkerColor(kBlue);
  westBG->SetLineColor(kBlue);
  westBG->Draw("PSAME");

  TF1 *westBGfit = new TF1("westBGfit","[0]",runNumber[0],runNumber[runNumber.size()-1]);
  westBGfit->SetParameter(0,0.1);
  westBGfit->SetLineColor(kBlue);
  westBGfit->SetLineStyle(6);
  westBGfit->SetLineWidth(3);
  westBG->Fit(westBGfit,"R");

  TLegend *westLeg = new TLegend(0.55,0.5,0.8,0.66);
  westLeg->AddEntry(westFG,"Foreground","ep");
  westLeg->AddEntry(westFGfit,TString::Format("Ave Rate = %0.2f Hz",westFGfit->GetParameter(0)),"l");
  westLeg->AddEntry(westBG,"Background","ep");
  westLeg->AddEntry(westBGfit,TString::Format("Ave Rate = %0.2f Hz",westBGfit->GetParameter(0)),"l");
  westLeg->Draw("SAME");

  gPad->SetLogy();
  cWest->Update();

  

  cEast->Print(TString::Format("integratedRatesByRun_%s_anaCh%s_%0.0f-%0.0f.pdf(",year.Data(),anaCh.Data(),binMin*10.,binMax*10.));
  cWest->Print(TString::Format("integratedRatesByRun_%s_anaCh%s_%0.0f-%0.0f.pdf)",year.Data(),anaCh.Data(),binMin*10.,binMax*10.));

}
