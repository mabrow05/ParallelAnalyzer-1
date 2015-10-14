#include <vector>
#include <cstdlib>

void plot_resid_vs_pos(string SRC) {
  
  UInt_t calPeriodLow = 1;
  UInt_t calPeriodHigh = 12;
  Char_t tempfile[500];

  string SRCE = SRC + "_East";
  string SRCW = SRC + "_West";
  
  //Setting vectors for side to hold the source name, runNumber, residual, x, and y positions

  vector <string> sourceE;
  vector <UInt_t> runE;
  vector <Double_t> residualE;
  vector <Double_t> xposE;
  vector <Double_t> yposE;
  vector <string> sourceW;
  vector <UInt_t> runW;
  vector <Double_t> residualW;
  vector <Double_t> xposW;
  vector <Double_t> yposW;

  // places to hold what we're reading in...
  string src;
  UInt_t rn;
  Double_t res;
  Double_t xE,yE,stddevE;
  Double_t xW,yW,stddevW;

  ifstream residualFile, positionFile, sourceFile;
  //East Residuals
  sprintf(tempfile, "../residuals/residuals_global_EvisPMTbyPMT_East_periods_%i-%i.dat", calPeriodLow, calPeriodHigh); 
  residualFile.open(tempfile);

  while (residualFile >> src >> rn >> res) {
    if (src==SRCE) {
      sourceE.push_back(src);
      runE.push_back(rn);
      residualE.push_back(res);
    }
  }
  residualFile.close();

  //West Residuals
  sprintf(tempfile, "../residuals/residuals_global_EvisPMTbyPMT_West_periods_%i-%i.dat", calPeriodLow, calPeriodHigh); 
  residualFile.open(tempfile);
  
  while (residualFile >> src >> rn >> res) {
    if (src==SRCW) {
      sourceW.push_back(src);
      runW.push_back(rn);
      residualW.push_back(res);
    }
  }
  residualFile.close();
  
  if (runE.size()!=runW.size()) {
    cout << "not same number of points on east and west" << endl;
    exit(0);
  }

  xposE.resize(runE.size(),0.);
  yposE.resize(runE.size(),0.); 
  xposW.resize(runW.size(),0.);
  yposW.resize(runW.size(),0.);

  UInt_t nsources = 0.;
  //Read in source positions
  for (UInt_t i=0; i<runE.size(); i++) {
    //determine number of sources in this run
    vector <string> srchold;
    sprintf(tempfile,"%s/source_list_%i.dat", getenv("SOURCE_LIST"), runE[i]);
    sourceFile.open(tempfile);
    sourceFile >> nsources;
    for (UInt_t n = 0; n<nsources; n++) {
      sourceFile >> src;
      srchold.push_back(src);
    }
    sourceFile.close();

    //read in source positions
    sprintf(tempfile,"%s/source_positions_%i.dat", getenv("SOURCE_POSITIONS"), runE[i]);
    positionFile.open(tempfile);
    
    if (nsources==1) {
      positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
    }
    if (nsources==2) {
      if (sourceE[i]=="Ce_East") {
	positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
      }
      if (sourceE[i]=="Sn_East") {
	if (srchold[0]=="Sn") {
	  positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
	}
	else if (srchold[1]=="Sn") {
	  positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
	  positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
	}
      }
      if (sourceE[i]=="Bi1_East" || sourceE[i]=="Bi2_East") {
	positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
	positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
      }
    }
    if (nsources==3) {
      if (sourceE[i]=="Ce_East") {
	positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
      }
      if (sourceE[i]=="Sn_East") {
	  positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
	  positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
      }
      if (sourceE[i]=="Bi1_East" || sourceE[i]=="Bi2_East") {
	positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
	positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
	positionFile >> xposE[i] >> yposE[i] >> stddevE >> xposW[i] >> yposW[i] >> stddevW;
      }
    }
    positionFile.close();
    cout << runE[i] << sourceE[i] << xposE[i] << endl;
  }

  //plot them
  const Int_t nn = 2;
  Double_t x[nn] = {-50., 50.};
  Double_t y[nn] = {0.0, 0.0};

  TGraph *gr0 = new TGraph(nn,x,y);
  //gr0->Draw("Same");
  gr0->SetLineWidth(2);
  gr0->SetLineColor(1);
  gr0->SetLineStyle(2);
  
  Char_t name[500];
  TCanvas *c1 = new TCanvas("c1", "c1", 1400, 1000);
  c1->Divide(1,2);
  c1->cd(1);
  TGraph *eastX = new TGraph(runE.size(),&xposE[0],&residualE[0]);
  eastX->SetMarkerStyle(20);
  eastX->SetMarkerColor(kBlue);
  sprintf(name, "East %s: x-position",SRC.c_str());
  eastX->SetTitle(name);
  eastX->GetXaxis()->SetTitle("x (mm)");
  eastX->GetYaxis()->SetTitle("residual (keV)");
  eastX->GetYaxis()->SetRangeUser(-50., 50.);
  //eastX->GetXaxis()->SetRangeUser(-50., 50.);
  eastX->GetXaxis()->SetLimits(-50., 50.);
  eastX->Draw("AP");
  gr0->Draw("SAME");
  
  //TCanvas *c2 = new TCanvas("c2", "c2", 1400, 500);
  c1->cd(2);
  TGraph *westX = new TGraph(runW.size(),&xposW[0],&residualW[0]);
  westX->SetMarkerStyle(20);
  westX->SetMarkerColor(kBlue);
  sprintf(name, "West %s: x-position",SRC.c_str());
  westX->SetTitle(name);
  westX->GetXaxis()->SetTitle("x (mm)");
  westX->GetYaxis()->SetRangeUser(-50., 50.);
  westX->GetXaxis()->SetLimits(-50., 50.);
  westX->GetYaxis()->SetTitle("residual (keV)");
  westX->Draw("AP");
  gr0->Draw("SAME");


  //Y-variation
  gr0->GetXaxis()->SetLimits(-20.,20.);

  TCanvas *c2 = new TCanvas("c2", "c2", 1400, 1000);
  c2->Divide(1,2);
  c2->cd(1);
  TGraph *eastY = new TGraph(runE.size(),&yposE[0],&residualE[0]);
  eastY->SetMarkerStyle(20);
  eastY->SetMarkerColor(kBlue);
  sprintf(name, "East %s: y-position",SRC.c_str());
  eastY->SetTitle(name);
  eastY->GetXaxis()->SetTitle("y (mm)");
  eastY->GetYaxis()->SetTitle("residual (keV)");
  eastY->GetYaxis()->SetRangeUser(-50., 50.);
  //eastY->GetXaxis()->SetRangeUser(-50., 50.);
  eastY->GetXaxis()->SetLimits(-20., 20.);
  eastY->Draw("AP");
  gr0->Draw("SAME");
  
  //TCanvas *c2 = new TCanvas("c2", "c2", 1400, 500);
  c2->cd(2);
  TGraph *westY = new TGraph(runW.size(),&yposW[0],&residualW[0]);
  westY->SetMarkerStyle(20);
  westY->SetMarkerColor(kBlue);
  sprintf(name, "West %s: y-position",SRC.c_str());
  westY->SetTitle(name);
  westY->GetXaxis()->SetTitle("y (mm)");
  westY->GetYaxis()->SetRangeUser(-50., 50.);
  westY->GetXaxis()->SetLimits(-20., 20.);
  westY->GetYaxis()->SetTitle("residual (keV)");
  westY->Draw("AP");
  gr0->Draw("SAME");

  TCanvas *c3 = new TCanvas("c3", "c3", 1400, 1000);
  c3->Divide(1,2);
  c3->cd(1);
  TGraph2D *east = new TGraph2D(runE.size(),&xposE[0],&yposE[0], &residualE[0]);
  east->SetMarkerStyle(20);
  east->SetMarkerColor(kBlue);
  sprintf(name, "East %s",SRC.c_str());
  east->SetTitle(name);
  east->GetXaxis()->SetTitle("x (mm)");
  east->GetYaxis()->SetTitle("y (mm)");
  east->GetZaxis()->SetTitle("residual (keV)");
  east->GetXaxis()->SetTitleOffset(1.8);
  east->GetYaxis()->SetTitleOffset(1.8);
  east->GetXaxis()->SetLimits(-50., 50.);
  east->GetYaxis()->SetLimits(-50., 50.);
  east->GetZaxis()->SetRangeUser(-50., 50.);
  //eastY->GetXaxis()->SetRangeUser(-50., 50.);
  //east->GetXaxis()->SetLimits(-20., 20.);
  east->GetXaxis()->CenterTitle();
  east->GetYaxis()->CenterTitle();
  east->GetZaxis()->CenterTitle();
  east->Draw("PCOLZ");
  //gr0->Draw("SAME");
  
  //TCanvas *c2 = new TCanvas("c2", "c2", 1400, 500);
  c3->cd(2);
  TGraph2D *west = new TGraph2D(runW.size(),&xposW[0],&yposW[0],&residualW[0]);
  west->SetMarkerStyle(20);
  west->SetMarkerColor(kBlue);
  sprintf(name, "West %s",SRC.c_str());
  west->SetTitle(name);
  west->GetXaxis()->SetTitle("x (mm)");
  west->GetYaxis()->SetTitle("y (mm))");
  west->GetZaxis()->SetTitle("residual (keV)");
  west->GetXaxis()->SetTitleOffset(1.8);
  west->GetYaxis()->SetTitleOffset(1.8);
  west->GetXaxis()->SetLimits(-50., 50.);
  west->GetYaxis()->SetLimits(-50., 50.);
  west->GetZaxis()->SetRangeUser(-50., 50.);
  west->GetXaxis()->CenterTitle();
  west->GetYaxis()->CenterTitle();
  west->GetZaxis()->CenterTitle();
  west->Draw("PCOLZ");
  //gr0->Draw("SAME");

  //TCanvas *c5 = new TCanvas("c5");
  //TH2D *h = (TH2D*)east->Project("x:y");
  //h->Draw("colz");
}
