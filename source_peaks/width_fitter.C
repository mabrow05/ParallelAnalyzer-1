
std::vector <Int_t> getPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"%s/residuals/PMT_runQuality_master.dat",getenv("ANALYSIS_CODE")); 
  ifstream pmt;
  std::cout << temp << std::endl;
  pmt.open(temp);
  Int_t run_hold;
  while (pmt >> run_hold >> pmtQuality[0] >> pmtQuality[1] >> pmtQuality[2]
	 >> pmtQuality[3] >> pmtQuality[4] >> pmtQuality[5]
	 >> pmtQuality[6] >> pmtQuality[7]) {
    if (run_hold==runNumber) break;
    if (pmt.fail()) break;
  }
  pmt.close();
  if (run_hold!=runNumber) {
    cout << "Run not found in PMT quality file!" << endl;
    exit(0);
  }
  return pmtQuality;
};


void width_fitter(Int_t calPeriod)
{

  //Read in sim and data widths

  std::vector < std::vector < Double_t > > simWidths(8,std::vector <Double_t> (500, 0.));
  std::vector < std::vector < Double_t > > dataWidths(8,std::vector <Double_t> (500,0.));
  std::vector < Double_t > old_k(8,0.);
  std::vector < Double_t > new_k(8,0.);
  std::vector < Double_t > slope(8,0.);

  Char_t tempfile[200];
  sprintf(tempfile, "%s/simulation_comparison/nPE_per_keV/nPE_per_keV_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);
  ifstream kfilein(tempfile);
 
  Int_t k=0;
  std::cout << "OLD k VALUES:\n";
  while (kfilein >> old_k[k]) {
    std::cout << old_k[k] << std::endl;
    k++;
  } 
  kfilein.close();


  sprintf(tempfile,"%s/residuals/source_runs_EvisWidth_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);
  //sprintf(tempfile,"../residuals/source_runs_EvisWidth_RunPeriod_%i.dat",calPeriod);
  
  ifstream dataFile(tempfile); 
  if (dataFile.is_open()) cout << tempfile << endl;
  
  sprintf(tempfile,"%s/residuals/SIM_source_runs_EvisWidth_RunPeriod_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);
  //sprintf(tempfile,"../residuals/SIM_source_runs_EvisWidth_runPeriod_%i.dat",calPeriod);
  
  ifstream simFile(tempfile);
  
  std::string srcNameData, srcNameSim;
  Int_t Run, simRun;

  Int_t i=0;
  std::vector < Double_t > dataWidths_hold(8,0.);
  std::vector < Double_t > simWidths_hold(8,0.);
  std::vector <Int_t> num(8,0);
  
  while (dataFile >> Run >> srcNameData
	 >> dataWidths_hold[0] >> dataWidths_hold[1] >> dataWidths_hold[2] >> dataWidths_hold[3] >>
	 dataWidths_hold[4] >> dataWidths_hold[5] >> dataWidths_hold[6] >> dataWidths_hold[7]) {

    //cout << Run << " " << srcNameData << endl;

    simFile >> simRun >> srcNameSim;
    if (simRun==Run && srcNameSim==srcNameData) {
      simFile >> simWidths_hold[0] >> simWidths_hold[1] >> simWidths_hold[2] >> simWidths_hold[3]
	      >> simWidths_hold[4] >> simWidths_hold[5] >> simWidths_hold[6] >> simWidths_hold[7];

      std::vector <Int_t> pmtQuality = getPMTQuality(Run);

      for (Int_t p=0; p<8; p++) {
	if (pmtQuality[p]) { 
	  
	  if (TMath::Abs(simWidths_hold[p]-dataWidths_hold[p])/simWidths_hold[p] < 0.2) { //Checking for outliers
	    dataWidths[p][num[p]] = dataWidths_hold[p];
	    simWidths[p][num[p]] = simWidths_hold[p];
	    if (srcNameData!="Cd") num[p]++; //Put peaks to exclude here
	    //num[p]++; 
	  }
	}
      }

      cout << simRun << " " << srcNameSim << endl;

      
    }
    else {
      cout << "Data and sim files don't match. Rerun MakeSourceCalibrationFiles!\n"; 
      exit(0);
    }
    if (dataFile.fail()) break;
    
  }
  dataFile.close();
  simFile.close();

  cout << i << endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 1400, 1000);

  TPad *p0 = new TPad("p0","East 1", 0.0, 0.5, 0.25, 1.0);
  TPad *p1 = new TPad("p1","East 2", 0.25, 0.5, 0.5, 1.0);
  TPad *p2 = new TPad("p2","East 3", 0.5, 0.5, 0.75, 1.0);
  TPad *p3 = new TPad("p3","East 4", 0.75, 0.5, 1.0, 1.0);
  TPad *p4 = new TPad("p4","West 1", 0.0, 0., 0.25, .5);
  TPad *p5 = new TPad("p5","West 2", 0.25, 0., 0.5, .5);
  TPad *p6 = new TPad("p6","West 3", 0.5, 0., 0.75, .5);
  TPad *p7 = new TPad("p7","West 4", 0.75, 0., 1.0, .5);
  p0->Draw();
  p1->Draw();
  p2->Draw();
  p3->Draw();
  p4->Draw(); 
  p5->Draw();
  p6->Draw(); 
  p7->Draw(); 


  TF1 *f1 = new TF1("f1","[0]*x",0., 170.);
  f1->SetParameter(0,1.);
  f1->SetParLimits(0,0.5,1.5);
  gStyle->SetOptFit();
  
  TString status = " ";
  Int_t nAttempts = 5;
  Int_t attempt = 0;

  
  p0->cd();

  if (num[0]>0) {

    TGraph *pmt0 = new TGraph(num[0], &simWidths[0][0], &dataWidths[0][0]);
    pmt0->SetMarkerColor(1);
    pmt0->SetLineColor(1);
    pmt0->SetMarkerStyle(20);
    pmt0->SetMarkerSize(0.75);
    pmt0->GetXaxis()->SetLimits(0.0,170.);
    pmt0->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt0->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt0->SetMinimum(0.0);
    pmt0->SetMaximum(170.);
    pmt0->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt0->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[0] = 1.;
    else slope[0] = f1->GetParameter(0);
    
  }
  
  else slope[0]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);

  p1->cd();

  if (num[1]>0) {
    TGraph *pmt1 = new TGraph(num[1], &simWidths[1][0], &dataWidths[1][0]);
    pmt1->SetMarkerColor(1);
    pmt1->SetLineColor(1);
    pmt1->SetMarkerStyle(20);
    pmt1->SetMarkerSize(0.75);
    pmt1->GetXaxis()->SetLimits(0.0,170.);
    pmt1->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt1->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt1->SetMinimum(0.0);
    pmt1->SetMaximum(170.);
    pmt1->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt1->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[1] = 1.;
    else slope[1] = f1->GetParameter(0);
  }

  else slope[1]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);

  p2->cd();


  if (num[2]>0) {
    TGraph *pmt2 = new TGraph(num[2], &simWidths[2][0], &dataWidths[2][0]);
    pmt2->SetMarkerColor(1);
    pmt2->SetLineColor(1);
    pmt2->SetMarkerStyle(20);
    pmt2->SetMarkerSize(0.75);
    pmt2->GetXaxis()->SetLimits(0.0,170.);
    pmt2->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt2->GetYaxis()->SetTitle("Actual Width (keV)");  
    pmt2->SetMinimum(0.0);
    pmt2->SetMaximum(170.);
    pmt2->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt2->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[2] = 1.;
    else slope[2] = f1->GetParameter(0);
  }

  else slope[2]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);


  p3->cd();

  if (num[3]>0) {

    TGraph *pmt3 = new TGraph(num[3], &simWidths[3][0], &dataWidths[3][0]);
    pmt3->SetMarkerColor(1);
    pmt3->SetLineColor(1);
    pmt3->SetMarkerStyle(20);
    pmt3->SetMarkerSize(0.75);
    pmt3->GetXaxis()->SetLimits(0.0,170.);
    pmt3->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt3->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt3->SetMinimum(0.0);
    pmt3->SetMaximum(170.);
    pmt3->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt3->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[3] = 1.;
    else slope[3] = f1->GetParameter(0);
  }
  
  else slope[3]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);


  p4->cd();

  if (num[4]>0) {
    
    TGraph *pmt4 = new TGraph(num[4], &simWidths[4][0], &dataWidths[4][0]);
    pmt4->SetMarkerColor(1);
    pmt4->SetLineColor(1);
    pmt4->SetMarkerStyle(20);
    pmt4->SetMarkerSize(0.75);
    pmt4->GetXaxis()->SetLimits(0.0,170.);
    pmt4->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt4->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt4->SetMinimum(0.0);
    pmt4->SetMaximum(170.);
    pmt4->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt4->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[4] = 1.;
    else slope[4] = f1->GetParameter(0);
  }

  else slope[4]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);


  p5->cd();

  if (num[5]>0) {

    TGraph *pmt5 = new TGraph(num[5], &simWidths[5][0], &dataWidths[5][0]);
    pmt5->SetMarkerColor(1);
    pmt5->SetLineColor(1);
    pmt5->SetMarkerStyle(20);
    pmt5->SetMarkerSize(0.75);
    pmt5->GetXaxis()->SetLimits(0.0,170.);
    pmt5->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt5->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt5->SetMinimum(0.0);
    pmt5->SetMaximum(170.);
    pmt5->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt5->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[5] = 1.;
    else slope[5] = f1->GetParameter(0);
  }

  else slope[5]=1.;

  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);


  p6->cd();

  if (num[6]>0) {

    TGraph *pmt6 = new TGraph(num[6], &simWidths[6][0], &dataWidths[6][0]);
    pmt6->SetMarkerColor(1);
    pmt6->SetLineColor(1);
    pmt6->SetMarkerStyle(20);
    pmt6->SetMarkerSize(0.75);
    pmt6->GetXaxis()->SetLimits(0.0,170.);
    pmt6->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt6->GetYaxis()->SetTitle("Actual Width (keV)");  
    pmt6->SetMinimum(0.0);
    pmt6->SetMaximum(170.);
    pmt6->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
    pmt6->Fit("f1");
    status = gMinuit->fCstatu;
    attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[6] = 1.;
    else slope[6] = f1->GetParameter(0);
  }
  
  else slope[6]=1.;
  
  status = " ";
  attempt = 0;
  f1->SetParameter(0,1.);


  p7->cd();

  if (num[7]>0) {
    
    TGraph *pmt7 = new TGraph(num[7], &simWidths[7][0], &dataWidths[7][0]);
    pmt7->SetMarkerColor(1);
    pmt7->SetLineColor(1);
    pmt7->SetMarkerStyle(20);
    pmt7->SetMarkerSize(0.75);
    pmt7->GetXaxis()->SetLimits(0.0,170.);
    pmt7->GetXaxis()->SetTitle("Simulated Width (keV)");
    pmt7->GetYaxis()->SetTitle("Actual Width (keV)");
    pmt7->SetMinimum(0.0);
    pmt7->SetMaximum(170.);
    pmt7->Draw("AP");
    
    while (status!=TString("CONVERGED ") && attempt!=nAttempts) {    
      pmt7->Fit("f1");
      status = gMinuit->fCstatu;
      attempt++;
    }
    if (status!=TString("CONVERGED ")) slope[7] = 1.;
    else slope[7] = f1->GetParameter(0);
  }
  
  else slope[7]=1.;


  TString pdffile = TString::Format("%s/simulation_comparison/nPE_per_keV/width_comp_%i.pdf",getenv("ANALYSIS_CODE"),calPeriod);
  c1->Print(pdffile);


  //Calculate new k 

  sprintf(tempfile,"%s/simulation_comparison/nPE_per_keV/nPE_per_keV_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);
  ofstream new_k_fileout(tempfile);

  sprintf(tempfile,"%s/simulation_comparison/nPE_per_keV/prev_nPE_per_keV_%i.dat",getenv("ANALYSIS_CODE"),calPeriod);
  ofstream old_k_fileout(tempfile);

  for (UInt_t j=0; j<8; j++) {
    if (j!=0) {
      new_k_fileout << std::endl;
      old_k_fileout << std::endl;
    }
    std::cout << old_k[j] << " " << slope[j] << std::endl;

    new_k[j] = old_k[j]/(slope[j]*slope[j]);
    new_k_fileout << new_k[j];
    old_k_fileout << old_k[j];
  }
  old_k_fileout.close();
  new_k_fileout.close();
  
}