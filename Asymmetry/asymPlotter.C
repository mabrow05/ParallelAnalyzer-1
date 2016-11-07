{

  std::vector <Double_t> oct_NCSU;
  std::vector <Double_t> asym_NCSU;
  std::vector <Double_t> err_NCSU;

  std::vector <Double_t> oct_UK;
  std::vector <Double_t> asym_UK;
  std::vector <Double_t> err_UK;
  
  ifstream infile("octval.dat");

  string hold[8];
  double octet, asym, err, pull;
  
  infile >> hold[0] >> hold[1] >> hold[2] >> hold[3] >> hold[4] >> hold[5] >> hold[6] >> hold[7];

  while (infile >> octet >> asym >> err >> pull) {
    oct_NCSU.push_back(octet);
    asym_NCSU.push_back(asym);
    err_NCSU.push_back(err);
  }

  infile.close();

  infile.open("octvalUK.dat");

  infile >> hold[0] >> hold[1] >> hold[2] >> hold[3];
    
  while (infile >> octet >> asym >> err >> pull) {
    oct_UK.push_back(octet);
    asym_UK.push_back(asym);
    err_UK.push_back(err);
  }

  std::vector <Double_t> diff;

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);

  c1->Divide(1,2);

  c1->cd(1);

  TGraphErrors *UK = new TGraphErrors(oct_UK.size(), &oct_UK[0], &asym_UK[0], 0, &err_UK[0]);
  UK->SetMarkerColor(kBlue);
  UK->SetMarkerStyle(20);
  UK->SetLineWidth(2);
  UK->SetLineColor(kBlue);
  UK->Draw("AP");

  TF1 *fit = new TF1("fit","[0]",oct_UK[0], oct_UK[oct_UK.size()-1]);
  //fit->SetLineColor(kRed);
  fit->SetLineWidth(3);
  fit->SetParameter(0,0.05);
  
  UK->Fit("fit","R");

  cout << "UK\t" << fit->GetParameter(0) << "\t+/-\t" << fit->GetParError(0) << std::endl;


  c1->cd(2);

  TGraphErrors *NCSU = new TGraphErrors(oct_NCSU.size(), &oct_NCSU[0], &asym_NCSU[0], 0, &err_NCSU[0]);
  NCSU->SetMarkerColor(kRed);
  NCSU->SetMarkerStyle(20);
  NCSU->SetLineWidth(2);
  NCSU->SetLineColor(kRed);
  NCSU->Draw("AP");

  TF1 *fit2 = new TF1("fit2","[0]",oct_NCSU[0], oct_NCSU[oct_NCSU.size()-1]);
  //fit->SetLineColor(kRed);
  fit2->SetLineWidth(3);
  fit2->SetParameter(0,0.05);
  
  NCSU->Fit("fit2","R");

  cout << "NCSU\t" << fit2->GetParameter(0) << "\t+/-\t" << fit2->GetParError(0) << std::endl;

  c1->Update();
}
  
		  
