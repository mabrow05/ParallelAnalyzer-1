

{

  TString year = "2011-2012";
  Int_t emin = 220;
  Int_t emax = 750;
    
  ifstream infile(TString::Format("fieldDipDistr_%s_%i-%i.txt",year.Data(),emin,emax));

  std::vector <Double_t> res(0);
  
  Double_t hold;
  
  while (infile >> hold) res.push_back(hold);

  
  //calculate RMS, mean, sigma
  TH1D *h = new TH1D("h","Field Non-Uniformity Correction",15,-0.025,0.025);
  Double_t RMS = 0., mean=0., sigma=0.;
  
  for (UInt_t i=0;i<res.size();++i) { 
    RMS+=(res[i]*res[i]);
    mean+=res[i];
    h->Fill(res[i]);
  }
  mean /= (Double_t)res.size();
  RMS = TMath::Sqrt(RMS/(Double_t)res.size());

  for (UInt_t i=0;i<res.size();++i) sigma+=TMath::Power(res[i]-mean,2.);
  sigma = TMath::Sqrt(sigma/((Double_t)res.size()-1.));


  std::cout << "RMS = " << RMS << std::endl;
  std::cout << "Mean = " << mean << std::endl;
  std::cout << "Sigma = " << sigma << std::endl;


  std::cout << std::endl;

  h->Fit("gaus");
  h->Draw();
}

  
