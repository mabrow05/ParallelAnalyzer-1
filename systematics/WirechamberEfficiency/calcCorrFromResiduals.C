

{

  gStyle->SetTitleSize(0.06,"t");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetTitleSize(0.06,"y");
  //gStyle->SetTitleOffset(0.85,"y");
  //gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  //gStyle->SetPadLeftMargin(0.12);
  gStyle->SetLabelSize(0.035,"xyz");

  //gStyle->SetOptFit(1111);
  gStyle->SetTitleX(0.5);
  //gStyle->SetStatX(0.75);
  //gStyle->SetStatY(0.80);

  TString year = "2012-2013";
  Int_t emin = 220;
  Int_t emax = 670;
    
  ifstream infile(TString::Format("EfficiencyDistr_%s_%i-%i.txt",year.Data(),emin,emax));

  std::vector <Double_t> res(0);
  
  Double_t hold;
  
  while (infile >> hold) res.push_back(hold);

  
  //calculate RMS, mean, sigma
  TH1D *h = new TH1D("h","Efficiency Correction",20,-0.01,0.01);
  h->GetXaxis()->SetTitle("#DeltaA/A");
  h->GetXaxis()->CenterTitle();
  h->SetLineWidth(3);
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

  h->Draw();
}

  
