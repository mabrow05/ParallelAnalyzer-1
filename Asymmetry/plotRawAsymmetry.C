{
  //#include <vector>
  gStyle->SetOptFit(1111);

  ifstream infile("rawAsymmetryByOctet_0-59_AnaCh7_50mm.dat");
  char in1[10], in2[10];
  vector <double> octet;
  vector <double> asym; 
  vector <double> errorY;
  vector <double> errorX;
  double Octet = 0.;

  while (infile >> in1 >> in2) {
    if (in1[0]=='0') {
      octet.push_back(Octet);
      asym.push_back(atof(in1));
      errorY.push_back(atof(in2));  
      errorX.push_back(0.);
    }
    Octet = Octet+1.;
  }

  for (int i=0; i<octet.size(); i++) {
    cout << octet[i] << " " << asym[i] << " " << errorY[i] << endl;
  }

  TGraphErrors *g = new TGraphErrors(octet.size(), &octet[0],&asym[0],&errorX[0], &errorY[0]);
  g->SetTitle("Raw Measured Asymmetry");
  g->SetMarkerStyle(20);
  g->SetLineWidth(2);
  g->GetXaxis()->SetLimits(-2., octet[octet.size()-1]+2.);
  
  TF1 *fit = new TF1("fit","[0]",octet[0], octet[octet.size()-1]);
  fit->SetLineColor(kRed);
  fit->SetLineWidth(3);

  g->Fit("fit","R");

  
  g->Draw("AP");

} 
