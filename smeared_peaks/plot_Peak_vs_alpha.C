{
#include <vector>

  vector <Double_t> alpha;
  vector <Double_t> peakDiff;

  ifstream infile;
  infile.open("alphaTrack.dat");
  
  Double_t alph=0., peak=0.;
  while (infile >> alph >> peak) {
    alpha.push_back(alph);
    peakDiff.push_back(peak);
  }
  
  UInt_t num = alpha.size();

  TGraph *gr = new TGraph(num-5,&alpha[5],&peakDiff[5]);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.75);
  gr->SetMarkerColor(kRed);
  gr->GetXaxis()->SetTitle("#alpha");
  gr->GetYaxis()->SetTitle("E_{KLM}^{weighted} - E_{true}^{smeared}");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->SetTitle("[E_{KLM}^{weighted} - E_{true}^{smeared}] vs. #alpha");
  gr->Draw("AP");
}
