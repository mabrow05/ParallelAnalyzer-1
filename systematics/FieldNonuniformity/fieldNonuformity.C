

void fieldNonuniformity(Double_t emin, Double_t emax) {

  std::vector <Double_t> A;

  TString basepath = "/extern/mabrow05/ucna/geant4work/output/bigSim/badField_2011-2012/";

  for (int i=0; i<10; ++i) {

    for (int j=0; j<1000; ++j) {

      TFile *f = new TFile(TString::Format("%s/B%i/analyzed_%i.root",basepath.Data(),i,j),"READ");
      TTree *t = (TTree*)f->Get("anaTree");

      Double_t rE = t->GetEntries(TString::Format("EdepSD[5]>0. && hitCountSD[15]==0 && primKE>%f && primKE<%f",emin, emax));
      Double_t rW = t->GetEntries(TString::Format("EdepSD[15]>0. && hitCountSD[5]==0 && primKE>%f && primKE<%f",emin, emax));

      A.push_back( (rE-rW)/(rE+rW) );

      delete f;

    }

  }

  TCanvas *c1 = new TCanvas("c1","c1");

  TH1D *h = new TH1D("h",TString::Format("Field Non-uniformity A_{raw} for %0.0f-%0.0f keV",emin,emax),100, -1.,1.);

  for (UInt_t b=0; b<A.size(); ++b) {
    h->Fill(A[i]);
  }

  f->Fit("gaus");



}
