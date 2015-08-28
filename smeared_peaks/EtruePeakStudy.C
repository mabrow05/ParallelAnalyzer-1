//Code to study the effect of changing alpha on the peaks in the 
// Bi spectrum
{
Char_t temp[500];
TChain *chain = new TChain("Evts");
for (int i=0; i<100; i++) {
  sprintf(temp,"/extern/mabrow05/ucna/geant4work/output/10mil_src_events/Bi207/Evts_%i.root",i);
  chain->AddFile(temp);
 }
Int_t PID; 
Double_t KE;
chain->SetBranchAddress("KE",&KE);
chain->SetBranchAddress("PID",&PID);
TRandom3 *rand = new TRandom3(0);
Int_t nevents = chain->GetEntries();
 cout << "events = " << nevents << endl;

for (double alpha = 0.01; alpha<.1; alpha+=0.01) {
  TH1D *h = new TH1D("h", "h", 400, 0., 1600.);  
  for (Int_t evt=0; evt<nevents; evt++) {
    chain->GetEvent(evt);
    if (PID==11 && KE>200.) {
      h->Fill(rand->Gaus(KE,sqrt(KE/alpha)));
      //cout << "filled event " << evt << endl;
    }
  }
  h->Draw();
  Int_t maxBin = h->GetMaximumBin();
  Double_t peak = h->GetXaxis()->GetBinCenter(maxBin);
  Double_t maxBinContent = h->GetBinContent(maxBin);
  Double_t high = 0., low=0.;

  for (int i=maxBin; i<400; i++) {
        if (h->GetBinContent(i+1) < 0.5*maxBinContent) {
          high= h->GetXaxis()->GetBinCenter(i+1);
          break;
        }
      }
  for (int i=maxBin; i<400; i--) {
        if (h->GetBinContent(i-1) < 0.5*maxBinContent) {
          low= h->GetXaxis()->GetBinCenter(i-1);
          break;
        }
      }
      

  TF1 *func = new TF1("func", "gaus", low, high);
  //h->SetParameter(1,peak);
  h->Fit("func", "LRQ");
  cout << alpha << " \t" << func->GetParameter(1) << endl;
 }
}
