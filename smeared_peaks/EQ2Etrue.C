#include <iostream>
#include <string>
#include <vector>

using namespace std;


void EQ2Etrue (int runPeriod) {

  //define output file
  Char_t filename[100];
  ofstream outfile; 
  sprintf(filename,"EQ2Etrue_srcPeriod_%i.dat",runPeriod);
  outfile.open(filename);
  //store the data from smeared files
  ifstream infile;
  vector < vector <double> > EQ_smeared;
  EQ_smeared.resize(8, vector <double> (4,0.));

  sprintf(filename,"smear_%i.dat",runPeriod);
  infile.open(filename);

  for (int i=0; i<8; i++) {
    for (int j=0; j<4; j++) {
      infile >> EQ_smeared[i][j];
      cout << EQ_smeared[i][j] << " ";
    } 
    cout << endl;
  }

  //Store the true peak values (weighted average of K,L,M shells (N if in data sheet)
  Double_t Etrue[4];
  Etrue[0] = 131.5956; //Ce139 
  Etrue[1] = 368.4938; //Sn113
  Etrue[2] = 501.420; // Bi207 lower peak
  Etrue[3] = 993.789; // Bi207 upper peak

  vector <TCanvas*> can (8);
  //TCanvas *c1 = new TCanvas("c1");
  vector <TGraph*> graph (8);
  Char_t title[100], name[100], canvas[100];
  TF1 *func = new TF1("func", "[0]+[1]*x+[2]*x*x", 0., 1000.);

  for (int i=0; i<8; i++) {
    sprintf(canvas, "c%i", i); 
    can[i] = new TCanvas(canvas);
    can[i]->cd();
    sprintf(name, "graph%i",i);
    sprintf(title, "PMT %i",i);
    graph[i] = new TGraph(4, &EQ_smeared[i][0], Etrue);
    graph[i]->SetMarkerStyle(20);
    graph[i]->SetMarkerSize(0.75);
    graph[i]->GetXaxis()->SetTitle("Smeared E_{Q} [keV]");
    graph[i]->GetYaxis()->SetTitle("E_{true} [keV]");
    graph[i]->SetTitle(title);
    graph[i]->Draw("AP");
    func->SetParameters(0., 1., 0.);
    
    graph[i]->Fit("func", "R");
    outfile << func->GetParameter(0) << " " 
	    << func->GetParameter(1) << " " 
	    << func->GetParameter(2) << endl;
  }
}
