#include <vector>
#include <iostream>
#include <cstdlib>

vector <vector <double> > returnSourcePosition (Int_t runNumber, string src) {
  Char_t temp[500];
  sprintf(temp,"%s/source_list_%i.dat",getenv("SOURCE_LIST"),runNumber);
  ifstream file(temp);
  cout << src << endl;
  int num = 0;
  file >> num;
  cout << num << endl;
  int srcNum = 0;
  string src_check;
  for (int i=0; i<num;srcNum++,i++) {
    file >> src_check;
    cout << src_check << endl;
    if (src_check==src) break;   
  }
  cout << "The source Number is: " << srcNum << endl;
  if (srcNum==num) {
    cout << "Didn't find source in that run\n"; exit(0);
  }
  file.close();
  
  sprintf(temp,"%s/source_positions_%i.dat",getenv("SOURCE_POSITIONS"),runNumber);
  file.open(temp);
  
  vector < vector < double > > srcPos;
  srcPos.resize(2,vector <double> (3,0.));
  
  for (int i=0; i<srcNum+1; i++) {
    for (int j=0; j<2; j++) {
      for (int jj=0; jj<3; jj++) {
	file >> srcPos[j][jj];
      }
    }
  }
  return srcPos;
}

void quick_plot_spectra(Int_t runNumber, string src) {
  Char_t temp[500];
  sprintf(temp,"/extern/UCNA/replay_pass4_MB/replay_pass4_%i.root",runNumber);
  TFile *data = new TFile(temp,"READ");
  TTree *TData = (TTree*)(data->Get("pass4"));
  sprintf(temp,"/extern/UCNA/reverse_cal_sim_MB/revCalSim_%i.root",runNumber);
  TFile *sim = new TFile(temp,"READ");
  TTree *TSim = (TTree*)(sim->Get("revCalSim"));
  cout << "Opened Files\n";

  vector < vector <double> > srcPos = returnSourcePosition(runNumber,src);
  
  for (int j=0; j<2; j++) {
      for (int jj=0; jj<3; jj++) {
	cout << srcPos[j][jj] << " ";
      }
      cout << endl;
    }
  
  cout << "Got source position\n";
  TH1D *data_E_all = new TH1D("data_E_all", "East All Evts", 400, 0., 1200.);
  TH1D *sim_E_all = new TH1D("sim_E_all", "East All Evts", 400, 0., 1200.);
  sim_E_all->SetLineColor(kRed);
  TH1D *data_W_all = new TH1D("data_W_all", "West All Evts", 400, 0., 1200.);
  TH1D *sim_W_all = new TH1D("sim_W_all", "West All Evts", 400, 0., 1200.);
  sim_W_all->SetLineColor(kRed);
  Char_t cuts[1000];
  
  //Type0

  TCanvas *c1 = new TCanvas("c1");

  sprintf(cuts,"PID_pass4==1 && type_pass4<3 && side_pass4==0 && EvisE>0. && xE_pass4>(%f-2.*%f) && xE_pass4<(%f+2.*%f) && yE_pass4>(%f-2.*%f) && yE_pass4<(%f+2.*%f)",srcPos[0][0],srcPos[0][2],srcPos[0][0],srcPos[0][2],srcPos[0][1],srcPos[0][2],srcPos[0][1],srcPos[0][2]);  
  TData->Draw("EvisE>>data_E_all", cuts);
  data_E_all->GetXaxis()->SetRangeUser(0.,250.);
  Double_t integralE = data_E_all->Integral();
  data_E_all->GetXaxis()->SetRangeUser(0.,1200.);

  sprintf(cuts,"PID==1 && side==0 && type<3 && EvisE>0.");
  TSim->Draw("EvisE>>sim_E_all",cuts,"SAME");
  sim_E_all->GetXaxis()->SetRangeUser(0.,250.);
  Double_t scaleValE = integralE/sim_E_all->Integral();
  sim_E_all->GetXaxis()->SetRangeUser(0.,1200.);
  sim_E_all->Scale(scaleValE);
  sim_E_all->Draw("SAME");

  TCanvas *c2 = new TCanvas("c2");
  
  sprintf(cuts,"PID_pass4==1 && type_pass4<3 && side_pass4==1 && EvisW>0. && xW_pass4>(%f-2.*%f) && xW_pass4<(%f+2.*%f) && yW_pass4>(%f-2.*%f) && yW_pass4<(%f+2.*%f)",srcPos[1][0],srcPos[1][2],srcPos[1][0],srcPos[1][2],srcPos[1][1],srcPos[1][2],srcPos[1][1],srcPos[1][2]);  
  TData->Draw("EvisW>>data_W_all", cuts);
  data_W_all->GetXaxis()->SetRangeUser(0.,250.);
  Double_t integralW = data_W_all->Integral();
  data_W_all->GetXaxis()->SetRangeUser(0.,1200.);

  sprintf(cuts,"PID==1 && side==1 && type<3 && EvisW>0.");
  TSim->Draw("EvisW>>sim_W_all",cuts,"SAME");
  sim_W_all->GetXaxis()->SetRangeUser(0.,250.);
  Double_t scaleValW = integralW/sim_W_all->Integral();
  sim_W_all->GetXaxis()->SetRangeUser(0.,1200.);
  sim_W_all->Scale(scaleValW);
  sim_W_all->Draw("SAME");

  
}
