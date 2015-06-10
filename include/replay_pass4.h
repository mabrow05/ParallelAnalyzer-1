double pmt_pass4[8]; // pedestal-subtracted PMT QADC

double xE_pass4, yE_pass4, xW_pass4, yW_pass4; // East and West MWPC (x,y) positions

int PID_pass4; // Event PID

int type_pass4; // Type 0, 1, 2/3

int side_pass4; // Trigger side: 0 = East, 1 = West

int posError_pass4; // MWPC Position Reconstruction Error

 //The reconstructed energy in each individual PMT and their weight
struct PMT_Evis {
  double Evis0;
  double Evis1; 
  double Evis2; 
  double Evis3; 
  double Evis4;
  double Evis5; 
  double Evis6; 
  double Evis7; 
  double weight0; 
  double weight1; 
  double weight2; 
  double weight3; 
  double weight4; 
  double weight5; 
  double weight6; 
  double weight7;
} pmt_Evis;

double EvisW, EvisE, EvisTot; // Weighted Visible energy as seen in the PMTs for each side
