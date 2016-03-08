
int PID, type, side; // basic analysis tags

double Eprim, AsymWeight; // initial energy from simulation

double Erecon; // smeared, weighted, and trigger func corrected energy with the EQ2Etrue conversion applied

struct Evis {
  double EvisE;
  double EvisW;
} evis;

struct Edep {
  double EdepE;
  double EdepW;
} edep;

struct EdepQ {
  double EdepQE;
  double EdepQW;
} edepQ;

struct time {
  double timeE;
  double timeW;
} Time;

struct MWPCEnergy {
  double MWPCEnergyE; 
  double MWPCEnergyW;
} mwpcE;

struct MWPCPos {
  double MWPCPosE[3];
  double MWPCPosW[3];
} mwpc_pos;

struct ScintPos {
  double ScintPosE[3];
  double ScintPosW[3];
} scint_pos;

struct ScintPosAdjusted {
  double ScintPosAdjE[3];
  double ScintPosAdjW[3];
} scint_pos_adj;

struct PMT_Evis {
  double Evis[8];
  double weight[8]; 
} pmt_Evis;
