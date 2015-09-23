
int PID, type, side;

double EvisTot; // smeared, weighted, and trigger func corrected energy summed over two sides scint and MWPC energy

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
  double MWPCPosE;
  double MWPCPosW;
} mwpc_pos;

struct ScintPos {
  double ScintPosE;
  double ScintPosW;
} scint_pos;

struct PMT_Evis {
  double Evis[8];
  double weight[8]; 
} pmt_Evis;
