#include "../include/DataTree.hh"


MWPC xE, xW, yE, yW;
MWPC old_xE, old_xW, old_yE, old_yW;
MWPC gaus_xE, gaus_xW, gaus_yE, gaus_yW;

int primaryID, PID, type, side; // basic analysis tags
int hitCountSD[24]; // holds how many times a SD was hit in simulation
double primKE, AsymWeight, primTheta; // initial energy from simulation

double primMom[3];

double Erecon, old_Erecon, gaus_Erecon; // smeared, weighted, and trigger func corrected energy summed over two sides scint and MWPC energy

int nClipped_EX, nClipped_EY, nClipped_WX, nClipped_WY; // Clipping numbers

float Cath_EX[16], Cath_EY[16], Cath_WX[16], Cath_WY[16];

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

struct cathRespPos {
  double cathRespPosE[3];
  double cathRespPosW[3];
} cathResp_pos;

struct PMT {
  double Evis[8];
  double etaEvis[8];
  double nPE[8];
  //double weight[8]; 
} pmt;




