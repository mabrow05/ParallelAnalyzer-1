double pmt[8]; // pedestal-subtracted PMT QADC

double cathodeEast[32]; // East pedestal-subtracted cathode MWPC PADC
double cathodeWest[32]; // West pedestal-subtracted cathode MWPC PADC

double xE, yE, xW, yW; // East and West MWPC (x,y) positions

int PID; // Event PID

int type; // Type 0, 1, 2/3

int side; // Trigger side: 0 = East, 1 = West

int posError; // MWPC Position Reconstruction Error

double timeE, timeW; // East/West blinded times [s]

double timeE_BB, timeW_BB; // East/West blinded times since last beam burst [s]

double UBtime; //UNBLINDED time since start of run [s]

double UBtime_BB; //UNBLINDED time since last beam burst [s]

double twoFoldE, twoFoldW; // East/West two-fold trigger TDC


 
