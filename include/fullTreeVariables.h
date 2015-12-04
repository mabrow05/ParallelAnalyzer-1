// Conversion factors
const static double tdcChannelToTime = 180./4096.; // TDC Channel to [ns]
const static double scalerCountsToTime = 1.0e-6;   // Scaler Counts to [s]

float Pdc30; // East MWPC Anode PADC
float Pdc34; // West MWPC Anode PADC

float Tdc016; // East Two-Fold Timing TDC
float Tdc017; // West Two-Fold Timing TDC
float Tdc00;  // East individual PMT trigger TDC
float Tdc01; 
float Tdc02; 
float Tdc03;
float Tdc08;  // West individual PMT trigger TDC
float Tdc09;
float Tdc010;
float Tdc011;


float Sis00; // Sis00 Input Register

float Qadc[8]; // East and West PMT QADC

float Pdc2[32]; // East MWPC Cathode PADC
float Padc[32]; // West MWPC Cathode PADC

float S83028; // UNBLINDED Clock in [us] Since Run Start
float S8200; // UNBLINDED Clock in [us] Since most recent beam burst

float Clk0; // East Blinded time
float Clk1; // West Blinded time
float Clk2; // East Blinded time since most recent beam-burst
float Clk3; // West Blinded time since most recent beam-burst

float Pdc38;  // Gate Valve UCN Monitor PADC
float Pdc39;  // Switcher UCN Monitor PADC
float Pdc310; // AFP Fe Foil UCN Monitor PADC
float Pdc311; // SCS UCN Monitor PADC

float Qadc9;  // East Top Veto QADC
float Tdc019; // East Top Veto TDC
float Pdc313; // East Drift Tube Veto TAC
float Pdc315; // West Drift Tube Veto TAC
float Qadc8;  // East Backing Veto QADC
float Tdc018; // East Backing Veto TDC
float Qadc10; // West Backing Veto QADC
float Tdc020; // West Backing Veto TDC

float Number; // event trigger number
float Delta0; // Time since last event 





