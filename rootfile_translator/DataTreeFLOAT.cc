#include "DataTreeFLOAT.hh"


DataTreeFLOAT::DataTreeFLOAT() : inputFile(NULL), outputFile(NULL), inputTree(NULL), outputTree(NULL) {
  
};

DataTreeFLOAT::~DataTreeFLOAT() {
  if (outputTree) delete outputTree;
  if (outputFile) {outputFile->Close();  delete outputFile;}
  if (inputTree) delete inputTree;
  if (inputFile) { inputFile->Close(); delete inputFile;}
};

void DataTreeFLOAT::makeOutputTree(std::string outputFileName, std::string outputTreeName) {
  outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputTree = new TTree(outputTreeName.c_str(),outputTreeName.c_str());

  outputTree->Branch("TriggerNum",&TriggerNum, "TriggerNum/I");
  outputTree->Branch("EvtN",&EvtN, "EvtN/I");
  outputTree->Branch("Sis00",&Sis00, "Sis00/I");
  outputTree->Branch("DeltaT",&DeltaT, "DeltaT/F");
  outputTree->Branch("Tof",&Tof, "Tof/F");   
  outputTree->Branch("TimeE",&TimeE, "TimeE/F");
  outputTree->Branch("TimeW",&TimeW, "TimeW/F");
  outputTree->Branch("TDCE",&TDCE,"TDCE/F");
  outputTree->Branch("TDCW",&TDCW,"TDCW/F");
  outputTree->Branch("TDCE1",&TDCE1,"TDCE1/F");
  outputTree->Branch("TDCE2",&TDCE2,"TDCE2/F");
  outputTree->Branch("TDCE3",&TDCE3,"TDCE3/F");
  outputTree->Branch("TDCE4",&TDCE4,"TDCE4/F");
  outputTree->Branch("TDCW1",&TDCW1,"TDCW1/F");
  outputTree->Branch("TDCW2",&TDCW2,"TDCW2/F");
  outputTree->Branch("TDCW3",&TDCW3,"TDCW3/F");
  outputTree->Branch("TDCW4",&TDCW4,"TDCW4/F");
  outputTree->Branch("xEmpm",&xE,"center/F:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/F:height");
  outputTree->Branch("yEmpm",&yE,"center/F:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/F:height");
  outputTree->Branch("xWmpm",&xW,"center/F:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/F:height");
  outputTree->Branch("yWmpm",&yW,"center/F:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/F:height");
  outputTree->Branch("Cathodes_Ex",&Cathodes_Ex,"Cathodes_Ex/F");
  outputTree->Branch("Cathodes_Ey",&Cathodes_Ey,"Cathodes_Ey/F");
  outputTree->Branch("Cathodes_Wx",&Cathodes_Wx,"Cathodes_Wx/F");
  outputTree->Branch("Cathodes_Wy",&Cathodes_Wy,"Cathodes_Wy/F");
  outputTree->Branch("ScintE", &ScintE, "q1/F:q2:q3:q4:e1:de1:e2:de2:e3:de3:e4:de4:energy:denergy:nPE1:nPE2:nPE3:nPE4");
  outputTree->Branch("ScintW", &ScintW, "q1/F:q2:q3:q4:e1:de1:e2:de2:e3:de3:e4:de4:energy:denergy:nPE1:nPE2:nPE3:nPE4");
  outputTree->Branch("EvisE",&EvisE,"EvisE/F");
  outputTree->Branch("EvisW",&EvisW,"EvisW/F");
  outputTree->Branch("CathSumE",&CathSumE,"CathSumE/F");
  outputTree->Branch("CathSumW",&CathSumW,"CathSumW/F"); 
  outputTree->Branch("CathMaxE",&CathMaxE,"CathMaxE/F");
  outputTree->Branch("CathMaxW",&CathMaxW,"CathMaxW/F"); 
  outputTree->Branch("EMWPC_E",&EMWPC_E,"EMWPC_E/F");
  outputTree->Branch("EMWPC_W",&EMWPC_W,"EMWPC_W/F"); 
  outputTree->Branch("AnodeE",&AnodeE,"AnodeE/F");
  outputTree->Branch("AnodeW",&AnodeW,"AnodeW/F"); 
  outputTree->Branch("PassedAnoE",&PassedAnoE,"PassedAnoE/O");
  outputTree->Branch("PassedAnoW",&PassedAnoW,"PassedAnoW/O");
  outputTree->Branch("PassedCathE",&PassedCathE,"PassedCathE/O");
  outputTree->Branch("PassedCathW",&PassedCathW,"PassedCathW/O");
  outputTree->Branch("TaggedBackE",&TaggedBackE,"TaggedBackE/O");
  outputTree->Branch("TaggedBackW",&TaggedBackW,"TaggedBackW/O");
  outputTree->Branch("TaggedTopE",&TaggedTopE,"TaggedTopE/O");
  outputTree->Branch("TaggedTopW",&TaggedTopW,"TaggedTopW/O");
  outputTree->Branch("TaggedDriftE",&TaggedDriftE,"TaggedDriftE/O");
  outputTree->Branch("TaggedDriftW",&TaggedDriftW,"TaggedDriftW/O");
  outputTree->Branch("EastBackADC",&EastBackADC,"EastBackADC/F");
  outputTree->Branch("WestBackADC",&WestBackADC,"WestBackADC/F");
  outputTree->Branch("EastBackTDC",&EastBackTDC,"EastBackTDC/F");
  outputTree->Branch("WestBackTDC",&WestBackTDC,"WestBackTDC/F");
  outputTree->Branch("EastDriftVetoADC",&EastDriftVetoADC,"EastDriftVetoADC/F");
  outputTree->Branch("WestDriftVetoADC",&WestDriftVetoADC,"WestDriftVetoADC/F");
  outputTree->Branch("EastTopVetoADC",&EastTopVetoADC,"EastTopVetoADC/F");
  outputTree->Branch("EastTopVetoTDC",&EastTopVetoTDC,"EastTopVetoTDC/F");
  outputTree->Branch("EvnbGood",&EvnbGood,"EvnbGood/O");
  outputTree->Branch("BkhfGood",&BkhfGood,"BkhfGood/O");

  outputTree->Branch("xeRC", &xeRC, "xeRC/I"); //x east response class. 
  outputTree->Branch("yeRC", &yeRC, "yeRC/I"); //y east response class... 
  outputTree->Branch("xwRC", &xwRC, "xwRC/I");
  outputTree->Branch("ywRC", &ywRC, "ywRC/I");

  outputTree->Branch("PID",&PID,"PID/I");
  outputTree->Branch("Side",&Side,"Side/I"); 
  outputTree->Branch("Type",&Type,"Type/I");
  outputTree->Branch("ProbIII",&ProbIII,"ProbIII/F");
  outputTree->Branch("Erecon",&Erecon,"Erecon/F");

  std::cout << "Created output tree " << outputTreeName << " in " << outputFileName << "...\n";
};

void DataTreeFLOAT::setupInputTree(std::string inputFileName, std::string inputTreeName) {
  inputFile = new TFile(inputFileName.c_str(),"READ");
  inputTree = (TTree*)(inputFile->Get(inputTreeName.c_str()));

  inputTree->SetBranchAddress("TriggerNum",&TriggerNum);
  inputTree->SetBranchAddress("EvtN",&EvtN);
  inputTree->SetBranchAddress("Sis00",&Sis00);
  inputTree->SetBranchAddress("DeltaT",&DeltaT);
  inputTree->SetBranchAddress("Tof",&Tof);   
  inputTree->SetBranchAddress("TimeE",&TimeE);
  inputTree->SetBranchAddress("TimeW",&TimeW);
  inputTree->SetBranchAddress("TDCE",&TDCE);
  inputTree->SetBranchAddress("TDCW",&TDCW);
  inputTree->SetBranchAddress("TDCE1",&TDCE1);
  inputTree->SetBranchAddress("TDCE2",&TDCE2);
  inputTree->SetBranchAddress("TDCE3",&TDCE3);
  inputTree->SetBranchAddress("TDCE4",&TDCE4);
  inputTree->SetBranchAddress("TDCW1",&TDCW1);
  inputTree->SetBranchAddress("TDCW2",&TDCW2);
  inputTree->SetBranchAddress("TDCW3",&TDCW3);
  inputTree->SetBranchAddress("TDCW4",&TDCW4);
  inputTree->SetBranchAddress("xE",&xE);
  inputTree->SetBranchAddress("yE",&yE);
  inputTree->SetBranchAddress("xW",&xW);
  inputTree->SetBranchAddress("yW",&yW);
  inputTree->SetBranchAddress("Cathodes_Ex",&Cathodes_Ex);
  inputTree->SetBranchAddress("Cathodes_Ey",&Cathodes_Ey);
  inputTree->SetBranchAddress("Cathodes_Wx",&Cathodes_Wx);
  inputTree->SetBranchAddress("Cathodes_Wy",&Cathodes_Wy);
  inputTree->SetBranchAddress("ScintE", &ScintE);
  inputTree->SetBranchAddress("ScintW", &ScintW);
  inputTree->SetBranchAddress("EvisE",&EvisE);
  inputTree->SetBranchAddress("EvisW",&EvisW);
  inputTree->SetBranchAddress("CathSumE",&CathSumE);
  inputTree->SetBranchAddress("CathSumW",&CathSumW); 
  inputTree->SetBranchAddress("CathMaxE",&CathMaxE);
  inputTree->SetBranchAddress("CathMaxW",&CathMaxW); 
  inputTree->SetBranchAddress("EMWPC_E",&EMWPC_E);
  inputTree->SetBranchAddress("EMWPC_W",&EMWPC_W); 
  inputTree->SetBranchAddress("AnodeE",&AnodeE);
  inputTree->SetBranchAddress("AnodeW",&AnodeW); 
  inputTree->SetBranchAddress("PassedAnoE",&PassedAnoE);
  inputTree->SetBranchAddress("PassedAnoW",&PassedAnoW);
  inputTree->SetBranchAddress("PassedCathE",&PassedCathE);
  inputTree->SetBranchAddress("PassedCathW",&PassedCathW);
  inputTree->SetBranchAddress("TaggedBackE",&TaggedBackE);
  inputTree->SetBranchAddress("TaggedBackW",&TaggedBackW);
  inputTree->SetBranchAddress("TaggedTopE",&TaggedTopE);
  inputTree->SetBranchAddress("TaggedTopW",&TaggedTopW);
  inputTree->SetBranchAddress("TaggedDriftE",&TaggedDriftE);
  inputTree->SetBranchAddress("TaggedDriftW",&TaggedDriftW);
  inputTree->SetBranchAddress("EastBackADC",&EastBackADC);
  inputTree->SetBranchAddress("WestBackADC",&WestBackADC);
  inputTree->SetBranchAddress("EastBackTDC",&EastBackTDC);
  inputTree->SetBranchAddress("WestBackTDC",&WestBackTDC);
  inputTree->SetBranchAddress("EastDriftVetoADC",&EastDriftVetoADC);
  inputTree->SetBranchAddress("WestDriftVetoADC",&WestDriftVetoADC);
  inputTree->SetBranchAddress("EastTopVetoADC",&EastTopVetoADC);
  inputTree->SetBranchAddress("EastTopVetoTDC",&EastTopVetoTDC);
  inputTree->SetBranchAddress("EvnbGood",&EvnbGood);
  inputTree->SetBranchAddress("BkhfGood",&BkhfGood);

  inputTree->SetBranchAddress("xeRC", &xeRC); //x east response class. 
  inputTree->SetBranchAddress("yeRC", &yeRC); //y east response class... 
  inputTree->SetBranchAddress("xwRC", &xwRC);
  inputTree->SetBranchAddress("ywRC", &ywRC);

  inputTree->SetBranchAddress("PID",&PID);
  inputTree->SetBranchAddress("Side",&Side); 
  inputTree->SetBranchAddress("Type",&Type);
  inputTree->SetBranchAddress("ProbIII",&ProbIII);
  inputTree->SetBranchAddress("Erecon",&Erecon);

  std::cout << "Prepared input tree " << inputTreeName << " in " << inputFileName << "...\n";
};


