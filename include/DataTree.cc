#include "makeOutputTree.hh"


DataTree::DataTree(std::string tree) : infile(NULL), outfile(NULL), inputTree(NULL), outputTree(NULL) {
  
};

DataTree::~DataTree() {
  if (outputTree) delete outputTree;
  if (outPutFile) {outputFile.Close();  delete outputFile;}
  if (inputTree) delete inputTree;
  if (inputFile) { inputFile.Close(); delete inputFile;}
};

void DataTree::makeOutputTree(std::string outputFile, std::string outputTreeName) {
  outputFile = new TFile(outputFile.c_str(),"RECREATE");
  outputTree = new TTree(outputTreeName.c_str());

  outputTree->Branch("TriggerNum",&TriggerNum, "TriggerNum/I");
  outputTree->Branch("EvtN",&EvtN, "EvtN/I");
  outputTree->Branch("Sis00",&Sis00, "Sis00/I");
  outputTree->Branch("DeltaT",&DeltaT, "DeltaT/D");
  outputTree->Branch("Tof",&Tof, "Tof/D");   
  outputTree->Branch("TimeE",&TimeE, "TimeE/D");
  outputTree->Branch("TimeW",&TimeW, "TimeW/D");
  outputTree->Branch("TDCE",&TDCE,"TDCE/D");
  outputTree->Branch("TDCW",&TDCW,"TDCW/D");
  outputTree->Branch("TDCE1",&TDCE1,"TDCE1/D");
  outputTree->Branch("TDCE2",&TDCE2,"TDCE2/D");
  outputTree->Branch("TDCE3",&TDCE3,"TDCE3/D");
  outputTree->Branch("TDCE4",&TDCE4,"TDCE4/D");
  outputTree->Branch("TDCW1",&TDCW1,"TDCW1/D");
  outputTree->Branch("TDCW2",&TDCW2,"TDCW2/D");
  outputTree->Branch("TDCW3",&TDCW3,"TDCW3/D");
  outputTree->Branch("TDCW4",&TDCW4,"TDCW4/D");
  outputTree->Branch("xE",&xE,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  outputTree->Branch("yE",&yE,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  outputTree->Branch("xW",&xW,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  outputTree->Branch("yW",&yW,"center/D:width:cathSum:maxValue:maxWire/I:mult:nClipped:err:rawCenter/D:height");
  outputTree->Branch("Cathodes_Ex",&Cathodes_Ex,"Cathodes_Ex/D");
  outputTree->Branch("Cathodes_Ey",&Cathodes_Ey,"Cathodes_Ey/D");
  outputTree->Branch("Cathodes_Wx",&Cathodes_Wx,"Cathodes_Wx/D");
  outputTree->Branch("Cathodes_Wy",&Cathodes_Wy,"Cathodes_Wy/D");
  outputTree->Branch("ScintE", &ScintE, "q1/D:q2:q3:q4:e1:de1:e2:de2:e3:de3:e4:de4:energy:denergy:nPE1:nPE2:nPE3:nPE4");
  outputTree->Branch("ScintW", &ScintW, "q1/D:q2:q3:q4:e1:de1:e2:de2:e3:de3:e4:de4:energy:denergy:nPE1:nPE2:nPE3:nPE4");
  outputTree->Branch("EvisE",&EvisE,"EvisE/D");
  outputTree->Branch("EvisW",&EvisW,"EvisW/D");
  outputTree->Branch("CathSumE",&CathSumE,"CathSumE/D");
  outputTree->Branch("CathSumW",&CathSumW,"CathSumW/D"); 
  outputTree->Branch("CathMaxE",&CathMaxE,"CathMaxE/D");
  outputTree->Branch("CathMaxW",&CathMaxW,"CathMaxW/D"); 
  outputTree->Branch("EMWPC_E",&EMWPC_E,"EMWPC_E/D");
  outputTree->Branch("EMWPC_W",&EMWPC_W,"EMWPC_W/D"); 
  outputTree->Branch("AnodeE",&AnodeE,"AnodeE/D");
  outputTree->Branch("AnodeW",&AnodeW,"AnodeW/D"); 
  outputTree->Branch("PassedAnoE",&PassedAnoE,"PassedAnoE/I");
  outputTree->Branch("PassedAnoW",&PassedAnoW,"PassedAnoW/I");
  outputTree->Branch("PassedCathE",&PassedCathE,"PassedCathE/I");
  outputTree->Branch("PassedCathW",&PassedCathW,"PassedCathW/I");
  outputTree->Branch("TaggedBackE",&TaggedBackE,"TaggedBackE/I");
  outputTree->Branch("TaggedBackW",&TaggedBackW,"TaggedBackW/I");
  outputTree->Branch("TaggedTopE",&TaggedTopE,"TaggedTopE/I");
  outputTree->Branch("TaggedTopW",&TaggedTopW,"TaggedTopW/I");
  outputTree->Branch("TaggedDriftE",&TaggedDriftE,"TaggedDriftE/I");
  outputTree->Branch("TaggedDriftW",&TaggedDriftW,"TaggedDriftW/I");
  outputTree->Branch("EastBackADC",&EastBackADC,"EastBackADC/D");
  outputTree->Branch("WestBackADC",&WestBackADC,"WestBackADC/D");
  outputTree->Branch("EastBackTDC",&EastBackTDC,"EastBackTDC/D");
  outputTree->Branch("WestBackTDC",&WestBackTDC,"WestBackTDC/D");
  outputTree->Branch("EastDriftVetoADC",&EastDriftVetoADC,"EastDriftVetoADC/D");
  outputTree->Branch("WestDriftVetoADC",&WestDriftVetoADC,"WestDriftVetoADC/D");
  outputTree->Branch("EastTopVetoADC",&EastTopVetoADC,"EastTopVetoADC/D");
  outputTree->Branch("EastTopVetoTDC",&EastTopVetoTDC,"EastTopVetoTDC/D");
  outputTree->Branch("EvnbGood",&EvnbGood,"EvnbGood/I");
  outputTree->Branch("BkhfGood",&BkhfGood,"BkhfGood/I");

  outputTree->Branch("PID",&PID,"PID/I");
  outputTree->Branch("Side",&Side,"Side/I"); 
  outputTree->Branch("Type",&Type,"Type/I");
  outputTree->Branch("ProbIII",&ProbIII,"ProbIII/D");
  outputTree->Branch("Erecon",&Erecon,"Erecon/D");

  std::cout << "Created output tree " << outputTreeName << " in " << outputFile << "...\n";
};

void DataTree::setupInputTree(std::string inputFile, std::string inputTreeName) {
  outputFile = new TFile(inputFile.c_str(),"READ");
  outputTree = (TTree*)(outputFile->Get(inputTreeName.c_str()));

  outputTree->SetBranchAddress("TriggerNum",&TriggerNum);
  outputTree->SetBranchAddress("EvtN",&EvtN);
  outputTree->SetBranchAddress("Sis00",&Sis00);
  outputTree->SetBranchAddress("DeltaT",&DeltaT);
  outputTree->SetBranchAddress("Tof",&Tof);   
  outputTree->SetBranchAddress("TimeE",&TimeE);
  outputTree->SetBranchAddress("TimeW",&TimeW);
  outputTree->SetBranchAddress("TDCE",&TDCE);
  outputTree->SetBranchAddress("TDCW",&TDCW);
  outputTree->SetBranchAddress("TDCE1",&TDCE1);
  outputTree->SetBranchAddress("TDCE2",&TDCE2);
  outputTree->SetBranchAddress("TDCE3",&TDCE3);
  outputTree->SetBranchAddress("TDCE4",&TDCE4);
  outputTree->SetBranchAddress("TDCW1",&TDCW1);
  outputTree->SetBranchAddress("TDCW2",&TDCW2);
  outputTree->SetBranchAddress("TDCW3",&TDCW3);
  outputTree->SetBranchAddress("TDCW4",&TDCW4);
  outputTree->SetBranchAddress("xE",&xE);
  outputTree->SetBranchAddress("yE",&yE);
  outputTree->SetBranchAddress("xW",&xW);
  outputTree->SetBranchAddress("yW",&yW);
  outputTree->SetBranchAddress("Cathodes_Ex",&Cathodes_Ex);
  outputTree->SetBranchAddress("Cathodes_Ey",&Cathodes_Ey);
  outputTree->SetBranchAddress("Cathodes_Wx",&Cathodes_Wx);
  outputTree->SetBranchAddress("Cathodes_Wy",&Cathodes_Wy);
  outputTree->SetBranchAddress("ScintE", &ScintE);
  outputTree->SetBranchAddress("ScintW", &ScintW);
  outputTree->SetBranchAddress("EvisE",&EvisE);
  outputTree->SetBranchAddress("EvisW",&EvisW);
  outputTree->SetBranchAddress("CathSumE",&CathSumE);
  outputTree->SetBranchAddress("CathSumW",&CathSumW); 
  outputTree->SetBranchAddress("CathMaxE",&CathMaxE);
  outputTree->SetBranchAddress("CathMaxW",&CathMaxW); 
  outputTree->SetBranchAddress("EMWPC_E",&EMWPC_E);
  outputTree->SetBranchAddress("EMWPC_W",&EMWPC_W); 
  outputTree->SetBranchAddress("AnodeE",&AnodeE);
  outputTree->SetBranchAddress("AnodeW",&AnodeW); 
  outputTree->SetBranchAddress("PassedAnoE",&PassedAnoE);
  outputTree->SetBranchAddress("PassedAnoW",&PassedAnoW);
  outputTree->SetBranchAddress("PassedCathE",&PassedCathE);
  outputTree->SetBranchAddress("PassedCathW",&PassedCathW);
  outputTree->SetBranchAddress("TaggedBackE",&TaggedBackE);
  outputTree->SetBranchAddress("TaggedBackW",&TaggedBackW);
  outputTree->SetBranchAddress("TaggedTopE",&TaggedTopE);
  outputTree->SetBranchAddress("TaggedTopW",&TaggedTopW);
  outputTree->SetBranchAddress("TaggedDriftE",&TaggedDriftE);
  outputTree->SetBranchAddress("TaggedDriftW",&TaggedDriftW);
  outputTree->SetBranchAddress("EastBackADC",&EastBackADC);
  outputTree->SetBranchAddress("WestBackADC",&WestBackADC);
  outputTree->SetBranchAddress("EastBackTDC",&EastBackTDC);
  outputTree->SetBranchAddress("WestBackTDC",&WestBackTDC);
  outputTree->SetBranchAddress("EastDriftVetoADC",&EastDriftVetoADC);
  outputTree->SetBranchAddress("WestDriftVetoADC",&WestDriftVetoADC);
  outputTree->SetBranchAddress("EastTopVetoADC",&EastTopVetoADC);
  outputTree->SetBranchAddress("EastTopVetoTDC",&EastTopVetoTDC);
  outputTree->SetBranchAddress("EvnbGood",&EvnbGood);
  outputTree->SetBranchAddress("BkhfGood",&BkhfGood);

  outputTree->SetBranchAddress("PID",&PID);
  outputTree->SetBranchAddress("Side",&Side); 
  outputTree->SetBranchAddress("Type",&Type);
  outputTree->SetBranchAddress("ProbIII",&ProbIII);
  outputTree->SetBranchAddress("Erecon",&Erecon);

  std::cout << "Created output tree " << outputTreeName << " in " << outputFile << "...\n";
};
