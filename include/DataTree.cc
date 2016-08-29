#include "DataTree.hh"


DataTree::DataTree() : inputFile(NULL), outputFile(NULL), inputTree(NULL), outputTree(NULL) {
  
};

DataTree::~DataTree() {
  //if (UCN_Mon_1_Rate) { std::cout << "trying to delete\n"; delete UCN_Mon_1_Rate; }
  //if (UCN_Mon_2_Rate) delete UCN_Mon_2_Rate;
  //if (UCN_Mon_3_Rate) delete UCN_Mon_3_Rate;
  //if (UCN_Mon_4_Rate) delete UCN_Mon_4_Rate;
  if (outputTree) delete outputTree;
  if (outputFile) {outputFile->Close();  delete outputFile;}
  if (inputTree) delete inputTree;
  if (inputFile) { inputFile->Close(); delete inputFile;}
};

void DataTree::makeOutputTree(std::string outputFileName, std::string outputTreeName) {
  outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputTree = new TTree(outputTreeName.c_str(),outputTreeName.c_str());

  outputTree->Branch("TriggerNum",&TriggerNum, "TriggerNum/I");
  outputTree->Branch("EvtN",&EvtN, "EvtN/I");
  outputTree->Branch("Sis00",&Sis00, "Sis00/I");
  outputTree->Branch("DeltaT",&DeltaT, "DeltaT/D");
  outputTree->Branch("Tof",&Tof, "Tof/D");   
  outputTree->Branch("TimeE",&TimeE, "TimeE/D");
  outputTree->Branch("TimeW",&TimeW, "TimeW/D");
  outputTree->Branch("Time",&Time, "Time/D");
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
  outputTree->Branch("EastBackADC",&EastBackADC,"EastBackADC/D");
  outputTree->Branch("WestBackADC",&WestBackADC,"WestBackADC/D");
  outputTree->Branch("EastBackTDC",&EastBackTDC,"EastBackTDC/D");
  outputTree->Branch("WestBackTDC",&WestBackTDC,"WestBackTDC/D");
  outputTree->Branch("EastDriftVetoADC",&EastDriftVetoADC,"EastDriftVetoADC/D");
  outputTree->Branch("WestDriftVetoADC",&WestDriftVetoADC,"WestDriftVetoADC/D");
  outputTree->Branch("EastTopVetoADC",&EastTopVetoADC,"EastTopVetoADC/D");
  outputTree->Branch("EastTopVetoTDC",&EastTopVetoTDC,"EastTopVetoTDC/D");
  outputTree->Branch("EvnbGood",&EvnbGood,"EvnbGood/O");
  outputTree->Branch("BkhfGood",&BkhfGood,"BkhfGood/O");

  outputTree->Branch("xeRC", &xeRC, "xeRC/I"); //x east response class. 
  outputTree->Branch("yeRC", &yeRC, "yeRC/I"); //y east response class... 
  outputTree->Branch("xwRC", &xwRC, "xwRC/I");
  outputTree->Branch("ywRC", &ywRC, "ywRC/I");

  outputTree->Branch("PID",&PID,"PID/I");
  outputTree->Branch("Side",&Side,"Side/I"); 
  outputTree->Branch("Type",&Type,"Type/I");
  outputTree->Branch("ProbIII",&ProbIII,"ProbIII/D");
  outputTree->Branch("Erecon",&Erecon,"Erecon/D");

  std::cout << "Created output tree " << outputTreeName << " in " << outputFileName << "...\n";
};

void DataTree::writeOutputFile() {
  
  if (outputFile) {
    outputFile->cd();
    UCN_Mon_1_Rate->Write("",TObject::kOverwrite);
    UCN_Mon_2_Rate->Write("",TObject::kOverwrite);
    UCN_Mon_3_Rate->Write("",TObject::kOverwrite);
    UCN_Mon_4_Rate->Write("",TObject::kOverwrite);
    outputTree->Write();
  }
};

void DataTree::setupInputTree(std::string inputFileName, std::string inputTreeName) {
  inputFile = new TFile(inputFileName.c_str(),"READ");
  inputTree = (TTree*)(inputFile->Get(inputTreeName.c_str()));
  //inputTree = (TTree*)(gROOT->FindObject(inputTreeName.c_str()));

  UCN_Mon_1_Rate = (TH1F*)(inputFile->Get("UCN_Mon_1_Rate"));
  UCN_Mon_2_Rate = (TH1F*)(inputFile->Get("UCN_Mon_2_Rate"));
  UCN_Mon_3_Rate = (TH1F*)(inputFile->Get("UCN_Mon_3_Rate"));
  UCN_Mon_4_Rate = (TH1F*)(inputFile->Get("UCN_Mon_4_Rate"));

  inputTree->SetBranchAddress("TriggerNum",&TriggerNum);
  inputTree->SetBranchAddress("EvtN",&EvtN);
  inputTree->SetBranchAddress("Sis00",&Sis00);
  inputTree->SetBranchAddress("DeltaT",&DeltaT);
  inputTree->SetBranchAddress("Tof",&Tof);   
  inputTree->SetBranchAddress("TimeE",&TimeE);
  inputTree->SetBranchAddress("TimeW",&TimeW);
  inputTree->SetBranchAddress("Time",&Time);
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


