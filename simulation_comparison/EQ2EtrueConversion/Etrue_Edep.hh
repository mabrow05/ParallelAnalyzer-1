/*-----------------------------------------------------------------------------
This is the header file associated with TreeExample3.cpp. In here we will define
all of our structs to hold data from the trees and also write all of our class 
declarations. 
-----------------------------------------------------------------------------*/


#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

using namespace std;

//First I want you to declare all the struct types necessary to read 
//in the branches. These will be used in the class declarations below

//I'll do the first for you. You can add any others you think you need ... 
//This is what will store the Wirechamber information

struct E_dep {
  Double_t EdepE, EdepW;
};

struct E_depQ {
  Double_t EdepQE, EdepQW;
};

struct MWPC_Energy {
  Double_t MWPCEnergyE, MWPCEnergyW;
};

struct MWPC_Pos {
  Double_t MWPCPosE[3], MWPCPosW[3]; 
};

struct Time {
  Double_t timeE, timeW; 
};


//Here is our class declaration. Now a class declaration is simply telling
//the compiler what it can expect in the classes you wish to write.

class DataTree {

//private variables/functions are those which are only accessible within the 
// class. So only the code which exists in the class can use these, and they
// can't be accessed by someone simply using an object of type DataTree

private:
  TFile *inputFile; 
  TTree *inputTree;

  
//Public variables/functions are those which can be accessed by the user
// when he creates an object of the type class. Imagine these as things you
// can use by typing (if your object was called t) t->setupInputTree(args)

public: 
  
  // You need to have a constructor and a destructor for your class. The
  //constructor is what is called when you instantiate an object of type
  //DataTree. Destructor is called when the object is being deleted..
  DataTree(); //Constructor
  ~DataTree(); //Destructor
  
  //These are some functions which will be useful. You need to define them
  //in the .cpp

  // This function is only declared. We don't see the code.
  void setupInputTree(string inputFileName, string inputTreeName); 

  // Note that this one is already DEFINED here. We see what it will do in the {}
  void getEvent(UInt_t N) {inputTree->GetEvent(N);} 
  Int_t getEntries() {return inputTree->GetEntriesFast();}

  //Here you need to define all the variables where you will store the tree
  //information. This is things like energy, position, etc. You don't have 
  //to store them all, but at least save several. Here are a few examples:

  Double_t primKE, primTheta;
  E_dep Edep;
  E_depQ EdepQ;
  MWPC_Pos MWPCPos;
  MWPC_Energy MWPCEnergy;
  Time time;

}; // End of class DataTree Declaration
