//This program gives the output root file with variables for Acp calculation
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "iostream"

void Dk3pi_cp_skim2 (TString input, TString output) {

  std::cout << "Starting" << std::endl;
  TChain* tc = new TChain("D0tree");
  tc->Add(input);

// TTree* tc = (TTree*) file1-> Get ("D0tree") ;//taking tree from the input root file 
  
  std::cout << "Got tree from file" << std::endl;

// Define variables that need to take from the root file
  Double_t eventRandom;
  Double_t isSignal;
  Double_t M;
  Double_t CMS_p;
  Double_t CMS_cosTheta;
  //Double_t K_kaonID;
  //Double_t K_pionID;
  Double_t K_binaryID;
  Double_t K_charge;
  Double_t K_p;    
  Double_t K_cosTheta;
  Double_t pi_1_charge;
  Double_t pi_2_charge;
  Double_t pi_3_charge;
  Double_t pi_1_p;
  Double_t pi_2_p;
  Double_t pi_3_p;
  Double_t pi_1_cosTheta;
  Double_t pi_2_cosTheta;
  Double_t pi_3_cosTheta;
  //Double_t pi_1_pionIDNN;
  //Double_t pi_2_pionIDNN;
  //Double_t pi_3_pionIDNN;
  
  std::cout << "Defined Varables to get from input root file" << std::endl;

  // Taking those variables from the root file to the tc tree
  tc->SetBranchAddress("eventRandom", &eventRandom);
  tc->SetBranchAddress("isSignal", &isSignal);
  tc->SetBranchAddress("M", &M);
  tc->SetBranchAddress("CMS_cosTheta", &CMS_cosTheta);
  tc->SetBranchAddress("CMS_p", &CMS_p);  
  tc->SetBranchAddress("K_charge", &K_charge);
  tc->SetBranchAddress("pi_1_charge", &pi_1_charge);
  tc->SetBranchAddress("pi_2_charge", &pi_2_charge);  
  tc->SetBranchAddress("pi_3_charge", &pi_3_charge);
  tc->SetBranchAddress("K_binaryID", &K_binaryID);
  tc->SetBranchAddress("K_p", &K_p);
  tc->SetBranchAddress("K_cosTheta", &K_cosTheta);
  tc->SetBranchAddress("pi_1_p", &pi_1_p);
  tc->SetBranchAddress("pi_2_p", &pi_2_p);
  tc->SetBranchAddress("pi_3_p", &pi_3_p);
  tc->SetBranchAddress("pi_1_cosTheta", &pi_1_cosTheta);
  tc->SetBranchAddress("pi_2_cosTheta", &pi_2_cosTheta);
  tc->SetBranchAddress("pi_3_cosTheta", &pi_3_cosTheta);
  //tc->SetBranchAddress("pi_1_pionIDNN", &pi_1_pionIDNN);
  //tc->SetBranchAddress("pi_2_pionIDNN", &pi_2_pionIDNN);
  //tc->SetBranchAddress("pi_3_pionIDNN", &pi_3_pionIDNN);

  std::cout << "Variables taken" << std::endl;

  // Taking the total number of entries from  D0tree
  Long64_t nEntries(tc->GetEntries());
  
  std::cout<< "Got Number of Entries in tree = "<< nEntries<< std::endl;
  
  // New File and Tree for storing 
  TFile *file2 = new TFile (output,"RECREATE");
  TTree* LowMod = new TTree("D0tree", "a new tree");//TTree *tree = new TTree(name, title) 
  std::cout << "New output root file defined" << std::endl;

  // Define what to right to the Tree as leafs
  LowMod->Branch("eventRandom", &eventRandom, "eventRandom/D");
  LowMod->Branch("isSignal", &isSignal, "isSignal/D");
  LowMod->Branch("M", &M, "M/D");
  LowMod->Branch("CMS_cosTheta", &CMS_cosTheta, "CMS_cosTheta/D");
  LowMod->Branch("CMS_p", &CMS_p, "CMS_p/D");  
  LowMod->Branch("K_charge", &K_charge, "K_charge/D");
  LowMod->Branch("pi_1_charge", &pi_1_charge, "pi_1_charge/D");
  LowMod->Branch("pi_2_charge", &pi_2_charge, "pi_2_charge/D");
  LowMod->Branch("pi_3_charge", &pi_3_charge, "pi_3_charge/D");
  //LowMod->Branch("K_kaonID", &K_kaonID, "K_kaonID/D");
  //LowMod->Branch("K_pionID", &K_pionID, "K_pionID/D");
  LowMod->Branch("K_binaryID", &K_binaryID, "K_binaryID/D");
  LowMod->Branch("K_p", &K_p, "K_p/D");
  LowMod->Branch("K_cosTheta", &K_cosTheta, "K_cosTheta/D");  
  LowMod->Branch("pi_1_p", &pi_1_p, "pi_1_p/D");
  LowMod->Branch("pi_2_p", &pi_2_p, "pi_2_p/D");
  LowMod->Branch("pi_3_p", &pi_3_p, "pi_3_p/D");
  LowMod->Branch("pi_1_cosTheta", &pi_1_cosTheta, "pi_1_cosTheta/D");  
  LowMod->Branch("pi_2_cosTheta", &pi_2_cosTheta, "pi_2_cosTheta/D");
  LowMod->Branch("pi_3_cosTheta", &pi_3_cosTheta, "pi_3_cosTheta/D");
  //LowMod->Branch("pi_1_pionIDNN", &pi_1_pionIDNN, "pi_1_pionIDNN/D");
  //LowMod->Branch("pi_2_pionIDNN", &pi_2_pionIDNN, "pi_2_pionIDNN/D");
  //LowMod->Branch("pi_3_pionIDNN", &pi_3_pionIDNN, "pi_3_pionIDNN/D");

  std::cout << "defined what to right to the Tree as leafs" << std::endl;

  std::cout << "Starting for loop" << std::endl;
  for (Long64_t iEntry(0); iEntry < nEntries; ++iEntry){
         tc->GetEntry(iEntry);//taking the entree corresponding to that number
	 if (eventRandom<0.25 && 1.80<M && M<1.92 && K_p<4.0 && pi_1_p<4.0 && pi_2_p<4.0 && pi_3_p<4.0){
	 LowMod->Fill();//filling the tree using the variables
            }
  }
                                                                                 
  std::cout << "End for loop" << std::endl;
  
  LowMod->Write(); //writing the tree TO THE root file
  std::cout << "Wrote the TREE to the root " << std::endl;
  file2->Close();
  
}


// root -l 'Dk3pi_cp_skim2.C("/b2diskb/dinura/dk3pi/MC/ntuples/k3pi_*.root","/b2diska/dinura/rootfiles/ntuples/MC/k3pi_complete.root")'





