//This program gives the output root file with variables for Acp calculation
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "iostream"

void lc_cp_skim2 (TString input, TString output) {

  std::cout << "Starting" << std::endl;
  TChain* tc = new TChain("lcp_pkk");
  tc->Add(input);

//taking tree from the input root file 
  
  std::cout << "Got tree from file" << std::endl;

// define variables that need to take from the root file
  Double_t eventRandom;
  Double_t Lambdac_isSignal;
  Double_t Lambdac_M;
  Double_t Lambdac_significanceOfDistance;
  Double_t Lambdac_flightDistance;
  Double_t Lambdac_p_cms;
  Double_t Lambdac_cosTheta_cms;
  Double_t Lambdac_p;
  Double_t Lambdac_E;
  //Double_t Lambdac_mcErrors;
  Double_t p_protonID;
  Double_t p_kaonID;
  Double_t p_pionID;
  Double_t p_trinaryID;
  Double_t p_p;
  Double_t p_E;
  Double_t p_dr;
  Double_t p_dz;
  //Double_t p_cosTheta_cms;
  Double_t p_cosTheta;
  Double_t p_charge;
  Double_t K1_kaonID;
  Double_t K1_pionID;
  Double_t K1_binaryID;
  Double_t K1_charge;
  Double_t K1_p;
  Double_t K1_E;
  Double_t K1_dr;
  Double_t K1_dz;
  Double_t K1_cosTheta;
  Double_t K2_kaonID;
  Double_t K2_pionID;
  Double_t K2_binaryID;
  Double_t K2_charge;
  Double_t K2_p;
  Double_t K2_E;
  Double_t K2_dr;
  Double_t K2_dz;
  Double_t K2_cosTheta;

  std::cout << "Defined Varables to get from input root file" << std::endl;

  // Taking those variables from the root file to the tc tree
  tc->SetBranchAddress("eventRandom", &eventRandom);
  tc->SetBranchAddress("Lambdac_isSignal", &Lambdac_isSignal);
  tc->SetBranchAddress("Lambdac_M", &Lambdac_M);
  tc->SetBranchAddress("Lambdac_cosTheta_cms", &Lambdac_cosTheta_cms);
  tc->SetBranchAddress("Lambdac_flightDistance", &Lambdac_flightDistance);
  tc->SetBranchAddress("Lambdac_significanceOfDistance", &Lambdac_significanceOfDistance);
  tc->SetBranchAddress("Lambdac_p_cms", &Lambdac_p_cms);
  tc->SetBranchAddress("Lambdac_p", &Lambdac_p);
  tc->SetBranchAddress("Lambdac_E", &Lambdac_E);
  //tc->SetBranchAddress("Lambdac_mcErrors", &Lambdac_mcErrors);
  tc->SetBranchAddress("p_charge", &p_charge);
  tc->SetBranchAddress("p_p", &p_p);
  tc->SetBranchAddress("p_E", &p_E);
  tc->SetBranchAddress("p_dr", &p_dr);
  tc->SetBranchAddress("p_dz", &p_dz);
  //tc->SetBranchAddress("p_cosTheta_cms", &p_cosTheta_cms);
  tc->SetBranchAddress("p_cosTheta", &p_cosTheta);
  tc->SetBranchAddress("p_protonID", &p_protonID);
  tc->SetBranchAddress("p_pionID", &p_pionID);
  tc->SetBranchAddress("p_kaonID", &p_kaonID);
  tc->SetBranchAddress("K1_kaonID", &K1_kaonID);
  tc->SetBranchAddress("K1_pionID", &K1_pionID);
  tc->SetBranchAddress("K1_charge", &K1_charge);
  tc->SetBranchAddress("K1_p", &K1_p);
  tc->SetBranchAddress("K1_E", &K1_E);
  tc->SetBranchAddress("K1_dr", &K1_dr);
  tc->SetBranchAddress("K1_dz", &K1_dz);
  tc->SetBranchAddress("K1_cosTheta", &K1_cosTheta);
  tc->SetBranchAddress("K2_kaonID", &K2_kaonID);
  tc->SetBranchAddress("K2_pionID", &K2_pionID);
  tc->SetBranchAddress("K2_charge", &K2_charge);
  tc->SetBranchAddress("K2_p", &K2_p);
  tc->SetBranchAddress("K2_E", &K2_E);
  tc->SetBranchAddress("K2_dr", &K2_dr);
  tc->SetBranchAddress("K2_dz", &K2_dz);
  tc->SetBranchAddress("K2_cosTheta", &K2_cosTheta);


  std::cout << "Variables taken" << std::endl;

  // taking the total number of entries from my Bplus tree
  Long64_t nEntries(tc->GetEntries());
  
  std::cout<< "Got Number of Entries in tree = "<< nEntries<< std::endl;
  
  //New File and Tree for storing 
  TFile *file2 = new TFile (output,"RECREATE");
  TTree* LowMod = new TTree("lcp_pkk", "a new tree");//TTree *tree = new TTree(name, title) 
  std::cout << "New output root file defined" << std::endl;

  //define what to right to  the Tree as leafs
  LowMod->Branch("eventRandom", &eventRandom, "eventRandom/D");
  LowMod->Branch("Lambdac_isSignal", &Lambdac_isSignal, "Lambdac_isSignal/D");
  LowMod->Branch("Lambdac_M", &Lambdac_M, "Lambdac_M/D");
  LowMod->Branch("Lambdac_cosTheta_cms", &Lambdac_cosTheta_cms, "Lambdac_cosTheta_cms/D");
  LowMod->Branch("Lambdac_flightDistance", &Lambdac_flightDistance, "Lambdac_flightDistance/D");
  LowMod->Branch("Lambdac_significanceOfDistance", &Lambdac_significanceOfDistance, "Lambdac_significanceOfDistance/D");
  LowMod->Branch("Lambdac_p_cms", &Lambdac_p_cms, "Lambdac_p_cms/D");
  LowMod->Branch("Lambdac_p", &Lambdac_p, "Lambdac_p/D");
  LowMod->Branch("Lambdac_E", &Lambdac_E, "Lambdac_E/D");
  //LowMod->Branch("Lambdac_mcErrors", &Lambdac_mcErrors, "Lambdac_mcErrors/D");
  LowMod->Branch("p_charge", &p_charge, "p_charge/D");
  //LowMod->Branch("p_cosTheta_cms", &p_cosTheta_cms, "p_cosTheta_cms/D");
  LowMod->Branch("p_cosTheta", &p_cosTheta, "p_cosTheta/D");
  LowMod->Branch("p_p", &p_p, "p_p/D");
  LowMod->Branch("p_E", &p_E, "p_E/D");
  LowMod->Branch("p_dr", &p_dr, "p_dr/D");
  LowMod->Branch("p_dz", &p_dz, "p_dz/D");
  //LowMod->Branch("p_protonID", &p_protonID, "p_protonID/D");
  //LowMod->Branch("p_kaonID", &p_kaonID, "p_kaonID/D");
  //LowMod->Branch("p_pionID", &p_pionID, "p_pionID/D");
  LowMod->Branch("p_trinaryID", &p_trinaryID, "p_trinaryID/D");
  //LowMod->Branch("K1_kaonID", &K1_kaonID, "K1_kaonID/D");
  //LowMod->Branch("K1_pionID", &K1_pionID, "K1_pionID/D");
  LowMod->Branch("K1_binaryID", &K1_binaryID, "K1_binaryID/D");
  LowMod->Branch("K1_charge", &K1_charge, "K1_charge/D");
  LowMod->Branch("K1_p", &K1_p, "K1_p/D");
  LowMod->Branch("K1_E", &K1_E, "K1_E/D");
  LowMod->Branch("K1_dr", &K1_dr, "K1_dr/D");
  LowMod->Branch("K1_dz", &K1_dz, "K1_dz/D");
  LowMod->Branch("K1_cosTheta", &K1_cosTheta, "K1_cosTheta/D");
  //LowMod->Branch("K2_kaonID", &K2_kaonID, "K2_kaonID/D");
  //LowMod->Branch("K2_pionID", &K2_pionID, "K2_pionID/D");
  LowMod->Branch("K2_binaryID", &K2_binaryID, "K2_binaryID/D");
  LowMod->Branch("K2_charge", &K2_charge, "K2_charge/D");
  LowMod->Branch("K2_p", &K2_p, "K2_p/D");
  LowMod->Branch("K2_E", &K2_E, "K2_E/D");
  LowMod->Branch("K2_dr", &K2_dr, "K2_dr/D");
  LowMod->Branch("K2_dz", &K2_dz, "K2_dz/D");
  LowMod->Branch("K2_cosTheta", &K2_cosTheta, "K2_cosTheta/D");

  std::cout << "defined what to right to  the Tree as leafs" << std::endl;
  
  std::cout << "Starting for loop" << std::endl;
  for (Long64_t iEntry(0); iEntry < nEntries; ++iEntry){
    
        tc->GetEntry(iEntry);//taking the entree corresponding to that number
	if ( 2.23<Lambdac_M && Lambdac_M<2.35){
		K1_binaryID = K1_kaonID/(K1_pionID + K1_kaonID);
	        K2_binaryID = K2_kaonID/(K2_pionID + K2_kaonID);
		p_trinaryID = p_protonID/(p_pionID + p_kaonID + p_protonID);
        LowMod->Fill();//filling the tree using the variables
	   } 
  }
  std::cout << "End for loop" << std::endl;
  
  LowMod->Write(); //writing the tree TO THE root file
  std::cout << "Wrote the TREE to the root " << std::endl;
  file2->Close();
  
}

// root -l 'lc_cp_skim2.C("/b2diska/janaka/XicLc_MC15/ntuples/*.root","/b2diska/dinura/rootfiles/ntuples/MC/pkk_complete.root")'

