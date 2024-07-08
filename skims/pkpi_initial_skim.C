//This program gives the output root file with variables for Acp calculation
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "iostream"

void lc_cp_skim3 (TString input, TString output) {

  std::cout << "Starting" << std::endl;
  TChain* tc = new TChain("lcp_pkpi");
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
  Double_t p_cosTheta;
  Double_t p_protonID;
  Double_t p_kaonID;
  Double_t p_pionID;
  Double_t p_trinaryID;
  //Double_t p_cosTheta_cms;
  Double_t p_p;
  Double_t p_E;
  Double_t p_dr;
  Double_t p_dz;
  Double_t p_charge;
  Double_t K_kaonID;
  Double_t K_pionID;
  Double_t K_binaryID;
  Double_t K_charge;
  Double_t K_p;
  Double_t K_E;
  Double_t K_dr;
  Double_t K_dz;
  Double_t K_cosTheta;
  Double_t pi_charge;
  Double_t pi_p;
  Double_t pi_E;
  Double_t pi_dz;
  Double_t pi_dr;
  Double_t pi_cosTheta;

  std::cout << "Defined Varables to get from input root file" << std::endl;

  // Taking those variables from the root file to the tc tree
  tc->SetBranchAddress("eventRandom", &eventRandom);
  tc->SetBranchAddress("Lambdac_isSignal", &Lambdac_isSignal);
  tc->SetBranchAddress("Lambdac_M", &Lambdac_M);
  tc->SetBranchAddress("Lambdac_cosTheta_cms", &Lambdac_cosTheta_cms);
  tc->SetBranchAddress("Lambdac_flightDistance", &Lambdac_flightDistance);
  tc->SetBranchAddress("Lambdac_significanceOfDistance", &Lambdac_significanceOfDistance);
  tc->SetBranchAddress("Lambdac_p_cms", &Lambdac_p_cms);
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
  tc->SetBranchAddress("K_kaonID", &K_kaonID);
  tc->SetBranchAddress("K_pionID", &K_pionID);
  tc->SetBranchAddress("K_charge", &K_charge);
  tc->SetBranchAddress("K_p", &K_p);
  tc->SetBranchAddress("K_E", &K_E);
  tc->SetBranchAddress("K_dr", &K_dr);
  tc->SetBranchAddress("K_dz", &K_dz);
  tc->SetBranchAddress("K_cosTheta", &K_cosTheta);
  tc->SetBranchAddress("pi_charge", &pi_charge);
  tc->SetBranchAddress("pi_p", &pi_p);
  tc->SetBranchAddress("pi_E", &pi_E);
  tc->SetBranchAddress("pi_dr", &pi_dr);
  tc->SetBranchAddress("pi_dz", &pi_dz);
  tc->SetBranchAddress("pi_cosTheta", &pi_cosTheta);


  std::cout << "Variables taken" << std::endl;

  // taking the total number of entries from my Bplus tree
  Long64_t nEntries(tc->GetEntries());
  
  std::cout<< "Got Number of Entries in tree = "<< nEntries<< std::endl;
  
  //New File and Tree for storing 
  TFile *file2 = new TFile (output,"RECREATE");
  TTree* LowMod = new TTree("lcp_pkpi", "a new tree");//TTree *tree = new TTree(name, title) 
  std::cout << "New output root file defined" << std::endl;

  //define what to right to  the Tree as leafs
  LowMod->Branch("eventRandom", &eventRandom, "eventRandom/D");
  LowMod->Branch("Lambdac_isSignal", &Lambdac_isSignal, "Lambdac_isSignal/D");
  LowMod->Branch("Lambdac_M", &Lambdac_M, "Lambdac_M/D");
  LowMod->Branch("Lambdac_cosTheta_cms", &Lambdac_cosTheta_cms, "Lambdac_cosTheta_cms/D");
  LowMod->Branch("Lambdac_flightDistance", &Lambdac_flightDistance, "Lambdac_flightDistance/D");
  LowMod->Branch("Lambdac_significanceOfDistance", &Lambdac_significanceOfDistance, "Lambdac_significanceOfDistance/D");
  LowMod->Branch("Lambdac_p_cms", &Lambdac_p_cms, "Lambdac_p_cms/D");
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
  //LowMod->Branch("K_kaonID", &K_kaonID, "K_kaonID/D");
  //LowMod->Branch("K_pionID", &K_pionID, "K_pionID/D");
  LowMod->Branch("K_binaryID", &K_binaryID, "K_binaryID/D");
  LowMod->Branch("K_charge", &K_charge, "K_charge/D");
  LowMod->Branch("K_p", &K_p, "K_p/D");
  LowMod->Branch("K_E", &K_E, "K_E/D");
  LowMod->Branch("K_dr", &K_dr, "K_dr/D");
  LowMod->Branch("K_dz", &K_dz, "K_dz/D");
  LowMod->Branch("K_cosTheta", &K_cosTheta, "K_cosTheta/D");
  LowMod->Branch("pi_charge", &pi_charge, "pi_charge/D");
  LowMod->Branch("pi_p", &pi_p, "pi_p/D");
  LowMod->Branch("pi_E", &pi_E, "pi_E/D");
  LowMod->Branch("pi_dr", &pi_dr, "pi_dr/D");
  LowMod->Branch("pi_dz", &pi_dz, "pi_dz/D");
  LowMod->Branch("pi_cosTheta", &pi_cosTheta, "pi_cosTheta/D");

  std::cout << "defined what to right to  the Tree as leafs" << std::endl;
  
  std::cout << "Starting for loop" << std::endl;
  for (Long64_t iEntry(0); iEntry < nEntries; ++iEntry){
    
        tc->GetEntry(iEntry);//taking the entree corresponding to that number
	if ( 2.23<Lambdac_M && Lambdac_M<2.35){
	K_binaryID = K_kaonID/(K_pionID + K_kaonID);
	p_trinaryID = p_protonID/(p_pionID + p_kaonID + p_protonID);
        LowMod->Fill();//filling the tree using the variables
	   } 
  } 
  std::cout << "End for loop" << std::endl;


  
  LowMod->Write(); //writing the tree TO THE root file
  std::cout << "Wrote the TREE to the root " << std::endl;
  file2->Close();
  
}


// root -l 'lc_cp_skim3.C("/b2diska/janaka/XicLc_MC15/ntuples/*.root","/b2diska/dinura/rootfiles/ntuples/MC/pkpi_complete.root")'
