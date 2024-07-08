//This program gives the output root file with variables for Acp calculation
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "iostream"

void lc_cp_skim1 (TString input, TString output) {

  std::cout << "Starting" << std::endl;
  TChain* tc = new TChain("lcp_ppipi");
  tc->Add(input);

// TTree* tc = (TTree*) file1-> Get ("lambdatree") ;//taking tree from the input root file 
  
  std::cout << "Got tree from file" << std::endl;

// define variables that need to take from the root file
  Double_t eventRandom;
  Double_t Signal;
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
  //Double_t p_flightDistance;
  //Double_t p_significanceOfDistance;
  //Double_t p_cosTheta_cms;
  Double_t p_cosTheta;
  Double_t p_charge;
  Double_t p_dz;
  Double_t p_dr;
  Double_t pi1_charge;
  Double_t pi2_charge;
  Double_t pi1_p;
  Double_t pi2_p;
  Double_t pi1_E;
  Double_t pi2_E;
  Double_t pi1_cosTheta;
  Double_t pi2_cosTheta;
  Double_t pi1_dz;
  Double_t pi2_dz;
  Double_t pi1_dr;
  Double_t pi2_dr;
  Double_t pi1_p_cms;
  Double_t pi2_p_cms;


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
  tc->SetBranchAddress("p_protonID", &p_protonID);
  tc->SetBranchAddress("p_pionID", &p_pionID);
  tc->SetBranchAddress("p_kaonID", &p_kaonID);
  //tc->SetBranchAddress("p_cosTheta_cms", &p_cosTheta_cms); 
  tc->SetBranchAddress("p_cosTheta", &p_cosTheta);
  tc->SetBranchAddress("p_p", &p_p);
  tc->SetBranchAddress("p_E", &p_E);
  //tc->SetBranchAddress("p_flightDistance", &p_flightDistance);
  //tc->SetBranchAddress("p_significanceOfDistance", &p_significanceOfDistance);  
  tc->SetBranchAddress("p_dz", &p_dz);
  tc->SetBranchAddress("p_dr", &p_dr);
  tc->SetBranchAddress("pi1_charge", &pi1_charge);
  tc->SetBranchAddress("pi2_charge", &pi2_charge);
  tc->SetBranchAddress("pi1_p", &pi1_p);
  tc->SetBranchAddress("pi2_p", &pi2_p);
  tc->SetBranchAddress("pi1_E", &pi1_E);
  tc->SetBranchAddress("pi2_E", &pi2_E);
  tc->SetBranchAddress("pi1_cosTheta", &pi1_cosTheta);
  tc->SetBranchAddress("pi2_cosTheta", &pi2_cosTheta);
  tc->SetBranchAddress("pi1_dr", &pi1_dr);
  tc->SetBranchAddress("pi2_dr", &pi2_dr);
  tc->SetBranchAddress("pi1_dz", &pi1_dz);
  tc->SetBranchAddress("pi2_dz", &pi2_dz);
  tc->SetBranchAddress("pi1_p_cms", &pi1_p_cms);
  tc->SetBranchAddress("pi2_p_cms", &pi2_p_cms);

  std::cout << "Variables taken" << std::endl;

  // taking the total number of entries from my Bplus tree
  Long64_t nEntries(tc->GetEntries());
  
  std::cout<< "Got Number of Entries in tree = "<< nEntries<< std::endl;
  
  //New File and Tree for storing 
  TFile *file2 = new TFile (output,"RECREATE");
  TTree* LowMod = new TTree("lcp_ppipi", "a new tree");//TTree *tree = new TTree(name, title) 
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
  //LowMod->Branch("p_protonID", &p_protonID, "p_protonID/D");
  //LowMod->Branch("p_kaonID", &p_kaonID, "p_kaonID/D");
  //LowMod->Branch("p_pionID", &p_pionID, "p_pionID/D");
  LowMod->Branch("p_trinaryID", &p_trinaryID, "p_trinaryID/D");
  //LowMod->Branch("p_cosTheta_cms", &p_cosTheta_cms, "p_cosTheta_cms/D");
  LowMod->Branch("p_cosTheta", &p_cosTheta, "p_cosTheta/D");
  LowMod->Branch("p_p", &p_p, "p_p/D");
  LowMod->Branch("p_E", &p_E, "p_E/D");
  //LowMod->Branch("p_flightDistance", &p_flightDistance, "p_flightDistance/D");
  //LowMod->Branch("p_significanceOfDistance", &p_significanceOfDistance, "p_significanceOfDistance/D");
  LowMod->Branch("p_dr", &p_dr, "p_dr/D");
  LowMod->Branch("p_dz", &p_dz, "p_dz/D");
  LowMod->Branch("pi1_charge", &pi1_charge, "pi1_charge/D");
  LowMod->Branch("pi2_charge", &pi2_charge, "pi2_charge/D");  
  LowMod->Branch("pi1_p", &pi1_p, "pi1_p/D");
  LowMod->Branch("pi2_p", &pi2_p, "pi2_p/D");
  LowMod->Branch("pi1_E", &pi1_E, "pi1_E/D");
  LowMod->Branch("pi2_E", &pi2_E, "pi2_E/D");
  LowMod->Branch("pi1_cosTheta", &pi1_cosTheta, "pi1_cosTheta/D");
  LowMod->Branch("pi2_cosTheta", &pi2_cosTheta, "pi2_cosTheta/D");
  LowMod->Branch("pi1_dr", &pi1_dr, "pi1_dr/D");
  LowMod->Branch("pi2_dr", &pi2_dr, "pi2_dr/D");
  LowMod->Branch("pi1_dz", &pi1_dz, "pi1_dz/D");
  LowMod->Branch("pi2_dz", &pi2_dz, "pi2_dz/D");
  LowMod->Branch("pi1_p_cms", &pi1_p_cms, "pi1_p_cms/D");
  LowMod->Branch("pi2_p_cms", &pi2_p_cms, "pi2_p_cms/D");

  std::cout << "defined what to right to  the Tree as leafs" << std::endl;
  
  std::cout << "Starting for loop" << std::endl;
  for (Long64_t iEntry(0); iEntry < nEntries; ++iEntry){
    
        tc->GetEntry(iEntry);//taking the entree corresponding to that number
	if ( 2.23<Lambdac_M && Lambdac_M<2.35){
	  p_trinaryID = p_protonID/(p_pionID + p_kaonID + p_protonID);
        LowMod->Fill();//filling the tree using the variables
	   } 
  }

  std::cout << "End for loop" << std::endl;


  
  LowMod->Write(); //writing the tree TO THE root file
  std::cout << "Wrote the TREE to the root " << std::endl;
  file2->Close();
  
}


// root -l 'lc_cp_skim1.C("/b2diska/janaka/XicLc_MC15/ntuples/*.root","/b2diska/dinura/rootfiles/ntuples/MC/ppipi_complete.root")'
