//This program is for best candidate selections by assign a new variable value for Ks
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "iostream"

void cosTheta_CMS_check (TString input) {

  std::cout << "Starting" << std::endl;
  TChain* T_Bplus1 = new TChain("D0tree");
  T_Bplus1->Add(input);

  std::cout << "Got tree from file" << std::endl;

// define variables that need to take from the root file

 Double_t isSignal;
 Double_t beamE;
 Double_t beamPx;
 Double_t beamPy;
 Double_t beamPz;
 Double_t cosTheta;
 Double_t CMS_cosTheta;
 Double_t M;
 Double_t P;

 Double_t mycosTheta;
 Double_t boost;
 
// Taking those variables from the root file to the T_Bplus tree
 std::cout << "Defined Varables to get from input root file" << std::endl;

 T_Bplus1->SetBranchAddress("isSignal", &isSignal);
 T_Bplus1->SetBranchAddress("beamE", &beamE);
 T_Bplus1->SetBranchAddress("beamPx", &beamPx);
 T_Bplus1->SetBranchAddress("beamPy", &beamPy);
 T_Bplus1->SetBranchAddress("beamPz", &beamPz);
 T_Bplus1->SetBranchAddress("cosTheta", &cosTheta);
 T_Bplus1->SetBranchAddress("CMS_cosTheta", &CMS_cosTheta);
 T_Bplus1->SetBranchAddress("M", &M);
 T_Bplus1->SetBranchAddress("P", &P);

  std::cout << "Variables taken" << std::endl;

// taking the total number of entries from my Bplus tree
  Long64_t nEntries1(T_Bplus1->GetEntries());
  
  std::cout<< "Got Number of Entries in tree = "<< nEntries1<< std::endl;
   
//New File and Tree for storing 
  TFile *file2 = new TFile ("/belle2work/dinurah/Acp_Lambdac/D0_MC/costheta_check.root","RECREATE");
  TTree* LowMod1 = new TTree("myD0tree", "a new tree");//TTree *tree = new TTree(name, title) 
 
  std::cout << "New output root file defined" << std::endl;

//define what to right to  the Tree as leafs

 LowMod1->Branch("isSignal", &isSignal,"isSignal");
 LowMod1->Branch("beamE", &beamE,"beamE/D");
 LowMod1->Branch("beamPx", &beamPx,"beamPx/D");
 LowMod1->Branch("beamPy", &beamPy,"beamPy/D");
 LowMod1->Branch("beamPz", &beamPz,"beamPz/D");
 LowMod1->Branch("cosTheta", &cosTheta,"cosTheta/D");
 LowMod1->Branch("CMS_cosTheta", &CMS_cosTheta,"CMS_cosTheta/D");
 LowMod1->Branch("M", &M,"M/D");
 LowMod1->Branch("P", &P,"P");

 LowMod1->Branch("mycosTheta", &mycosTheta,"mycosTheta/D"); // cosTheta_CMS distribution made by myself
 LowMod1->Branch("boost", &boost,"boost");
 
 std::cout << "defined what to right to  the Tree as leafs" << std::endl;
  
  std::cout << "Starting for loop" << std::endl;
  for (Long64_t iEntry(0); iEntry < nEntries1; ++iEntry){
    
        T_Bplus1->GetEntry(iEntry);//taking the entree corresponding to that number
		if ( 1.75 < M && M < 1.95){
		  boost = beamPz/(sqrt(beamPx*beamPx + beamPy*beamPy + beamPz*beamPz));
		mycosTheta = cosTheta * boost;
		}
		LowMod1->Fill();//filling the tree using the variables
	  	 } 

    std::cout << "End for loop" << std::endl;
  
  LowMod1->Write(); //writing the tree TO THE root file
 
  std::cout << "Wrote the TREE to the root " << std::endl;
  file2->Close();
  
}

// root -l 'cosTheta_CMS_check.C("/belle2work/dinurah/Acp_Lambdac/D0_MC/out_CFT.root")'
