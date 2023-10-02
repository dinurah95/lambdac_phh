#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"

void skimLcSubset(TString infile, TString outfile){

  auto newfile = new TFile(outfile,"RECREATE");

  TChain* tc = new TChain("lcp_pkk");
  tc->Add(infile);

  // disable all branches
  tc->SetBranchStatus("*",0);

  // turn on the ones we want
  tc->SetBranchStatus("eventRandom");

  tc->SetBranchStatus("Lambdac_isSignal");
  tc->SetBranchStatus("Lambdac_M");
  tc->SetBranchStatus("Lambdac_cosTheta_cms");
  tc->SetBranchStatus("Lambdac_flightDistance");
  tc->SetBranchStatus("Lambdac_significanceOfDistance");
  tc->SetBranchStatus("Lambdac_dcosTheta");
  tc->SetBranchStatus("Lambdac_cosAngleBetweenMomentumAndVertexVector");
  tc->SetBranchStatus("Lambdac_p_cms");
  tc->SetBranchStatus("Lambdac_vtxChi2");

  //new variables
  tc->SetBranchStatus("Lambdac_px");
  tc->SetBranchStatus("Lambdac_py");
  tc->SetBranchStatus("Lambdac_pz");
  tc->SetBranchStatus("Lambdac_pt");
  tc->SetBranchStatus("Lambdac_p");
  tc->SetBranchStatus("Lambdac_dr");
  tc->SetBranchStatus("Lambdac_dcosTheta");
  tc->SetBranchStatus("Lambdac_mcDecayTime");

  tc->SetBranchStatus("p_nCDCHits");
  tc->SetBranchStatus("p_firstPXDLayer");
  tc->SetBranchStatus("p_firstSVDLayer");
  tc->SetBranchStatus("p_charge");
  tc->SetBranchStatus("p_protonID_noSVD");
  tc->SetBranchStatus("p_kaonID_noSVD");
  tc->SetBranchStatus("p_pionID_noSVD");

  //new variables
  tc->SetBranchStatus("p_px");
  tc->SetBranchStatus("p_py");
  tc->SetBranchStatus("p_pz");
  tc->SetBranchStatus("p_pt");
  tc->SetBranchStatus("p_p");
  tc->SetBranchStatus("p_dr");
  tc->SetBranchStatus("p_nPXDHits");
  tc->SetBranchStatus("p_nSVDHits");

  tc->SetBranchStatus("K1_nCDCHits");
  tc->SetBranchStatus("K1_firstPXDLayer");
  tc->SetBranchStatus("K1_firstSVDLayer");
  tc->SetBranchStatus("K1_kaonID_noSVD");
  tc->SetBranchStatus("K1_pionID_noSVD");

  //new variables
  tc->SetBranchStatus("K1_px");
  tc->SetBranchStatus("K1_py");
  tc->SetBranchStatus("K1_pz");
  tc->SetBranchStatus("K1_pt");
  tc->SetBranchStatus("K1_p");
  tc->SetBranchStatus("K1_dr");
  tc->SetBranchStatus("K1_nPXDHits");
  tc->SetBranchStatus("K1_nSVDHits");

  tc->SetBranchStatus("K2_nCDCHits");
  tc->SetBranchStatus("K2_firstPXDLayer");
  tc->SetBranchStatus("K2_firstSVDLayer");
  tc->SetBranchStatus("K2_kaonID_noSVD");
  tc->SetBranchStatus("K2_pionID_noSVD");

  //new variables
  tc->SetBranchStatus("K2_px");
  tc->SetBranchStatus("K2_py");
  tc->SetBranchStatus("K2_pz");
  tc->SetBranchStatus("K2_pt");
  tc->SetBranchStatus("K2_p");
  tc->SetBranchStatus("K2_dr");
  tc->SetBranchStatus("K2_nPXDHits");
  tc->SetBranchStatus("K2_nSVDHits");

  auto newtree = tc->CloneTree(0);
  newtree->CopyEntries(tc);
  delete tc;
}

  /*TChain* tc2 = new TChain("lcp_ppipi");
  tc2->Add(infile);

  // disable all branches
  tc2->SetBranchStatus("*",0);

  // turn on the ones we want
  tc2->SetBranchStatus("eventRandom");

  tc2->SetBranchStatus("Lambdac_isSignal");
  tc2->SetBranchStatus("Lambdac_M");
  tc2->SetBranchStatus("Lambdac_cosTheta_cms");
  tc2->SetBranchStatus("Lambdac_flightDistance");
  tc2->SetBranchStatus("Lambdac_flightDistanceErr");
  tc2->SetBranchStatus("Lambdac_significanceOfDistance");
  tc2->SetBranchStatus("Lambdac_dcosTheta");
  tc2->SetBranchStatus("Lambdac_cosAngleBetweenMomentumAndVertexVector");
  tc2->SetBranchStatus("Lambdac_p_cms");
  tc2->SetBranchStatus("Lambdac_vtxChi2");

  //new variables
  tc2->SetBranchStatus("Lambdac_px");
  tc2->SetBranchStatus("Lambdac_py");
  tc2->SetBranchStatus("Lambdac_pz");
  tc2->SetBranchStatus("Lambdac_pt");
  tc2->SetBranchStatus("Lambdac_p");
  tc2->SetBranchStatus("Lambdac_dr");
  tc2->SetBranchStatus("Lambdac_dcosTheta");
  tc2->SetBranchStatus("Lambdac_mcDecayTime");

  tc2->SetBranchStatus("p_nCDCHits");
  tc2->SetBranchStatus("p_firstPXDLayer");
  tc2->SetBranchStatus("p_firstSVDLayer");
  tc2->SetBranchStatus("p_charge");
  tc2->SetBranchStatus("p_protonID_noSVD");
  tc2->SetBranchStatus("p_kaonID_noSVD");
  tc2->SetBranchStatus("p_pionID_noSVD");

  //new variables
  tc2->SetBranchStatus("p_px");
  tc2->SetBranchStatus("p_py");
  tc2->SetBranchStatus("p_pz");
  tc2->SetBranchStatus("p_pt");
  tc2->SetBranchStatus("p_p");
  tc2->SetBranchStatus("p_dr");
  tc2->SetBranchStatus("p_nPXDHits");
  tc2->SetBranchStatus("p_nSVDHits");

  tc2->SetBranchStatus("pi1_nCDCHits");
  tc2->SetBranchStatus("pi1_firstPXDLayer");
  tc2->SetBranchStatus("pi1_firstSVDLayer");
  tc2->SetBranchStatus("pi1_kaonID_noSVD");
  tc2->SetBranchStatus("pi1_pionID_noSVD");
  
  //new variables
  tc2->SetBranchStatus("pi1_px");
  tc2->SetBranchStatus("pi1_py");
  tc2->SetBranchStatus("pi1_pz");
  tc2->SetBranchStatus("pi1_pt");
  tc2->SetBranchStatus("pi1_p");
  tc2->SetBranchStatus("pi1_dr");
  tc2->SetBranchStatus("pi1_nPXDHits");
  tc2->SetBranchStatus("pi1_nSVDHits");

  tc2->SetBranchStatus("pi2_nCDCHits");
  tc2->SetBranchStatus("pi2_firstPXDLayer");
  tc2->SetBranchStatus("pi2_firstSVDLayer");
  tc2->SetBranchStatus("pi2_kaonID_noSVD");
  tc2->SetBranchStatus("pi2_pionID_noSVD");

  //new variables
  tc2->SetBranchStatus("pi2_px");
  tc2->SetBranchStatus("pi2_py");
  tc2->SetBranchStatus("pi2_pz");
  tc2->SetBranchStatus("pi2_pt");
  tc2->SetBranchStatus("pi2_p");
  tc2->SetBranchStatus("pi2_dr");
  tc2->SetBranchStatus("pi2_nPXDHits");
  tc2->SetBranchStatus("pi2_nSVDHits");

  auto newtree2 = tc2->CloneTree(0);
  newtree2->CopyEntries(tc2);
  delete tc2;


  TChain* tc3 = new TChain("lcp_pkpi");
  tc3->Add(infile);

  // disable all branches
  tc3->SetBranchStatus("*",0);

  // turn on the ones we want
  tc3->SetBranchStatus("eventRandom");

  tc3->SetBranchStatus("Lambdac_isSignal");
  tc3->SetBranchStatus("Lambdac_M");
  tc3->SetBranchStatus("Lambdac_cosTheta_cms");
  tc3->SetBranchStatus("Lambdac_flightDistance");
  tc3->SetBranchStatus("Lambdac_flightDistanceErr");
  tc3->SetBranchStatus("Lambdac_significanceOfDistance");
  tc3->SetBranchStatus("Lambdac_dcosTheta");
  tc3->SetBranchStatus("Lambdac_cosAngleBetweenMomentumAndVertexVector");
  tc3->SetBranchStatus("Lambdac_p_cms");
  tc3->SetBranchStatus("Lambdac_vtxChi2");

  //new variables
  tc3->SetBranchStatus("Lambdac_px");
  tc3->SetBranchStatus("Lambdac_py");
  tc3->SetBranchStatus("Lambdac_pz");
  tc3->SetBranchStatus("Lambdac_pt");
  tc3->SetBranchStatus("Lambdac_p");
  tc3->SetBranchStatus("Lambdac_dr");
  tc3->SetBranchStatus("Lambdac_dcosTheta");
  tc3->SetBranchStatus("Lambdac_mcDecayTime");

  tc3->SetBranchStatus("p_nCDCHits");
  tc3->SetBranchStatus("p_firstPXDLayer");
  tc3->SetBranchStatus("p_firstSVDLayer");
  tc3->SetBranchStatus("p_charge");
  tc3->SetBranchStatus("p_protonID_noSVD");
  tc3->SetBranchStatus("p_kaonID_noSVD");
  tc3->SetBranchStatus("p_pionID_noSVD");

  //new variables
  tc3->SetBranchStatus("p_px");
  tc3->SetBranchStatus("p_py");
  tc3->SetBranchStatus("p_pz");
  tc3->SetBranchStatus("p_pt");
  tc3->SetBranchStatus("p_p");
  tc3->SetBranchStatus("p_dr");
  tc3->SetBranchStatus("p_nPXDHits");
  tc3->SetBranchStatus("p_nSVDHits");

  tc3->SetBranchStatus("K_nCDCHits");
  tc3->SetBranchStatus("K_firstPXDLayer");
  tc3->SetBranchStatus("K_firstSVDLayer");
  tc3->SetBranchStatus("K_kaonID_noSVD");
  tc3->SetBranchStatus("K_pionID_noSVD");

  //new variables
  tc3->SetBranchStatus("K_px");
  tc3->SetBranchStatus("K_py");
  tc3->SetBranchStatus("K_pz");
  tc3->SetBranchStatus("K_pt");
  tc3->SetBranchStatus("K_p");
  tc3->SetBranchStatus("K_dr");
  tc3->SetBranchStatus("K_nPXDHits");
  tc3->SetBranchStatus("K_nSVDHits");

  tc3->SetBranchStatus("pi_nCDCHits");
  tc3->SetBranchStatus("pi_firstPXDLayer");
  tc3->SetBranchStatus("pi_firstSVDLayer");
  tc3->SetBranchStatus("pi_kaonID_noSVD");
  tc3->SetBranchStatus("pi_pionID_noSVD");

  //new variables
  tc3->SetBranchStatus("pi_px");
  tc3->SetBranchStatus("pi_py");
  tc3->SetBranchStatus("pi_pz");
  tc3->SetBranchStatus("pi_pt");
  tc3->SetBranchStatus("pi_p");
  tc3->SetBranchStatus("pi_dr");
  tc3->SetBranchStatus("pi_nPXDHits");
  tc3->SetBranchStatus("pi_nSVDHits");

  auto newtree3 = tc3->CloneTree(0);
  newtree3->CopyEntries(tc3);
  delete tc3;

  newfile->Write();
  delete newfile;
  }*/

/*void skimSigmaSubset(TString infile, TString outfile){

  TChain* tc4 = new TChain("lcp_sigkpi");
  tc4->Add(infile);

  // disable all branches
  tc4->SetBranchStatus("*",0);

  // turn on the ones we want
  tc4->SetBranchStatus("eventRandom");

  tc4->SetBranchStatus("Lambdac_isSignal");
  tc4->SetBranchStatus("Lambdac_M");
  tc4->SetBranchStatus("Lambdac_cosTheta_cms");
  tc4->SetBranchStatus("Lambdac_flightDistance");
  tc4->SetBranchStatus("Lambdac_flightDistanceErr");
  tc4->SetBranchStatus("Lambdac_significanceOfDistance");
  tc4->SetBranchStatus("Lambdac_dcosTheta");
  tc4->SetBranchStatus("Lambdac_cosAngleBetweenMomentumAndVertexVector");
  tc4->SetBranchStatus("Lambdac_p_cms");
  tc4->SetBranchStatus("Lambdac_vtxChi2");

  tc4->SetBranchStatus("Sigma_p_nCDCHits");
  tc4->SetBranchStatus("Sigma_p_firstPXDLayer");
  tc4->SetBranchStatus("Sigma_p_firstSVDLayer");
  tc4->SetBranchStatus("Sigma_p_charge");
  tc4->SetBranchStatus("Sigma_p_protonID_noSVD");
  tc4->SetBranchStatus("Sigma_p_kaonID_noSVD");
  tc4->SetBranchStatus("Sigma_p_pionID_noSVD");

  tc4->SetBranchStatus("K_nCDCHits");
  tc4->SetBranchStatus("K_firstPXDLayer");
  tc4->SetBranchStatus("K_firstSVDLayer");
  tc4->SetBranchStatus("K_kaonID_noSVD");
  tc4->SetBranchStatus("K_pionID_noSVD");
  tc4->SetBranchStatus("pi_nCDCHits");
  tc4->SetBranchStatus("pi_firstPXDLayer");
  tc4->SetBranchStatus("pi_firstSVDLayer");
  tc4->SetBranchStatus("pi_kaonID_noSVD");
  tc4->SetBranchStatus("pi_pionID_noSVD");

  auto newtree4 = tc4->CloneTree(0);
  newtree4->CopyEntries(tc4);
  delete tc4;

}

void skimTruth(TString infile, TString outfile){

  TChain* tc = new TChain("xicp");
  tc->Add(infile);

  TFile* newfile = new TFile(outfile,"RECREATE");
  TTree* newtree = tc->CloneTree(0);

  double issig;
  tc->SetBranchAddress("Xic_isSignal",&issig);

  for( int i = 0; i < tc->GetEntries(); ++i ){
    tc->GetEntry(i);

    if( issig!=1 ) newtree->Fill();
  }

  newtree->AutoSave();
  delete tc;
  delete newfile;
}

void skimFullXicSample(TString infile, TString treename, TString outfile, bool tight=false){

  TChain* tc = new TChain(treename);
  tc->Add(infile);

  TFile* newfile = new TFile(outfile,"RECREATE");
  TTree* newtree = tc->CloneTree(0);

  double xicm, lcm, chiprob, pcms;
  double p_piid, p_kid, p_pid;
  double k_piid, k_kid;
  double p_nPXD, k_nPXD, pi_nPXD, p_nSVD, k_nSVD, pi_nSVD;
  double p_SVD, k_SVD, pi_SVD;
  tc->SetBranchAddress("Xic_M",&xicm);
  tc->SetBranchAddress("Xic_Lambdac_M",&lcm);
  tc->SetBranchAddress("Xic_chiProb",&chiprob);
  tc->SetBranchAddress("Xic_p_CMS",&pcms);
  tc->SetBranchAddress("Xic_p_pionID",&p_piid);
  tc->SetBranchAddress("Xic_p_kaonID",&p_kid);
  tc->SetBranchAddress("Xic_p_protonID",&p_pid);
  tc->SetBranchAddress("Xic_K_pionID",&k_piid);
  tc->SetBranchAddress("Xic_K_kaonID",&k_kid);
  tc->SetBranchAddress("Xic_p_nPXDHits",&p_nPXD);
  tc->SetBranchAddress("Xic_K_nPXDHits",&k_nPXD);
  tc->SetBranchAddress("Xic_pi_nPXDHits",&pi_nPXD);
  tc->SetBranchAddress("Xic_p_nSVDHits",&p_nSVD);
  tc->SetBranchAddress("Xic_K_nSVDHits",&k_nSVD);
  tc->SetBranchAddress("Xic_pi_nSVDHits",&pi_nSVD);
  tc->SetBranchAddress("Xic_p_firstSVDLayer",&p_SVD);
  tc->SetBranchAddress("Xic_K_firstSVDLayer",&k_SVD);
  tc->SetBranchAddress("Xic_pi_firstSVDLayer",&pi_SVD);

  double p_px, p_py, p_pz, p_E;
  double pi_px, pi_py, pi_pz, pi_E;
  double K_px, K_py, K_pz, K_E;

  tc->SetBranchAddress("Xic_p_px",&p_px);
  tc->SetBranchAddress("Xic_p_py",&p_py);
  tc->SetBranchAddress("Xic_p_pz",&p_pz);
  tc->SetBranchAddress("Xic_p_E",&p_E);

  tc->SetBranchAddress("Xic_pi_px",&pi_px);
  tc->SetBranchAddress("Xic_pi_py",&pi_py);
  tc->SetBranchAddress("Xic_pi_pz",&pi_pz);
  tc->SetBranchAddress("Xic_pi_E",&pi_E);

  tc->SetBranchAddress("Xic_K_px",&K_px);
  tc->SetBranchAddress("Xic_K_py",&K_py);
  tc->SetBranchAddress("Xic_K_pz",&K_pz);
  tc->SetBranchAddress("Xic_K_E",&K_E);

  double g1E9E21, g2E9E21, g1MVA, g2MVA;

  tc->SetBranchAddress("g1_clusterE9E21",&g1E9E21);
  tc->SetBranchAddress("g2_clusterE9E21",&g2E9E21);
  tc->SetBranchAddress("g1_clusterPulseShapeDiscriminationMVA",&g1MVA);
  tc->SetBranchAddress("g2_clusterPulseShapeDiscriminationMVA",&g2MVA);

  for( int i = 0; i < tc->GetEntries(); ++i ){
    tc->GetEntry(i);
    float cpe = sqrt(fabs(p_px*p_px+p_py*p_py+p_pz*p_pz+0.13957*0.13957));
    float cmass = sqrt(fabs((cpe+pi_E+K_E)*(cpe+pi_E+K_E)-(p_px+pi_px+K_px)*(p_px+pi_px+K_px)-(p_py+pi_py+K_py)*(p_py+pi_py+K_py)-(p_pz+pi_pz+K_pz)*(p_pz+pi_pz+K_pz)));

    float tri = p_pid/(p_pid+p_kid+p_piid);

    // loose
    if( !tight && chiprob > 0.001 && pcms > 2.5 ) newtree->Fill();

    // mod
    if( tight && (cmass < 1.858 || (cmass > 1.881 && cmass < 2.0) || cmass > 2.02) && chiprob > 0.01 && pcms > 2.5 && tri > 0.9 && k_kid > 0.6) newtree->Fill();

    // tight
    //    if( tight && (cmass < 1.858 || (cmass > 1.881 && cmass < 2.0) || cmass > 2.02) && chiprob > 0.01 && pcms > 2.5 && tri > 0.9 && k_kid > 0.6 && lcm>2.276 && lcm<2.298 && (xicm-lcm)<0.4 && g1E9E21>0.9 && g2E9E21>0.9 && g1MVA>0.03 && g2MVA>0.03) newtree->Fill();

  }

  newtree->AutoSave();
  delete tc;
  delete newfile;
}

void skimFullLcSample(TString infile, TString treename, TString outfile, bool tight=false){

  TChain* tc = new TChain(treename);
  tc->Add(infile);

  TFile* newfile = new TFile(outfile,"RECREATE");
  TTree* newtree = tc->CloneTree(0);

  double lambdacm, lcm, chiprob, pcms;
  double p_piid, p_kid, p_pid;
  double k1_piid, k1_kid, k2_piid, k2_kid;
  double p_nPXD, k_nPXD, pi_nPXD, p_nSVD, k_nSVD, pi_nSVD;
  double p_SVD, k_SVD, pi_SVD;
  tc->SetBranchAddress("Lambdac_M",&lambdacm);
  tc->SetBranchAddress("Lambdac_chiProb",&chiprob);
  tc->SetBranchAddress("Lambdac_p_cms",&pcms);
  tc->SetBranchAddress("p_pionID",&p_piid);
  tc->SetBranchAddress("p_kaonID",&p_kid);
  tc->SetBranchAddress("p_protonID",&p_pid);

  for( int i = 0; i < tc->GetEntries(); ++i ){
    tc->GetEntry(i);

    float tri = p_pid/(p_pid+p_kid+p_piid);

    // loose
    if( !tight && chiprob > 0.001 ) newtree->Fill();

    // tight
    //    if( tight && chiprob > 0.01 && pcms > 2.5 && tri > 0.9 && k1_kid > 0.6 && k2_kid > 0.6) newtree->Fill();
  }

  newtree->AutoSave();
  delete tc;
  delete newfile;
}


void skimXicSample(TString infile, TString treename, TString outfile, bool tight=false){

  TChain* tc = new TChain(treename);
  tc->Add(infile);

  TFile* newfile = new TFile(outfile,"RECREATE");
  TTree* newtree = tc->CloneTree(0);

  double xicm, lcm, chiprob, pcms;
  double p_piid, p_kid, p_pid;
  double k1_piid, k1_kid, k2_piid, k2_kid;
  double p_nPXD, k_nPXD, pi_nPXD, p_nSVD, k_nSVD, pi_nSVD;
  double p_SVD, k_SVD, pi_SVD;
  tc->SetBranchAddress("Xic_M",&xicm);
  tc->SetBranchAddress("Xic_chiProb",&chiprob);
  tc->SetBranchAddress("Xic_p_cms",&pcms);
  tc->SetBranchAddress("p_pionID",&p_piid);
  tc->SetBranchAddress("p_kaonID",&p_kid);
  tc->SetBranchAddress("p_protonID",&p_pid);

  for( int i = 0; i < tc->GetEntries(); ++i ){
    tc->GetEntry(i);

    float tri = p_pid/(p_pid+p_kid+p_piid);

    // loose
    if( !tight && chiprob > 0.001 ) newtree->Fill();

    // tight
    //    if( tight && chiprob > 0.01 && pcms > 2.5 && tri > 0.9 && k1_kid > 0.6 && k2_kid > 0.6) newtree->Fill();
  }

  newtree->AutoSave();
  delete tc;
  delete newfile;
  }*/

void skim(){
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_uu_230613/sub00/*.root","lcp.uu_00.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_uu_230613/sub01/*.root","lcp.uu_01.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_uu_230613/sub02/*.root","lcp.uu_02.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_uu_230613/sub03/*.root","lcp.uu_03.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_uu_230613/sub04/*.root","lcp.uu_04.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_uu_230613/sub05/*.root","lcp.uu_05.root");
  

  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_dd_230613/sub00/*.root","lcp.dd_00.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_dd_230613/sub01/*.root","lcp.dd_01.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_dd_230613/sub02/*.root","lcp.dd_02.root");
  
  
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_ss_230613/sub00/*.root","lcp.ss_00.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_ss_230613/sub01/*.root","lcp.ss_01.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_ss_230613/sub02/*.root","lcp.ss_02.root");
  
  
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_cc_230613/sub00/*.root","lcp.cc_00.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_cc_230613/sub01/*.root","lcp.cc_01.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_cc_230613/sub02/*.root","lcp.cc_02.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_cc_230613/sub03/*.root","lcp.cc_03.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_cc_230613/sub04/*.root","lcp.cc_04.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_cc_230613/sub05/*.root","lcp.cc_05.root");

  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_tp_230613/sub00/*.root","lcp.tp_00.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_tp_230613/sub01/*.root","lcp.tp_01.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_tp_230613/sub02/*.root","lcp.tp_02.root");

  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_mi_230613/sub00/*.root","lcp.mi_00.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_mi_230613/sub01/*.root","lcp.mi_01.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_mi_230613/sub02/*.root","lcp.mi_02.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_mi_230613/sub03/*.root","lcp.mi_03.root");

  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_ch_230613/sub00/*.root","lcp.ch_00.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_ch_230613/sub01/*.root","lcp.ch_01.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_ch_230613/sub02/*.root","lcp.ch_02.root");
  skimLcSubset("/belle2work/BelleII/Xic2Sigmahh/all_ntuples/XicAcp_prompt_ch_230613/sub03/*.root","lcp.ch_03.root");
  
}
