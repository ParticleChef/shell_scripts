//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 17 20:40:59 2015 by ROOT version 5.34/19
// from TTree t/t
// found on file: SingleEle_ntuple_1.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
#include <iostream>
#include <fstream>
#include <string>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include "/u/user/particlechef/rate/withEM_Tower.h"
#include "/u/user/particlechef/rate/withoutEM_SingleCrys.h"
#include "/u/user/particlechef/rate/RegionOfInterest.h"

using namespace std;
class test {
private:

map<TString, TH1*> maphist;

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           bunchN;
   vector<int>     *pileup;
   Float_t         beamSpotX0;
   Float_t         beamSpotY0;
   Float_t         beamSpotZ0;
   Float_t         beamSpotX0Error;
   Float_t         beamSpotY0Error;
   Float_t         beamSpotZ0Error;
   Float_t         beamWidthX;
   Float_t         beamWidthY;
   Float_t         beamSigmaZ;
   Float_t         beamWidthXError;
   Float_t         beamWidthYError;
   Float_t         beamSigmaZError;
   Int_t           genPartN;
   vector<float>   *genPartE;
   vector<float>   *genPartPt;
   vector<float>   *genPartEta;
   vector<float>   *genPartPhi;
   vector<int>     *genPartCharge;
   vector<int>     *genPartId;
   vector<float>   *propgenPartE;
   vector<float>   *propgenPartPt;
   vector<float>   *propgenPartEta;
   vector<float>   *propgenPartPhi;
   vector<int>     *propgenPartCharge;
   vector<int>     *propgenPartId;
   Int_t           simTrkN;
   vector<float>   *simTrkPt;
   vector<float>   *simTrkEta;
   vector<float>   *simTrkPhi;
   vector<int>     *simTrkId;
   vector<int>     *simTrkType;
   vector<float>   *simTrkVx;
   vector<float>   *simTrkVy;
   vector<float>   *simTrkVz;
   vector<float>   *simVx;
   vector<float>   *simVy;
   vector<float>   *simVz;
   Float_t         gammaBrem;
   Int_t           hardGamma;
   Int_t           egCrysN;
   vector<float>   *egCrysE;
   vector<float>   *egCrysEt;
   vector<float>   *egCrysEta;
   vector<float>   *egCrysPhi;
   vector<float>   *egCrysGx;
   vector<float>   *egCrysGy;
   vector<float>   *egCrysGz;
   vector<int>     *egCrysCharge;
   Int_t           l1tkegN;
   vector<float>   *l1tkegEt;
   vector<float>   *l1tkegEta;
   vector<float>   *l1tkegPhi;
   Int_t           l1tkegIsoN;
   vector<float>   *l1tkegIsoEt;
   vector<float>   *l1tkegIsoEta;
   vector<float>   *l1tkegIsoPhi;
   Int_t           egN;
   vector<float>   *egE;
   vector<float>   *egEt;
   vector<float>   *egEta;
   vector<float>   *egPhi;
   vector<float>   *egGx;
   vector<float>   *egGy;
   vector<float>   *egGz;
   vector<float>   *egclusterBS;
   vector<float>   *egclusterPt;
   vector<float>   *egclusterPtCore;
   vector<float>   *egclusterEta;
   vector<float>   *egclusterPhi;
   vector<double>  *egHoE;
   vector<double>  *clusterHoE;
   vector<double>  *seedHoE;
   vector<bool>    *clusterIsEG;
   vector<bool>    *seedIsEG;
   vector<bool>    *TriggerTower_IsBarrel;
   vector<double>  *seedHadEt;
   vector<double>  *seedEmEt;
   vector<double>  *HadEt;
   vector<double>  *EmEt;
   Int_t           fHitN;
   vector<int>     *fHitDisk;
   vector<int>     *fHitBlade;
   vector<int>     *fHitSide;
   vector<float>   *fHitGx;
   vector<float>   *fHitGy;
   vector<float>   *fHitGz;
   vector<int>     *fClSize;
   vector<int>     *fClSizeX;
   vector<int>     *fClSizeY;
   Int_t           bHitN;
   vector<int>     *bHitLayer;
   vector<int>     *bHitLadder;
   vector<float>   *bHitGx;
   vector<float>   *bHitGy;
   vector<float>   *bHitGz;
   vector<int>     *bClSize;
   vector<int>     *bClSizeX;
   vector<int>     *bClSizeY;
   

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bunchN;   //!
   TBranch        *b_pileup;   //!
   TBranch        *b_beamSpotX0;   //!
   TBranch        *b_beamSpotY0;   //!
   TBranch        *b_beamSpotZ0;   //!
   TBranch        *b_beamSpotX0Error;   //!
   TBranch        *b_beamSpotY0Error;   //!
   TBranch        *b_beamSpotZ0Error;   //!
   TBranch        *b_beamWidthX;   //!
   TBranch        *b_beamWidthY;   //!
   TBranch        *b_beamSigmaZ;   //!
   TBranch        *b_beamWidthXError;   //!
   TBranch        *b_beamWidthYError;   //!
   TBranch        *b_beamSigmaZError;   //!
   TBranch        *b_genPartN;   //!
   TBranch        *b_genPartE;   //!
   TBranch        *b_genPartPt;   //!
   TBranch        *b_genPartEta;   //!
   TBranch        *b_genPartPhi;   //!
   TBranch        *b_genPartCharge;   //!
   TBranch        *b_genPartId;   //!
   TBranch        *b_propgenPartE;   //!
   TBranch        *b_propgenPartPt;   //!
   TBranch        *b_propgenPartEta;   //!
   TBranch        *b_propgenPartPhi;   //!
   TBranch        *b_propgenPartCharge;   //!
   TBranch        *b_propgenPartId;   //!
   TBranch        *b_simTrkN;   //!
   TBranch        *b_simTrkPt;   //!
   TBranch        *b_simTrkEta;   //!
   TBranch        *b_simTrkPhi;   //!
   TBranch        *b_simTrkId;   //!
   TBranch        *b_simTrkType;   //!
   TBranch        *b_simTrkVx;   //!
   TBranch        *b_simTrkVy;   //!
   TBranch        *b_simTrkVz;   //!
   TBranch        *b_simVx;   //!
   TBranch        *b_simVy;   //!
   TBranch        *b_simVz;   //!
   TBranch        *b_gammaBrem;   //!
   TBranch        *b_hardGamma;   //!
   TBranch        *b_egCrysN;   //!
   TBranch        *b_egCrysE;   //!
   TBranch        *b_egCrysEt;   //!
   TBranch        *b_egCrysEta;   //!
   TBranch        *b_egCrysPhi;   //!
   TBranch        *b_egCrysGx;   //!
   TBranch        *b_egCrysGy;   //!
   TBranch        *b_egCrysGz;   //!
   TBranch        *b_egCrysCharge;   //!
   TBranch        *b_l1tkegN;   //!
   TBranch        *b_l1tkegEt;   //!
   TBranch        *b_l1tkegEta;   //!
   TBranch        *b_l1tkegPhi;   //!
   TBranch        *b_l1tkegIsoN;   //!
   TBranch        *b_l1tkegIsoEt;   //!
   TBranch        *b_l1tkegIsoEta;   //!
   TBranch        *b_l1tkegIsoPhi;   //!
   TBranch        *b_egN;   //!
   TBranch        *b_egE;   //!
   TBranch        *b_egEt;   //!
   TBranch        *b_egEta;   //!
   TBranch        *b_egPhi;   //!
   TBranch        *b_egGx;   //!
   TBranch        *b_egGy;   //!
   TBranch        *b_egGz;   //!
   TBranch        *b_egclusterBS;   //!
   TBranch        *b_egclusterPt;   //!
   TBranch        *b_egclusterPtCore;   //!
   TBranch        *b_egclusterEta;   //!
   TBranch        *b_egclusterPhi;   //!
   TBranch        *b_egHoE;   //!
   TBranch        *b_clusterHoE;   //!
   TBranch        *b_seedHoE;   //!
   TBranch        *b_clusterIsEG;   //!
   TBranch        *b_seedIsEG;   //!
   TBranch        *b_TriggerTower_IsBarrel;   //!
   TBranch        *b_seedHadEt;   //!
   TBranch        *b_seedEmEt;   //!
   TBranch        *b_HadEt;   //!
   TBranch        *b_EmEt;   //!
   TBranch        *b_fHitN;   //!
   TBranch        *b_fHitDisk;   //!
   TBranch        *b_fHitBlade;   //!
   TBranch        *b_fHitSide;   //!
   TBranch        *b_fHitGx;   //!
   TBranch        *b_fHitGy;   //!
   TBranch        *b_fHitGz;   //!
   TBranch        *b_fClSize;   //!
   TBranch        *b_fClSizeX;   //!
   TBranch        *b_fClSizeY;   //!
   TBranch        *b_bHitN;   //!
   TBranch        *b_bHitLayer;   //!
   TBranch        *b_bHitLadder;   //!
   TBranch        *b_bHitGx;   //!
   TBranch        *b_bHitGy;   //!
   TBranch        *b_bHitGz;   //!
   TBranch        *b_bClSize;   //!
   TBranch        *b_bClSizeX;   //!
   TBranch        *b_bClSizeY;   //!



   int Ele, Pos;
   int skip;
   float EgN;
   int eta_region;

   int withoutEM_count_Ele, withEM_count_Ele;
   double pass_count, pass_count_EleorPos, pass_count_Ele, pass_count_Pos;
   int pass_count_wo4thPix, pass_count_wo3thPix, pass_count_wo2thPix, pass_count_wo1thPix;
   int woEM_pass_count_wo4thPix, woEM_pass_count_wo3thPix, woEM_pass_count_wo2thPix, woEM_pass_count_wo1thPix;
   int wEM_pass_count_wo4thPix, wEM_pass_count_wo3thPix, wEM_pass_count_wo2thPix, wEM_pass_count_wo1thPix;



   double all_cut_pass_eg;
   int all_cut_pass_Ele, withoutEM_pass_Ele, withEM_pass_Ele, withoutEM_pass_Pos, withEM_pass_Pos;
   int all_cut_pass_Pos;
   int fourth_layer_missing;
   int third_layer_missing;
   int second_layer_missing;
   int first_layer_missing;

   double  L1_Dphi_cut1, L1_Dphi_cut2;
   double  L2_Dphi_cut1, L2_Dphi_cut2;
   double  L3_Dphi_cut1, L3_Dphi_cut2;
   double  L4_Dphi_cut1, L4_Dphi_cut2;
   double  D1_Dphi_cut1, D1_Dphi_cut2;
   double  D2_Dphi_cut1, D2_Dphi_cut2;
   double  D3_Dphi_cut1, D3_Dphi_cut2;

   double dPhi012;
   double dPhi013;
   double dPhi014;
   double dPhi023;
   double dPhi024;
   double dPhi034;

   double  L012_DPhi_cut1, L012_DPhi_cut2;
   double  L012_DEta_cut1, L012_DEta_cut2;
   double  L012_DR_cut1, L012_DR_cut2;
   
   double  L013_DPhi_cut1, L013_DPhi_cut2;
   double  L013_DEta_cut1, L013_DEta_cut2;
   double  L013_DR_cut1, L013_DR_cut2;
   
   double  L014_DPhi_cut1, L014_DPhi_cut2;
   double  L014_DEta_cut1, L014_DEta_cut2;
   double  L014_DR_cut1, L014_DR_cut2;
   
   double  L023_DPhi_cut1, L023_DPhi_cut2;
   double  L023_DEta_cut1, L023_DEta_cut2;
   double  L023_DR_cut1, L023_DR_cut2;
   
   double  L024_DPhi_cut1, L024_DPhi_cut2;
   double  L024_DEta_cut1, L024_DEta_cut2;
   double  L024_DR_cut1, L024_DR_cut2;
   
   double  L034_DPhi_cut1, L034_DPhi_cut2;
   double  L034_DEta_cut1, L034_DEta_cut2;
   double  L034_DR_cut1, L034_DR_cut2;
   
   double  L123_DPhi_cut1, L123_DPhi_cut2;
   double  L123_DEta_cut1, L123_DEta_cut2;
   double  L123_DR_cut1, L123_DR_cut2;
   
   double  L124_DPhi_cut1, L124_DPhi_cut2;
   double  L124_DEta_cut1, L124_DEta_cut2;
   double  L124_DR_cut1, L124_DR_cut2;
   
   double  L134_DPhi_cut1, L134_DPhi_cut2;
   double  L134_DEta_cut1, L134_DEta_cut2;
   double  L134_DR_cut1, L134_DR_cut2;
   
   double  L234_DPhi_cut1, L234_DPhi_cut2;
   double  L234_DEta_cut1, L234_DEta_cut2;
   double  L234_DR_cut1, L234_DR_cut2;

   double  L12_eta_upper, L13_eta_upper, L14_eta_upper, L23_eta_upper, L24_eta_upper, L34_eta_upper;
   double  L12_phi_upper, L13_phi_upper, L14_phi_upper, L23_phi_upper, L24_phi_upper, L34_phi_upper;
   double  L12_R_upper, L13_R_upper, L14_R_upper, L23_R_upper, L24_R_upper, L34_R_upper;
   double  L12_eta_bellow, L13_eta_bellow, L14_eta_bellow, L23_eta_bellow, L24_eta_bellow, L34_eta_bellow;
   double  L12_phi_bellow, L13_phi_bellow, L14_phi_bellow, L23_phi_bellow, L24_phi_bellow, L34_phi_bellow;
   double  L12_R_bellow, L13_R_bellow, L14_R_bellow, L23_R_bellow, L24_R_bellow, L34_R_bellow;

   double dPhi;
   double dEta;
   double dR;
   double dPhi_1, dPhi_2, dPhi_3;
   double dEta_1, dEta_2, dEta_3;
   double dR_1, dR_2, dR_3;
   TVector3 first_temp, second_temp;
   int _pass_Ele, _pass_Pos;

   int L012_pass_Ele, L012_pass_Pos; 
   int L013_pass_Ele, L013_pass_Pos;
   int L014_pass_Ele, L014_pass_Pos;
   int L023_pass_Ele, L023_pass_Pos;
   int L024_pass_Ele, L024_pass_Pos;
   int L034_pass_Ele, L034_pass_Pos;
   int L123_pass_Ele, L123_pass_Pos;
   int L124_pass_Ele, L124_pass_Pos;
   int L134_pass_Ele, L134_pass_Pos;
   int L234_pass_Ele, L234_pass_Pos;

   int L12_EM_Ele, L12_EM_Pos; 
   int L13_EM_Ele, L13_EM_Pos;
   int L14_EM_Ele, L14_EM_Pos;
   int L23_EM_Ele, L23_EM_Pos;
   int L24_EM_Ele, L24_EM_Pos;
   int L34_EM_Ele, L34_EM_Pos;

   std::vector<TVector3> first_layer_hits;
   std::vector<TVector3> second_layer_hits;
   std::vector<TVector3> third_layer_hits;
   std::vector<TVector3> fourth_layer_hits;

   double r; // r for radius of pixel tracker layer
   int layers[5];  // initialize as 0, layers contain # of hits on each pixel layer

   TVector3 emvector;
   float EgEt;
   float EgEta;
   float EgPhi;

   std::vector<int> first_layer_hits_Ele_or_Pos;
   std::vector<int> second_layer_hits_Ele_or_Pos;
   std::vector<int> third_layer_hits_Ele_or_Pos;
   std::vector<int> fourth_layer_hits_Ele_or_Pos;
   std::vector<int> hitted_layers;

   void MakeHistograms(TString hname, int nbins, float xmin, float xmax);
   TH1* GetHist(TString hname);
   void FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins);

   void StorePixelHit( int region);
   double StandaloneDPhi( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit );
   double StandaloneDEta( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit );
   double EMmatchingDEta(TVector3& first_layer, TVector3& second_layer, TVector3& egvector);
   double EMmatchingDPhi(TVector3& first_layer, TVector3& second_layer, TVector3& egvector);
   int Signal_window_check( double upper, double value, double lower, int Ele_Pos);
   void FillCutFlow(TString cut, float weight);
   void SetROI(int region);
   void SetSingalBoundary( int region);
   void TriggeringWith_1st2ndPixel(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_1st3rdPixel(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_2nd3rdPixel(int nthFirstHit, int nthSecondHit);
   void TriggeringWithout_4thPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
   void TriggeringWithout_3rdPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
   void TriggeringWithout_2ndPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
   void TriggeringWithout_1stPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);

   inline float deltaPhi(float phi1, float phi2){
     float result = phi1 - phi2;
     while(result > float(M_PI)) result -= float(2*M_PI);
     while(result <= -float(M_PI)) result += float(2*M_PI);
    return result;
   }

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TFile *file3;
   TTree* pixtrk_tree;

   int useDR;

   int count_Entry; 
   int pass_egobjects_check;
   int ntnEg2; 
   int event_denominator; 
   int event_nominator; 

   vector<float> ntEgEt; 
   vector<float> ntL1TkEleEt, ntL1TkEleIsoEt; 
   vector<float> ntEgEta; 
   vector<bool> ntclusterIsEG; 
   vector<float> ntL1TkEleEta, ntL1TkEleIsoEta; 
   vector<float> ntL1TkElePhi, ntL1TkEleIsoPhi;
   vector<bool> ntCl_match; 
   vector<bool> withoutEM_match; 
   vector<bool> withEM_match; 

   vector<bool> ntCl_match_wo4thPix, ntCl_match_wo3thPix, ntCl_match_wo2thPix, ntCl_match_wo1thPix;
   vector<int> Npass_woEM_wo4thPix, Npass_woEM_wo3thPix, Npass_woEM_wo2thPix, Npass_woEM_wo1thPix;
   vector<int> Npass_wEM_wo4thPix, Npass_wEM_wo3thPix, Npass_wEM_wo2thPix, Npass_wEM_wo1thPix;

   vector<float> dphi_L12EM_wo4thPix, dphi_L12EM_wo3thPix, dphi_L12EM_wo2thPix, dphi_L12EM_wo1thPix;
   vector<float> deta_L12EM_wo4thPix, deta_L12EM_wo3thPix, deta_L12EM_wo2thPix, deta_L12EM_wo1thPix;

   vector<int> pass_Ele; 
   vector<int> pass_Pos; 
   vector<int> pass_ElePos; 

   vector<int> ntfirstPix;
   vector<int> ntsecondPix;
   vector<int> ntthirdPix;
   vector<int> ntfourthPix;

   int L1TkEleN;
   int L1TkEleIsoN;

   vector<float>  bHitGx1;
   vector<float>  bHitGy1;
   vector<float>  bHitGz1;
   vector<float>  bHitGx2;
   vector<float>  bHitGy2;
   vector<float>  bHitGz2;
   vector<float>  bHitGx3;
   vector<float>  bHitGy3;
   vector<float>  bHitGz3;
   vector<float>  bHitGx4;
   vector<float>  bHitGy4;
   vector<float>  bHitGz4;

   vector<float> Eget;
   vector<float> Egeta;
   vector<float> Egphi;
   vector<float> Egx;
   vector<float> Egy;
   vector<float> Egz;
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
#ifdef SINGLE_TREE

       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/u/user/moon/L1Pixel_KNU/MC/PhaseII_SE140_CMSSW620_SLHC12_patch1/");
       if(!f || !f->IsOpen()){
	   f = new TFile("/u/user/moon/L1Pixel_KNU/MC/PhaseII_SE140_CMSSW620_SLHC12_patch1/");
       }
       f->GetObject("demo/l1pixel_tree","");
#else
       TChain *chain = new TChain("demo/l1pixel_tree","");

//       chain->Add("/u/user/moon/L1Pixel_KNU/MC/PhaseII_SE140_CMSSW620_SLHC12_patch1/SingleElectron_140PU_1.root/demo/l1pixel_tree");
  //     chain->Add("/u/user/moon/L1Pixel_KNU/MC/PhaseII_SE140_CMSSW620_SLHC12_patch1/SingleElectron_140PU_2.root/demo/l1pixel_tree");
    //   chain->Add("/u/user/moon/L1Pixel_KNU/MC/PhaseII_SE140_CMSSW620_SLHC12_patch1/SingleElectron_140PU_3.root/demo/l1pixel_tree");
      // chain->Add("/u/user/moon/L1Pixel_KNU/MC/PhaseII_SE140_CMSSW620_SLHC12_patch1/SingleElectron_140PU_4.root/demo/l1pixel_tree");
//       chain->Add("/u/user/moon/L1Pixel_KNU/MC/PhaseII_SE140_CMSSW620_SLHC12_patch1/SingleElectron_140PU_5.root/demo/l1pixel_tree");
       chain->Add("/u/user/moon/L1Pixel_KNU/MC/PhaseII_SE140_CMSSW620_SLHC12_patch1/*.root/demo/l1pixel_tree");
       

      tree = chain;
#endif
   }
   Init(tree);
 
   Ele = 1, Pos = 2;
   skip = 0;
/*
   file3 = new TFile("This_is_Test_Tree_output","recreate");
   pixtrk_tree = new TTree("t","t");

   count_Entry = 1;
   pixtrk_tree->Branch("totalEvent", &count_Entry, "count_Entry/I");   pixtrk_tree->Branch("totalEgN", &EgN, "EgN/F");
   pixtrk_tree->Branch("ntnEg2", &ntnEg2, "ntnEg2/I");
   pixtrk_tree->Branch("event_with_selection_passed", &event_denominator, "event_denominator/I");
   pixtrk_tree->Branch("event_with_trigger_passed", &event_nominator, "event_nominator/I");

   pixtrk_tree->Branch("L1TkEleN", &L1TkEleN, "L1TkEleN/I");
   pixtrk_tree->Branch("L1TkEleIsoN", &L1TkEleIsoN, "L1TkEleIsoN/I");
   pixtrk_tree->Branch("ntEgEt",&ntEgEt);
   pixtrk_tree->Branch("ntL1TkEleEt",&ntL1TkEleEt);
   pixtrk_tree->Branch("ntL1TkEleIsoEt",&ntL1TkEleIsoEt);
   pixtrk_tree->Branch("ntEgEta",&ntEgEta);
   pixtrk_tree->Branch("ntclusterIsEG",&ntclusterIsEG);
   pixtrk_tree->Branch("ntL1TkEleEta",&ntL1TkEleEta);
   pixtrk_tree->Branch("ntL1TkEleIsoEta",&ntL1TkEleIsoEta);
   pixtrk_tree->Branch("ntL1TkElePhi",&ntL1TkElePhi);
   pixtrk_tree->Branch("ntL1TkEleIsoPhi",&ntL1TkEleIsoPhi);
   pixtrk_tree->Branch("ntCl_match",&ntCl_match);
   pixtrk_tree->Branch("withoutEM_match",&withoutEM_match);
   pixtrk_tree->Branch("withEM_match",&withEM_match);

   pixtrk_tree->Branch("ntCl_match_wo4thPix",&ntCl_match_wo4thPix);
   pixtrk_tree->Branch("ntCl_match_wo3thPix",&ntCl_match_wo3thPix);
   pixtrk_tree->Branch("ntCl_match_wo2thPix",&ntCl_match_wo2thPix);
   pixtrk_tree->Branch("ntCl_match_wo1thPix",&ntCl_match_wo1thPix);

   pixtrk_tree->Branch("Npass_woEM_wo4thPix", &Npass_woEM_wo4thPix);
   pixtrk_tree->Branch("Npass_woEM_wo3thPix", &Npass_woEM_wo3thPix);
   pixtrk_tree->Branch("Npass_woEM_wo2thPix", &Npass_woEM_wo2thPix);
   pixtrk_tree->Branch("Npass_woEM_wo1thPix", &Npass_woEM_wo1thPix);

   pixtrk_tree->Branch("Npass_wEM_wo4thPix", &Npass_wEM_wo4thPix);
   pixtrk_tree->Branch("Npass_wEM_wo3thPix", &Npass_wEM_wo3thPix);
   pixtrk_tree->Branch("Npass_wEM_wo2thPix", &Npass_wEM_wo2thPix);
   pixtrk_tree->Branch("Npass_wEM_wo1thPix", &Npass_wEM_wo1thPix);

   pixtrk_tree->Branch("deta_L12EM_wo4thPix", &deta_L12EM_wo4thPix);
   pixtrk_tree->Branch("deta_L12EM_wo3thPix", &deta_L12EM_wo3thPix);
   pixtrk_tree->Branch("deta_L12EM_wo2thPix", &deta_L12EM_wo2thPix);
   pixtrk_tree->Branch("deta_L12EM_wo1thPix", &deta_L12EM_wo1thPix);

   pixtrk_tree->Branch("dphi_L12EM_wo4thPix", &dphi_L12EM_wo4thPix);
   pixtrk_tree->Branch("dphi_L12EM_wo3thPix", &dphi_L12EM_wo3thPix);
   pixtrk_tree->Branch("dphi_L12EM_wo2thPix", &dphi_L12EM_wo2thPix);
   pixtrk_tree->Branch("dphi_L12EM_wo1thPix", &dphi_L12EM_wo1thPix);

   pixtrk_tree->Branch("pass_Ele", &pass_Ele);
   pixtrk_tree->Branch("pass_Pos", &pass_Pos);
   pixtrk_tree->Branch("pass_ElePos", &pass_ElePos);

   pixtrk_tree->Branch("ntfirstPix",&ntfirstPix);
   pixtrk_tree->Branch("ntsecondPix",&ntsecondPix);
   pixtrk_tree->Branch("ntthirdPix",&ntthirdPix);
   pixtrk_tree->Branch("ntfourthPix",&ntfourthPix);

*/

}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pileup = 0;
   genPartE = 0;
   genPartPt = 0;
   genPartEta = 0;
   genPartPhi = 0;
   genPartCharge = 0;
   genPartId = 0;
   propgenPartE = 0;
   propgenPartPt = 0;
   propgenPartEta = 0;
   propgenPartPhi = 0;
   propgenPartCharge = 0;
   propgenPartId = 0;
   simTrkPt = 0;
   simTrkEta = 0;
   simTrkPhi = 0;
   simTrkId = 0;
   simTrkType = 0;
   simTrkVx = 0;
   simTrkVy = 0;
   simTrkVz = 0;
   simVx = 0;
   simVy = 0;
   simVz = 0;
   egCrysE = 0;
   egCrysEt = 0;
   egCrysEta = 0;
   egCrysPhi = 0;
   egCrysGx = 0;
   egCrysGy = 0;
   egCrysGz = 0;
   egCrysCharge = 0;
   egE = 0;
   egEt = 0;
   egEta = 0;
   egPhi = 0;
   egGx = 0;
   egGy = 0;
   egGz = 0;
   egclusterBS = 0;
   egclusterPt = 0;
   egclusterPtCore = 0;
   egclusterEta = 0;
   egclusterPhi = 0;
   egHoE = 0;
   clusterHoE = 0;
   seedHoE = 0;
   clusterIsEG = 0;
   seedIsEG = 0;
   TriggerTower_IsBarrel = 0;
   seedHadEt = 0;
   seedEmEt = 0;
   HadEt = 0;
   EmEt = 0;
   fHitDisk = 0;
   fHitBlade = 0;
   fHitSide = 0;
   fHitGx = 0;
   fHitGy = 0;
   fHitGz = 0;
   fClSize = 0;
   fClSizeX = 0;
   fClSizeY = 0;
   bHitLayer = 0;
   bHitLadder = 0;
   bHitGx = 0;
   bHitGy = 0;
   bHitGz = 0;
   bClSize = 0;
   bClSizeX = 0;
   bClSizeY = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
//   fChain->SetBranchAddress("bunchN", &bunchN, &b_bunchN);
   fChain->SetBranchAddress("pileup", &pileup, &b_pileup);
   fChain->SetBranchAddress("beamSpotX0", &beamSpotX0, &b_beamSpotX0);
   fChain->SetBranchAddress("beamSpotY0", &beamSpotY0, &b_beamSpotY0);
   fChain->SetBranchAddress("beamSpotZ0", &beamSpotZ0, &b_beamSpotZ0);
   fChain->SetBranchAddress("beamSpotX0Error", &beamSpotX0Error, &b_beamSpotX0Error);
   fChain->SetBranchAddress("beamSpotY0Error", &beamSpotY0Error, &b_beamSpotY0Error);
   fChain->SetBranchAddress("beamSpotZ0Error", &beamSpotZ0Error, &b_beamSpotZ0Error);
   fChain->SetBranchAddress("beamWidthX", &beamWidthX, &b_beamWidthX);
   fChain->SetBranchAddress("beamWidthY", &beamWidthY, &b_beamWidthY);
   fChain->SetBranchAddress("beamSigmaZ", &beamSigmaZ, &b_beamSigmaZ);
   fChain->SetBranchAddress("beamWidthXError", &beamWidthXError, &b_beamWidthXError);
   fChain->SetBranchAddress("beamWidthYError", &beamWidthYError, &b_beamWidthYError);
   fChain->SetBranchAddress("beamSigmaZError", &beamSigmaZError, &b_beamSigmaZError);
//   fChain->SetBranchAddress("genPartN", &genPartN, &b_genPartN);
   fChain->SetBranchAddress("genPartE", &genPartE, &b_genPartE);
   fChain->SetBranchAddress("genPartPt", &genPartPt, &b_genPartPt);
   fChain->SetBranchAddress("genPartEta", &genPartEta, &b_genPartEta);
   fChain->SetBranchAddress("genPartPhi", &genPartPhi, &b_genPartPhi);
   fChain->SetBranchAddress("genPartCharge", &genPartCharge, &b_genPartCharge);
   fChain->SetBranchAddress("genPartId", &genPartId, &b_genPartId);
   fChain->SetBranchAddress("propgenPartE", &propgenPartE, &b_propgenPartE);
   fChain->SetBranchAddress("propgenPartPt", &propgenPartPt, &b_propgenPartPt);
   fChain->SetBranchAddress("propgenPartEta", &propgenPartEta, &b_propgenPartEta);
   fChain->SetBranchAddress("propgenPartPhi", &propgenPartPhi, &b_propgenPartPhi);
   fChain->SetBranchAddress("propgenPartCharge", &propgenPartCharge, &b_propgenPartCharge);
   fChain->SetBranchAddress("propgenPartId", &propgenPartId, &b_propgenPartId);
   fChain->SetBranchAddress("simTrkN", &simTrkN, &b_simTrkN);
   fChain->SetBranchAddress("simTrkPt", &simTrkPt, &b_simTrkPt);
   fChain->SetBranchAddress("simTrkEta", &simTrkEta, &b_simTrkEta);
   fChain->SetBranchAddress("simTrkPhi", &simTrkPhi, &b_simTrkPhi);
   fChain->SetBranchAddress("simTrkId", &simTrkId, &b_simTrkId);
   fChain->SetBranchAddress("simTrkType", &simTrkType, &b_simTrkType);
   fChain->SetBranchAddress("simTrkVx", &simTrkVx, &b_simTrkVx);
   fChain->SetBranchAddress("simTrkVy", &simTrkVy, &b_simTrkVy);
   fChain->SetBranchAddress("simTrkVz", &simTrkVz, &b_simTrkVz);
   fChain->SetBranchAddress("simVx", &simVx, &b_simVx);
   fChain->SetBranchAddress("simVy", &simVy, &b_simVy);
   fChain->SetBranchAddress("simVz", &simVz, &b_simVz);
//   fChain->SetBranchAddress("gammaBrem", &gammaBrem, &b_gammaBrem);
//   fChain->SetBranchAddress("hardGamma", &hardGamma, &b_hardGamma);
   fChain->SetBranchAddress("egCrysN", &egCrysN, &b_egCrysN);
   fChain->SetBranchAddress("egCrysE", &egCrysE, &b_egCrysE);
   fChain->SetBranchAddress("egCrysEt", &egCrysEt, &b_egCrysEt);
   fChain->SetBranchAddress("egCrysEta", &egCrysEta, &b_egCrysEta);
   fChain->SetBranchAddress("egCrysPhi", &egCrysPhi, &b_egCrysPhi);
   fChain->SetBranchAddress("egCrysGx", &egCrysGx, &b_egCrysGx);
   fChain->SetBranchAddress("egCrysGy", &egCrysGy, &b_egCrysGy);
   fChain->SetBranchAddress("egCrysGz", &egCrysGz, &b_egCrysGz);
//   fChain->SetBranchAddress("egCrysCharge", &egCrysCharge, &b_egCrysCharge);
   fChain->SetBranchAddress("egHoE", &egHoE, &b_egHoE);
   fChain->SetBranchAddress("clusterHoE", &clusterHoE, &b_clusterHoE);
   fChain->SetBranchAddress("seedHoE", &seedHoE, &b_seedHoE);
   fChain->SetBranchAddress("clusterIsEG", &clusterIsEG, &b_clusterIsEG);
   fChain->SetBranchAddress("seedIsEG", &seedIsEG, &b_seedIsEG);
   fChain->SetBranchAddress("TriggerTower_IsBarrel", &TriggerTower_IsBarrel, &b_TriggerTower_IsBarrel);
   fChain->SetBranchAddress("seedHadEt", &seedHadEt, &b_seedHadEt);
   fChain->SetBranchAddress("seedEmEt", &seedEmEt, &b_seedEmEt);
   fChain->SetBranchAddress("HadEt", &HadEt, &b_HadEt);
   fChain->SetBranchAddress("EmEt", &EmEt, &b_EmEt);
   fChain->SetBranchAddress("l1tkegN", &l1tkegN, &b_l1tkegN);
   fChain->SetBranchAddress("l1tkegEta", &l1tkegEta, &b_l1tkegEta);
   fChain->SetBranchAddress("l1tkegEt", &l1tkegEt, &b_l1tkegEt);
   fChain->SetBranchAddress("l1tkegPhi", &l1tkegPhi, &b_l1tkegPhi);
   fChain->SetBranchAddress("l1tkegIsoN", &l1tkegIsoN, &b_l1tkegIsoN);
   fChain->SetBranchAddress("l1tkegIsoEta", &l1tkegIsoEta, &b_l1tkegIsoEta);
   fChain->SetBranchAddress("l1tkegIsoEt", &l1tkegIsoEt, &b_l1tkegIsoEt);
   fChain->SetBranchAddress("l1tkegIsoPhi", &l1tkegIsoPhi, &b_l1tkegIsoPhi);
   fChain->SetBranchAddress("egN", &egN, &b_egN);
   fChain->SetBranchAddress("egE", &egE, &b_egE);
   fChain->SetBranchAddress("egEt", &egEt, &b_egEt);
   fChain->SetBranchAddress("egEta", &egEta, &b_egEta);
   fChain->SetBranchAddress("egPhi", &egPhi, &b_egPhi);
   fChain->SetBranchAddress("egGx", &egGx, &b_egGx);
   fChain->SetBranchAddress("egGy", &egGy, &b_egGy);
   fChain->SetBranchAddress("egGz", &egGz, &b_egGz);
//   fChain->SetBranchAddress("egclusterBS", &egclusterBS, &b_egclusterBS);
//   fChain->SetBranchAddress("egclusterPt", &egclusterPt, &b_egclusterPt);
//   fChain->SetBranchAddress("egclusterPtCore", &egclusterPtCore, &b_egclusterPtCore);
//   fChain->SetBranchAddress("egclusterEta", &egclusterEta, &b_egclusterEta);
//   fChain->SetBranchAddress("egclusterPhi", &egclusterPhi, &b_egclusterPhi);
   fChain->SetBranchAddress("fHitN", &fHitN, &b_fHitN);
   fChain->SetBranchAddress("fHitDisk", &fHitDisk, &b_fHitDisk);
   fChain->SetBranchAddress("fHitBlade", &fHitBlade, &b_fHitBlade);
   fChain->SetBranchAddress("fHitSide", &fHitSide, &b_fHitSide);
   fChain->SetBranchAddress("fHitGx", &fHitGx, &b_fHitGx);
   fChain->SetBranchAddress("fHitGy", &fHitGy, &b_fHitGy);
   fChain->SetBranchAddress("fHitGz", &fHitGz, &b_fHitGz);
   fChain->SetBranchAddress("fClSize", &fClSize, &b_fClSize);
   fChain->SetBranchAddress("fClSizeX", &fClSizeX, &b_fClSizeX);
   fChain->SetBranchAddress("fClSizeY", &fClSizeY, &b_fClSizeY);
   fChain->SetBranchAddress("bHitN", &bHitN, &b_bHitN);
   fChain->SetBranchAddress("bHitLayer", &bHitLayer, &b_bHitLayer);
   fChain->SetBranchAddress("bHitLadder", &bHitLadder, &b_bHitLadder);
   fChain->SetBranchAddress("bHitGx", &bHitGx, &b_bHitGx);
   fChain->SetBranchAddress("bHitGy", &bHitGy, &b_bHitGy);
   fChain->SetBranchAddress("bHitGz", &bHitGz, &b_bHitGz);
   fChain->SetBranchAddress("bClSize", &bClSize, &b_bClSize);
   fChain->SetBranchAddress("bClSizeX", &bClSizeX, &b_bClSizeX);
   fChain->SetBranchAddress("bClSizeY", &bClSizeY, &b_bClSizeY);

   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void test::MakeHistograms(TString hname, int nbins, float xmin, float xmax){

 maphist[hname] =  new TH1F(hname.Data(),hname.Data(),nbins,xmin,xmax);

}

TH1* test::GetHist(TString hname){

 TH1* h = NULL;
 std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
 if(mapit != maphist.end()) return mapit->second;

 return h;

}

void test::FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins){

 if(GetHist(histname)) GetHist(histname)->Fill(value, w);
 else{
//     cout << "Making histogram..." << endl;
     MakeHistograms(histname, nbins, xmin, xmax);
     if(GetHist(histname)) GetHist(histname)->Fill(value, w);
 }
}


double test::StandaloneDPhi( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit ){
  if( first_hit == 0 ){

     TVector3 temp;

     if( second_hit == 1 ) temp.SetXYZ( first_layer_hits[which_second_hit].X(), first_layer_hits[which_second_hit].Y(), first_layer_hits[which_second_hit].Z() );
     if( second_hit == 2 ) temp.SetXYZ( second_layer_hits[which_second_hit].X(), second_layer_hits[which_second_hit].Y(), second_layer_hits[which_second_hit].Z() );
     if( second_hit == 3 ) temp.SetXYZ( third_layer_hits[which_second_hit].X(), third_layer_hits[which_second_hit].Y(), third_layer_hits[which_second_hit].Z() );


     if( third_hit == 2 ) return deltaPhi( (second_layer_hits[which_third_hit] - temp).Phi(), temp.Phi());
     if( third_hit == 3 ) return deltaPhi( (third_layer_hits[which_third_hit] - temp).Phi(), temp.Phi());
     if( third_hit == 4 ) return deltaPhi( (fourth_layer_hits[which_third_hit] - temp).Phi(), temp.Phi());

  }
  if( first_hit != 0 ){
    TVector3 temp_first_layer;
    TVector3 temp_second_layer;
    TVector3 temp_third_layer;

    if( first_hit == 1 ) temp_first_layer = first_layer_hits[which_first_hit];
    if( first_hit == 2 ) temp_first_layer = second_layer_hits[which_first_hit];

    if( second_hit == 2 ) temp_second_layer = second_layer_hits[which_second_hit];
    if( second_hit == 3 ) temp_second_layer = third_layer_hits[which_second_hit];


    if( third_hit == 3 ) temp_third_layer = third_layer_hits[which_third_hit];
    if( third_hit == 4 ) temp_third_layer = fourth_layer_hits[which_third_hit];

    return deltaPhi( (temp_third_layer - temp_second_layer).Phi(), (temp_second_layer - temp_first_layer).Phi());
  }
  return 0.;
}

double test::StandaloneDEta( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit ){

    TVector3 temp_first_layer;
    TVector3 temp_second_layer;
    TVector3 temp_third_layer;

    if( first_hit == 1 ) temp_first_layer = first_layer_hits[which_first_hit];
    if( first_hit == 2 ) temp_first_layer = second_layer_hits[which_first_hit];

    if( second_hit == 2 ) temp_second_layer = second_layer_hits[which_second_hit];
    if( second_hit == 3 ) temp_second_layer = third_layer_hits[which_second_hit];


    if( third_hit == 3 ) temp_third_layer = third_layer_hits[which_third_hit];
    if( third_hit == 4 ) temp_third_layer = fourth_layer_hits[which_third_hit];

    return (temp_third_layer - temp_second_layer).PseudoRapidity() - (temp_second_layer - temp_first_layer).PseudoRapidity();
}

double test::EMmatchingDEta(TVector3& first_layer, TVector3& second_layer, TVector3& egvector){

    TVector3 pixelVector = second_layer - first_layer;
    TVector3 EM_pixelVector = egvector - second_layer;
    return EM_pixelVector.Eta() - pixelVector.Eta();

}

double test::EMmatchingDPhi(TVector3& first_layer, TVector3& second_layer, TVector3& egvector){


    TVector3 pixelVector = second_layer - first_layer;
    TVector3 EM_pixelVector = egvector - second_layer;
    return deltaPhi( EM_pixelVector.Phi(), pixelVector.Phi());


}

int test::Signal_window_check( double upper, double value, double lower, int Ele_Pos){

  if( Ele_Pos == 1 ){ // 1 is Electron
    if( value <= upper && value >= lower){
      return true;
    }
    else
     return false;
  }
  if( Ele_Pos == 2 ){ // 2 is Positron
    if( value >= -upper && value <= -lower){
      return true;
    }
    else
     return false;
  }
 return 0;
}

void test::FillCutFlow(TString cut, float weight){


  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);

  }
  else{
    test::MakeHistograms("cutflow", 6,0.,6.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"MinEtCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"EtaCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"EvtCut");


  }
}
void test::StorePixelHit( int region){

        for(int a=0; a<bHitN; a++){
           int Dphi_Ele_pass = 0;
           int Dphi_Pos_pass = 0;
           double Dphi = 0.;
           int el_or_po = 0; // electron = 1, positron = 2, both = 3
           TVector3 current_hit;
           current_hit.SetXYZ( bHitGx->at(a), bHitGy->at(a), bHitGz->at(a) );
           Dphi = deltaPhi( current_hit.Phi(), EgPhi);

           double r = sqrt( pow(bHitGx->at(a),2) + pow(bHitGy->at(a),2) );

           if( r < 5.5 ){ // First layer
             if( Dphi < L1_Dphi_cut1 && Dphi > L1_Dphi_cut2 ){
                Dphi_Ele_pass = 1; el_or_po = 1;
             }
             if( Dphi > -L1_Dphi_cut1 && Dphi < -L1_Dphi_cut2 ){
                Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
             }
             if( Dphi_Ele_pass || Dphi_Pos_pass ){
                layers[1]++;
                first_layer_hits.push_back( TVector3(bHitGx->at(a), bHitGy->at(a), bHitGz->at(a)));
                first_layer_hits_Ele_or_Pos.push_back(el_or_po);
             }
           }

           if( region == 1 || region == 2 || region == 3 ){
              if( r > 5.5 && r < 8.5 ){ // Second layer
                if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
                   Dphi_Ele_pass = 1; el_or_po = 1;
                }
                if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
                   Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
                }
                if( Dphi_Ele_pass || Dphi_Pos_pass ){
                  layers[2]++;
                  second_layer_hits.push_back( TVector3(bHitGx->at(a), bHitGy->at(a), bHitGz->at(a)));
                  second_layer_hits_Ele_or_Pos.push_back(el_or_po);
                }
              }
           }
           if( region == 1 || region == 2){
             if( r > 8.5 && r < 13 ){ // Third layer 
               if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[3]++;
                 third_layer_hits.push_back( TVector3(bHitGx->at(a), bHitGy->at(a), bHitGz->at(a)));
                 third_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }
           if( region == 1){
              if( r > 13 && r < 18 ){ // Fourth layer
                if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
                   Dphi_Ele_pass = 1; el_or_po = 1;
                }
                if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
                   Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
                }
                if( Dphi_Ele_pass || Dphi_Pos_pass ){
                  layers[4]++;
                  fourth_layer_hits.push_back( TVector3(bHitGx->at(a), bHitGy->at(a), bHitGz->at(a)));
                  fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
                }
              }
           }
       }
/*
	for(int a=0; a<fHitN; a++){
           int Dphi_Ele_pass = 0;
           int Dphi_Pos_pass = 0;
           double Dphi = 0.;
           int el_or_po = 0; // electron = 1, positron = 2, both = 3
           TVector3 current_hit;
           current_hit.SetXYZ( fHitGx->at(a), fHitGy->at(a), fHitGz->at(a) );
           Dphi = deltaPhi(current_hit.Phi(), EgPhi);

           double disk_z = fabs(fHitGz->at(a));

           if( region == 4 || region ==5 ){
             if( disk_z > 28 && disk_z < 36 ){ // first disk
               if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[2]++;
                 second_layer_hits.push_back( TVector3(fHitGx->at(a), fHitGy->at(a), fHitGz->at(a)));
                 second_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }
           if( region == 4 || region == 5 ){
             if( disk_z > 36 && disk_z < 44 ){ // second disk
               if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[3]++;
                 third_layer_hits.push_back( TVector3(fHitGx->at(a), fHitGy->at(a), fHitGz->at(a)));
                 third_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }
           if( region == 3 ){
             if( disk_z > 28 && disk_z < 36 ){ // first disk
               if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[3]++;
                 third_layer_hits.push_back( TVector3(fHitGx->at(a), fHitGy->at(a), fHitGz->at(a)));
                 third_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }
           if( region == 4 || region == 5 ){
             if( disk_z > 44 && disk_z < 54 ){ // third disk
               if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[4]++;
                 fourth_layer_hits.push_back( TVector3(fHitGx->at(a), fHitGy->at(a), fHitGz->at(a)));
                 fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }

           if( region == 2 ){
             if( disk_z > 28 && disk_z < 36 ){ // first disk
               if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[4]++;
                 fourth_layer_hits.push_back( TVector3(fHitGx->at(a), fHitGy->at(a), fHitGz->at(a)));
                 fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }
           if( region == 3){
             if( disk_z > 36 && disk_z < 44 ){ // second disk
               if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[4]++;
                 fourth_layer_hits.push_back( TVector3(fHitGx->at(a), fHitGy->at(a), fHitGz->at(a)));
                 fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }
        }
*/	
      ntfirstPix.push_back(layers[1]);
      ntsecondPix.push_back(layers[2]);
      ntthirdPix.push_back(layers[3]);
      ntfourthPix.push_back(layers[4]);

}

void test::SetROI(int region){

    L1_Dphi_cut1 = ROI_func( EgEt, region, 0, 0);
    L1_Dphi_cut2 = ROI_func( EgEt, region, 0, 1);

    L2_Dphi_cut1 = ROI_func( EgEt, region, 1, 0);
    L2_Dphi_cut2 = ROI_func( EgEt, region, 1, 1);

    L3_Dphi_cut1 = ROI_func( EgEt, region, 2, 0);
    L3_Dphi_cut2 = ROI_func( EgEt, region, 2, 1);

    L4_Dphi_cut1 = ROI_func( EgEt, region, 3, 0);
    L4_Dphi_cut2 = ROI_func( EgEt, region, 3, 1);

}

void test::SetSingalBoundary( int region){

      L12_phi_upper = SW_func2( 1,0, EgEt, region,0);
      L12_phi_bellow = SW_func2( 1,0, EgEt, region,1);

      L13_phi_upper = SW_func2( 1,1, EgEt, region,0);
      L13_phi_bellow = SW_func2( 1,1, EgEt, region,1);

      L14_phi_upper = SW_func2( 1,2, EgEt, region,0);
      L14_phi_bellow = SW_func2( 1,2, EgEt, region,1);

      L23_phi_upper = SW_func2( 1,3, EgEt, region,0);
      L23_phi_bellow = SW_func2( 1,3, EgEt, region,1);

      L24_phi_upper = SW_func2( 1,4, EgEt, region,0);
      L24_phi_bellow = SW_func2( 1,4, EgEt, region,1);

      L34_phi_upper = SW_func2( 1,5, EgEt, region,0);
      L34_phi_bellow = SW_func2( 1,5, EgEt, region,1);

      L12_eta_upper = SW_func2( 2,0, EgEt, region,0);
      L12_eta_bellow = SW_func2( 2,0, EgEt, region,1);

      L13_eta_upper = SW_func2( 2,1, EgEt, region,0);
      L13_eta_bellow = SW_func2( 2,1, EgEt, region,1);

      L14_eta_upper = SW_func2( 2,2, EgEt, region,0);
      L14_eta_bellow = SW_func2( 2,2, EgEt, region,1);

      L23_eta_upper = SW_func2( 2,3, EgEt, region,0);
      L23_eta_bellow = SW_func2( 2,3, EgEt, region,1);

      L24_eta_upper = SW_func2( 2,4, EgEt, region,0);
      L24_eta_bellow = SW_func2( 2,4, EgEt, region,1);

      L34_eta_upper = SW_func2( 2,5, EgEt, region,0);
      L34_eta_bellow = SW_func2( 2,5, EgEt, region,1);

      L12_R_upper = SW_func2( 3,0, EgEt, region,0);
      L12_R_bellow = SW_func2( 3,0, EgEt, region,1);

      L13_R_upper = SW_func2( 3,1, EgEt, region,0);
      L13_R_bellow = SW_func2( 3,1, EgEt, region,1);

      L14_R_upper = SW_func2( 3,2, EgEt, region,0);
      L14_R_bellow = SW_func2( 3,2, EgEt, region,1);
      
      L23_R_upper = SW_func2( 3,3, EgEt, region,0);
      L23_R_bellow = SW_func2( 3,3, EgEt, region,1);

      L24_R_upper = SW_func2( 3,4, EgEt, region,0);
      L24_R_bellow = SW_func2( 3,4, EgEt, region,1);

      L34_R_upper = SW_func2( 3,5, EgEt, region,0);
      L34_R_bellow = SW_func2( 3,5, EgEt, region,1);

      L012_DPhi_cut1 = SW_func1( 1, 0, EgEt, region, 0);
      L012_DPhi_cut2 = SW_func1( 1, 0, EgEt, region, 1);

      L013_DPhi_cut1 = SW_func1( 1, 1, EgEt, region, 0);
      L013_DPhi_cut2 = SW_func1( 1, 1, EgEt, region, 1);

      L014_DPhi_cut1 = SW_func1( 1, 2, EgEt, region, 0);
      L014_DPhi_cut2 = SW_func1( 1, 2, EgEt, region, 1);

      L023_DPhi_cut1 = SW_func1( 1, 3, EgEt, region, 0);
      L023_DPhi_cut2 = SW_func1( 1, 3, EgEt, region, 1);

      L024_DPhi_cut1 = SW_func1( 1, 4, EgEt, region, 0);
      L024_DPhi_cut2 = SW_func1( 1, 4, EgEt, region, 1);

      L034_DPhi_cut1 = SW_func1( 1, 5, EgEt, region, 0);
      L034_DPhi_cut2 = SW_func1( 1, 5, EgEt, region, 1);

      L123_DPhi_cut1 = SW_func1( 1, 6, EgEt, region, 0);
      L123_DPhi_cut2 = SW_func1( 1, 6, EgEt, region, 1);

      L124_DPhi_cut1 = SW_func1( 1, 7, EgEt, region, 0);
      L124_DPhi_cut2 = SW_func1( 1, 7 ,EgEt, region, 1);

      L134_DPhi_cut1 = SW_func1( 1, 8, EgEt, region, 0);
      L134_DPhi_cut2 = SW_func1( 1, 8, EgEt, region, 1);

      L234_DPhi_cut1 = SW_func1( 1, 9, EgEt, region, 0);
      L234_DPhi_cut2 = SW_func1( 1, 9, EgEt, region, 1);

      L123_DEta_cut1 = SW_func1( 2, 0, EgEt, region, 0);
      L123_DEta_cut2 = SW_func1( 2, 0, EgEt, region, 1);

      L124_DEta_cut1 = SW_func1( 2, 1, EgEt, region, 0);
      L124_DEta_cut2 = SW_func1( 2, 1 ,EgEt, region, 1);

      L134_DEta_cut1 = SW_func1( 2, 2, EgEt, region, 0);
      L134_DEta_cut2 = SW_func1( 2, 2, EgEt, region, 1);

      L234_DEta_cut1 = SW_func1( 2, 3, EgEt, region, 0);
      L234_DEta_cut2 = SW_func1( 2, 3, EgEt, region, 1);

      L123_DR_cut1 = SW_func1( 3, 0, EgEt, region, 0);
      L123_DR_cut2 = SW_func1( 3, 0, EgEt, region, 1);

      L124_DR_cut1 = SW_func1( 3, 1, EgEt, region, 0);
      L124_DR_cut2 = SW_func1( 3, 1 ,EgEt, region, 1);

      L134_DR_cut1 = SW_func1( 3, 2, EgEt, region, 0);
      L134_DR_cut2 = SW_func1( 3, 2, EgEt, region, 1);

      L234_DR_cut1 = SW_func1( 3, 3, EgEt, eta_region, 0);
      L234_DR_cut2 = SW_func1( 3, 3, EgEt, eta_region, 1);


}
void test::TriggeringWithout_4thPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit);
dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthThirdHit );
dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthSecondHit, nthThirdHit);

dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
dR_1   = sqrt( pow(dPhi_1,2) + pow(dEta_1,2) );

dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], third_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], third_layer_hits[nthThirdHit], emvector);
dR_2   = sqrt( pow(dPhi_2,2) + pow(dEta_2,2) );

dPhi_3 = EMmatchingDPhi(second_layer_hits[nthSecondHit], third_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(second_layer_hits[nthSecondHit], third_layer_hits[nthThirdHit], emvector);
dR_3   = sqrt( pow(dPhi_3,2) + pow(dEta_3,2) );

L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );
L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );
L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );

if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele) )
   L12_EM_Ele = 1;
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele) )
   L12_EM_Pos = 1;

if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele) )
   L13_EM_Ele = 1;
if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele) )
   L13_EM_Pos = 1;

if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele) )
   L23_EM_Ele = 1;
if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele) )
   L23_EM_Pos = 1;

if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Ele ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele ) )
   L123_pass_Ele = 1;
if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Pos ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele ) ) 
   L123_pass_Pos = 1;

}

void test::TriggeringWithout_3rdPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit);
dPhi014 = StandaloneDPhi( 0, 1, 4, 0, nthFirstHit, nthThirdHit );
dPhi024 = StandaloneDPhi( 0, 2, 4, 0, nthSecondHit, nthThirdHit);

dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
dR_1   = sqrt( pow(dPhi_1,2) + pow(dEta_1,2) );

dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dR_2   = sqrt( pow(dPhi_2,2) + pow(dEta_2,2) );

dPhi_3 = EMmatchingDPhi(second_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(second_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dR_3   = sqrt( pow(dPhi_3,2) + pow(dEta_3,2) );

L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );
L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );
L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );


if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele) ) 
   L12_EM_Ele = 1;
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele) )
   L12_EM_Pos = 1;

if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele) )
   L14_EM_Ele = 1;
if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele) )
   L14_EM_Pos = 1;

if( Signal_window_check(L24_phi_upper, dPhi_3, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_3, L24_eta_bellow, Ele) ) 
   L24_EM_Ele = 1;
if( Signal_window_check(L24_phi_upper, dPhi_3, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_3, L24_eta_bellow, Ele) )
   L24_EM_Pos = 1;

if( Signal_window_check(L124_DPhi_cut1, dPhi, L124_DPhi_cut2, Ele ) && Signal_window_check(L124_DEta_cut1, dEta, L124_DEta_cut2, Ele ) )
   L124_pass_Ele = 1;
if( Signal_window_check(L124_DPhi_cut1, dPhi, L124_DPhi_cut2, Pos ) && Signal_window_check(L124_DEta_cut1, dEta, L124_DEta_cut2, Ele ) )
   L124_pass_Pos = 1;
}

void test::TriggeringWithout_2ndPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthSecondHit);
dPhi014 = StandaloneDPhi( 0, 1, 4, 0, nthFirstHit, nthThirdHit );
dPhi034 = StandaloneDPhi( 0, 3, 4, 0, nthSecondHit, nthThirdHit);


dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
dR_1   = sqrt( pow(dPhi_1,2) + pow(dEta_1,2) );

dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dR_2   = sqrt( pow(dPhi_2,2) + pow(dEta_2,2) );

dPhi_3 = EMmatchingDPhi(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dR_3   = sqrt( pow(dPhi_3,2) + pow(dEta_3,2) );

L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );
L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );
L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );

if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele) ) 
   L13_EM_Ele = 1;
if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele) )
   L13_EM_Pos = 1;

if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele) )
   L14_EM_Ele = 1;
if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele) )
   L14_EM_Pos = 1;

if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele) )
   L34_EM_Ele = 1;
if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele) ) 
   L34_EM_Pos = 1;

if( Signal_window_check(L134_DPhi_cut1, dPhi, L134_DPhi_cut2, Ele ) && Signal_window_check(L134_DEta_cut1, dEta, L134_DEta_cut2, Ele ) )
  L134_pass_Ele = 1;
if( Signal_window_check(L134_DPhi_cut1, dPhi, L134_DPhi_cut2, Pos ) && Signal_window_check(L134_DEta_cut1, dEta, L134_DEta_cut2, Ele ) ) 
  L134_pass_Pos = 1;
}

void test::TriggeringWithout_1stPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthFirstHit, nthSecondHit);
dPhi024 = StandaloneDPhi( 0, 2, 4, 0, nthFirstHit, nthThirdHit );
dPhi034 = StandaloneDPhi( 0, 3, 4, 0, nthSecondHit, nthThirdHit);

dPhi_1 = EMmatchingDPhi(second_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(second_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
dR_1   = sqrt( pow(dPhi_1,2) + pow(dEta_1,2) );

dPhi_2 = EMmatchingDPhi(second_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(second_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dR_2   = sqrt( pow(dPhi_2,2) + pow(dEta_2,2) );

dPhi_3 = EMmatchingDPhi(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dR_3   = sqrt( pow(dPhi_3,2) + pow(dEta_3,2) );

L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );
L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );
L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );

if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele) )
   L23_EM_Ele = 1;
if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele) )
   L23_EM_Pos = 1;

if( Signal_window_check(L24_phi_upper, dPhi_2, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_2, L24_eta_bellow, Ele) )
   L24_EM_Ele = 1;
if( Signal_window_check(L24_phi_upper, dPhi_2, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_2, L24_eta_bellow, Ele) )
   L24_EM_Pos = 1;

if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele) )
   L34_EM_Ele = 1;
if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele) )
   L34_EM_Pos = 1;

if( Signal_window_check(L234_DPhi_cut1, dPhi, L234_DPhi_cut2, Ele ) && Signal_window_check(L234_DEta_cut1, dEta, L234_DEta_cut2, Ele ) )
   L234_pass_Ele = 1;
if( Signal_window_check(L234_DPhi_cut1, dPhi, L234_DPhi_cut2, Pos ) && Signal_window_check(L234_DEta_cut1, dEta, L234_DEta_cut2, Ele ) )
   L234_pass_Pos = 1;
}
void test::TriggeringWith_1st2ndPixel(int nthFirstHit, int nthSecondHit){
dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );

if( L012_pass_Ele &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (second_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || second_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L012_pass_Pos &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (second_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || second_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_1st3rdPixel(int nthFirstHit, int nthSecondHit){
dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthSecondHit );
L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );

 if( L013_pass_Ele &&
     (first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
 if( L013_pass_Pos &&
     (first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;

}
void test::TriggeringWith_2nd3rdPixel(int nthFirstHit, int nthSecondHit){
dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthFirstHit, nthSecondHit );
L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );

 if( L023_pass_Ele &&
     (second_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || second_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
 if( L023_pass_Pos &&
     (second_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || second_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

#endif // #ifdef test_cxx

