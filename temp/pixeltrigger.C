#define test_cxx
#include "test_180314.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <TLorentzVector.h>

using namespace std;

void test::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   useDR == 1;


   TFile *test_file = new TFile("Test_File_140PU_180315.root","recreate");
//   TFile *test_file = new TFile("temp.root","recreate");

   TTree *test_sw = new TTree("sw","sw");

   test_sw->Branch("bHitGx1",&bHitGx1);
   test_sw->Branch("bHitGy1",&bHitGy1);
   test_sw->Branch("bHitGz1",&bHitGz1);

   test_sw->Branch("bHitGx2",&bHitGx2);
   test_sw->Branch("bHitGy2",&bHitGy2);
   test_sw->Branch("bHitGz2",&bHitGz2);

   test_sw->Branch("bHitGx3",&bHitGx3);
   test_sw->Branch("bHitGy3",&bHitGy3);
   test_sw->Branch("bHitGz3",&bHitGz3);

   test_sw->Branch("bHitGx4",&bHitGx4);
   test_sw->Branch("bHitGy4",&bHitGy4);
   test_sw->Branch("bHitGz4",&bHitGz4);

   test_sw->Branch("Eget",&Eget);
   test_sw->Branch("Egeta",&Egeta);
   test_sw->Branch("Egphi",&Egphi);
   test_sw->Branch("Egx",&Egx);
   test_sw->Branch("Egy",&Egy);
   test_sw->Branch("Egz",&Egz);

/*
   count_Entry = 1;
   test_sw->Branch("totalEvent", &count_Entry, "count_Entry/I");   pixtrk_tree->Branch("totalEgN", &EgN, "EgN/F");
   test_sw->Branch("ntnEg2", &ntnEg2, "ntnEg2/I");
   test_sw->Branch("event_with_selection_passed", &event_denominator, "event_denominator/I");
   test_sw->Branch("event_with_trigger_passed", &event_nominator, "event_nominator/I");
*/
   test_sw->Branch("L1TkEleN", &L1TkEleN, "L1TkEleN/I");
   test_sw->Branch("L1TkEleIsoN", &L1TkEleIsoN, "L1TkEleIsoN/I");
   test_sw->Branch("ntEgEt",&ntEgEt);
   test_sw->Branch("ntL1TkEleEt",&ntL1TkEleEt);
   test_sw->Branch("ntL1TkEleIsoEt",&ntL1TkEleIsoEt);
   test_sw->Branch("ntEgEta",&ntEgEta);
   test_sw->Branch("ntclusterIsEG",&ntclusterIsEG);
   test_sw->Branch("ntL1TkEleEta",&ntL1TkEleEta);
   test_sw->Branch("ntL1TkEleIsoEta",&ntL1TkEleIsoEta);
   test_sw->Branch("ntL1TkElePhi",&ntL1TkElePhi);
   test_sw->Branch("ntL1TkEleIsoPhi",&ntL1TkEleIsoPhi);
   test_sw->Branch("ntCl_match",&ntCl_match);
   test_sw->Branch("withoutEM_match",&withoutEM_match);
   test_sw->Branch("withEM_match",&withEM_match);

   test_sw->Branch("ntCl_match_wo4thPix",&ntCl_match_wo4thPix);
   test_sw->Branch("ntCl_match_wo3thPix",&ntCl_match_wo3thPix);
   test_sw->Branch("ntCl_match_wo2thPix",&ntCl_match_wo2thPix);
   test_sw->Branch("ntCl_match_wo1thPix",&ntCl_match_wo1thPix);

   test_sw->Branch("Npass_woEM_wo4thPix", &Npass_woEM_wo4thPix);
   test_sw->Branch("Npass_woEM_wo3thPix", &Npass_woEM_wo3thPix);
   test_sw->Branch("Npass_woEM_wo2thPix", &Npass_woEM_wo2thPix);
   test_sw->Branch("Npass_woEM_wo1thPix", &Npass_woEM_wo1thPix);

   test_sw->Branch("Npass_wEM_wo4thPix", &Npass_wEM_wo4thPix);
   test_sw->Branch("Npass_wEM_wo3thPix", &Npass_wEM_wo3thPix);
   test_sw->Branch("Npass_wEM_wo2thPix", &Npass_wEM_wo2thPix);
   test_sw->Branch("Npass_wEM_wo1thPix", &Npass_wEM_wo1thPix);

   test_sw->Branch("deta_L12EM_wo4thPix", &deta_L12EM_wo4thPix);
   test_sw->Branch("deta_L12EM_wo3thPix", &deta_L12EM_wo3thPix);
   test_sw->Branch("deta_L12EM_wo2thPix", &deta_L12EM_wo2thPix);
   test_sw->Branch("deta_L12EM_wo1thPix", &deta_L12EM_wo1thPix);

   test_sw->Branch("dphi_L12EM_wo4thPix", &dphi_L12EM_wo4thPix);
   test_sw->Branch("dphi_L12EM_wo3thPix", &dphi_L12EM_wo3thPix);
   test_sw->Branch("dphi_L12EM_wo2thPix", &dphi_L12EM_wo2thPix);
   test_sw->Branch("dphi_L12EM_wo1thPix", &dphi_L12EM_wo1thPix);

   test_sw->Branch("pass_Ele", &pass_Ele);
   test_sw->Branch("pass_Pos", &pass_Pos);
   test_sw->Branch("pass_ElePos", &pass_ElePos);


  //nentries = 100;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { //nentries
//   for (Long64_t jentry=0; jentry<30;jentry++) { //nentries
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //
      int passOrnot = 0;
      vector<int> flag;
      vector<int> flag1;
      vector<int> flag2;
      vector<int> flag3;
      vector<int> flag4;

      if (!(jentry%1) ) cout << "Processing entry " << jentry << "/" << nentries << endl;
      FillCutFlow("NoCut", 1.);
        
      EgN=egN;

      pass_egobjects_check = 0;
      ntnEg2 = 0;
      ntEgEt.clear();
      ntEgEta.clear();
      ntclusterIsEG.clear();
      ntCl_match.clear();
      withoutEM_match.clear();
      withEM_match.clear();
      pass_Ele.clear();
      pass_Pos.clear();
      pass_ElePos.clear();

      ntCl_match_wo4thPix.clear();
      ntCl_match_wo3thPix.clear();
      ntCl_match_wo2thPix.clear();
      ntCl_match_wo1thPix.clear();

      Npass_woEM_wo4thPix.clear();
      Npass_woEM_wo3thPix.clear();
      Npass_woEM_wo2thPix.clear();
      Npass_woEM_wo1thPix.clear();

      Npass_wEM_wo4thPix.clear();
      Npass_wEM_wo3thPix.clear();
      Npass_wEM_wo2thPix.clear();
      Npass_wEM_wo1thPix.clear();

      ntfirstPix.clear();
      ntsecondPix.clear();
      ntthirdPix.clear();
      ntfourthPix.clear();

      ntL1TkEleEt.clear();
      ntL1TkEleEta.clear();
      ntL1TkElePhi.clear();

      ntL1TkEleIsoEt.clear();
      ntL1TkEleIsoEta.clear();
      ntL1TkEleIsoPhi.clear();
    
      dphi_L12EM_wo4thPix.clear();
      dphi_L12EM_wo3thPix.clear();
      dphi_L12EM_wo2thPix.clear();
      dphi_L12EM_wo1thPix.clear();

      deta_L12EM_wo4thPix.clear();
      deta_L12EM_wo3thPix.clear();
      deta_L12EM_wo2thPix.clear();
      deta_L12EM_wo1thPix.clear();

      all_cut_pass_eg = 0;
      event_nominator = 0;
      event_denominator = 0;
   
      bHitGx1.clear();
      bHitGy1.clear();
      bHitGz1.clear();
      bHitGx2.clear();
      bHitGy2.clear();
      bHitGz2.clear();
      bHitGx3.clear();
      bHitGy3.clear();
      bHitGz3.clear();
      bHitGx4.clear();
      bHitGy4.clear();
      bHitGz4.clear();

      Eget.clear();
      Egeta.clear();
      Egphi.clear();
      Egx.clear();
      Egy.clear();
      Egz.clear();
	  
      
      // cout flow
      float EgEtCut  = 0;
      for(int k=0; k<EgN; k++) {
        EgEt =egEt ->at(k);
        if(EgEt > 9.5 ) EgEtCut = 1;
      }
      if(EgEtCut == 0) continue;
      FillCutFlow("MinEtCut", 1.);

      int EtaCutFlow = 0;
      
//      for(int h=0; h < bHitN; h++){
//cout<<"h = "<<h<<endl;
      for(int k=0; k<EgN; k++) {
        EgEta=egEta->at(k);
        EgEt =egEt ->at(k);
        //if(fabs(EgEta) < 1.3 && EgEt > 10) EtaCutFlow = 1; // for first η region
        //if(fabs(EgEta) > 1.3 && fabs(EgEta) < 1.6 && EgEt > 10) EtaCutFlow = 1;
        //if(fabs(EgEta) > 1.6 && fabs(EgEta) < 1.9 && EgEt > 10) EtaCutFlow = 1;
        //if(fabs(EgEta) > 1.9 && fabs(EgEta) < 2.5 && EgEt > 10) EtaCutFlow = 1;
//        if(fabs(EgEta) < 2.5 && EgEt > 9.5) EtaCutFlow = 1;
        if(fabs(EgEta) < 1.3 && EgEt > 9.5) EtaCutFlow = 1;
      }
      cout<<"EtaCutFlow = "<<EtaCutFlow<<endl;
      if(EtaCutFlow == 0) continue;
      FillCutFlow("EtaCut", 1.);

      //save L1TkElectron
      L1TkEleN = l1tkegN;
      for(int i = 0; i < l1tkegN; i++){
         ntL1TkEleEt.push_back(l1tkegEt->at(i));
	 ntL1TkEleEta.push_back(l1tkegEta->at(i));
         ntL1TkElePhi.push_back(l1tkegPhi->at(i));
      }

      //save L1TkElectron_Iso
      L1TkEleIsoN = l1tkegIsoN;
      for(int i = 0; i < l1tkegIsoN; i++){
         ntL1TkEleIsoEt.push_back(l1tkegIsoEt->at(i));
	 ntL1TkEleIsoEta.push_back(l1tkegIsoEta->at(i));
         ntL1TkEleIsoPhi.push_back(l1tkegIsoPhi->at(i));
      }

     // find egamma objects passing pixtrk signal windows
     for( int q=0; q<EgN; q++){ 
      EgEt =egEt ->at(q);
      EgEta=egEta->at(q);
      EgPhi=egPhi->at(q);

      float EgGx = egGx->at(q);
      float EgGy = egGy->at(q);
      float EgGz = egGz->at(q);
      emvector.SetXYZ(EgGx,EgGy,EgGz);

      if(EgEt < 9.5) continue;

      eta_region = 0; // initialize variable 
      if( fabs(EgEta) < 1.3 ) eta_region =1;
      if( fabs(EgEta) < 1.6 && fabs(EgEta) > 1.3 ) eta_region =2;
      if( fabs(EgEta) < 1.9 && fabs(EgEta) > 1.6 ) eta_region =3;
      if( fabs(EgEta) < 2.5 && fabs(EgEta) > 1.9 ) eta_region =4;
      //if( fabs(EgEta) < 2.8 && fabs(EgEta) > 2.5 ) eta_region =5;
      //if( eta_region != 4 ) continue;
      if( fabs(EgEta) > 2.5 ) continue;


      pass_egobjects_check = 1;
      ntnEg2++;
      ntEgEt.push_back(EgEt);
      ntEgEta.push_back(EgEta);

      //ntclusterIsEG.push_back(clusterIsEG->at(q));
//      ntclusterIsEG.push_back(seedIsEG->at(q));
      
      // set regin of interest
      // for η < 1.3, Δφ < 0.05 
      if(eta_region==1)SetROI(2); 
      else SetROI(eta_region);

      // initialize pixel hit variables
      first_layer_hits.clear();
      second_layer_hits.clear();
      third_layer_hits.clear();
      fourth_layer_hits.clear();

      first_layer_hits_Ele_or_Pos.clear();
      second_layer_hits_Ele_or_Pos.clear();
      third_layer_hits_Ele_or_Pos.clear();
      fourth_layer_hits_Ele_or_Pos.clear();
      hitted_layers.clear();

      
      layers[0] = 1; // beam spot
      layers[1] = 0; layers[2] = 0; layers[3] = 0; layers[4] = 0;
      r = 0;
      StorePixelHit(eta_region); // save pixel hits in Region of Interest for the given eta region

      // check which pixel has hits
       for( int i=1; i < 5; i++){ 
          if( layers[i] != 0 ){ 
            hitted_layers.push_back(i); 
          }
       }
   
       // set pixtrk signal boundary
       if(eta_region==1)SetSingalBoundary(2);
       else SetSingalBoundary(eta_region);

       // PixTRK algorithm 
       pass_count = 0;
       pass_count_wo4thPix = 0, pass_count_wo3thPix = 0, pass_count_wo2thPix = 0, pass_count_wo1thPix = 0;
       woEM_pass_count_wo4thPix = 0, woEM_pass_count_wo3thPix = 0, woEM_pass_count_wo2thPix = 0, woEM_pass_count_wo1thPix = 0;
       wEM_pass_count_wo4thPix = 0, wEM_pass_count_wo3thPix = 0, wEM_pass_count_wo2thPix = 0, wEM_pass_count_wo1thPix = 0;
       withoutEM_count_Ele = 0, withEM_count_Ele = 0;

       fourth_layer_missing = 0;
       third_layer_missing = 0;
       second_layer_missing = 0;
       first_layer_missing = 0;

       // loop over every 3 out of 4 pixel combination 
       for( std::vector<int>::iterator first_hit = hitted_layers.begin(); first_hit != hitted_layers.end(); first_hit++){
          for ( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
              for ( std::vector<int>::iterator third_hit = second_hit+1; third_hit != hitted_layers.end(); third_hit++){
                 pass_count_EleorPos = 0;
                 pass_count_Ele = 0;
                 pass_count_Pos = 0;              
                 
                 // loop over every pixel hits in the given pixel combination
                 for( int k=0; k < layers[*first_hit]; k++){
                    for( int i=0; i < layers[*second_hit]; i++){
                        _pass_Ele = 0, _pass_Pos = 0;
                        L012_pass_Ele = 0, L012_pass_Pos = 0;
                        L013_pass_Ele = 0, L013_pass_Pos = 0;
                        L023_pass_Ele = 0, L023_pass_Pos = 0;

                        if( *first_hit == 1 && *second_hit == 2 )
                          TriggeringWith_1st2ndPixel(k,i);

                        if( *first_hit == 1 && *second_hit == 3 )
                          TriggeringWith_1st3rdPixel(k,i);

                        if( *first_hit == 2 && *second_hit == 3 )
                          TriggeringWith_2nd3rdPixel(k,i);

                        // skip only if both _pass_Ele and _pass_Pos are 0 i.e., both electron and positron signal window are not satisfied
                        if( !_pass_Ele && !_pass_Pos ) continue;

                        for( int j=0; j < layers[*third_hit]; j++){
                            all_cut_pass_Ele = 0, all_cut_pass_Pos = 0;
                            withoutEM_pass_Ele = 0, withEM_pass_Ele = 0;
                            withoutEM_pass_Pos = 0, withEM_pass_Pos = 0;

                            L012_pass_Ele = 0, L012_pass_Pos = 0;
                            L013_pass_Ele = 0, L013_pass_Pos = 0;
                            L014_pass_Ele = 0, L014_pass_Pos = 0;
                            L023_pass_Ele = 0, L023_pass_Pos = 0;
                            L024_pass_Ele = 0, L024_pass_Pos = 0;
                            L034_pass_Ele = 0, L034_pass_Pos = 0;
                            L123_pass_Ele = 0, L123_pass_Pos = 0;
                            L124_pass_Ele = 0, L124_pass_Pos = 0;
                            L134_pass_Ele = 0, L134_pass_Pos = 0;
                            L234_pass_Ele = 0, L234_pass_Pos = 0;

                            L12_EM_Ele = 0, L12_EM_Pos = 0;
                            L13_EM_Ele = 0, L13_EM_Pos = 0;
                            L14_EM_Ele = 0, L14_EM_Pos = 0;
                            L23_EM_Ele = 0, L23_EM_Pos = 0;
                            L24_EM_Ele = 0, L24_EM_Pos = 0;
                            L34_EM_Ele = 0, L34_EM_Pos = 0;

            	            dPhi = StandaloneDPhi( *first_hit, *second_hit, *third_hit, k, i, j );
                            dEta = StandaloneDEta( *first_hit, *second_hit, *third_hit, k, i, j );
                            dR = sqrt( pow(dEta,2) + pow(dPhi,2) );

                              if( *first_hit == 1 && *second_hit == 2 && *third_hit == 3 ){ // for efficiency counting  !!caution of the position of this codition
                                // This is for the case that the first hit is in the first pixel layer and the second hit is in the second pixel layer and the third hit is in the third layer. 
                                TriggeringWithout_4thPixel(k, i, j);

                                if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[j] == 1 || third_layer_hits_Ele_or_Pos[j] ==3)){
                                      if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele && L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                         all_cut_pass_Ele = 1; 
                                      if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele)
                                         withoutEM_pass_Ele = 1;
                                      if(L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                         withEM_pass_Ele = 1;
                                }
 
                                if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[j] == 2 || third_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos && L12_EM_Pos && L13_EM_Pos && L23_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L12_EM_Pos && L13_EM_Pos && L23_EM_Pos) withEM_pass_Pos = 1;
                                 }

				if( withEM_pass_Ele == 1 ){
				    float r1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
				    float r2 = sqrt(pow(second_layer_hits[i].X(),2)+pow(second_layer_hits[i].Y(),2));
				    float r3 = sqrt(pow(third_layer_hits[j].X(),2)+pow(third_layer_hits[j].Y(),2));
				    if( r1 < 5.5 ){
				    bHitGx1.push_back(first_layer_hits[k].X());
				    bHitGy1.push_back(first_layer_hits[k].Y());
				    bHitGz1.push_back(first_layer_hits[k].Z());
				    }
				    if( r2 > 5.5 && r2 < 8.5 ){
				    bHitGx2.push_back(second_layer_hits[i].X());
				    bHitGy2.push_back(second_layer_hits[i].Y());
				    bHitGz2.push_back(second_layer_hits[i].Z());
				    }
				    if( r3 > 8.5 && r3 < 13 ){
				    bHitGx3.push_back(third_layer_hits[j].X());
				    bHitGy3.push_back(third_layer_hits[j].Y());
				    bHitGz3.push_back(third_layer_hits[j].Z());
				    }

				    Eget.push_back(EgEt);
				    Egeta.push_back(EgEta);
				    Egphi.push_back(EgPhi);
				    Egx.push_back(egGx->at(q));
				    Egy.push_back(egGy->at(q));
				    Egz.push_back(egGz->at(q));

				    cout<<"event : "<<jentry<<endl;
				    cout<<"EgEta : "<< EgEta <<endl;
				    cout<<"EgPhi : "<< EgPhi <<endl;
				    cout<<"======================================="<<endl;
				}

                               if( all_cut_pass_Ele || all_cut_pass_Pos){
                                  pass_count_wo4thPix = 1;
                                }

                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo4thPix++;
                                }
                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo4thPix++;
                                }
				//
			      

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit loop
                                  cout << "without fourth layer: " << (k+1) * (i+1) * ( j+1) << endl;
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
                                  fourth_layer_missing = 1;
                   		 }
                               }
                              }//if(1,2,3layer)

                              if( *first_hit == 1 && *second_hit == 2 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                // This is for the case that the first hit is in the first pixel layer and the second is in the second layer and the third hit is in the fourth layer.
                                TriggeringWithout_3rdPixel(k, i, j);

                                if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele && L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) all_cut_pass_Ele = 1;
                                    if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) withEM_pass_Ele = 1;
                                } 

                                if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos && L12_EM_Pos && L14_EM_Pos && L24_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L12_EM_Pos && L14_EM_Pos && L24_EM_Pos) withEM_pass_Pos = 1;
                                }

				if( withEM_pass_Ele == 1 ){
				    float r1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
				    float r2 = sqrt(pow(second_layer_hits[i].X(),2)+pow(second_layer_hits[i].Y(),2));
				    float r4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
				    if( r1 < 5.5 ){
				    bHitGx1.push_back(first_layer_hits[k].X());
				    bHitGy1.push_back(first_layer_hits[k].Y());
				    bHitGz1.push_back(first_layer_hits[k].Z());
				    }
				    if( r2 > 5.5 && r2 < 8.5 ){
				    bHitGx2.push_back(second_layer_hits[i].X());
				    bHitGy2.push_back(second_layer_hits[i].Y());
				    bHitGz2.push_back(second_layer_hits[i].Z());
				    }
				    if( r4 > 13 && r4 < 18 ){
				    bHitGx4.push_back(fourth_layer_hits[j].X());
				    bHitGy4.push_back(fourth_layer_hits[j].Y());
				    bHitGz4.push_back(fourth_layer_hits[j].Z());
				    }
				    Eget.push_back(EgEt);
				    Egeta.push_back(EgEta);
				    Egphi.push_back(EgPhi);
				    Egx.push_back(egGx->at(q));
				    Egy.push_back(egGy->at(q));
				    Egz.push_back(egGz->at(q));
				    cout<<"event : "<<jentry<<endl;
				    cout<<"EgEta : "<< EgEta <<endl;
				    cout<<"EgPhi : "<< EgPhi <<endl;
				    cout<<"======================================="<<endl;

				}

                                if( all_cut_pass_Ele || all_cut_pass_Pos ){
                                  pass_count_wo3thPix = 1;
                                }
                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo3thPix++;
                                }
                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo3thPix++;
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                  cout << "without third layer: " << (k+1) * (i+1) * ( j+1) << endl;
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
            	     	          third_layer_missing = 1;
                   		 }
                               }
                              }//if 1,2,4
		    
                              if( *first_hit == 1 && *second_hit == 3 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                TriggeringWithout_2ndPixel(k, i, j);

                                if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 1 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele && L13_EM_Ele && L14_EM_Ele && L34_EM_Ele)all_cut_pass_Ele = 1; 
                                    if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L13_EM_Ele && L14_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                }

                                if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos && L13_EM_Pos && L14_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L13_EM_Pos && L14_EM_Pos && L34_EM_Pos) withEM_pass_Pos = 1;
                                 }

				if( withEM_pass_Ele == 1 ){
				    float r1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
				    float r3 = sqrt(pow(third_layer_hits[i].X(),2)+pow(third_layer_hits[i].Y(),2));
				    float r4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
				    if( r1 < 5.5 ){
				    bHitGx1.push_back(first_layer_hits[k].X());
				    bHitGy1.push_back(first_layer_hits[k].Y());
				    bHitGz1.push_back(first_layer_hits[k].Z());
				    }
				    if( r3 > 8.5 && r3 < 13 ){
				    bHitGx3.push_back(third_layer_hits[i].X());
				    bHitGy3.push_back(third_layer_hits[i].Y());
				    bHitGz3.push_back(third_layer_hits[i].Z());
				    }
				    if( r4 > 13 && r4 < 18 ){
				    bHitGx4.push_back(fourth_layer_hits[j].X());
				    bHitGy4.push_back(fourth_layer_hits[j].Y());
				    bHitGz4.push_back(fourth_layer_hits[j].Z());
				    }
				    Eget.push_back(EgEt);
				    Egeta.push_back(EgEta);
				    Egphi.push_back(EgPhi);
				    Egx.push_back(egGx->at(q));
				    Egy.push_back(egGy->at(q));
				    Egz.push_back(egGz->at(q));
				    cout<<"event : "<<jentry<<endl;
				    cout<<"EgEta : "<< EgEta <<endl;
				    cout<<"EgPhi : "<< EgPhi <<endl;
				    cout<<"======================================="<<endl;
				}

                                if( all_cut_pass_Ele || all_cut_pass_Pos){
                                  pass_count_wo2thPix = 1;
                                }
                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo2thPix++;
                                }

                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo2thPix++;
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                  cout << "without second layer: " << (k+1) * (i+1) * ( j+1) << endl;
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
                                  second_layer_missing = 1; 
                   		 }
                               }
                              }//if 1,3,4

                              if( *first_hit == 2 && *second_hit == 3 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                TriggeringWithout_1stPixel(k, i, j);

                                if( (second_layer_hits_Ele_or_Pos[k] == 1 || second_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 1 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele && L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) all_cut_pass_Ele = 1;
                                    if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                } 

                                if( (second_layer_hits_Ele_or_Pos[k] == 2 || second_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos && L23_EM_Pos && L24_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L23_EM_Pos && L24_EM_Pos && L34_EM_Pos) withEM_pass_Pos = 1;
                                }

				if( withEM_pass_Ele == 1 ){
				    float r2 = sqrt(pow(second_layer_hits[k].X(),2)+pow(second_layer_hits[k].Y(),2));
				    float r3 = sqrt(pow(third_layer_hits[i].X(),2)+pow(third_layer_hits[i].Y(),2));
				    float r4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
				    if( r2 > 5.5 && r2 < 8.5 ){
				    bHitGx2.push_back(second_layer_hits[k].X());
				    bHitGy2.push_back(second_layer_hits[k].Y());
				    bHitGz2.push_back(second_layer_hits[k].Z());
				    }
				    if( r3 > 8.5 && r3 < 13 ){
				    bHitGx3.push_back(third_layer_hits[i].X());
				    bHitGy3.push_back(third_layer_hits[i].Y());
				    bHitGz3.push_back(third_layer_hits[i].Z());
				    }
				    if( r4 > 13 && r4 < 18 ){
				    bHitGx4.push_back(fourth_layer_hits[j].X());
				    bHitGy4.push_back(fourth_layer_hits[j].Y());
				    bHitGz4.push_back(fourth_layer_hits[j].Z());
				    }
				    Eget.push_back(EgEt);
				    Egeta.push_back(EgEta);
				    Egphi.push_back(EgPhi);
				    Egx.push_back(egGx->at(q));
				    Egy.push_back(egGy->at(q));
				    Egz.push_back(egGz->at(q));
				    cout<<"event : "<<jentry<<endl;
				    cout<<"EgEta : "<< EgEta <<endl;
				    cout<<"EgPhi : "<< EgPhi <<endl;
				    cout<<"======================================="<<endl;
				}
                                if( all_cut_pass_Ele || all_cut_pass_Pos){
                                  pass_count_wo1thPix = 1;
                                }
                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo1thPix++;
                                }
                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo1thPix++;
                                }
                    
                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                  cout << "without first layer: " << (k+1) * (i+1) * ( j+1) << endl;
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
                                  first_layer_missing = 1;
                   		 } 
                               }
                              }

                         if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1) {pass_count_EleorPos++; pass_count = 1;}
                         if( all_cut_pass_Ele == 1 ) {pass_count_Ele++; }
                         if( all_cut_pass_Pos == 1 ) pass_count_Pos++;
                         if( withoutEM_pass_Ele == 1 ) withoutEM_count_Ele = 1;
                         if( withEM_pass_Ele == 1 ) withEM_count_Ele = 1;
                       } // loop for third layer hits
                   } // loop for second layer hits       
                 } // loop for first layer hits

                 pass_Ele.push_back(pass_count_Ele);
                 pass_Pos.push_back(pass_count_Pos);
                 pass_ElePos.push_back(pass_count_EleorPos);
          }          
        }
      }
      if( pass_count ){ 
          ntCl_match.push_back(true);
          all_cut_pass_eg = 1; 
      }
      else ntCl_match.push_back(false);

      if( pass_count_wo4thPix ){
          ntCl_match_wo4thPix.push_back(true);
      }
      else ntCl_match_wo4thPix.push_back(false);

      Npass_woEM_wo4thPix.push_back(woEM_pass_count_wo4thPix);
      Npass_wEM_wo4thPix.push_back(wEM_pass_count_wo4thPix);

      if( pass_count_wo3thPix ){
          ntCl_match_wo3thPix.push_back(true);
      }
      else ntCl_match_wo3thPix.push_back(false);

      Npass_woEM_wo3thPix.push_back(woEM_pass_count_wo3thPix);
      Npass_wEM_wo3thPix.push_back(wEM_pass_count_wo3thPix);

      if( pass_count_wo2thPix ){
          ntCl_match_wo2thPix.push_back(true);
      }
      else ntCl_match_wo2thPix.push_back(false);

      Npass_woEM_wo2thPix.push_back(woEM_pass_count_wo2thPix);
      Npass_wEM_wo2thPix.push_back(wEM_pass_count_wo2thPix);

      if( pass_count_wo1thPix ){
          ntCl_match_wo1thPix.push_back(true);
      }
      else ntCl_match_wo1thPix.push_back(false);

      Npass_woEM_wo1thPix.push_back(woEM_pass_count_wo1thPix);
      Npass_wEM_wo1thPix.push_back(wEM_pass_count_wo1thPix);

      if(withoutEM_count_Ele){
         withoutEM_match.push_back(true);
      }
      else withoutEM_match.push_back(false);

      if(withEM_count_Ele){
         withEM_match.push_back(true);
      }
      else withEM_match.push_back(false);
         /////////////////////////////////
  } // end of egamma loop 
//      } //end of bHit loop
  if(pass_egobjects_check){ event_denominator = 1; FillCutFlow("EvtCut", 1.);}
  if(all_cut_pass_eg) event_nominator = 1; 
     test_sw->Fill();
 } // end of entries loop 
   test_file->Write();
}

