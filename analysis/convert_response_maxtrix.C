
#include "TPad.h"
//response matrix binning
#include <TApplication.h> 
#include <TSystem.h> 
#include <TGraph.h>
#include <TROOT.h> 
#include <TFile.h>  
#include <TTree.h>  
#include <TCanvas.h>  
#include <TF1.h>  
#include <TH1F.h> 
 #include <TH2F.h> 

#include <iostream>
#include <string>    
using namespace std; 
const double XBINS[3]={1000, 0,  250 };
const double YBINS[3]={1000, 0,  250 };

void convert(TString fname, TString fname_out="drm.root", double fluence=1){

   TFile *f = new TFile(fname);
   TTree *events=(TTree*)f->Get("events");

   Double_t        edep[384];
   Double_t        collected[384];
   Double_t        charge[384];
   Double_t        charge2[384];
   Int_t           eventID;
   Double_t        gunPos[3];
   Double_t        gunVec[3];
   Double_t        E0;
   Int_t           numTracks;
   Double_t        hitx[30];
   Double_t        hity[30];
   Double_t        hitz[30];
   Int_t           parent[30];
   Int_t           pdg[30];
   Int_t           pixel[30];
   Double_t        energy[30];
   Double_t        time[30];

   // Set branch addresses.
   events->SetBranchAddress("edep",edep);
   events->SetBranchAddress("collected",collected);
   events->SetBranchAddress("charge",charge);
   events->SetBranchAddress("charge2",charge2);
   events->SetBranchAddress("eventID",&eventID);
   events->SetBranchAddress("gunPos",gunPos);
   events->SetBranchAddress("gunVec",gunVec);
   events->SetBranchAddress("E0",&E0);
   events->SetBranchAddress("numTracks",&numTracks);
   events->SetBranchAddress("hitx",hitx);
   events->SetBranchAddress("hity",hity);
   events->SetBranchAddress("hitz",hitz);
   events->SetBranchAddress("parent",parent);
   events->SetBranchAddress("pdg",pdg);
   events->SetBranchAddress("pixel",pixel);
   events->SetBranchAddress("energy",energy);
   events->SetBranchAddress("time",time);

   Long64_t nentries = events->GetEntries();
   cout<<"entries:"<<nentries<<endl;
   Long64_t nbytes = 0;
   TH2F *hresp[32];
   TH2F *hresp_mean;
   TH2F *hresp_fine;
   Long64_t i;
   hresp_mean=new TH2F(Form("hresp_mean"),Form("response matrix for detectors except BKG, CFL and fine grids; Energy (keV); Deposited energy (keV)"), 
			   XBINS[0], XBINS[1], XBINS[2],
			   YBINS[0], YBINS[1], YBINS[2]);
   hresp_fine=new TH2F(Form("hresp_fine"),Form("response matrix for detectors with fine grids; Energy (keV); Deposited energy (keV)"), 
			   XBINS[0], XBINS[1], XBINS[2],
			   YBINS[0], YBINS[1], YBINS[2]);
   for(i=0;i<32;i++)
   {
	   hresp[i]=new TH2F(Form("hresp_%ld",i),Form("hresp_%ld; Energy (keV); Deposited energy (keV)",i), 
			   XBINS[0], XBINS[1], XBINS[2],
			   YBINS[0], YBINS[1], YBINS[2]);
   }
   int find_grids[]={11,12,13, 17,18,19};
   int j=0;
   int k=0;
   for (i=0; i<nentries;i++) {
      nbytes += events->GetEntry(i);
	  for(j=0;j<384;j++){
		  int det=j%32;
		   bool is_fine_grid=false;
		  if(edep[j]>0)hresp[det]->Fill(E0, edep[j]);
		  if(det ==8 ||det==9)continue;
		  for(k=0;k<6;k++){
			  if(det == find_grids[k]){
				  is_fine_grid=true;
			  }
		  }
		  if(is_fine_grid){
			  hresp_fine->Fill(E0, edep[j]);
		  }
		  else{
			  hresp_mean->Fill(E0, edep[j]);
		  }
		  
	  }
   }
   TFile fout(fname_out,"recreate");
   fout.cd();
   for(j=0;j<32;j++){
	   hresp[j]->Write();
   }
   hresp_fine->Write();
   hresp_mean->Write();
   fout.Close();
   TCanvas c1;
   for(j=0;j<32;j++){
	   hresp[j]->Draw("colz");
	   gPad->SetLogz();
    gPad->WaitPrimitive();
   }

}
