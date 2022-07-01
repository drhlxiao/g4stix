
const double XBINS[3]={1000, 0,  250 };
const double YBINS[3]={1000, 0,  250 };
//response matrix binning

void convert(TString fname, TString fname_out="drm.root", double fluence=1){

   TFile *f = new TFile(fname);
   TTree *tree=(TTree*)f->Get("events");

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
   Long64_t nbytes = 0;
   TH2F *hresp[32];
   Long64_t i;
   for(i=0;i<32;i++)
   {
	   hresp[i]=new TH2F(Form("hresp_%d",i),Form("hresp_%d; Energy (keV); Deposited energy (keV)",i), 
			   XBINS[0], XBINS[1], XBINS[2],
			   YBINS[0], YBINS[1], YBINS[2]);
   }
   int j=0;
   for (i=0; i<nentries;i++) {
      nbytes += events->GetEntry(i);
	  for(j=0;j<384;j++){
		  int det=j%32;
		  if(edep[j]>0)hresp[det]->Fill(E0, edep[j]);
	  }
   }
   TFile fout(fname_out,"recreate");
   for(j=0;j<32;j++){
	   hresp[j]->Write();
   }
   fout->Close();
   TCanvas c1;
   for(j=0;j<32;j++){
	   hresp[j]->Draw("colz");
	   gPad->SetLogZ();
    gPad->WaitPrimitive();
   }

}
