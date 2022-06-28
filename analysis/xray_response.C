{
	//////////////////////////////////////////////////////////
	//   This file has been automatically generated 
	//     (Thu Jun  3 22:44:24 2021 by ROOT version6.22/00)
	//   from TTree events/events
	//   found on file: xray_response.root
	//////////////////////////////////////////////////////////


	//Reset ROOT and connect tree file
	gROOT->Reset();
	TFile *f = new TFile("xray_response.root");
	TTree *events=(TTree*)f->Get("events");

	//Declaration of leaves types
	Double_t        edep[384];
	Double_t        collect[384];
	Double_t        charge[384];
	Double_t        charge2[384];
	Int_t           eventID;
	Double_t        particle_position[3];
	Double_t        particle_direction[3];
	Double_t        particle_energy;
	Int_t           itrack;
	Int_t           parent[30];
	Int_t           pdg[30];
	Int_t           hitPixelID[30];
	Double_t        energy[30];
	Double_t        time[30];

	// Set branch addresses.
	events->SetBranchAddress("edep",edep);
	events->SetBranchAddress("collect",collect);
	events->SetBranchAddress("charge",charge);
	events->SetBranchAddress("charge2",charge2);
	events->SetBranchAddress("eventID",&eventID);
	events->SetBranchAddress("particle_position",particle_position);
	events->SetBranchAddress("particle_direction",particle_direction);
	events->SetBranchAddress("particle_energy",&particle_energy);
	events->SetBranchAddress("num_tracks",&itrack);
	events->SetBranchAddress("parent",parent);
	events->SetBranchAddress("pdg",pdg);
	events->SetBranchAddress("hitPixelID",hitPixelID);
	events->SetBranchAddress("energy",energy);
	events->SetBranchAddress("time",time);

	//     This is the loop skeleton
	//       To read only selected branches, Insert statements like:
	// events->SetBranchStatus("*",0);  // disable all branches
	// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

	Long64_t nentries = events->GetEntries();
	TH1F h("h","h",249,1,250);

	Long64_t nbytes = 0;
	   for (Long64_t i=0; i<nentries;i++) {
	      nbytes += events->GetEntry(i);
		  double edepSum=0;
		  for(int j=0;j<384;j++){
			  edepSum+=edep[j];
		  }
		  if(edepSum>2){
			  h.Fill(particle_energy/1e3);
		  }
	   }
	   h.Draw("ep");
}
