{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Wed Jun  8 17:38:00 2022 by ROOT version6.22/08)
//   from TTree grids/grids
//   found on file: stix_grid_parameters.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("stix_grid_parameters.root");
   if (!f) {
      f = new TFile("stix_grid_parameters.root");
   }
    f->GetObject("grids",tree);

//Declaration of leaves types
   Int_t           is_front;
   Int_t           det_idx;
   Int_t           i_strip;
   Float_t         polygon_x[5];
   Float_t         polygon_y[5];
   Int_t           is_nominal_parameters;

   // Set branch addresses.
   grids->SetBranchAddress("is_front",&is_front);
   grids->SetBranchAddress("det_idx",&det_idx);
   grids->SetBranchAddress("i_strip",&i_strip);
   grids->SetBranchAddress("polygon_x",polygon_x);
   grids->SetBranchAddress("polygon_y",polygon_y);
   grids->SetBranchAddress("is_nominal_parameters",&is_nominal_parameters);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// grids->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = grids->GetEntries();

   Long64_t nbytes = 0;
//   for (Long64_t i=0; i<nentries;i++) {
//      nbytes += grids->GetEntry(i);
//   }
}
