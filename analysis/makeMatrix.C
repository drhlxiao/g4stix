#include <iostream>
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
#include <string>    
using namespace std; 
int getScienceBin(Double_t energy){

	double energyRanges[]={0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150, 1e9};
	if(energy<4)return 0;
	else if(energy>=150)return 31;
	else{
		for(int i=0;i<31;i++)
		{
			if(energy>=energyRanges[i] && energy<energyRanges[i+1])return i;
		}
	}
	return 31;

}

void makeMatrix(TString filein,  TString fout,
		double eStep=0.1,
		double maxEnergy=250 , Long64_t entries=0, double flux=1)
{
	TFile *  f = new TFile(filein);
	//f->GetObject("events",tree);
	TTree*	events=(TTree*)f->Get("events");

	//Declaration of leaves types
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

	//     This is the loop skeleton
	//       To read only selected branches, Insert statements like:
	// events->SetBranchStatus("*",0);  // disable all branches
	// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname


	int num=(int)(maxEnergy/eStep);
	cout<<"Number of bins:"<<num<<endl;
	cout<<"Max Energy :"<<maxEnergy<<endl;
	cout<<"Output:"<<fout<<endl;
	TFile fo(fout,"recreate");

	TH2F *hresp[32];
	TH2F *hcoll[32];
	TH2F *hreal[32];

	TH2F *hresp_stix[32];
	TH2F *hcoll_stix[32];
	TH2F *hreal_stix[32];

	int k=0;
	cout<<"Creating histograms"<<endl;
	for(k=0;k<32;k++) {
		hresp[k]=new TH2F(Form("hresp_edep_%d", k),Form("Energy depositions (Det%d) - flux: %.2f ph/(cm2*keV); Energy (keV); Deposited energy (keV)",k, flux), num, 0, maxEnergy, num, 0, maxEnergy);
		hcoll[k]=new TH2F(Form("hresp_coll_%d", k),Form("collected energies (Det%d) - flux: %.2f ph/(cm2*keV); Energy (keV); Collected Energy (keV)",k, flux), num, 0, maxEnergy, num, 0, maxEnergy);
		hreal[k]=new TH2F(Form("hresp_real_%d", k),Form("Recorded energies (Det%d) - flux: %.2f ph/(cm2*keV); Energy (keV); Recorded energy (keV) ",k, flux), num, 0, maxEnergy, num,0, maxEnergy);

		hresp_stix[k]=new TH2F(Form("hresp_edep_stix_%d", k),Form("Energy depositions (Det%d); Energy (keV); Deposited energy (science bin)",k), num, 0, maxEnergy, 33, 0, 33);
		hcoll_stix[k]=new TH2F(Form("hresp_coll_stix_%d", k),Form("Collected energies (Det%d); Energy (keV); Collected Energy (science bin)",k), num, 0, maxEnergy, 33, 0, 33);
		hreal_stix[k]=new TH2F(Form("hresp_real_stix_%d", k),Form("Recorded energies (Det%d); Energy (keV); Recorded energy (science bin) ",k), num, 0, maxEnergy, 33,0, 33);

	}

	Long64_t nentries = events->GetEntries();
	cout<<"Number of enetries:"<<nentries<<endl;
	if(entries>0)nentries=entries;

	Long64_t nbytes = 0;
	double weight=1./flux;
	for (Long64_t i=0; i<nentries;i++) {
		nbytes += events->GetEntry(i);
		if(i%10000==0)cout<<100*i/(nentries+0.0)<<endl;
		for(k=0;k<384;k++){
			if(edep[k]<=0)continue;
			int det=(int)(k/12);
			hresp[det]->Fill(E0, edep[k], weight);
			hcoll[det]->Fill(E0, collected[k], weight);
			hreal[det]->Fill(E0, charge2[k], weight);
			hresp_stix[det]->Fill(E0, getScienceBin(edep[k]));
			hcoll_stix[det]->Fill(E0, getScienceBin(collected[k]));
			hreal_stix[det]->Fill(E0, getScienceBin(charge2[k]));
		}
	}
	fo.cd();
	cout<<"Writing histograms"<<endl;
	for(k=0;k<32;k++) {
		hresp[k]->Write();
		hcoll[k]->Write();
		hreal[k]->Write();
		hresp_stix[k]->Write();
		hcoll_stix[k]->Write();
		hreal_stix[k]->Write();
	}
	fo.Close();
	cout<<"Output file:"<<fout<<endl;

}

void Help(){
	cout<<"make_matrix  -i input.root -o output.root  -b <energy step 0.1 >  -m <max energy, 250> -f  <flux, in units of photons/(cm2*keV)>"<<endl; 
}
int main(int argc, char *argv[])
{
	TString outputFilename = "defaultOuput.root";
	TString inputFilename;

	if (argc == 1)  Help();
	int s = 0;
	TString sel;

	double eStep=0.1;
	double eMax=250;
	Long64_t entries=0;
	double flux=1;
	
	while (s < argc - 1) {

		sel = argv[++s];
		if (sel == "-h" || sel == "--help" || sel == "--h") {
			Help();
			return 0;

		} else if (sel == "-o") {
			outputFilename = argv[++s];
			if (!outputFilename.Contains(".root")) {
				Help();
				return 0;
			}
		} else if (sel == "-i") {
			inputFilename= argv[++s];
			if (!inputFilename.Contains(".root")) {
				Help();
				return 0;
			}
		} else if (sel == "-b") 
		{
			eStep= atof(argv[++s]);
		}
		else if (sel == "-m") {
			eMax= atof(argv[++s]);
		}
		else if (sel == "-n") {
			entries= atoi(argv[++s]);
		}
		else if (sel == "-f") {
			flux = atof(argv[++s]);
		}

	}
	makeMatrix(inputFilename, outputFilename, eStep, eMax, entries, flux);

	return 0;
} 
