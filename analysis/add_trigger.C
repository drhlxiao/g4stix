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
const int N = 37;
using namespace std; 

void updateTree(TString filein)
{
	TFile *  f = new TFile(filein, "update");
	TTree*	events=(TTree*)f->Get("events");
	Double_t        edep[384];
	Int_t nHits[32];
	TBranch *tb=events->Branch("nHits",nHits, "nHits[32]/I");
	Long64_t nentries = events->GetEntries();
	int det;
	int j=0;
	for (Long64_t i=0; i<nentries;i++) {
		nbytes += events->GetEntry(i);

		if(i%10000==0)cout<<100*i/(nentries+0.0)<<endl;
		for(j=0;j<32;j++)nHits[j]=0;

		for(k=0;k<384;k++){
			if(charge2[k]>=4){
				det=(int)(k/12);
				nHits[det]++;
			}
		};
		tb->Fill();
	}
	events->Write();
	f.Close();
	cout<<"Output file:"<<filein<<endl;

}

void Help(){
	cout<<"make_matrix  -i input.root updateTree"<<endl; 
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
	updateTree(inputFilename);

	return 0;
} 
