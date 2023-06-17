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
const int N = 38;
const float threshold=4;
using namespace std; 
int getScienceBin(Double_t energy){
	double energyRanges[]={0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 25, 
		28, 32, 36, 40, 45, 50, 56, 63, 70, 76, 84, 100, 120, 150, 1e9};
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
bool isFineGrid(int det){
	if( (10 <= det && det <=12)|| (16 <=det&& det<=18))return true;
	else return false;
}
void normMatrix(TH2F *hrep, TH1F *hspec){
	//normalized response matrix by input spectrum
	//not neccessary a flat spectrum
	int nbins=hspec->GetXaxis()->GetNbins();
	int nbinsX=hrep->GetXaxis()->GetNbins();
	int nbinsY=hrep->GetYaxis()->GetNbins();
	for(int i=1;i<nbinsX+1;i++){
		for(int j=1;i<nbinsY+1;j++){
			double nEvents=hspec->GetBinContent(i);
			double counts=hrep->GetBinContent(i,j);
			if(nEvents>0)hrep->SetBinContent(i,j, counts/nEvents);
		}
	}


}

void makeMatrix(TString filein,  TString fout,
		double eStep=0.1,
		double maxEnergy=250 , Long64_t entries=0, double radius=1, bool excludeDoubleHits=false)
{
	int num=(int)(maxEnergy/eStep);
	cout<<"Number of bins:"<<num<<endl;
	cout<<"Max Energy :"<<maxEnergy<<endl;

	TFile *  f = new TFile(filein);
	//f->GetObject("events",tree);
	TTree*	source=(TTree*)f->Get("source");
	Double_t        E0;
	Int_t           eventID;
	source->SetBranchAddress("E0",&E0);
	source->SetBranchAddress("eventID",&eventID);
	Int_t totalEntries=source->GetEntries();
	Int_t totalEvents;
	TH1F *hsource=new TH1F("hSource", Form("Source energy spectrum(Events: %d, radius: %f); Energy (keV); Counts", totalEvents, radius), num, 0, maxEnergy);
	for(int i=0;i<totalEntries; i++){
		source->GetEntry(i);
		totalEvents=eventID;
		hsource->Fill(E0);
	}
	// photons/cm^2*keV
	cout<<"Output:"<<fout<<endl;
	TFile fo(fout,"recreate");
	hsource->Write();

	Double_t flux=totalEvents/(3.1415*radius*radius*maxEnergy);

	TTree*	events=(TTree*)f->Get("events");

	//Declaration of leaves types
	Double_t        edep[384];
	Double_t        collected[384];
	Double_t        charge[384];
	Double_t        charge2[384];
	Double_t        gunPos[3];
	Double_t        gunVec[3];
	Int_t           numTracks;
	Double_t        hitx[30];
	Double_t        hity[30];
	Double_t        hitz[30];
	Int_t           parent[30];
	Int_t           pdg[30];
	Int_t           pixel[30];
	Double_t        energy[30];
	Double_t        time[30];
	Int_t nHits[32];

	// Set branch addresses.
	events->SetBranchStatus("*",0);
	events->SetBranchStatus("edep",1);  // activate branchname
	events->SetBranchStatus("collected",1);  // activate branchname
	events->SetBranchStatus("charge",1);  // activate branchname
	events->SetBranchStatus("charge2",1);  
	events->SetBranchStatus("E0",1);  
	events->SetBranchStatus("nHits",1);  

	events->SetBranchAddress("edep",edep);
	events->SetBranchAddress("collected",collected);
	events->SetBranchAddress("charge",charge);
	events->SetBranchAddress("charge2",charge2);
	events->SetBranchAddress("eventID",&eventID);
	events->SetBranchAddress("E0",&E0);
	events->SetBranchAddress("nHits",nHits);

	//     This is the loop skeleton
	//       To read only selected branches, Insert statements like:
	// events->SetBranchStatus("*",0);  // disable all branches
	// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname



	TH2F *hresp[N];
	TH2F *hcoll[N];
	TH2F *hreal[N];

	TH1F *hspc_resp[N];
	TH1F *hspc_coll[N];
	TH1F *hspc_real[N];
	TH1F *hspc_stix[N];

	TH2F *h_pix_edep[12];
	TH2F *h_pix_real[12];
	int k=0;
	const char *dh=excludeDoubleHits? "double hits excluded":"double hits included";

	for(k=0;k<12;k++) {
		h_pix_edep[k]=new TH2F(Form("h_pix_resp_edep_%d", k),Form("Energy depositions (pix %d) - flux: %.2f ph/(cm2*keV) - %s; Energy (keV); Deposited energy (keV)",k, flux, dh), num, 0, maxEnergy, num, 0, maxEnergy);
		h_pix_real[k]=new TH2F(Form("h_pix_resp_real_%d", k),Form("Energy depositions (pix %d) - flux: %.2f ph/(cm2*keV) - %s ; Energy (keV); Deposited energy (keV)",k, flux,dh), num, 0, maxEnergy, num, 0, maxEnergy);
	}

	TH2F *hresp_stix[N];
	TH2F *hcoll_stix[N];
	TH2F *hreal_stix[N];

	cout<<"Creating histograms"<<endl;
	for(k=0;k<N;k++) {

		hspc_resp[k]=new TH1F(Form("hspc_resp%d", k),
				Form("Energy depositions (Det%d) - flux: %.2f ph/(cm2*keV) - %s ; Energy (keV); counts",k, flux,dh),
				num, 0, maxEnergy);


		hresp[k]=new TH2F(Form("hresp_edep_%d", k),
				Form("Energy depositions (Det%d) - flux: %.2f ph/(cm2*keV) - %s ; Energy (keV); Deposited energy (keV)",k, flux,dh),
				num, 0, maxEnergy, num, 0, maxEnergy);
		hcoll[k]=new TH2F(Form("hresp_coll_%d", k),
				Form("collected energies (Det%d) - flux: %.2f ph/(cm2*keV) - %s ; Energy (keV); Collected Energy (keV)",k, flux,dh), 
				num, 0, maxEnergy, num, 0, maxEnergy);
		hreal[k]=new TH2F(Form("hresp_real_%d", k),
				Form("Recorded energies (Det%d) - flux: %.2f ph/(cm2*keV) - %s; Energy (keV); Recorded energy (keV) ",k, flux,dh),
				num, 0, maxEnergy, num,0, maxEnergy);
		hresp_stix[k]=new TH2F(Form("hresp_edep_stix_%d", k),
				Form("Energy depositions (Det%d) - %s ; Energy (keV); Deposited energy (science bin)",k,dh),
				num, 0, maxEnergy, 33, 0, 33);
		hcoll_stix[k]=new TH2F(Form("hresp_coll_stix_%d", k),
				Form("Collected energies (Det%d) -%s ; Energy (keV); Collected Energy (science bin)",k,dh), 
				num, 0, maxEnergy, 33, 0, 33);
		hreal_stix[k]=new TH2F(Form("hresp_real_stix_%d", k),
				Form("Recorded energies (Det%d) -%s ; Energy (keV); Recorded energy (science bin) ",k,dh),
				num, 0, maxEnergy, 33,0, 33);
	}
	TString titles[]={
		"All detector summed response", //32
		"Fine grid detector summed response", //33
		"Detectors summed response (except fine grids, CFL and BKG)", //34
		"CFL response ",//35
		"BKG response ",
		"All big pixel summed (except CFL and BKG)",
	}; //36

	for(k=32;k<N;k++){
		TString title=titles[k-32];
		hresp[k]->SetTitle(title);
		hcoll[k]->SetTitle(title);
		hreal[k]->SetTitle(title);
		hresp_stix[k]->SetTitle(title);
		hcoll_stix[k]->SetTitle(title);
		hreal_stix[k]->SetTitle(title);
	}


	Long64_t nentries = events->GetEntries();
	cout<<"Number of enetries:"<<nentries<<endl;
	if(entries>0)nentries=entries;

	Long64_t nbytes = 0;
	int det=0;
	double weight=1./flux;
	cout<<"Weight: "<<weight<<endl;
	int detHits=0;
	int kk=0;
	for (Long64_t i=0; i<nentries;i++) {
		nbytes += events->GetEntry(i);
		if(i%10000==0)cout<<100*i/(nentries+0.0)<<endl;


		for(k=0;k<384;k++){
			if(edep[k]<=4)continue;
			if(charge2[k]<=0)continue;

			det=(int)(k/12);


			if(excludeDoubleHits){
				detHits=0;
				for(kk=0;kk<12;kk++){
					if( charge2[det*12+kk] >threshold)detHits++;
				}
				if(detHits>=2)continue;
				//not triggered
			}


			//sum response
			hresp[det]->Fill(E0, edep[k]);
			hcoll[det]->Fill(E0, collected[k]);
			hreal[det]->Fill(E0, charge2[k]);
			hresp_stix[det]->Fill(E0, getScienceBin(edep[k]));
			hcoll_stix[det]->Fill(E0, getScienceBin(collected[k]) );
			hreal_stix[det]->Fill(E0, getScienceBin(charge2[k]) );

			int pix=k%12;

			h_pix_edep[pix]->Fill(E0, edep[k] );
			h_pix_real[pix]->Fill(E0, charge2[k] );

			int sumDetID=32;
			//fill all detectors

			//cout<<det<<","<<pix<<"," << E0 <<" , "<<edep[k]<<", "<<<<endl;
			hresp[sumDetID]->Fill(E0, edep[k] );
			hcoll[sumDetID]->Fill(E0, collected[k] );
			hreal[sumDetID]->Fill(E0, charge2[k] );
			hresp_stix[sumDetID]->Fill(E0, getScienceBin(edep[k]) );
			hcoll_stix[sumDetID]->Fill(E0, getScienceBin(collected[k]) );
			hreal_stix[sumDetID]->Fill(E0, getScienceBin(charge2[k]) );
			//sum response

			if(isFineGrid(det)){
				sumDetID=33;
			}
			else if(det==8){
				sumDetID=35;
			}else if(det==9){
				sumDetID=36;
			}
			else{
				sumDetID=34;//except fine grids, cfl and bkg
			}


			hresp[sumDetID]->Fill(E0, edep[k] );
			hcoll[sumDetID]->Fill(E0, collected[k] );
			hreal[sumDetID]->Fill(E0, charge2[k] );
			hresp_stix[sumDetID]->Fill(E0, getScienceBin(edep[k]) );
			hcoll_stix[sumDetID]->Fill(E0, getScienceBin(collected[k]) );
			hreal_stix[sumDetID]->Fill(E0, getScienceBin(charge2[k]) );
			if(sumDetID==34 && k%12<8 ){
				//big pixels except cfl and bkg
				sumDetID=37;
				hresp[sumDetID]->Fill(E0, edep[k] );
				hcoll[sumDetID]->Fill(E0, collected[k] );
				hreal[sumDetID]->Fill(E0, charge2[k] );
				hresp_stix[sumDetID]->Fill(E0, getScienceBin(edep[k]) );
				hcoll_stix[sumDetID]->Fill(E0, getScienceBin(collected[k]) );
				hreal_stix[sumDetID]->Fill(E0, getScienceBin(charge2[k]) );

			}


		}
	}
	fo.cd();

	fo.mkdir("singleDetector");
	fo.cd("singleDetector");
	cout<<"Writing histograms"<<endl;
	for(k=0;k<N;k++) {
		hresp[k]->Scale(weight);
		hcoll[k]->Scale(weight);
		hreal[k]->Scale(weight);
		hresp_stix[k]->Scale(weight);
		hcoll_stix[k]->Scale(weight);
		hreal_stix[k]->Scale(weight);

		normMatrix(hresp[k], hsource);
		normMatrix(hcoll[k], hsource);
		normMatrix(hreal[k], hsource);

		normMatrix(hresp_stix[k], hsource);
		normMatrix(hcoll_stix[k], hsource);
		normMatrix(hreal_stix[k], hsource);
		//normalization
	}
	for(k=0;k<32;k++) {

		hresp[k]->Write();
		hcoll[k]->Write();
		hreal[k]->Write();
		hresp_stix[k]->Write();
		hcoll_stix[k]->Write();
		hreal_stix[k]->Write();
	}

	fo.cd();
	fo.mkdir("detectorSum");
	fo.cd("detectorSum");
	for(k=32;k<N;k++)
	{
		hresp[k]->Write();
		hcoll[k]->Write();
		hreal[k]->Write();
		hresp_stix[k]->Write();
		hcoll_stix[k]->Write();
		hreal_stix[k]->Write();
	}


	fo.cd();
	fo.mkdir("pixelSum");
	fo.cd("pixelSum");
	for(k=0;k<12;k++) {
		h_pix_real[k]->Scale(weight);
		h_pix_edep[k]->Scale(weight);
		h_pix_real[k]->Write();
		h_pix_edep[k]->Write();
	}
	fo.Close();
	cout<<"Output file:"<<fout<<endl;

}

void Help(){
	cout<<"make_matrix  -i input.root -o output.root  -b <energy step 0.1 >  -m <max energy, 250> -r  <radius of the source plane> -x <exclude double-hits if 1> "<<endl; 
}
int main(int argc, char *argv[])
{
	TString outputFilename = "defaultOuput.root";
	TString inputFilename;

	if (argc == 1)  Help();
	int s = 0;
	TString sel;

	double eStep=0.1;
	double eMax=150;
	Long64_t entries=0;
	double radius=9;
	bool excludeDoubleHits=false;

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
		else if (sel == "-r") {
			radius = atof(argv[++s]);
		}
		//		else if (sel == "-x") {
		//			excludeDoubleHits = atoi(argv[++s]);
		//		}

	}
	if(inputFilename=="" ){
		Help();
		return 0;
	}
	//	if(excludeDoubleHits){
	outputFilename="resp_"+inputFilename;
	outputFilename="single_hit_"+outputFilename;
	cout<<"Creating response matrix for hits excluding hits..."<<endl;
	excludeDoubleHits=false;
	makeMatrix(inputFilename, outputFilename, eStep, eMax, entries, radius, excludeDoubleHits);
	outputFilename="all_hits_"+outputFilename;
	cout<<"Creating response matrix for all hits..."<<endl;
	excludeDoubleHits=true;
	makeMatrix(inputFilename, outputFilename, eStep, eMax, entries, radius, excludeDoubleHits);

	return 0;
} 
