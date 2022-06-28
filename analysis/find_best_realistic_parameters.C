{
	ifstream infile("exp_data.csv"); 
	double energy[250]={0};
	double counts[250]={0};
	int n=0;


	 if(!infile.good())
	 { 
		  cout<<"can not open file..."<<endl; 
		   return 0;
	 } 
	 int maxExpValue;
	 while(!infile.eof()){
		 infile>>energy[n]>>counts[n];
		 n++;
	 }
	 TCanvas c;
	 c.Divide(2);
	 TGraph gr(n,energy, counts);
	 gr.Draw("ALP");


	 infile.close(); 
	 TFile f("ba133.root");
	 TTree *evt=(TTree*)f.Get("events");
	 evt->Draw("edep>>hedep(400,0,100)","edep>0");
	 TH1F *h=(TH1F*)gDirectory->Get("hedep");
	double maxValue=h->GetBinContent(h->GetMaximumBin());
	h->Scale(1/maxValue);
	 h->Draw("same");
}
