{
	TFile *_file1 = TFile::Open("../roots/tungsten_grid_det30_E_pow_m5_track_killed.root");
	TFile *_file0 = TFile::Open("../roots/tungsten_grid_det30_E_pow_m5.root");
	auto *ev2=(TTree*)_file1->Get("events");
	auto *ev1=(TTree*)_file0->Get("events");
	ev2->Draw("edep>>h2(200,0,100)","edep>4");
	ev1->Draw("edep>>h1(200,0,100)","edep>4");
	TCanvas c1;
	c1.cd();
	h1->Draw("hist+e");
	h2->Draw("hist+e+same");

}
