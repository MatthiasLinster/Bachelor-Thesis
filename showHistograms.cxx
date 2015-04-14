void showHistograms() {
	TFile* file = new TFile("hist.root");
	TH1D* transverseHist = (TH1D*)file->GetObjectChecked("transverseMomentum","TH1D");
	TH1D* rapidityHist = (TH1D*)file->GetObjectChecked("rapidity","TH1D");
	
	TCanvas* canvas = new TCanvas("canvas","Distributions",1600,640);
	canvas->Divide(2,1);
	
	canvas->cd(1);
	transverseHist->DrawClone();
	canvas->cd(2);
	rapidityHist->DrawClone();
}
