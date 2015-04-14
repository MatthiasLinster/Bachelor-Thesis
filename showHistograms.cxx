void showHistograms() {
	// Lade Histogramme
	TFile* file = new TFile("hist.root");
	TH1D* transverseHist = (TH1D*)file->GetObjectChecked("transverseMomentum","TH1D");
	TH1D* rapidityHist = (TH1D*)file->GetObjectChecked("rapidity","TH1D");
	
	// Skalierung mit f^2/Lambda^4 (2.9e-10) und Umrechnung in pb (3.9e8)
	transverseHist->Scale(2.9e-10*3.9e8);
	rapidityHist->Scale(2.9e-10*3.9e8);
	
	// Überschrift + Achsenbeschriftung
	transverseHist->SetTitle("p_{T}-Verteilung bei #sqrt{s} = xx TeV;Transversalimpuls p_{T} / GeV;diff. WQ #frac{d#sigma}{dp_{T}} / #frac{pb}{GeV}");
	rapidityHist->SetTitle("Rapiditaetsverteilung bei #sqrt{s} = xx TeV;Rapiditaet y / 1;diff. WQ #frac{d#sigma}{dy} / pb");
	
	// Stats-Panel deaktivieren
	transverseHist->SetStats(0);
	rapidityHist->SetStats(0);
	
	// Canvas zum Zeichnen erstellen
	TCanvas* canvas = new TCanvas("canvas","Distributions",1600,640);
	canvas->Divide(2,1);
	
	// Histogramme zeichnen
	canvas->cd(1);
	transverseHist->DrawClone();
	canvas->cd(2);
	rapidityHist->DrawClone();
	
	// Check des Integrals
	cout << transverseHist->Integral() << "," << rapidityHist->Integral() << endl;
}
