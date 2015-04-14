#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>

#include "LHAPDF/LHAPDF.h"
#include "TH1D.h"
#include "TFile.h"

#include "MonteCarlo.h"
#include "MonteCarloFunction.h"
#include "PhaseSpaceMonteCarloFunction.h"
#include "TwoTwoOneFinalMassiveHadronScattering.h"

// Standard Anzahl an Schritten, wenn dies vom Benutzer nicht festgelegt wird
#define DEFAULT_STEP_COUNT 1000000

// Parameter
const double S_PP = 14000; /* Proton-Proton-Schwerpunktsenergie sqrt(s) in GeV */
const double M_H = 125.1; /* Higgs-Masse in GeV */
const double COS2_THETA_W = 1.0 - 0.23155; /* cos^2(theta_W) = 1 - sin^2(theta_W) */

// Funktion für den differentiellen Wirkungsquerschnitt (ausgerechnet bis auf nicht analytisch bestimmbare Parton-Parton-Prozesse)
inline double diffCrossSection(double* x,size_t numDim,void* params) {
	// Theta(x1*x2*s-m_h^2)
	if (x[0]*x[1]*S_PP*S_PP<M_H*M_H)
		return 0.0;
	
	// LHAPDF aus Parameter abfragen
	LHAPDF::PDF* pdf = (LHAPDF::PDF*)params;
	
	// Ergebnisvariable initialisieren
	double result = 0.0;
	
	// Alle möglichen Parton-Parton-Reaktionen durchgehen und zum Ergebnis summieren
	for (int i=1;i<7;i++) {
		result += pdf->xfxQ(i,x[0],M_H)/x[0]*pdf->xfxQ(-i,x[1],M_H)/x[1] 										/* f_q(x1)*f_qbar(x2) */
				  * COS2_THETA_W/144./M_PI													/* cos^2(theta_W)/(144 PI) */
				  * (x[0]*x[1]*S_PP*S_PP-M_H*M_H)*(x[0]*x[1]*S_PP*S_PP-M_H*M_H)*(x[0]*x[1]*S_PP*S_PP-M_H*M_H) 	/* (x1*x2*s-m_h^2)^3 */
				  / (x[0]*x[1]*S_PP*S_PP)/(x[0]*x[1]*S_PP*S_PP); 												/* (x1*x2*s)^2 */
				  
		result += pdf->xfxQ(-i,x[0],M_H)/x[0]*pdf->xfxQ(i,x[1],M_H)/x[1] 										/* f_qbar(x1)*f_q(x2) */
				  * COS2_THETA_W/144./M_PI													/* cos^2(theta_W)/(144 PI) */
				  * (x[0]*x[1]*S_PP*S_PP-M_H*M_H)*(x[0]*x[1]*S_PP*S_PP-M_H*M_H)*(x[0]*x[1]*S_PP*S_PP-M_H*M_H) 	/* (x1*x2*s-m_h^2)^3 */
				  / (x[0]*x[1]*S_PP*S_PP)/(x[0]*x[1]*S_PP*S_PP); 												/* (x1*x2*s)^2 */
	}
	
	// Ergebnis zurückgeben
	return result;
}

// Matrix-Element des Parton-Parton-Streuprozesses
inline double matrixElement(FourVector* momentum,void* params) {
	return 4./3. * 2. * COS2_THETA_W * (momentum[0]*momentum[3]) * (momentum[1]*momentum[3]);
}

// Liefert den Transversalimpuls eines Ereignisses (hier: des massiven Teilchens)
inline double transverseMomentum(FourVector* momentum,void* params) {
	return sqrt(momentum[0].x[1]*momentum[0].x[1]+momentum[0].x[2]*momentum[0].x[2]);
}

// Liefert die Rapidität eines Ereignisses (hier: des massiven Teilchens)
inline double rapidity(FourVector* momentum,void* params) {
	return 1./2.*log((momentum[0].x[0]+momentum[0].x[3])/(momentum[0].x[0]-momentum[0].x[3]));
}

int main(int argc,char** argv) {
	// Anzahl an Schritten
	unsigned long numSteps = DEFAULT_STEP_COUNT;

	// Überprüfen, ob die Anzahl an Schritten als Argument übergeben wurde, wenn ja, Variable updaten
	if (argc!=2) 
		printf("Keine Schrittanzahl übergeben, verwende Standardwert von 1.000.000 ...\n");
	else {
		numSteps = strtoul(argv[1],NULL,0);
		printf("Anzahl an Integrationsschritten: %lu\n",numSteps);
	}
	
	// Integrationsbereich
	double xl[] = { 0.,0. };
	double xu[] = { 1.,1. };
	
	// LHAPDF initialisieren
	LHAPDF::PDF* pdf = LHAPDF::mkPDF("MSTW2008lo90cl",0);
	
	// Monte-Carlo anlegen
	MonteCarlo mc(2,numSteps);
	
	// Funktion setzeeen
	MonteCarloFunction mcFunction(diffCrossSection,2,xl,xu,pdf);
	mc.setIntegrationFunction(&mcFunction);
	
	// Integrieren
	double result,error;
	//mc.integrate(&result,&error);
	
	// Ergebnis + Fehler ausgeben
	printf("cross section: %f +/- %f\n",result,error);
	
	// Alternativer Weg über Phasenraum und Matrixelement, dazu Streuprozess-Objekt erstellen
	TwoTwoOneFinalMassiveHadronScattering scattering(M_H,S_PP,M_H,pdf);
	
	// Dem Prozess beitragende einzelne Parton-Prozess registrieren (jeweils Quark-Antiquark)
	for (int i=1;i<7;i++) {
		scattering.addPartonProcess(i,-i,matrixElement);
		scattering.addPartonProcess(-i,i,matrixElement);
	}
	
	// Histogramme für die Verteilungen nach Rapidität und Transversalimpuls anlegen
	PhaseSpaceMonteCarloFunction transverseTransformFunction(scattering.getPhaseSpace(),transverseMomentum,NULL);
	scattering.addHistogram("transverseMomentum","Transverse Momentum",50,0.,1000.,&transverseTransformFunction);
	
	PhaseSpaceMonteCarloFunction rapidityTransformFunction(scattering.getPhaseSpace(),rapidity,NULL);
	scattering.addHistogram("rapidity","rapidity",50,-5.,5.,&rapidityTransformFunction);
	
	// Integrieren und Ergebnis ausgeben
	scattering.getCrossSection(&result,&error,numSteps);
	printf("cross section: %f +/- %f\n",result,error);
	
	// Histogramme zu Root-Histogrammen konvertieren
	Histogram* transverseMomentumHistogram = scattering.getHistogram("transverseMomentum");
	Histogram* rapidityHistogram = scattering.getHistogram("rapidity");
	TH1D* transverseRootHist = transverseMomentumHistogram->getRootTH1D();
	TH1D* rapidityRootHist = rapidityHistogram->getRootTH1D();
	
	// Histogramm in Root-File abspeichern
	TFile* file = new TFile("hist.root","RECREATE");
	transverseRootHist->Write();
	rapidityRootHist->Write();
	file->Close();
	
	// Speicher freigeben
	delete transverseRootHist;
	delete rapidityRootHist;
	delete file;
	
	return 0;
}
