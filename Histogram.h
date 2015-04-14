#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "TH1D.h"

#include "MonteCarloFunction.h"

/**
 * Struktur für Histogramme
 * @author Matthias Linster
 * @date 03.09.2014
 */
struct Histogram {
	// eindeutiger Name des Histogramms
	const char* name;
	// Titel des Histogramms
	const char* title;
	// Anzahl Bin des Histogramms
	unsigned int numBins;
	// Linke Grenze des Histogramm-Bereichs
	double xl;
	// Rechte Grenze des Histogramm-Bereichs
	double xu;
	// Breite der Bins
	double binWidth;
	// Funktion, die das zu füllende Bin bestimmt
	MonteCarloFunction* pTransformFunction;
	
	// Integral des Histogramms zur Normierung
	double histIntegral;
	
	// Zeiger auf die Arrays für den Wert der Bins und der Anzahl an Einträgen
	double* binIntegral = NULL;
	unsigned long* numBinEntries = NULL;
	
	/**
	 * Konstruktor, erstellt ein Histogramm-Objekt
	 * @param name Eindeutiger Name des Histogramms
	 * @param title Titel des Histogramms
	 * @param numBins Anzahl Bins des Histogramms
	 * @param xl Linke Grenze des Histogramm-Bereichs
	 * @param xu Rechte Grenze des Histogramm-Bereichs
	 * @param pTransformFunction Funktion, die das zu füllende Bin bestimmt
	 */
	Histogram(const char* name,const char* title,const unsigned int numBins,double xl,double xu,MonteCarloFunction* pTransformFunction) :
			name(name),
			title(title),
			numBins(numBins),
			xl(xl),
			xu(xu),
			pTransformFunction(pTransformFunction) {
		// Bin-Breite ausrechnen
		binWidth = (xu-xl)/numBins;
		
		// Variable für das Histogramm-Integral auf 0 setzen
		histIntegral = 0.;
		
		// Arrays für die Bins erzeugen und auf 0 setzen
		binIntegral = new double[numBins];
		numBinEntries = new unsigned long[numBins];
		
		for (unsigned int i=0;i<numBins;i++) {
			binIntegral[i] = 0.;
			numBinEntries[i] = 0;
		}
	}
	
	/**
	 * Destruktor
	 */
	~Histogram() {
		// Dynamisch angelegte Arrays löschen
		if (binIntegral!=NULL) {
			delete[] binIntegral;
			binIntegral = NULL;
		}
		
		if (numBinEntries!=NULL) {
			delete[] numBinEntries;
			numBinEntries = NULL;
		}
	}
	
	/**
	 * Gibt ein ROOT TH1D-Objekt zurück, das aus den Histogramm-Daten erzeugt wurde
	 * @return TH1D-Objekt mit dem Histogramm
	 */
	TH1D* getRootTH1D() const {
		// Histogramm erzeugen
		TH1D* rootHist = new TH1D(name,title,numBins,xl,xu);
		
		// Bins füllen
		for (unsigned int i=0;i<numBins;i++) {
			rootHist->SetBinContent(i+1,binIntegral[i]);
		}
		
		return rootHist;
	}
};

#endif
