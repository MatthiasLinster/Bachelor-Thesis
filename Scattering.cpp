#include "Scattering.h"

Scattering::Scattering() {
}

Scattering::~Scattering() {
	// Speicher für dynamisch erzeugte Variablen freigeben
	if (m_pXl!=NULL)  {
		delete[] m_pXl;
		m_pXl = NULL;
	}
	
	if (m_pXu!=NULL) {
		delete[] m_pXu;
		m_pXu = NULL;
	}
	
	for (Histogram* histogram : m_histograms) {
		delete histogram;
	}
}

void Scattering::addHistogram(const char* name,const char* title,const unsigned int numBins,const double xl,const double xu,MonteCarloFunction* pTransformFunction) {
	// Histogramm mit passenden Arrays erstellen
	Histogram* histogram = new Histogram(name, title, numBins, xl, xu, pTransformFunction);
	
	// Der Liste hinzufügen
	m_histograms.push_back(histogram);
}

Histogram* Scattering::getHistogram(const char* name) const {
	// Histogramm suchen und zurückgeben
	for (Histogram* histogram : m_histograms) {
		if (strcmp(name,histogram->name)==0) {
			return histogram;
		}
	}
	
	// Falls nicht gefunden, NULL zurückgeben
	return NULL;
}

PhaseSpace* Scattering::getPhaseSpace() const {
	return m_pPhaseSpace;
}

void Scattering::getCrossSection(double* crossSection,double* statError,const unsigned long numMonteCarloSteps) {
	// Falls genaueres Ergebnis bekannt, dieses zurückgeben, sonst ausrechnen lassen
	if (m_numLastMonteCarloSteps>=numMonteCarloSteps) {
		*crossSection = m_crossSection;
		*statError = m_crossSectionStatError;
	}
	else 
		calcCrossSection(crossSection,statError,numMonteCarloSteps);
}

void Scattering::initScattering(const size_t numDimensions,double* pXl,double* pXu) {
	// Speichern, dabei übergebene Arrays kopieren
	m_numDimensions = numDimensions;
	m_pXl = new double[numDimensions];
	m_pXu = new double[numDimensions];
	
	for (size_t i=0;i<numDimensions;i++) {
		m_pXl[i] = pXl[i];
		m_pXu[i] = pXu[i];
	}
}

void Scattering::calcCrossSection(double* crossSection,double* statError,const unsigned long numMonteCarloSteps) {
	// Passendes Monte Carlo anlegen
	MonteCarlo monteCarlo(m_numDimensions,numMonteCarloSteps);
	
	// Funktion für das Monte-Carlo erzeugen und setzen, Funktionswerte werden durch den differentiellen Wirkungsquerschnitt, der vom Streu-Objekt behandelt wird, gegeben
	ScatteringMonteCarloFunction monteCarloFunction(this,m_numDimensions,m_pXl,m_pXu);
	monteCarlo.setIntegrationFunction(&monteCarloFunction);
	
	// Histogramme hinzufügen
	for (Histogram* histogram : m_histograms) {
		monteCarlo.addHistogram(histogram);
	}
	
	// Integrieren
	monteCarlo.integrate(&m_crossSection,&m_crossSectionStatError);
	m_numLastMonteCarloSteps = numMonteCarloSteps;
	
	// Ergebnisvariablen updaten
	*crossSection = m_crossSection;
	*statError = m_crossSectionStatError;
}
