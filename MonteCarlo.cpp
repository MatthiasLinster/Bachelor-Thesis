#include "MonteCarlo.h"

MonteCarlo::MonteCarlo(const size_t numDimensions,const unsigned long numIntegrationSteps) :
		m_numDimensions(numDimensions),
		m_numIntegrationSteps(numIntegrationSteps) {
}

MonteCarlo::~MonteCarlo() {
}

void MonteCarlo::setIntegrationFunction(MonteCarloFunction* pFunction) {
	// Werte speichern
	m_pFunction = pFunction;
}

void MonteCarlo::addHistogram(Histogram* histogram) {
	// Histogramm der Liste hinzufügen
	m_histograms.push_back(histogram);
}

Histogram* MonteCarlo::getHistogram(const char* name) const {
	// Histogramm suchen und zurückgeben
	for (Histogram* histogram : m_histograms) {
		if (strcmp(name,histogram->name)==0) {
			return histogram;
		}
	}
	
	// Im Fehlerfall NULL-Zeiger zurückgeben
	return NULL;
}

void MonteCarlo::integrate(double* result,double* error) {
	// Hilfvariablen (Ergebnissummen, Array für den erzeugten Punkt, Funktionswert, x-Wert des Bins sowie Bin-Nummer)
	double resultSum = 0;
	double errorSum = 0;
	double samplingPoint[m_numDimensions];
	double fx;
	double binX;
	unsigned int binNumber;
	
	// Zufallszahlengenerator initialisieren
	gsl_qrng* q = gsl_qrng_alloc(gsl_qrng_sobol,m_numDimensions);
	
	// Monte Carlo durchführen
	for (unsigned long i=0;i<m_numIntegrationSteps;i++) {
		// Punkt erzeugen
		gsl_qrng_get(q,samplingPoint);
		
		// Summen aktualisieren
		fx = m_pFunction->getFunctionValue(samplingPoint);
		resultSum += fx;	
		errorSum += fx*fx;
		
		// Histogramme aktualisieren
		for (Histogram* histogram : m_histograms) {
			// Bin bestimmen
			binX = histogram->pTransformFunction->getFunctionValue(samplingPoint);
			
			// Sicherstellen, dass der Punkt im Histogramm-Bereich liegt
			if (binX>=histogram->xl && binX<=histogram->xu) {
				// Bin-Nr. bestimmen
				binNumber = (unsigned int)((binX-histogram->xl)/histogram->binWidth);
				
				// Integral des Histogramm-Bereichs updaten
				histogram->histIntegral += fx;
				
				// Bin-Integral updaten
				histogram->binIntegral[binNumber] += fx;
				histogram->numBinEntries[binNumber]++;
			}
		}
	}
	
	// Volumen des Integrationsbereiches
	double V = m_pFunction->getIntegrationVolume();
	
	// Ergebnis mitteln
	*result = V/((double)m_numIntegrationSteps)*resultSum;
	
	// Fehler
	double var = (V*V/m_numIntegrationSteps*errorSum-(*result)*(*result))/(m_numIntegrationSteps-1);
	*error = sqrt(var);
	
	// Histogramme berechnen
	for (Histogram* histogram : m_histograms) {
		// Integral in allen Bins ausrechnen
		for (unsigned int i=0;i<histogram->numBins;i++) {	
			// Monte-Carlo-Integral bestimmen
			histogram->binIntegral[i] = V/m_numIntegrationSteps*histogram->binIntegral[i];
		}
	}
		
	// Zufallsgenerator freigeben
	gsl_qrng_free(q);
}
