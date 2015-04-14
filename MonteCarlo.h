#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <gsl/gsl_qrng.h>
#include <math.h>
#include <string.h>
#include <list>

#include "MonteCarloFunction.h"
#include "FourVector.h"
#include "Histogram.h"

/**
 * Klasse zur Monte-Carlo-Integration
 * @author Matthias Linster
 * @date 03.09.2014
 */
class MonteCarlo {
	
	public:
		/** 
		 * Konstruktor, erstellt ein Monte Carlo mit der angegebenen Dimension und Anzahl Integrationsschritten
		 * @param numDimensions Anzahl Dimensionen
		 * @param numIntegrationSteps Anzahl an Integrationsschritten
		 */
		MonteCarlo(const size_t numDimensions,const unsigned long numIntegrationSteps);
		/**
		 * Destruktor
		 */
		virtual ~MonteCarlo();
		
		/**
		 * Setzt die zu integrierende Funktion
		 * @param pFunction zu integrierende Funktion
		 */
		virtual void setIntegrationFunction(MonteCarloFunction* pFunction);
		
		/**
		 * Fügt dem Monte Carlo ein anzulegendes Histogramm hinzu
		 * @param histogram zu erzeugendes Histogramm
		 */
		virtual void addHistogram(Histogram* histogram);
		
		/**
		 * Liefert ein erzeugtes Histogramm zurück
		 * @param name Name des gewünschten Histogramms
		 * @return Erzeugtes Histogramm
		 */
		virtual Histogram* getHistogram(const char* name) const;
		
		/**
		 * Berechnet das Integral mit der bereits übergebenen Anzahl an Schritten
		 * @param result Zeiger auf die Ergebnisvariable
		 * @param error Zeiger auf die Variable für den statistischen Fehler der Monte-Carlo-Integration
		 */
		virtual void integrate(double* result,double* error);
	
	protected:
		// Anzahl an Dimensionen
		size_t m_numDimensions;
		// Anzahl an Integrationsschritten
		unsigned long m_numIntegrationSteps;
		
		// zu integrierende Funktion
		MonteCarloFunction* m_pFunction;
		// zu erzeugende Verteilungen
		std::list<Histogram*> m_histograms;
	
};

#endif
