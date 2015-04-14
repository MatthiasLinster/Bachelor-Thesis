#ifndef SCATTERING_H
#define SCATTERING_H

#include <string.h>
#include <list>

#include "PhaseSpace.h"
#include "MonteCarlo.h"
#include "MonteCarloFunction.h"
#include "ScatteringMonteCarloFunction.h"

/**
 * Basis-Klasse für Streuprozesse
 * @author Matthias Linster
 * @date 03.09.2014
 */
class Scattering {
	
	public:
		/**
		 * Default-Konstruktor und Destruktor
		 */
		Scattering();
		virtual ~Scattering();
		
		/**
		 * Fügt dem Streuprozess eine zu histogrammierende Variable hinzu
		 * @param name Eindeutiger Name des Histogramms
		 * @param title Titel des Histogramms
		 * @param numBins Anzahl Bin des Histogramms
		 * @param xl Untere Grenze des Histogramm-Bereichs
		 * @param xu Obere Grenze des Histogramm-Bereichs
		 * @param pFunction Formel der zu histogrammierenden Variable
		 */
		virtual void addHistogram(const char* name,const char* title,const unsigned int numBins,const double xl,const double xu,MonteCarloFunction* pFunction);
		
		/**
		 * Liefert den Zeiger auf ein erzeugtes Histogramm
		 * @param name Eindeutiger Name des Histogramms, unter dem es angelegt wurde
		 * @return Zeiger auf das Histogramm
		 */
		virtual Histogram* getHistogram(const char* name) const;
		
		/**
		 * Liefert den Phasenraum des aktuellen Streuprozesses
		 * @return Zeiger auf den Phasenraum des Streuprozesses
		 */
		virtual PhaseSpace* getPhaseSpace() const;
		
		/**
		 * Gibt den Wirkungsquerschnitt für diesen Streuprozess zurück.
		 * Falls dieser noch nicht berechnet wurde, wird die Berechnung durchgeführt. Wurde die Berechnung bereits mit einer höheren Schrittanzahl durchgeführt, so wird dieser in der Regel genauere Wert zurückgegeben.
		 * @param crossSection Zeiger auf die Variable, in der der Wirkungsquerschnitt gespeichert werden soll
		 * @param statError Zeiger auf die Variable, in der der statistische Fehler des Monte-Carlo-Laufs gespeichert werden soll
		 * @param numMonteCarloSteps Anzahl Schritte der Monte-Carlo-Integration, Standardwert 1.000.000
		 */
		virtual void getCrossSection(double* crossSection,double* statError,const unsigned long numMonteCarloSteps = 1000000);
		
		/**
		 * Liefert das Volumen des Integrationsbereiches inkl. Phasenraumvolumen
		 * @return Volumen des kompletten Integrationsbereiches
		 */
		virtual double getIntegrationVolume() const = 0;
		
		/**
		 * Liefert den Funktionswert des differentiellen Wirkungsquerschnitts an der übergebenen Stelle
		 * @param pSamplingPoint Zeiger auf eine Variable, die den Punkt, an dem die Funktion ausgewertet werden soll, enthält
		 * @return Funktionswert an der angegebenen Stelle
		 */
		virtual double getFunctionValue(double* pSamplingPoint) = 0;
		
	protected:
		// Phasenraum des Streuprozesses
		PhaseSpace* m_pPhaseSpace;
		// Anzahl Dimensionen der Integration
		size_t m_numDimensions;
		
		// Grenzen für die Integration
		double* m_pXl = NULL;
		double* m_pXu = NULL;
		
		// Liste mit den zu erstellenden Diagrammen
		std::list<Histogram*> m_histograms;
		
		// Speichervariablen für den Wirkungsquerschnitt, den Fehler sowie die Anzahl der Monte-Carlo-Schritte
		double m_crossSection = 0.;
		double m_crossSectionStatError = 0.;
		unsigned long m_numLastMonteCarloSteps = 0;
		
		/**
		 * Initialisiert das Streu-Objekt
		 * @param numDimensions Anzahl Dimensionen des Integrationsbereiches (inkl. Phasenraum)
		 * @param pXl Zeiger auf ein Array mit den unteren Grenzen der Integration
		 * @param pXu Zeiger auf ein Array mit den oberen Grenzen der Integration
		 */
		virtual void initScattering(const size_t numDimensions,double* pXl,double* pXu);
		
		/**
		 * Berechnet den Wirkungsquerschnitt.
		 * @param crossSection Zeiger auf die Ergebnisvariable
		 * @param statError Zeiger auf die Variable für den statistischen Fehler
		 * @param numMonteCarloSteps Anzahl Schritte des Monte Carlo
		 */
		virtual void calcCrossSection(double* crossSection,double* statError,const unsigned long numMonteCarloSteps);
		
};

#endif
