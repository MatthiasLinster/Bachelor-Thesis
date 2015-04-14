#ifndef DECAY_H
#define DECAY_H

#include <list>

#include "PhaseSpace.h"
#include "MonteCarlo.h"
#include "DecayMonteCarloFunction.h"
#include "FunctionTypes.h"

/**
 * Hilfsstruktur für die Speicherung der Zerfallskanäle
 */
struct DecayChannel {
	matrix_element_ptr pMatrixElement;
	void* pParams;
};

/**
 * Basis-Klasse für Teilchenzerfälle. Berechnet die Zerfallsbreite im Ruhesystem des zerfallenden Teilchens.
 * @author Matthias Linster
 * @date 16.09.2014
 */
class Decay {

	public:
		/**
		 * Konstruktor, erzeugt den Zerfall eines Teilchens mit Masse mass
		 * @param mass Masse des zerfallenden Teilchens in dessen Ruhesystem
		 */
		Decay(const double mass);
		/**
		 * Destruktor
		 */
		virtual ~Decay();
		
		/**
		 * Fügt dem Zerfall einen Zerfallskanal hinzu.
		 * @param pMatrixElement Matrix-Element des Zerfalls
		 * @param pParams Parameter des Matrix-Elements
		 */
		virtual void addDecayChannel(matrix_element_ptr pMatrixElement,void* pParams);
		
		/**
		 * Gibt die Zerfallsbreite für diesen Prozess zurück
		 * Falls diese noch nicht berechnet wurde, wird die Berechnung durchgeführt. Wurde die Berechnung bereits mit einer höheren Schrittanzahl durchgeführt, so wird dieser in der Regel genauere Wert zurückgegeben.
		 * @param decayWidth Zeiger auf die Variable, in der die Zerfallsbreite gespeichert werden soll
		 * @param statError Zeiger auf die Variable, in der der statistische Fehler des Monte-Carlo-Laufs gespeichert werden soll
		 * @param numMonteCarloSteps Anzahl Schritte der Monte-Carlo-Integration, Standardwert 1.000.000
		 */
		virtual void getDecayWidth(double* decayWidth,double* statError,const unsigned long numMonteCarloSteps = 1000000);
		
		/**
		 * Liefert das Volumen des Integrationsbereiches inkl. Phasenraumvolumen
		 * @return Volumen des kompletten Integrationsbereiches
		 */
		virtual double getIntegrationVolume() const;
		
		/**
		 * Liefert den Funktionswert der differentiellen Zerfallsbreite an der übergebenen Stelle
		 * @param pSamplingPoint Zeiger auf eine Variable, die den Punkt, an dem die Funktion ausgewertet werden soll, enthält
		 * @return Funktionswert an der angegebenen Stelle
		 */
		virtual double getFunctionValue(double* pSamplingPoint);
	
	protected:
		/**
		 * Berechnet die Zerfallsbreite.
		 * @param decayWidth Zeiger auf die Ergebnisvariable
		 * @param statError Zeiger auf die Variable für den statistischen Fehler
		 * @param numMonteCarloSteps Anzahl Schritte des Monte Carlo
		 */
		virtual void calcDecayWidth(double* decayWidth,double* statError,const unsigned long numMonteCarloSteps);
	
		// Phasenraum des Zerfalls
		PhaseSpace* m_pPhaseSpace;
		// Masse des zerfallenden Teilchens in dessen Ruhesystem
		double m_mass;
		
		// Zerfallskanäle
		std::list<DecayChannel*> m_decayChannels;
		
		// Monte-Carlo-Ergebnisse
		unsigned long m_numLastMonteCarloSteps = 0;
		double m_decayWidth;
		double m_statError;

};

#endif
