#ifndef DECAY_MONTE_CARLO_FUNCTION_H
#define DECAY_MONTE_CARLO_FUNCTION_H

#include "MonteCarloFunction.h"
#include "Decay.h"

// Forward-Deklaration für die Scattering-Klasse
class Decay;

/**
 * Klasse, die dem Monte Carlo Zugriff auf ein Zerfallsobjekt gibt, das die zu integrierende Funktion beinhaltet (Zerfallsbreite)
 * @author Matthias Linster
 * @date 14.09.2014
 */
class DecayMonteCarloFunction : public MonteCarloFunction {
	
	public:
		/**
		 * Konstruktor, erstellt das Funktionsobjekt
		 * @param pDecay Zerfalls-Objekt, das betrachtet wird
		 * @param numDimensions Anzahl Dimensionen des Integrationsbereiches
		 * @param pXl Array mit den unteren Grenzen des Integrationsbereichs
		 * @param pXu Array mit den oberen Grenzen des Integrationsbereichs
		 * @param pParams Zeiger auf die Parameter der Funktion
		 */
		DecayMonteCarloFunction(Decay* pDecay,size_t numDimensions,double* pXl,double* pXu,void* pParams = NULL);
		
		/**
		 * Destruktor
		 */
		virtual ~DecayMonteCarloFunction();
	
		/**
		 * Liefert das Volumen des Integrationsbereichs
		 * @return Volumen des Integrationsbereichs
		 */
		virtual double getIntegrationVolume() const;
		
		/**
		 * Liefert den Funktionswert an einer definierten Stelle
		 * @param pSamplingPoint Punkt, an dem die Funktion ausgewertet werden soll
		 * @return Funktionswert an der entsprechenden Stelle
		 */
		virtual double getFunctionValue(double* pSamplingPoint);
	
	protected:
		// Zeiger auf das Streuobjekt
		Decay* m_pDecay;
	
};

#endif
