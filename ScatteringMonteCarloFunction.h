#ifndef SCATTERING_MONTE_CARLO_FUNCTION_H
#define SCATTERING_MONTE_CARLO_FUNCTION_H

#include "MonteCarloFunction.h"
#include "Scattering.h"

// Forward-Deklaration für die Scattering-Klasse
class Scattering;

/**
 * Klasse, die dem Monte Carlo Zugriff auf ein Streuobjekt gibt, das die zu integrierende Funktion beinhaltet (Wirkungsquerschnitt)
 * @author Matthias Linster
 * @date 03.09.2014
 */
class ScatteringMonteCarloFunction : public MonteCarloFunction {
	
	public:
		/**
		 * Konstruktor, erstellt das Funktionsobjekt
		 * @param pScattering Streu-Objekt, das betrachtet wird
		 * @param numDimensions Anzahl Dimensionen des Integrationsbereiches
		 * @param pXl Array mit den unteren Grenzen des Integrationsbereichs
		 * @param pXu Array mit den oberen Grenzen des Integrationsbereichs
		 * @param pParams Zeiger auf die Parameter der Funktion
		 */
		ScatteringMonteCarloFunction(Scattering* pScattering,size_t numDimensions,double* pXl,double* pXu,void* pParams = NULL);
		
		/**
		 * Destruktor
		 */
		virtual ~ScatteringMonteCarloFunction();
	
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
		Scattering* m_pScattering;
	
};

#endif
