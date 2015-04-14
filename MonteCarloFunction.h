#ifndef MONTE_CARLO_FUNCTION_H
#define MONTE_CARLO_FUNCTION_H

#include <stdlib.h>

#include "FunctionTypes.h"

/**
 * Klasse für die Funktion, die per Monte Carlo integriert werden soll
 * @author Matthias Linster
 * @date 03.09.2014
 */
class MonteCarloFunction {
	
	public:
		/**
		 * Konstruktor, erstellt die gewünschte Funktion
		 * @param pFunction zu integrierende Funktion
		 * @param numDimensions Anzahl Dimensionen des Integrationsbereiches
		 * @param pXl Zeiger auf das Array mit den unteren Grenzen des Integrationsbereichs
		 * @param pXu Zeiger auf das Array mit den oberen Grenzen des Integrationsbereichs
		 * @param pParams Zeiger auf die Parameter der Funktion
		 */
		MonteCarloFunction(int_function_ptr pFunction,size_t numDimensions,double* pXl,double* pXu,void* pParams = NULL);
		
		/**
		 * Destruktor
		 */
		virtual ~MonteCarloFunction();
		
		/**
		 * Liefert das Volumen des Integrationsbereichs der Funktion
		 * @return Volumen des Integrationsbereiches der Funktion
		 */
		virtual double getIntegrationVolume() const;
		
		/**
		 * Liefert den Funktionswert der Funktion an einer bestimmten Stelle
		 * @param pSamplingPoint Punkt, an dem die Funktion ausgewertet werden soll
		 * @return Funktionswert an der entsprechenden Stelle
		 */
		virtual double getFunctionValue(double* pSamplingPoint);
		
	protected:
		// Zeiger auf die Funktionsdefinition
		int_function_ptr m_pFunction;
		// Anzahl Dimensionen des Integrationsbereiches
		size_t m_numDimensions;
		// Array mit den unteren Grenzen des Integrationsbereichs
		double* m_pXl;
		// Array mit den oberen Grenzen des Integrationsbereichs
		double* m_pXu;
		// Parameter der Funktion
		void* m_pParams;
	
};

#endif
