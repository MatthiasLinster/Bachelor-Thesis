#ifndef PHASE_SPACE_MONTE_CARLO_FUNCTION_H
#define PHASE_SPACE_MONTE_CARLO_FUNCTION_H

#include "MonteCarloFunction.h"
#include "FourVector.h"
#include "PhaseSpace.h"
#include "FunctionTypes.h"

/**
 * Funktion für die Monte-Carlo-Integration, die vom Phasenraum abhängt
 * @author Matthias Linster
 * @date 03.09.2014
 */
class PhaseSpaceMonteCarloFunction : public MonteCarloFunction {
	
	public:
		/**
		 * Kontruktor, erzeugt das Funktionsobjekt
		 * @param pPhaseSpace Zeiger auf den Phasenraum des Prozesses
		 * @param pFunction Funktion
		 * @param pParams Zeiger auf die Parameter der Funktion
		 */
		PhaseSpaceMonteCarloFunction(PhaseSpace* pPhaseSpace,phase_space_function_ptr pFunction,void* pParams);
		
		/**
		 * Destruktor
		 */
		virtual ~PhaseSpaceMonteCarloFunction();
		
		/**
		 * Liefert den Funktionswert an einer definierten Stelle.
		 * @param pSamplingPoint Stelle, an der die Funktion ausgewertet werden soll
		 * @return Funktionswert an dieser Stelle
		 */
		virtual double getFunctionValue(double* pSamplingPoint);
	
	protected:
		// Funktion
		phase_space_function_ptr m_pPhaseSpaceFunction;
		// Phasenraum
		PhaseSpace* m_pPhaseSpace;
	
};

#endif
