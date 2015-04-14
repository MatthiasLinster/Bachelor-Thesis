#include "PhaseSpaceMonteCarloFunction.h"

PhaseSpaceMonteCarloFunction::PhaseSpaceMonteCarloFunction(PhaseSpace* pPhaseSpace,phase_space_function_ptr pFunction,void* pParams) :
		MonteCarloFunction(NULL,0,NULL,NULL,pParams) {
	// Werte speichern
	m_pPhaseSpace = pPhaseSpace;
	m_pPhaseSpaceFunction = pFunction;
}

PhaseSpaceMonteCarloFunction::~PhaseSpaceMonteCarloFunction() {
}

double PhaseSpaceMonteCarloFunction::getFunctionValue(double* pSamplingPoint) {
	// Letzten Phasenraumpunkt abfragen
	FourVector* particleMomentum = m_pPhaseSpace->getLastPhaseSpacePoint();
	
	return m_pPhaseSpaceFunction(particleMomentum,m_pParams);
}
