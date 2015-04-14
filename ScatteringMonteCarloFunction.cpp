#include "ScatteringMonteCarloFunction.h"

ScatteringMonteCarloFunction::ScatteringMonteCarloFunction(Scattering* pScattering,size_t numDimensions,double* pXl,double* pXu,void* pParams) : 
		MonteCarloFunction(NULL,numDimensions,pXl,pXu,pParams) {
	// Werte speichern
	m_pScattering = pScattering;
}

ScatteringMonteCarloFunction::~ScatteringMonteCarloFunction() {
}

double ScatteringMonteCarloFunction::getIntegrationVolume() const {
	// An Scattering-Objekt weiterleiten
	return m_pScattering->getIntegrationVolume();
}

double ScatteringMonteCarloFunction::getFunctionValue(double* pX) {
	// An Scattering-Objekt weiterleiten
	return m_pScattering->getFunctionValue(pX);
}
