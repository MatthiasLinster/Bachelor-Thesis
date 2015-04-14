#include "DecayMonteCarloFunction.h"

DecayMonteCarloFunction::DecayMonteCarloFunction(Decay* pDecay,size_t numDimensions,double* pXl,double* pXu,void* pParams) : 
		MonteCarloFunction(NULL,numDimensions,pXl,pXu,pParams) {
	// Werte speichern
	m_pDecay = pDecay;
}

DecayMonteCarloFunction::~DecayMonteCarloFunction() {
}

double DecayMonteCarloFunction::getIntegrationVolume() const {
	// An Decay-Objekt weiterleiten
	return m_pDecay->getIntegrationVolume();
}

double DecayMonteCarloFunction::getFunctionValue(double* pX) {
	// An Decay-Objekt weiterleiten
	return m_pDecay->getFunctionValue(pX);
}
