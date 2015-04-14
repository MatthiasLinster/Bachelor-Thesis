#include "MonteCarloFunction.h"

MonteCarloFunction::MonteCarloFunction(int_function_ptr pFunction,size_t numDimensions,double* pXl,double* pXu,void* pParams) :
	m_pFunction(pFunction),
	m_numDimensions(numDimensions),
	m_pXl(pXl),
	m_pXu(pXu),
	m_pParams(pParams) {
}

MonteCarloFunction::~MonteCarloFunction() {
}

double MonteCarloFunction::getIntegrationVolume() const {
	// Volumen des Integrationsbereiches berechnen
	double V = 1;
	
	for (size_t i=0;i<m_numDimensions;i++) {
		V *= (m_pXu[i]-m_pXl[i]);
	}
	
	return V;
}

double MonteCarloFunction::getFunctionValue(double* pSamplingPoint) {
	// Zufallspunkt auf Integrationsbereich ausdehnen
	for (size_t i=0;i<m_numDimensions;i++) {
		pSamplingPoint[i] = m_pXl[i] + pSamplingPoint[i] * (m_pXu[i] - m_pXl[i]);
	}
	
	// Funktionswert ausrechnen und zurückgeben
	return m_pFunction(pSamplingPoint,m_numDimensions,m_pParams);
}
