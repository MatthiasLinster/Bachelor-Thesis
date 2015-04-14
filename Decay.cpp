#include "Decay.h"

Decay::Decay(const double mass) :
		m_mass(mass) {
}

Decay::~Decay() {
	// Zerfallskanäle freigeben
	for (DecayChannel* decayChannel : m_decayChannels) {
		delete decayChannel;
	}
	
	m_decayChannels.clear();
}

void Decay::addDecayChannel(matrix_element_ptr pMatrixElement,void* pParams) {
	// Zerfallskanal der Liste hinzufügen
	DecayChannel* decayChannel = new DecayChannel { pMatrixElement, pParams };
	m_decayChannels.push_back(decayChannel);
}

void Decay::getDecayWidth(double* decayWidth,double* statError,const unsigned long numMonteCarloSteps) {
	// Falls genaueres Ergebnis bekannt, dieses zurückgeben, sonst Ergebnis berechnen lassen
	if (m_numLastMonteCarloSteps>numMonteCarloSteps) {
		*decayWidth = m_decayWidth;
		*statError = m_statError;
	}
	else
		calcDecayWidth(decayWidth,statError,numMonteCarloSteps);
}

void Decay::calcDecayWidth(double* decayWidth,double* statError,const unsigned long numMonteCarloSteps) {
	// Passendes Monte Carlo anlegen
	MonteCarlo monteCarlo(m_pPhaseSpace->getNumDimensions(),numMonteCarloSteps);
	
	// Grenzen der Integration
	double* pXl = new double[m_pPhaseSpace->getNumDimensions()];
	double* pXu = new double[m_pPhaseSpace->getNumDimensions()];
	
	for (unsigned int i=0;i<m_pPhaseSpace->getNumDimensions();i++) {
		pXl[i] = 0.;
		pXu[i] = 1.;
	}
	
	// Funktion für das Monte-Carlo erzeugen und setzen, Funktionswerte werden durch die differentielle Zerfallsbreite, die vom Zerfallsobjekt behandelt wird, gegeben
	DecayMonteCarloFunction monteCarloFunction(this,m_pPhaseSpace->getNumDimensions(),pXl,pXu);
	monteCarlo.setIntegrationFunction(&monteCarloFunction);
	
	// Integrieren
	monteCarlo.integrate(&m_decayWidth,&m_statError);
	m_numLastMonteCarloSteps = numMonteCarloSteps;
	
	// Ergebnisvariablen updaten
	*decayWidth = m_decayWidth;
	*statError = m_statError;
}

double Decay::getIntegrationVolume() const {
	// Integrationsvolumen entspricht in der Regel dem Phasenraumvolumen
	return m_pPhaseSpace->getPhaseSpaceVolume();
}

double Decay::getFunctionValue(double* pSamplingPoint) {
	// Impulse der Teilchen aus Phasenraum generieren, zerfallendes Teilchen dabei in dessen Ruhesystem betrachten
	FourVector momentum[4];
	momentum[0] = FourVector(m_mass,0.,0.,0.);
	double phaseSpaceFactor = 0.;
	
	m_pPhaseSpace->getNextPhaseSpacePoint(pSamplingPoint,momentum[0],&momentum[1],&phaseSpaceFactor);
	
	// Matrix-Elemente der einzelnen Zerfallskanäle aus den Phasenraumpunkten berechnen
	double matrixElementFactor = 0.;
	
	for (DecayChannel* decayChannel : m_decayChannels) {
		matrixElementFactor += decayChannel->pMatrixElement(momentum,decayChannel->pParams);
	}
	
	// differentielle Zerfallsbreite im Ruhesystem des zerfallenden Teilchens
	return 1./2./m_mass*matrixElementFactor*phaseSpaceFactor;
}
