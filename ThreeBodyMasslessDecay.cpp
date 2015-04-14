#include "ThreeBodyMasslessDecay.h"

ThreeBodyMasslessDecay::ThreeBodyMasslessDecay(double massDecayingParticle,matrix_element_ptr pMatrixElement,void* pParams) :
		Decay(massDecayingParticle) {
	// Werte speichern
	if (pMatrixElement!=NULL)
		addDecayChannel(pMatrixElement,pParams);
	
	// Phasenraum erzeugen
	m_pPhaseSpace = new ThreeBodyMasslessPhaseSpace();
}

ThreeBodyMasslessDecay::~ThreeBodyMasslessDecay() {
	// Phasenraum löschen
	delete m_pPhaseSpace;
}
