#include "ThreeBodyMassiveDecay.h"

ThreeBodyMassiveDecay::ThreeBodyMassiveDecay(double massDecayingParticle,const double massParticle1,const double massParticle2,const double massParticle3,matrix_element_ptr pMatrixElement,void* pParams) :
		Decay(massDecayingParticle) {
	// Werte speichern
	if (pMatrixElement!=NULL)
		addDecayChannel(pMatrixElement,pParams);
	
	// Phasenraum erzeugen
	m_pPhaseSpace = new ThreeBodyMassivePhaseSpace(massParticle1,massParticle2,massParticle3);
}

ThreeBodyMassiveDecay::~ThreeBodyMassiveDecay() {
	// Phasenraum löschen
	delete m_pPhaseSpace;
}
