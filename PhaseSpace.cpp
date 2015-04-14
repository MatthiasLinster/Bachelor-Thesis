#include "PhaseSpace.h"

PhaseSpace::PhaseSpace() {
}

PhaseSpace::~PhaseSpace() {
	// Speicher des letzten Phasenraum-Punktes freigeben, falls vorhanden
	if (m_lastPhaseSpacePoint!=NULL) {
		delete[] m_lastPhaseSpacePoint;
		m_lastPhaseSpacePoint = NULL;
	}
}

size_t PhaseSpace::getNumDimensions() const {
	return m_numDimensions;
}

double PhaseSpace::getPhaseSpaceVolume() const {
	return 1.0;
}

FourVector* PhaseSpace::getLastPhaseSpacePoint() const {
	return m_lastPhaseSpacePoint;
}
