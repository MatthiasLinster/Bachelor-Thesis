#ifndef THREE_BODY_MASSLESS_PHASE_SPACE_H
#define THREE_BODY_MASSLESS_PHASE_SPACE_H

#include "PhaseSpace.h"
#include "TwoBodyOneFinalMassivePhaseSpace.h"

class ThreeBodyMasslessPhaseSpace : public PhaseSpace {
	
	public:
		ThreeBodyMasslessPhaseSpace();
		~ThreeBodyMasslessPhaseSpace();
		
		// Volumen des Phasenraums überschreiben
		virtual double getPhaseSpaceVolume() const;
		// Implementierung der Phasenraum-Punkt-Funktion
		virtual void getNextPhaseSpacePoint(const double* pSamplingPoint,FourVector const& pInTotalMomentum,FourVector* pOutMomentum,double* pPhaseSpaceFactor);
		
	private:
		// Zwei Phasenraum-Objekte zur Generierung der Teilchen
		TwoBodyOneFinalMassivePhaseSpace* m_ps1,*m_ps2;

};

#endif
