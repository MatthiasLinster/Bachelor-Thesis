#ifndef THREE_BODY_MASSIVE_PHASE_SPACE_H
#define THREE_BODY_MASSIVE_PHASE_SPACE_H

#include "PhaseSpace.h"
#include "TwoBodyMassivePhaseSpace.h"

#include <iostream>

class ThreeBodyMassivePhaseSpace : public PhaseSpace {
	
	public:
		ThreeBodyMassivePhaseSpace(const double massParticle1,const double massParticle2,const double massParticle3);
		~ThreeBodyMassivePhaseSpace();
		
		// Volumen des Phasenraums überschreiben
		virtual double getPhaseSpaceVolume() const;
		// Implementierung der Phasenraum-Punkt-Funktion
		virtual void getNextPhaseSpacePoint(const double* pSamplingPoint,FourVector const& pInTotalMomentum,FourVector* pOutMomentum,double* pPhaseSpaceFactor);
		
	protected:
		// Källen-Funktion
		double kaellenFunction(const double x1,const double x2,const double x3) const;
	  
		// Massen der Teilchen
		double m_massParticles[3];
	
		// Zwei Phasenraum-Objekte zur Generierung der Teilchen
		TwoBodyMassivePhaseSpace* m_ps1,*m_ps2;

};

#endif
