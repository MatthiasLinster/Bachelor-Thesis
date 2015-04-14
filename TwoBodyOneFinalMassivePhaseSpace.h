#ifndef TWO_BODY_ONE_FINAL_MASSIVE_PHASE_SPACE_H
#define TWO_BODY_ONE_FINAL_MASSIVE_PHASE_SPACE_H

#include <math.h>

#include "PhaseSpace.h"

/**
 * Klasse für den Zwei-Körper-Phasenraum mit einem massiven Teilchen im Endzustand, alle anderen Teilchen sind masselos.
 * Das massive Teilchen entspricht jeweils dem ersten Teilchen.
 * @author Matthias Linster
 * @date 03.09.2014
 */
class TwoBodyOneFinalMassivePhaseSpace : public PhaseSpace {
	
	public:
		/**
		 * Konstruktor, erzeugt den Phasenraum.
		 * Übergeben werden muss die Masse des massiven Teilchens.
		 * @param mass Masse des 1. Teilchens (= massiv)
		 */
		TwoBodyOneFinalMassivePhaseSpace(double mass);
		~TwoBodyOneFinalMassivePhaseSpace();
		
		void setMass(double mass);
		// Volumen des Phasenraums überschreiben
		virtual double getPhaseSpaceVolume() const;
		// Implementierung der Phasenraum-Punkt-Funktion
		virtual void getNextPhaseSpacePoint(const double* pSamplingPoint,FourVector const& pInTotalMomentum,FourVector* pOutMomentum,double* pPhaseSpaceFactor);
		
	protected:
		// Masse des 1. Teilchens
		double m_massParticle1;
	
};

#endif
