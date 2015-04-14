#ifndef TWO_BODY_MASSIVE_PHASE_SPACE_H
#define TWO_BODY_MASSIVE_PHASE_SPACE_H

#include <math.h>

#include "PhaseSpace.h"

/**
 * Klasse für den Zwei-Körper-Phasenraum mit zwei massiven Teilchen im Endzustand.
 * @author Matthias Linster
 * @date 03.09.2014
 */
class TwoBodyMassivePhaseSpace : public PhaseSpace {
	
	public:
		/**
		 * Konstruktor, erzeugt den Phasenraum.
		 * Übergeben werden muss die Masse des massiven Teilchens.
		 * @param mass Masse des 1. Teilchens (= massiv)
		 */
		TwoBodyMassivePhaseSpace(const double massParticle1,const double massParticle2);
		~TwoBodyMassivePhaseSpace();
		
		void setMass(const unsigned int particle,const double mass);
		// Volumen des Phasenraums überschreiben
		virtual double getPhaseSpaceVolume() const;
		// Implementierung der Phasenraum-Punkt-Funktion
		virtual void getNextPhaseSpacePoint(const double* pSamplingPoint,FourVector const& pInTotalMomentum,FourVector* pOutMomentum,double* pPhaseSpaceFactor);
		
	protected:
		// Masse der Teilchen
		double m_massParticles[2];
		
		/**
		 * Berechnet die Wurzel der Kaellen-Funktion sqrt(x1^2 + x2^2 + x3^2 - 2*x1*x2 - 2*x1*x3 . 2*x2*x3)
		 * @param x1 Erster Parameter
		 * @param x2 Zweiter Parameter
		 * @param x3 Dritter Parameter
		 * @return Wurzel der Kaellen-Funktion sqrt(x1^2 + x2^2 + x3^2 - 2*x1*x2 - 2*x1*x3 . 2*x2*x3)
		 */
		double sqrtKaellenFunction(const double x1,const double x2,const double x3) const;
	
};

#endif
