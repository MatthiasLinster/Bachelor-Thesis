#ifndef THREE_BODY_MASSIVE_DECAY_H
#define THREE_BODY_MASSIVE_DECAY_H

#include "Decay.h"
#include "ThreeBodyMassivePhaseSpace.h"
#include "FunctionTypes.h"

/**
 * Klasse für den Zerfall eines Teilchens in drei massive Teilchen im Endzustand
 * @author Matthias Linster
 * @date 16.09.2014
 */
class ThreeBodyMassiveDecay : public Decay {
	
	public:
		/**
		 * Konstruktor, falls kein Matrix-Element angegeben wird, muss dieses über die addDecayChannel-Methode hinzugefügt werden!
		 * @param massDecayingParticle Masse des zerfallenden Teilchens in dessen Ruhesystem
		 * @param pMatrixElement Matrix-Element des Zerfalls
		 * @param pParams Parameter des Matrix-Elements
		 */
		ThreeBodyMassiveDecay(const double massDecayingParticle,const double massParticle1,const double massParticle2,const double massParticle3,matrix_element_ptr pMatrixElement = NULL,void* pParams = NULL);
		/**
		 * Destruktor
		 */
		~ThreeBodyMassiveDecay();
	
};

#endif
