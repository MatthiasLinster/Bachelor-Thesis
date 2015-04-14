#ifndef THREE_BODY_MASSLESS_DECAY_H
#define THREE_BODY_MASSLESS_DECAY_H

#include "Decay.h"
#include "ThreeBodyMasslessPhaseSpace.h"
#include "FunctionTypes.h"

/**
 * Klasse für den Zerfall eines Teilchens in drei masselose Teilchen im Endzustand
 * @author Matthias Linster
 * @date 16.09.2014
 */
class ThreeBodyMasslessDecay : public Decay {
	
	public:
		/**
		 * Konstruktor, falls kein Matrix-Element angegeben wird, muss dieses über die addDecayChannel-Methode hinzugefügt werden!
		 * @param massDecayingParticle Masse des zerfallenden Teilchens in dessen Ruhesystem
		 * @param pMatrixElement Matrix-Element des Zerfalls
		 * @param pParams Parameter des Matrix-Elements
		 */
		ThreeBodyMasslessDecay(const double massDecayingParticle,matrix_element_ptr pMatrixElement = NULL,void* pParams = NULL);
		/**
		 * Destruktor
		 */
		~ThreeBodyMasslessDecay();
	
};

#endif
