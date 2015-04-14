#ifndef PHASE_SPACE_H
#define PHASE_SPACE_H

#include <stdlib.h>

#include "FourVector.h"

/**
 * Klasse PhaseSpace als Mutterklasse für alle Phasenraum-Klassen.
 * Stellt grundlegende Methoden bereit, die alle Phasenraum-Klassen haben sollten.
 * @author Matthias Linster
 * @date 03.09.2014
 */
class PhaseSpace {

	public:
		/**
		 * Default-Konstruktor und Destruktor
		 */
		PhaseSpace();
		virtual ~PhaseSpace();
		
		/**
		 * Liefert das Volumen des Phasenraums z.B. für die Monte-Carlo-Simulation
		 * @return Volumen des Phasenraums
		 */
		virtual double getPhaseSpaceVolume() const;
		
		/**
		 * Liefert die Anzahl der Dimensionen des Phasenraums
		 * @return Anzahl der Phasenraumdimensionen
		 */
		virtual size_t getNumDimensions() const;
		
		/**
		 * Liefert einen Phasenraum-Punkt basierend auf den übergebenen Zufallszahlen.
		 * Benötigt genauso viele Zufallszahlen wie der Phasenraum Dimensionen hat.
		 * @param pSamplingPoint Array aus Zufallszahlen, es müssen so viele übergeben werden, wie der Phasenraum Dimensionen hat
		 * @param pInTotalMomentum Gesamt-4-Impuls der einlaufenden Teilchen
		 * @param pOutMomentum Array der 4-Impulse der beiden auslaufenden Teilchen
		 * @param pPhaseSpaceFactor Zeiger auf die Variablen, in dem der Faktor, den der Phasenraum zum Integral beiträgt, abgelegt werden soll
		 */
		virtual void getNextPhaseSpacePoint(const double* pSamplingPoint,FourVector const& pInTotalMomentum,FourVector* pOutMomentum,double* pPhaseSpaceFactor) = 0;
		
		/**
		 * Liefert den zuletzt generierten Punkt des Phasenraums
		 * @return Array von 4-Impulsen für alle Teilchen (einlaufend in den ersten Komponenten sowie auslaufend in den hinteren)
		 */
		virtual FourVector* getLastPhaseSpacePoint() const;
		
	protected:
		// Hilfsvariable für den letzten Phasenraum-Punkt
		FourVector* m_lastPhaseSpacePoint = NULL;
		// Anzahl Dimensionen des Phasenraums, sollte vom Konstruktor der abgeleiteten Klasse gesetzt werden
		size_t m_numDimensions = 1;

};

#endif
