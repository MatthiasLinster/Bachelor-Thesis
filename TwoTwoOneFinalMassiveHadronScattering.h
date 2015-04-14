#ifndef TWO_TWO_ONE_FINAL_MASSIVE_HADRON_SCATTERING_H
#define TWO_TWO_ONE_FINAL_MASSIVE_HADRON_SCATTERING_H

#include <stdlib.h>
#include <list>

#include "LHAPDF/LHAPDF.h"

#include "FourVector.h"
#include "Scattering.h"
#include "TwoBodyOneFinalMassivePhaseSpace.h"
#include "MonteCarlo.h"
#include "ScatteringMonteCarloFunction.h"
#include "FunctionTypes.h"

// Struktur für Parton-Parton-Subprozesse
struct PartonProcess {
	// Beteiligte Partonen in beiden Hadronen
	int parton1,parton2;
	// Matrixelement des aktuellen Prozesses
	matrix_element_ptr pMatrixElement;
	// Parameter des Matrix-Elements
	void* pParams;
};

/**
 * Klasse für den harten Streuprozess zwischen zwei Hadronen mit einem massiven Teilchen im Endzustand
 * @author Matthias Linster
 * @date 03.09.2014
 */
class TwoTwoOneFinalMassiveHadronScattering : public Scattering {
	
	public:
		/**
		 * Kontruktor zur Erzeugung eines Objektes mit einer definierten Schwerpunktsenergie und Parton-Verteilungsfunktion
		 * @param mass Masse des massiven Teilchens im Endzustand
		 * @param comEnergyCollider Schwerpunktsenergie der kollidierenden Hadronen
		 * @param renormalizationScale Skala für die Partondichte-Verteilungen
		 * @param pPdf Zeiger auf das LHAPDF-Objekt mit der Parton-Verteilungsfunktion
		 */
		TwoTwoOneFinalMassiveHadronScattering(const double mass,const double comEnergyCollider,const double factorizationScale,LHAPDF::PDF* pPdf);
		
		/**
		 * Destruktor
		 */
		~TwoTwoOneFinalMassiveHadronScattering();
		
		/**
		 * Fügt einen Prozess zwischen zwei kollidierenden Partonen der beiden Hadronen hinzu.
		 * @param parton1 PDG-ID des 1. Partons
		 * @param parton2 PDG-ID des 2. Partons
		 * @param pMatrixElement Zeiger auf das Matrixelement des Prozesses
		 * @param pParams Zeiger auf die Parameter des Matrixelements
		 */
		void addPartonProcess(const int parton1,const int parton2,matrix_element_ptr pMatrixElement,void* pParams = NULL);
		
		/**
		 * Liefert das Integrationsvolumen des Streuprozesses
		 * @return Volumen des Integrationsbereiches
		 */
		virtual double getIntegrationVolume() const;
		
		/**
		 * Liefert den Funktionswert für die Integration, d.h. den differentiellen Wirkungsquerschnitt inkl. Matrix-Element, Flussfaktor und Phasenraum
		 * @param pSamplingPoint Punkt, an dem die Funktion ausgewertet werden soll
		 * @return Funktionswert an der entsprechenden Stelle
		 */
		virtual double getFunctionValue(double* pSamplingPoint);
		
	protected:
		// Masse des massiven Teilchens im Endzustand
		double m_massFinal;
		// Schwerpunktsenergie des Hadron-Hadron-Systems
		double m_comEnergyCollider;
		// Energie-Skala für LHAPDF
		double m_factorizationScale;
		// Zeiger auf das LHAPDF-Objekt mit der Parton-Verteilungsfunktion
		LHAPDF::PDF* m_pPdf;
		// Liste aller beitragenden Prozesse
		std::list<PartonProcess> m_partonProcesses;
	
};

#endif
