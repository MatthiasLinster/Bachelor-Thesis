#include "ThreeBodyMasslessPhaseSpace.h"

ThreeBodyMasslessPhaseSpace::ThreeBodyMasslessPhaseSpace() {	
	// Zwei 2-Körper-Phasenräume für Rekursion
	m_ps1 = new TwoBodyOneFinalMassivePhaseSpace(0);
	m_ps2 = new TwoBodyOneFinalMassivePhaseSpace(0);
	
	// Insgesamt 5 Dimensionen: 3 (Teilchen) * 3 (Dimensionen) - 4 (Delta-Funktion) = 5
	m_numDimensions = 5;
}

ThreeBodyMasslessPhaseSpace::~ThreeBodyMasslessPhaseSpace() {
	// Angelegte Phasenräume löschen
	delete m_ps1;
	delete m_ps2;
}

double ThreeBodyMasslessPhaseSpace::getPhaseSpaceVolume() const {
	// Phasenraum-Volumen aus den beiden einzelnen Phasenräumen
	return m_ps1->getPhaseSpaceVolume()*m_ps2->getPhaseSpaceVolume();
}

void ThreeBodyMasslessPhaseSpace::getNextPhaseSpacePoint(const double* pSamplingPoint,FourVector const& pInTotalMomentum,FourVector* pOutMomentum,double* pPhaseSpaceFactor) {
	// Invariante Masse des einlaufenden Teilchens
	double invMassSqIn = pInTotalMomentum*pInTotalMomentum;
	
	// Invariante Masse des "mittleren" Teilchens erzeugen
	double qMassSq = pSamplingPoint[0]*invMassSqIn;
	
	// 1. Phasenraum das erste auslaufende Teilchen sowie das "mittlere" Teilchen generieren lassen
	m_ps1->setMass(sqrt(qMassSq));
	FourVector momPS1[2];
	double factorPS1 = 0.;
	m_ps1->getNextPhaseSpacePoint(&pSamplingPoint[1],pInTotalMomentum,momPS1,&factorPS1);
	
	// Aus dem "mittleren" Teilchen den zweiten Phasenraum die beiden anderen äußeren Teilchen generieren lassen
	FourVector momPS2[2];
	double factorPS2 = 0.;
	m_ps2->getNextPhaseSpacePoint(&pSamplingPoint[3],momPS1[0],momPS2,&factorPS2);
	
	// Alle erzeugten Impulse in das Ergebnis-Array kopieren
	pOutMomentum[0] = FourVector(momPS1[1].x[0],momPS1[1].x[1],momPS1[1].x[2],momPS1[1].x[3]);
	pOutMomentum[1] = FourVector(momPS2[0].x[0],momPS2[0].x[1],momPS2[0].x[2],momPS2[0].x[3]);
	pOutMomentum[2] = FourVector(momPS2[1].x[0],momPS2[1].x[1],momPS2[1].x[2],momPS2[1].x[3]);
	
	// Phasenraum-Faktor bestimmen
	*pPhaseSpaceFactor = 1./2./M_PI*factorPS1*factorPS2*invMassSqIn;
}
