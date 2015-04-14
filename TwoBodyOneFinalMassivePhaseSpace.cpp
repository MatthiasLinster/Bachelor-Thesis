#include "TwoBodyOneFinalMassivePhaseSpace.h"

TwoBodyOneFinalMassivePhaseSpace::TwoBodyOneFinalMassivePhaseSpace(double mass) {
	// Masse speichern, 2 Dimensionen (Winkel theta und phi)
	m_massParticle1 = mass;
	m_numDimensions = 2;
}

TwoBodyOneFinalMassivePhaseSpace::~TwoBodyOneFinalMassivePhaseSpace() {
}

void TwoBodyOneFinalMassivePhaseSpace::setMass(double mass) {
	m_massParticle1 = mass;
}

double TwoBodyOneFinalMassivePhaseSpace::getPhaseSpaceVolume() const {
	// Integration in theta -> Volumen PI-0 = PI und phi -> Volumen 2*PI
	return M_PI*2.*M_PI;
}

void TwoBodyOneFinalMassivePhaseSpace::getNextPhaseSpacePoint(const double* pSamplingPoint,FourVector const& pInTotalMomentum,FourVector* pOutMomentum,double* pPhaseSpaceFactor) {
	// Winkel zwischen 0 und PI bzw. 0 und 2*PI erzeugen
	double theta = pSamplingPoint[0]*M_PI;
	double phi = pSamplingPoint[1]*2.*M_PI;
	
	// Schwerpunktsenergie berechnen
	double s = pInTotalMomentum*pInTotalMomentum;
	double comEnergy = sqrt(s);
	
	// Labor-Boost-Faktoren bestimmen
	double gamma,beta[3];
	gamma = pInTotalMomentum.x[0]/comEnergy;
	double gamma2[3] = { gamma*gamma/(1+gamma), gamma*gamma/(1+gamma), gamma*gamma/(1+gamma) };
	
	for (int i=1;i<4;i++) {
		beta[i-1] = -pInTotalMomentum.x[i]/pInTotalMomentum.x[0];
		gamma2[i-1] *= beta[i-1];
	}

	// evtl. letzten Phasenraum-Punkt löschen
	if (m_lastPhaseSpacePoint!=NULL) 
		delete[] m_lastPhaseSpacePoint;
	
	// Neuen Phasenraum-Punkt erzeugen	
	m_lastPhaseSpacePoint = new FourVector[2];
	
	// Betrag des Impulses der beiden äußeren Teilchen
	double absMomentum = (s-m_massParticle1*m_massParticle1)/2./comEnergy;
	
	// Impulskomponenten des 3. massiven Teilchens im CMS
	double p3_0 = sqrt(m_massParticle1*m_massParticle1+absMomentum*absMomentum);
	double p3_1 = absMomentum*cos(phi)*sin(theta);
	double p3_2 = absMomentum*sin(phi)*sin(theta);
	double p3_3 = absMomentum*cos(theta);
	
	// Impulskomponenten des 4. masselosen Teilchens im CMS
	double p4_0 = absMomentum;
	double p4_1 = -p3_1;
	double p4_2 = -p3_2;
	double p4_3 = -p3_3;
	
	// Impuls des 3. und 4. Teilchens boosten
	m_lastPhaseSpacePoint[0] = FourVector(p3_0*gamma - p3_1*gamma*beta[0] - p3_2*gamma*beta[1] - p3_3*gamma*beta[2],
										  -p3_0*gamma*beta[0] + p3_1 + p3_1*gamma2[0]*beta[0] + p3_2*gamma2[1]*beta[0] + p3_3*gamma2[2]*beta[0],
										  -p3_0*gamma*beta[1] + p3_1*gamma2[0]*beta[1] + p3_2 + p3_2*gamma2[1]*beta[1] + p3_3*gamma2[2]*beta[1],
										  -p3_0*gamma*beta[2] + p3_1*gamma2[0]*beta[2] + p3_2*gamma2[1]*beta[2] + p3_3 + p3_3*gamma2[2]*beta[2]);
	m_lastPhaseSpacePoint[1] = FourVector(p4_0*gamma - p4_1*gamma*beta[0] - p4_2*gamma*beta[1] - p4_3*gamma*beta[2],
										  -p4_0*gamma*beta[0] + p4_1 + p4_1*gamma2[0]*beta[0] + p4_2*gamma2[1]*beta[0] + p4_3*gamma2[2]*beta[0],
										  -p4_0*gamma*beta[1] + p4_1*gamma2[0]*beta[1] + p4_2 + p4_2*gamma2[1]*beta[1] + p4_3*gamma2[2]*beta[1],
										  -p4_0*gamma*beta[2] + p4_1*gamma2[0]*beta[2] + p4_2*gamma2[1]*beta[2] + p4_3 + p4_3*gamma2[2]*beta[2]);
	
	// Ins Ergebnis-Array schreiben
	pOutMomentum[0] = m_lastPhaseSpacePoint[0];
	pOutMomentum[1] = m_lastPhaseSpacePoint[1];
	
	// Phasenraum-Faktor (Theta-Funktion dabei berücksichtigen)
	if (comEnergy>m_massParticle1)
		*pPhaseSpaceFactor = 1./32./M_PI/M_PI*(s-m_massParticle1*m_massParticle1)/s*sin(theta);
	else
		*pPhaseSpaceFactor = 0.;
}
