#include "TwoBodyMassivePhaseSpace.h"

TwoBodyMassivePhaseSpace::TwoBodyMassivePhaseSpace(const double massParticle1,const double massParticle2) {
	// Masse speichern, 2 Dimensionen (Winkel theta und phi)
	m_massParticles[0] = massParticle1;
	m_massParticles[1] = massParticle2;
	m_numDimensions = 2;
}

TwoBodyMassivePhaseSpace::~TwoBodyMassivePhaseSpace() {
}

void TwoBodyMassivePhaseSpace::setMass(const unsigned int particle,const double mass) {
	m_massParticles[particle-1] = mass;
}

double TwoBodyMassivePhaseSpace::getPhaseSpaceVolume() const {
	// Integration in theta -> Volumen PI-0 = PI und phi -> Volumen 2*PI
	return M_PI*2.*M_PI;
}

void TwoBodyMassivePhaseSpace::getNextPhaseSpacePoint(const double* pSamplingPoint,FourVector const& pInTotalMomentum,FourVector* pOutMomentum,double* pPhaseSpaceFactor) {
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
	
	// Theta-Funktion
	if (comEnergy>(m_massParticles[0]+m_massParticles[1])) {
		// Phasenraum-Faktor
		*pPhaseSpaceFactor = 1./32./M_PI/M_PI/s*sqrtKaellenFunction(s,m_massParticles[0]*m_massParticles[0],m_massParticles[1]*m_massParticles[1])*sin(theta);
		
		// Betrag des Impulses der beiden äußeren Teilchen
		double absMomentum = sqrtKaellenFunction(s,m_massParticles[0]*m_massParticles[0],m_massParticles[1]*m_massParticles[1])/2./comEnergy;
		
		// Impulskomponenten des 3. massiven Teilchens im CMS
		double p1_0 = sqrt(m_massParticles[0]*m_massParticles[0]+absMomentum*absMomentum);
		double p1_1 = absMomentum*cos(phi)*sin(theta);
		double p1_2 = absMomentum*sin(phi)*sin(theta);
		double p1_3 = absMomentum*cos(theta);
		
		// Impulskomponenten des 4. masselosen Teilchens im CMS
		double p2_0 = sqrt(m_massParticles[1]*m_massParticles[1]+absMomentum*absMomentum);
		double p2_1 = -p1_1;
		double p2_2 = -p1_2;
		double p2_3 = -p1_3;
		
		// Impuls des 3. und 4. Teilchens boosten
		m_lastPhaseSpacePoint[0] = FourVector(p1_0*gamma - p1_1*gamma*beta[0] - p1_2*gamma*beta[1] - p1_3*gamma*beta[2],
											  -p1_0*gamma*beta[0] + p1_1 + p1_1*gamma2[0]*beta[0] + p1_2*gamma2[1]*beta[0] + p1_3*gamma2[2]*beta[0],
											  -p1_0*gamma*beta[1] + p1_1*gamma2[0]*beta[1] + p1_2 + p1_2*gamma2[1]*beta[1] + p1_3*gamma2[2]*beta[1],
											  -p1_0*gamma*beta[2] + p1_1*gamma2[0]*beta[2] + p1_2*gamma2[1]*beta[2] + p1_3 + p1_3*gamma2[2]*beta[2]);
		m_lastPhaseSpacePoint[1] = FourVector(p2_0*gamma - p2_1*gamma*beta[0] - p2_2*gamma*beta[1] - p2_3*gamma*beta[2],
											  -p2_0*gamma*beta[0] + p2_1 + p2_1*gamma2[0]*beta[0] + p2_2*gamma2[1]*beta[0] + p2_3*gamma2[2]*beta[0],
											  -p2_0*gamma*beta[1] + p2_1*gamma2[0]*beta[1] + p2_2 + p2_2*gamma2[1]*beta[1] + p2_3*gamma2[2]*beta[1],
											  -p2_0*gamma*beta[2] + p2_1*gamma2[0]*beta[2] + p2_2*gamma2[1]*beta[2] + p2_3 + p2_3*gamma2[2]*beta[2]);
	}
	else {
		// Phasenraum-Faktor
		*pPhaseSpaceFactor = 0.;
		
		// Dummy-Phasenraum-Punkt erzeugen (um NaN-Fehler zu verhindern)
		m_lastPhaseSpacePoint[0] = FourVector(0,0,0,0);
		m_lastPhaseSpacePoint[1] = FourVector(0,0,0,0);
	}
	
	// Ins Ergebnis-Array schreiben
	pOutMomentum[0] = m_lastPhaseSpacePoint[0];
	pOutMomentum[1] = m_lastPhaseSpacePoint[1];
}

double TwoBodyMassivePhaseSpace::sqrtKaellenFunction(const double x1,const double x2,const double x3) const {
	return sqrt(x1*x1 + x2*x2 + x3*x3 - 2*x1*x2 -2*x1*x3 - 2*x2*x3);
}
