#include "ThreeBodyMassivePhaseSpace.h"

ThreeBodyMassivePhaseSpace::ThreeBodyMassivePhaseSpace(const double massParticle1,const double massParticle2,const double massParticle3) {	
	// Massen speichern
	m_massParticles[0] = massParticle1;
	m_massParticles[1] = massParticle2;
	m_massParticles[2] = massParticle3;
	
	// Insgesamt 5 Dimensionen: 3 (Teilchen) * 3 (Dimensionen) - 4 (Delta-Funktion) = 5
	m_numDimensions = 5;
}

ThreeBodyMassivePhaseSpace::~ThreeBodyMassivePhaseSpace() {
}

double ThreeBodyMassivePhaseSpace::getPhaseSpaceVolume() const {
	// Phasenraum-Volumen aus den beiden einzelnen Phasenräumen
	return 1.0;
}

void ThreeBodyMassivePhaseSpace::getNextPhaseSpacePoint(const double* pSamplingPoint,FourVector const& pInTotalMomentum,FourVector* pOutMomentum,double* pPhaseSpaceFactor) {
	// Schwerpunktsenergie berechnen
	double s = pInTotalMomentum*pInTotalMomentum;
	double comEnergy = sqrt(s);
	
	// Phasenraumvariablen generieren
	double M2sqMin = (m_massParticles[0]+m_massParticles[1])*(m_massParticles[0]+m_massParticles[1]);
	double M2sqMax = (comEnergy-m_massParticles[2])*(comEnergy-m_massParticles[2]);
	double M2sq = M2sqMin + pSamplingPoint[0]*(M2sqMax-M2sqMin);
	double M2 = sqrt(M2sq);
	
	//double M2 = pSamplingPoint[0]*(comEnergy - m_massParticles[2] - m_massParticles[1] - m_massParticles[0]) + m_massParticles[1] + m_massParticles[0];
	//double M2sq = M2*M2;
	double theta2 = pSamplingPoint[1]*M_PI;
	double phi2 = pSamplingPoint[2]*2.*M_PI;
	double theta3 = pSamplingPoint[3]*M_PI;
	double phi3 = pSamplingPoint[4]*2.*M_PI;
	
	// Generiere p1 und p2 in deren Schwerpunktssystem
	double absP2 = sqrt(kaellenFunction(M2sq,m_massParticles[0]*m_massParticles[0],m_massParticles[1]*m_massParticles[1])/4./M2sq);
	double p2_0 = sqrt(m_massParticles[1]*m_massParticles[1] + absP2*absP2);
	double p2_1 = absP2 * cos(phi2) * sin(theta2);
	double p2_2 = absP2 * sin(phi2) * sin(theta2);
	double p2_3 = absP2 * cos(theta2);
	
	double p1_0 = sqrt(m_massParticles[0]*m_massParticles[0] + absP2*absP2);
	double p1_1 = -p2_1;
	double p1_2 = -p2_2;
	double p1_3 = -p2_3;
	
	// Generiere p3 im Schwerpunktssystem aller 3 Teilchen
	double absP3 = sqrt(kaellenFunction(s,M2sq,m_massParticles[2]*m_massParticles[2])/4./s);
	double p3_0 = sqrt(m_massParticles[2]*m_massParticles[2]+absP3*absP3);
	double p3_1 = absP3 * cos(phi3) * sin(theta3);
	double p3_2 = absP3 * sin(phi3) * sin(theta3);
	double p3_3 = absP3 * cos(theta3);
	
	// Boost von p1 & p2 ins Schwerpunktssystem von p1+p2+p3
	double gamma12,beta12[3];
	gamma12 = (comEnergy - p3_0)/M2;
	double gamma212[3] = { gamma12*gamma12/(1+gamma12), gamma12*gamma12/(1+gamma12), gamma12*gamma12/(1+gamma12) };
	
	beta12[0] = p3_1/(comEnergy-p3_0);
	beta12[1] = p3_2/(comEnergy-p3_0);
	beta12[2] = p3_3/(comEnergy-p3_0);
	
	for (int i=0;i<3;i++) {
		gamma212[i] *= beta12[i];
	}
	
	double p112_0 = p1_0*gamma12 - p1_1*gamma12*beta12[0] - p1_2*gamma12*beta12[1] - p1_3*gamma12*beta12[2];
	double p112_1 = -p1_0*gamma12*beta12[0] + p1_1 + p1_1 * gamma212[0] * beta12[0] + p1_2 * gamma212[1] * beta12[0] + p1_3 * gamma212[2] * beta12[0];
	double p112_2 = -p1_0*gamma12*beta12[1] + p1_1 * gamma212[0] * beta12[1] + p1_2 + p1_2 * gamma212[1] * beta12[1] + p1_3 * gamma212[2] * beta12[1];
	double p112_3 = -p1_0*gamma12*beta12[2] + p1_1 * gamma212[0] * beta12[2] + p1_2 * gamma212[1] * beta12[2] + p1_3 + p1_3 * gamma212[2] * beta12[2];
	
	double p212_0 = p2_0*gamma12 - p2_1*gamma12*beta12[0] - p2_2*gamma12*beta12[1] - p2_3*gamma12*beta12[2];
	double p212_1 = -p2_0*gamma12*beta12[0] + p2_1 + p2_1 * gamma212[0] * beta12[0] + p2_2 * gamma212[1] * beta12[0] + p2_3 * gamma212[2] * beta12[0];
	double p212_2 = -p2_0*gamma12*beta12[1] + p2_1 * gamma212[0] * beta12[1] + p2_2 + p2_2 * gamma212[1] * beta12[1] + p2_3 * gamma212[2] * beta12[1];
	double p212_3 = -p2_0*gamma12*beta12[2] + p2_1 * gamma212[0] * beta12[2] + p2_2 * gamma212[1] * beta12[2] + p2_3 + p2_3 * gamma212[2] * beta12[2];
	
	// Boost ins Laborsystem, d.h. in das System, in dem pInTotalMomentum gegeben ist
	double gammaLab,betaLab[3];
	gammaLab = pInTotalMomentum.x[0]/comEnergy;
	double gamma2Lab[3] = { gammaLab*gammaLab/(1+gammaLab), gammaLab*gammaLab/(1+gammaLab), gammaLab*gammaLab/(1+gammaLab) };
	
	for (int i=1;i<4;i++) {
		betaLab[i-1] = -pInTotalMomentum.x[i]/pInTotalMomentum.x[0];
		gamma2Lab[i-1] *= betaLab[i-1];
	}
	
	// evtl. letzten Phasenraum-Punkt löschen
	if (m_lastPhaseSpacePoint!=NULL) 
		delete[] m_lastPhaseSpacePoint;
	
	// Neuen Phasenraum-Punkt erzeugen	
	m_lastPhaseSpacePoint = new FourVector[3];
	
	// Boost durchführen
	m_lastPhaseSpacePoint[0] = FourVector(p112_0*gammaLab - p112_1*gammaLab*betaLab[0] - p112_2*gammaLab*betaLab[1] - p112_3*gammaLab*betaLab[2],
											  -p112_0*gammaLab*betaLab[0] + p112_1 + p112_1*gamma2Lab[0]*betaLab[0] + p112_2*gamma2Lab[1]*betaLab[0] + p112_3*gamma2Lab[2]*betaLab[0],
											  -p112_0*gammaLab*betaLab[1] + p112_1*gamma2Lab[0]*betaLab[1] + p112_2 + p112_2*gamma2Lab[1]*betaLab[1] + p112_3*gamma2Lab[2]*betaLab[1],
											  -p112_0*gammaLab*betaLab[2] + p112_1*gamma2Lab[0]*betaLab[2] + p112_2*gamma2Lab[1]*betaLab[2] + p112_3 + p112_3*gamma2Lab[2]*betaLab[2]);
	m_lastPhaseSpacePoint[1] = FourVector(p212_0*gammaLab - p212_1*gammaLab*betaLab[0] - p212_2*gammaLab*betaLab[1] - p212_3*gammaLab*betaLab[2],
											  -p212_0*gammaLab*betaLab[0] + p212_1 + p212_1*gamma2Lab[0]*betaLab[0] + p212_2*gamma2Lab[1]*betaLab[0] + p212_3*gamma2Lab[2]*betaLab[0],
											  -p212_0*gammaLab*betaLab[1] + p212_1*gamma2Lab[0]*betaLab[1] + p212_2 + p212_2*gamma2Lab[1]*betaLab[1] + p212_3*gamma2Lab[2]*betaLab[1],
											  -p212_0*gammaLab*betaLab[2] + p212_1*gamma2Lab[0]*betaLab[2] + p212_2*gamma2Lab[1]*betaLab[2] + p212_3 + p212_3*gamma2Lab[2]*betaLab[2]);
	m_lastPhaseSpacePoint[2] = FourVector(p3_0*gammaLab - p3_1*gammaLab*betaLab[0] - p3_2*gammaLab*betaLab[1] - p3_3*gammaLab*betaLab[2],
											  -p3_0*gammaLab*betaLab[0] + p3_1 + p3_1*gamma2Lab[0]*betaLab[0] + p3_2*gamma2Lab[1]*betaLab[0] + p3_3*gamma2Lab[2]*betaLab[0],
											  -p3_0*gammaLab*betaLab[1] + p3_1*gamma2Lab[0]*betaLab[1] + p3_2 + p3_2*gamma2Lab[1]*betaLab[1] + p3_3*gamma2Lab[2]*betaLab[1],
											  -p3_0*gammaLab*betaLab[2] + p3_1*gamma2Lab[0]*betaLab[2] + p3_2*gamma2Lab[1]*betaLab[2] + p3_3 + p3_3*gamma2Lab[2]*betaLab[2]);
	
	// Phasenraumgewicht
	double V = ((comEnergy-m_massParticles[2])*(comEnergy-m_massParticles[2]) - (m_massParticles[0]+m_massParticles[1])*(m_massParticles[0]+m_massParticles[1]))*2.*M_PI*M_PI*2.*M_PI*M_PI;
	*pPhaseSpaceFactor = sqrt(kaellenFunction(s,M2sq,m_massParticles[2]*m_massParticles[2]))*sqrt(kaellenFunction(M2sq,m_massParticles[0]*m_massParticles[0],m_massParticles[1]*m_massParticles[1]))/64./s/M2sq*sin(theta2)*sin(theta3)*V/2./M_PI/2./M_PI/2./M_PI/2./M_PI/2./M_PI;
	
	// Ins Ergebnis-Array schreiben
	pOutMomentum[0] = m_lastPhaseSpacePoint[0];
	pOutMomentum[1] = m_lastPhaseSpacePoint[1];
	pOutMomentum[2] = m_lastPhaseSpacePoint[2];
	
}

double ThreeBodyMassivePhaseSpace::kaellenFunction(const double x1,const double x2,const double x3) const {
	return x1*x1 + x2*x2 + x3*x3 - 2*x1*x2 -2*x1*x3 - 2*x2*x3;
}
