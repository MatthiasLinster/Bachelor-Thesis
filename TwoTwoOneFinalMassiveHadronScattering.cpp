#include "TwoTwoOneFinalMassiveHadronScattering.h"

TwoTwoOneFinalMassiveHadronScattering::TwoTwoOneFinalMassiveHadronScattering(const double massFinal,const double comEnergyCollider,const double factorizationScale,LHAPDF::PDF* pPdf) {
	// Werte speichern
	m_massFinal = massFinal;
	m_comEnergyCollider = comEnergyCollider;
	m_factorizationScale = factorizationScale;
	m_pPdf = pPdf;
	
	// Phasenraum erzeugen
	m_pPhaseSpace = new TwoBodyOneFinalMassivePhaseSpace(massFinal);
	
	// Objekt initialisieren
	double* pXl = new double[4] { 0.,0.,0.,0. };
	double* pXu = new double[4] { 1.,1.,1.,1. };
	
	initScattering(4,pXl,pXu);
	
	delete[] pXl;
	delete[] pXu;
}

TwoTwoOneFinalMassiveHadronScattering::~TwoTwoOneFinalMassiveHadronScattering() {
	// Phasenraum löschen
	delete m_pPhaseSpace;
}

void TwoTwoOneFinalMassiveHadronScattering::addPartonProcess(int parton1,int parton2,matrix_element_ptr pMatrixElement,void* pParams) {
	// Prozess der Liste hinzufügen
	PartonProcess partonProcess = { parton1, parton2, pMatrixElement, pParams };
	m_partonProcesses.push_back(partonProcess);
}

double TwoTwoOneFinalMassiveHadronScattering::getIntegrationVolume() const {
	// Phasenraumvolumen zurückgeben
	return m_pPhaseSpace->getPhaseSpaceVolume();
}

double TwoTwoOneFinalMassiveHadronScattering::getFunctionValue(double* pSamplingPoint) {
	// Funktionswert noch 0
	double fx = 0.;
	
	// Array für die 4-Impulse aller beteiligten Teilchen (Anfangs- & Endzustand)
	FourVector particleMomentum[4];
	
	// 4-Impulse der einlaufenden Partonen im Laborsystem erzeugen
	particleMomentum[0] = FourVector(1.,0.,0.,1.);
	particleMomentum[0] *= pSamplingPoint[0] * m_comEnergyCollider/2.;
		
	particleMomentum[1] = FourVector(1.,0.,0.,-1.);
	particleMomentum[1] *= pSamplingPoint[1] * m_comEnergyCollider/2.;
		
	// Gesamtimpuls und Schwerpunktsenergie berechnen
	FourVector sum = particleMomentum[0] + particleMomentum[1];
	double s = sum*sum;
		
	// Phasenraum-Punkte generieren lassen
	double phaseSpaceFactor;
	m_pPhaseSpace->getNextPhaseSpacePoint(&pSamplingPoint[2],sum,&particleMomentum[2],&phaseSpaceFactor);
		
	// Flussfaktor
	double fluxFactor = 1./2./s;
	
	// Prozesse addieren
	for (PartonProcess const& process : m_partonProcesses) {
		fx += m_pPdf->xfxQ(process.parton1,pSamplingPoint[0],m_factorizationScale)/pSamplingPoint[0]*m_pPdf->xfxQ(process.parton2,pSamplingPoint[1],m_factorizationScale)/pSamplingPoint[1] * process.pMatrixElement(particleMomentum,process.pParams);
	}
	
	// Flussfaktor, Phasenraumfaktor berücksichtigen
	fx *= fluxFactor * phaseSpaceFactor;
	
	return fx;
}
