#include <stdio.h>

#include "ThreeBodyMassiveDecay.h"

// Standard Anzahl an Schritten, wenn dies vom Benutzer nicht festgelegt wird
#define DEFAULT_STEP_COUNT 1000000

double matrixElement(FourVector* pMomentum,void* pParams) {
	return 48. * 2. * 0.76845 * (pMomentum[3]*pMomentum[1]) * (pMomentum[3]*pMomentum[2]);
}

int main(int argc,char** argv) {
	// Anzahl an Schritten
	unsigned long numSteps = DEFAULT_STEP_COUNT;

	// Überprüfen, ob die Anzahl an Schritten als Argument übergeben wurde, wenn ja, Variable updaten
	if (argc!=2) 
		printf("Keine Schrittanzahl übergeben, verwende Standardwert von 1.000.000 ...\n");
	else {
		numSteps = strtoul(argv[1],NULL,0);
		printf("Anzahl an Integrationsschritten: %lu\n",numSteps);
	}
	
	// Zerfallsobjekt erstellen
	ThreeBodyMassiveDecay decay(125.1,4.2,4.2,0.,matrixElement,NULL);
	//ThreeBodyMasslessDecay decay(125.1,matrixElement,NULL);
	
	// Zerfallsbreite berechnen
	double result,error;
	decay.getDecayWidth(&result,&error,numSteps);
	
	// Ergebnis + Fehler ausgeben
	printf("decay width: %f +/- %f\n",result,error);

	return 0;
}
