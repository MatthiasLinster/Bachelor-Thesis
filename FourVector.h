#ifndef FOUR_VECTOR_H
#define FOUR_VECTOR_H

/**
 * Struktur für einen Vierer-Vektor
 * @author Matthias Linster
 * @date 03.09.2014
 */
struct FourVector {
	// Komponenten
	double x[4];
	
	/**
	 * Default-Konstruktur, erstellt den Nullvektor
	 */
	FourVector() {
		x[0] = 0.;
		x[1] = 0.;
		x[2] = 0.;
		x[3] = 0.;
	}
	
	/**
	 * Erstellt einen Vierervektor mit den entsprechenden Komponenten
	 * @param x0 0. Komponente
	 * @param x1 1. Komponente
	 * @param x2 2. Komponente
	 * @param x3 3. Komponente
	 */
	FourVector(double x0,double x1,double x2,double x3) {
		x[0] = x0;
		x[1] = x1;
		x[2] = x2;
		x[3] = x3;
	}
	
	/**
	 * +-Operator zur Addition zweier Vierer-Vektoren
	 */
	inline FourVector operator+(FourVector v) const {
		return FourVector(x[0]+v.x[0],x[1]+v.x[1],x[2]+v.x[2],x[3]+v.x[3]);
	}
	
	/**
	 * *=-Operator als Multiplikation eines Vierer-Vektors mit einer reellen Zahl
	 */
	inline FourVector& operator*=(double a) {
		x[0] *= a;
		x[1] *= a;
		x[2] *= a;
		x[3] *= a;
		
		return *this;
	}
	
	/**
	 * *-Operator für das Minkowski-Skalarprodukt
	 */
	inline double operator*(FourVector v) const {
		return x[0]*v.x[0]-x[1]*v.x[1]-x[2]*v.x[2]-x[3]*v.x[3];
	}
};

#endif
