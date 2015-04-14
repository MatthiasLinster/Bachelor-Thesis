#ifndef FUNCTION_TYPES_H
#define FUNCTION_TYPES_H

#include "FourVector.h"

// Struktur einer Monte-Carlo-Funktion double f(double* x,size_t numDimensions,void* params), vgl. GSL
typedef double (*int_function_ptr)(double*,size_t,void*);

// Art einer vom Phasenraum abhängigen Funktion double f(FourVector* 4-Impulse,void* Parameter)
typedef double (*phase_space_function_ptr)(FourVector*,void*);

// Funktionsstruktur für Matrixelemente: double f(FourVector*,void*), wobei sich FourVector* auf die 4-Impulse der beteiligten Teilchen und void* auf die Parameter bezieht
typedef double (*matrix_element_ptr)(FourVector*,void*);

#endif
