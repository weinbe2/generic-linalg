// Wed Aug 24 15:46:50 EDT 2016
// Evan S Weinberg
// Include file for Gaussian Elimination

#ifndef ESW_INVERTER_GELIM
#define ESW_INVERTER_GELIM

#include <string>
#include <complex>
using std::complex;

#include "inverter_struct.h"
#include "verbosity.h"

// Performs gaussian elimination.
int gaussian_elimination(double* x, double* b, double** matrix, int size);
int gaussian_elimination(complex<double>* x, complex<double>* b, complex<double>** matrix, int size);

#endif