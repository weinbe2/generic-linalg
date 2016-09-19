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

// Performs gaussian elimination against multirhs.
int gaussian_elimination_multi_rhs(double** x, double** b, double** matrix, int n_rhs, int size);
int gaussian_elimination_multi_rhs(complex<double>** x, complex<double>** b, complex<double>** matrix, int n_rhs, int size);

// Compute a matrix inverse using gaussian elimination (uses multirhs under the hood). minv and matrix can be the same pointer.
int gaussian_elimination_matrix_inverse(double** minv, double** matrix, int size);
int gaussian_elimination_matrix_inverse(complex<double>** minv, complex<double>** matrix, int size);

#endif