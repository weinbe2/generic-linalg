// Tue May  3 16:49:07 EDT 2016
// Evan S Weinberg
// Include file for eigenvalue routines. 

#ifndef ESW_EIGENVALUE
#define ESW_EIGENVALUE

// Uncomment to print warnings. 
//#define VERBOSE_WARN

#include <complex>

using std::complex; 

// Struct that contains information that
// all matrix functions return.
// Return for various matrix functions.
struct eigenvalue_info
{
  double relative_diff; // relative difference on last step.
  int iter; // number of iterations.
  bool success; // did we reach residual?
  std::string name; // name of algorithm.
};

// Performs power iterations to find the largest eigenvalue. 
eigenvalue_info eig_vector_poweriter(double* eig, double* phi0, int size, int max_iter, double relres, void (*matrix_vector)(double*,double*,void*), void* extra_info);

// To add: general power iterator? I really should in a Ritz algorithm, but to start, something that
// just explicitly re-orthogonalizes. 

#endif // define ESW_EIGENVALUE
