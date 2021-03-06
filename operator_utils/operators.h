#ifndef U1_OPERATORS
#define U1_OPERATORS

#include <complex>
using namespace std;

struct staggered_u1_op
{
  complex<double> *lattice;
  double mass;
  int x_fine;
  int y_fine; 
  int Nc; // only relevant for square laplace. 
  double wilson_coeff; // only relevant for staggered + 2 link laplace. 
};

enum op_type
{
  STAGGERED = 0,
  LAPLACE = 1,
  LAPLACE_NC2 = 2,
  G5_STAGGERED = 3,
  STAGGERED_NORMAL = 4,
  STAGGERED_INDEX = 5,
};

// Get the stencil size (maximum number of hops) as a function of the operator.
int get_stencil_size(op_type opt);

// Square lattice.
// Kinetic term for a 2D laplace w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
void square_laplace(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square laplacian w/ u1. 
void square_laplace_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square lattice.
// Kinetic term for a 2D staggered w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" doesn't include anything.
// Apsi[x][y] = m psi[x,y] - U[y][x,x+1] 
void square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square lattice.
// Kinetic term for a 2D staggered w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" is a cast to a complex gauge_field[N*N*2], 
//    loaded by the function read_lattice_u1. 
// Apsi[x][y] = m psi[x,y] - U[y][x,x+1] 
void square_staggered_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// \gamma_5
void gamma_5(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square \gamma_5 staggered 2d operator w/out u1 function.
void square_staggered_gamma5(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square \gamma_5 staggered 2d operator w/ u1 function.
void square_staggered_gamma5_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square staggered dagger 2d operator w/ u1 function.
void square_staggered_dagger_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Staggered normal equations.
void square_staggered_normal_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square staggered 2d operator w/ u1 function, D_{eo} only.
void square_staggered_deo_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square staggered 2d operator w/ u1 function, D_{oe} only.
void square_staggered_doe_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Prepare an even rhs for an even/odd preconditioned solve.
// Takes in rhs_orig, returns rhs_e. 
void square_staggered_eoprec_prepare(complex<double>* rhs_e, complex<double>* rhs_orig, void* extra_data);

// Square staggered 2d operator w/ u1 function, m^2 - D_{eo} D_{oe} [zeroes odd explicitly]
void square_staggered_m2mdeodoe_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Reconstruct the full solution for an even/odd preconditioned solve.
// Takes in lhs_e, rhs_o, returns lhs_full (copying over the even part from lhs_e)
void square_staggered_eoprec_reconstruct(complex<double>* lhs_full, complex<double>* lhs_e, complex<double>* rhs_o, void* extra_data);


// Square lattice.
// Kinetic term for a 2D staggered w/ period bc, plus 2-link laplace.
// Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" is a cast to a complex gauge_field[N*N*2], 
//    loaded by the function read_lattice_u1. 
// Apsi[x][y] = m psi[x,y] - U[y][x,x+1] 
void square_staggered_2linklaplace_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Operators for symmetric shifts.
void staggered_symmshift_x(complex<double>* lhs, complex<double>* rhs, void* extra_data);
void staggered_symmshift_y(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

// Staggered index operator. (See arXiv 1410.5733, 1203.2560)
void staggered_index_operator(complex<double>* lhs, complex<double>* rhs, void* extra_data);




#endif

