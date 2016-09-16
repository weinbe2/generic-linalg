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
};

enum op_type
{
    STAGGERED = 0,
    LAPLACE = 1,
    LAPLACE_NC2 = 2,
    G5_STAGGERED = 3,
    STAGGERED_NORMAL = 4,
    STAGGERED_INDEX = 5
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

// Square \gamma_5 staggered 2d operator w/ u1 function.
void square_staggered_gamma5_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Staggered normal equations.
void square_staggered_normal_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Operators for symmetric shifts.
void staggered_symmshift_x(complex<double>* lhs, complex<double>* rhs, void* extra_data);
void staggered_symmshift_y(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

// Staggered index operator. (See arXiv 1410.5733, 1203.2560)
void staggered_index_operator(complex<double>* lhs, complex<double>* rhs, void* extra_data);




#endif

