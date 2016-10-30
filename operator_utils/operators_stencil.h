// Build and return the stencil (look in the stencil directory) for a U(1) operator.
// Depends on the lattice class and the staggered_u1_op struct.

// Also defines some stencil functions unique to our operators, for ex for even/odd preconditioning
// for the staggered operator. 

#ifndef U1_STENCIL
#define U1_STENCIL

#include "operators.h"
#include "lattice.h"
#include "coarse_stencil.h"


// Functions for generating a stencil.

// Given an allocated but otherwise empty stencil, appropriately fill it with the 2d staggered links.
void get_square_staggered_u1_stencil(stencil_2d* stenc, staggered_u1_op* stagif);

// Given an allocated but otherwise empty stencil, appropriately fill it with the 2d g5_staggered links.
void get_square_staggered_gamma5_u1_stencil(stencil_2d* stenc, staggered_u1_op* stagif);

// Given an allocated but otherwise empty stencil, appropriately fill it with staggered dagger links.
void get_square_staggered_dagger_u1_stencil(stencil_2d* stenc, staggered_u1_op* stagif);


// Functions for applying a specialized stencil.

// Prepare an even rhs for an even/odd preconditioned solve.
// Takes in rhs_orig, returns rhs_e. 
void apply_square_staggered_eoprec_prepare_stencil(complex<double>* rhs_e, complex<double>* rhs_orig, stencil_2d* stenc);

// Square staggered 2d operator w/ u1 function, m^2 - D_{eo} D_{oe} [zeroes odd explicitly]
void apply_square_staggered_m2mdeodoe_stencil(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Reconstruct the full solution for an even/odd preconditioned solve.
// Takes in lhs_e, rhs_o, returns lhs_full (copying over the even part from lhs_e)
void apply_square_staggered_eoprec_reconstruct_stencil(complex<double>* lhs_full, complex<double>* lhs_e, complex<double>* rhs_o, stencil_2d* stenc);

#endif // U1_STENCIL
