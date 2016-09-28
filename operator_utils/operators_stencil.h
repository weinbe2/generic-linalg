// Build and return the stencil (look in the stencil directory) for a U(1) operator.
// Depends on the lattice class and the staggered_u1_op struct.

#ifndef U1_STENCIL
#define U1_STENCIL

#include "operators.h"
#include "lattice.h"
#include "coarse_stencil.h"

// Given an allocated but otherwise empty stencil, appropriately fill it with the 2d staggered links.
void get_square_staggered_u1_stencil(stencil_2d* stenc, staggered_u1_op* stagif);

// Given an allocated but otherwise empty stencil, appropriately fill it with the 2d g5_staggered links.
void get_square_staggered_gamma5_u1_stencil(stencil_2d* stenc, staggered_u1_op* stagif);

// Given an allocated but otherwise empty stencil, appropriately fill it with staggered dagger links.
void get_square_staggered_dagger_u1_stencil(stencil_2d* stenc, staggered_u1_op* stagif);

#endif // U1_STENCIL
