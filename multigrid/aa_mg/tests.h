// Header file for various test routines. 

#include <complex>

#include "mg.h"
#include "mg_complex.h"
#include "generic_vector.h"
#include "coarse_stencil.h"
#include "null_gen.h"
#include "operators.h"

// Test eigenvalue overlaps. 
void test_eigenvalue_overlap(mg_operator_struct_complex* mgstruct, staggered_u1_op* stagif, op_type opt, int set_eigen, int set_cv);

// Test constructing the coarse operator, comparing to using prolong/restrict of fine.
void test_stencil_construct(mg_operator_struct_complex* mgstruct, int level, int stencil_size);

// Test constructing the coarse operator, comparing applying the full coarse stencil to the piece-by-piece stencil.
void test_stencil_piece(mg_operator_struct_complex* mgstruct, int level, int stencil_size);

// Test efficient construct of coarse operator, compared to brute force construct.
void test_stencil_efficient_construct(mg_operator_struct_complex* mgstruct, int level, int stencil_size); 