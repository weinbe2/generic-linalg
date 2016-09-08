// Header file for various test routines. 

#include <complex>

#include "mg.h"
#include "mg_complex.h"
#include "generic_vector.h"
#include "coarse_stencil.h"
#include "null_gen.h"

// Test constructing the coarse operator, comparing to using prolong/restrict of fine.
void test_stencil_construct(mg_operator_struct_complex* mgstruct, int level, int stencil_size);