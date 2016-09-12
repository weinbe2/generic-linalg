// A header for routines related to null vector generation.

#ifndef MG_NULL_GEN
#define MG_NULL_GEN

#include "mg_complex.h"

// How should we generate null vectors?
enum mg_null_gen_type
{
    NULL_GCR = 0,                    // Generate null vectors with GCR
    NULL_BICGSTAB = 1,               // Generate null vectors with BiCGStab
    NULL_CG = 2,                    // Generate null vectors with CG
    NULL_MINRES = 3,                // Generate null vectors with MinRes
    NULL_ARPACK = 4,                // Generate null vectors as low eigenvectors with Arpack. 
    NULL_BICGSTAB_L = 5,            // Generate null vectors with BiCGStab-l
};

// What structure are we preserving when we block? None, E/O, Corners?
enum blocking_strategy
{
    BLOCK_NONE = 0,                   // Block fully.
    BLOCK_EO = 1,                     // Even/odd
    BLOCK_CORNER = 2,                 // Corners
    BLOCK_TOPO = 3                    // Chirality defined by taste singlet
};

// Function to partition null vectors on the top level.
// Supports partitioning one null vector at a time.
void null_partition_staggered(mg_operator_struct_complex* mgstruct, int num_null_vec, blocking_strategy bstrat, Lattice* Lat);

// Function to partition null vectors on the coarse level. Should update "Lattice" object for all levels.
void null_partition_coarse(mg_operator_struct_complex* mgstruct, int num_null_vec, blocking_strategy bstrat);


#endif // MG_NULL_GEN