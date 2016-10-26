// A header for routines related to null vector generation.

#ifndef MG_NULL_GEN
#define MG_NULL_GEN

#include <vector>

#include "generic_inverters.h"
#include "inverter_struct.h"
#include "mg_complex.h"
#include "lattice.h"

// What structure are we preserving when we block? None, E/O, Corners?
enum blocking_strategy
{
    BLOCK_NONE = 0,                   // Block fully.
    BLOCK_EO = 1,                     // Even/odd
    BLOCK_CORNER = 2,                 // Corners
    BLOCK_TOPO = 3                    // Chirality defined by taste singlet
};

struct null_vector_params
{
    op_type opt_null;
    int n_null_vector; 
    minv_inverter null_gen; 
    std::vector<double> null_precisions; 
    std::vector<int> null_max_iters; 
    bool null_restart; 
    int null_restart_freq; 
    int null_bicgstab_l; 
    double null_mass; 
    int null_partitions; 
    blocking_strategy bstrat;
    bool do_global_ortho_conj; 
    bool do_ortho_eo; 
};

// Function to partition null vectors on the top level.
// Supports partitioning one null vector at a time.
void null_partition_staggered(mg_operator_struct_complex* mgstruct, int num_null_vec, blocking_strategy bstrat, Lattice* Lat);

// Function to partition null vectors on the coarse level. Should update "Lattice" object for all levels.
void null_partition_coarse(mg_operator_struct_complex* mgstruct, int num_null_vec, blocking_strategy bstrat);

// Function to generate free field null vectors (up to a gauge transform).
void null_generate_random_smooth(mg_operator_struct_complex* mgstruct, null_vector_params* nvec_params, bool do_gauge_transform, std::complex<double>* gauge_trans);

// Function to generate null vectors by throwing a random source and smoothing.
void null_generate_random_smooth(mg_operator_struct_complex* mgstruct, null_vector_params* nvec_params, inversion_verbose_struct* verb, std::mt19937* generator);

#endif // MG_NULL_GEN