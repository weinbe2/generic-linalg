
#ifndef MG_COMPLEX
#define MG_COMPLEX

#include <complex>
using std::complex;

#include "mg.h"
#include "lattice.h"
#include "verbosity.h"
#include "coarse_stencil.h"
#include "operators.h" // needed for 'dagger' functions, since different ops have different daggered versions. 
#include "generic_inverters.h"

// For multilevel: smooth then down, or recursive Krylov?
enum mg_multilevel_type
{
    MLEVEL_SMOOTH = 0,
    MLEVEL_RECURSIVE = 1
};


// Apply the current coarse operator. 
void coarse_square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

// Apply the current fine operator.
void fine_square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

// Apply the hermitian conjugate of the current coarse operator. 
void coarse_square_staggered_dagger(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

// Apply the hermitian conjugate of the current fine operator.
void fine_square_staggered_dagger(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

// Apply the normal coarse equation (uses the two above defined coarse opts).
void coarse_square_staggered_normal(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

// Apply the normal fine equation (uses the two above defined fine opts).
void fine_square_staggered_normal(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

// Useful mg functions.
struct mg_operator_struct_complex;

// Normalize the projectors based on the coarse blocksize. 
void block_normalize(mg_operator_struct_complex* mgstruct);

// Orthonormalize the projectors based on the coarse blocksize.
void block_orthonormalize(mg_operator_struct_complex* mgstruct);

// Prolong from the coarse lattice to the fine lattice. 
void prolong(complex<double>* x_fine, complex<double>* x_coarse, mg_operator_struct_complex* mgstruct);

// Restrict a fine vector to a coarse vector using the info in mgstruct.
void restrict(complex<double>* x_coarse, complex<double>* x_fine, mg_operator_struct_complex* mgstruct);

// Push the operator struct down a level, updating curr_* variables.
void level_down(mg_operator_struct_complex* mgstruct);

// Pull the operator struct up a level, updating curr_* variables.
void level_up(mg_operator_struct_complex* mgstruct);

// Generate a coarse stencil from a fine stencil.
// ignore_shifts = true -> set to zero before building stencil. Otherwise, build shifts directly into new stencil.
void generate_coarse_from_fine_stencil(stencil_2d* stenc_coarse, stencil_2d* stenc_fine, mg_operator_struct_complex* mgstruct, bool ignore_shifts);
    
// Assumes ignore_shifts = false.
void generate_coarse_from_fine_stencil(stencil_2d* stenc_coarse, stencil_2d* stenc_fine, mg_operator_struct_complex* mgstruct);

// Dslash tracker
struct dslash_tracker
{
    int n_refine;
    int* krylov; // dslashes from the Krylov
    int* presmooth; // dslashes from the presmooth
    int* postsmooth; // dslashes from the postsmooth
    int* residual; // dslashes from the residual equation in multigrid.
    int* nullvectors; // dslashes from null vector generation/
    
    dslash_tracker(int refine) : n_refine(refine)
    {
        krylov = new int[refine+1];
        presmooth = new int[refine+1];
        postsmooth = new int[refine+1];
        residual = new int[refine+1];
        nullvectors = new int[refine+1];
        
        for (int i = 0; i < refine+1; i++)
        {
            krylov[i] = presmooth[i] = postsmooth[i] = residual[i] = nullvectors[i] = 0;
        }
    }
    
    ~dslash_tracker()
    {
        delete[] krylov;
        delete[] presmooth;
        delete[] postsmooth;
        delete[] residual;
        delete[] nullvectors; 
    }
            
};

// Multigrid operator struct.
struct mg_operator_struct_complex
{
    int x_fine; // Fine x dimension.
    int y_fine; // Fine y dimension. 
    int n_refine; // How many refinements? 1 = two level, 2 = three level, etc. 
    
    int* blocksize_x; // How much to block in x direction.
    int* blocksize_y; // How much to block in y direction. 
    unsigned int Nc; // What's Nc on the top level? Square laplace only.
    
    Lattice** latt; // Array of pointers to lattice classes for each level.
    
    stencil_2d** stencils; // Array of pointers to stencils for each level.
    
    bool have_dagger_stencil; // Do we have dagger stencils?
    stencil_2d** dagger_stencils; // Array of pointers to stencils of daggered operator for each level.
                                  // Only created if needed (eg, for normal smoother).
    
    int n_vector; // Number of vectors. 
    complex<double>*** null_vectors; // Holds the null vectors. First index level, second n_vector, third size.
    void (*matrix_vector)(complex<double>*, complex<double>*, void*);
    void (*matrix_vector_dagger)(complex<double>*, complex<double>*, void*);
    void* matrix_extra_data; 
    
    // Track current state.
    int curr_level;
    
    int curr_dof_fine; // degrees of freedom per site. 
    int curr_x_fine;
    int curr_y_fine; 
    int curr_fine_size; 
    
    int curr_dof_coarse; 
    int curr_x_coarse;
    int curr_y_coarse; 
    int curr_coarse_size; 
    
    // Track dslashes
    dslash_tracker* dslash_count; 
};

// Preconditioning struct
struct mg_precond_struct_complex
{
    // What inner smoother?
    minv_inverter in_smooth_type; // MR or GCR
    
    // Set the relaxation parameter (MR only.)
    double omega_smooth;
    
    // How many pre-smooth steps?
    int* n_pre_smooth;
    
    // How many post-smooth steps?
    int* n_post_smooth;
    
    // Are we solving with the coarse normal equations?
    bool normal_eqn_mg; 
    
    // Are we smoothing with the normal equations?
    bool normal_eqn_smooth; 
    
    // How do we do recursive MG?
    mg_multilevel_type mlevel_type; // SMOOTH, RECURSIVE
    
    // What inner solver should we use?
    inner_solver in_solve_type; // MR, CG, CR, or GCR
    int n_max; // Max steps for inner solver?
    int n_restart; // restart frequency (relevant for CG, CR, GCR).
    double* rel_res; // Rel_res for inner solver?
    
    // What's the mg_info?
    mg_operator_struct_complex* mgstruct; 
    
    // What's matrix function are we dealing with?
    // This is the fine function.
    void (*fine_matrix_vector)(complex<double>*, complex<double>*, void*);
    
    // This is the coarse function. 
    void (*coarse_matrix_vector)(complex<double>*, complex<double>*, void*);
    
    // If we need the dagger (thus far only for 'normal_eqn_smooth'), what are the dagger functions?
    // This is the fine function.
    void (*fine_matrix_vector_dagger)(complex<double>*, complex<double>*, void*);
    
    // This is the coarse function. 
    void (*coarse_matrix_vector_dagger)(complex<double>*, complex<double>*, void*);
    
    // If we need the normal operator, what is it?
    void (*fine_matrix_vector_normal)(complex<double>*, complex<double>*, void*);
    
    // This is the coarse function. 
    void (*coarse_matrix_vector_normal)(complex<double>*, complex<double>*, void*);
    
    
    void* matrix_extra_data; 
};

// MG preconditioner!
void mg_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data, inversion_verbose_struct* verb = 0);

#endif // MG_REAL
