
#include "mg.h"
#include "coordinate.h"

#ifndef MG_COMPLEX
#define MG_COMPLEX

// General multigrid projector function!
void coarse_square_laplace(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

// Apply the current coarse operator. 
void coarse_square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

// Apply the current fine operator.
void fine_square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

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

// Multigrid operator struct.
struct mg_operator_struct_complex
{
    int x_fine; // Fine x dimension.
    int y_fine; // Fine y dimension. 
    int n_refine; // How many refinements? 1 = two level, 2 = three level, etc. 
    
    int* blocksize_x; // How much to block in x direction.
    int* blocksize_y; // How much to block in y direction. 
    int Nc; // What's Nc on the top level? Square laplace only.
    
    int eo; // 0 if no even-odd blocking, 1 if we use even-odd blocking. 
    int n_vector; // Number of vectors. 
    complex<double>*** null_vectors; // Holds the null vectors. First index level, second n_vector, third size.
    void (*matrix_vector)(complex<double>*, complex<double>*, void*);
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
};

// Preconditioning struct
struct mg_precond_struct_complex
{
    // What inner smoother?
    inner_solver in_smooth_type; // MR or GCR
    
    // Set the relaxation parameter (MR only.)
    double omega_smooth;
    
    // How many pre-smooth steps?
    int n_pre_smooth;
    
    // How many post-smooth steps?
    int n_post_smooth;
    
    // What inner solver should we use?
    inner_solver in_solve_type; // MR, CG, CR, or GCR
    int n_max; // Max steps for inner solver?
    int n_restart; // restart frequency (relevant for CG, CR, GCR).
    double rel_res; // Rel_res for inner solver?
    
    // What's the mg_info?
    mg_operator_struct_complex* mgstruct; 
    
    // What's matrix function are we dealing with?
    // This is the fine function.
    void (*fine_matrix_vector)(complex<double>*, complex<double>*, void*);
    
    // This is the coarse function. 
    void (*coarse_matrix_vector)(complex<double>*, complex<double>*, void*);
    void* matrix_extra_data; 
};

// MG preconditioner!
void mg_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data, inversion_verbose_struct* verb = 0);

#endif // MG_REAL
