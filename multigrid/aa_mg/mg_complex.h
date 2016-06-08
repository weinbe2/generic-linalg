
#include "mg.h"

#ifndef MG_COMPLEX
#define MG_COMPLEX

// General multigrid projector function!
void coarse_square_laplace(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

// General multigrid projector function!
void coarse_square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data); 

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

// Multigrid operator struct.
struct mg_operator_struct_complex
{
    int x_fine; // Fine x dimension.
    int y_fine; // Fine y dimension. 
    int blocksize_x; // How much to block in x direction.
    int blocksize_y; // How much to block in y direction. 
    int n_vector; // Number of vectors. 
    complex<double>** projectors; // Holds the projectors. First index n_vector, second size.
    void (*matrix_vector)(complex<double>*, complex<double>*, void*);
    void* matrix_extra_data; 
};

// Preconditioning struct
struct mg_precond_struct_complex
{
    // What inner smoother?
    inner_solver in_smooth_type; // MINRES or GCR
    
    // How many MinRes pre-smooth steps?
    int n_pre_smooth;
    
    // How many MinRes post-smooth steps?
    int n_post_smooth;
    
    // What inner solver should we use?
    inner_solver in_solve_type; // MINRES, CG, or GCR
    int n_step; // Max steps for inner solver?
    double rel_res; // Rel_res for inner solver?
    
    // What's the mg_info?
    mg_operator_struct_complex* mgstruct; 
    
    // What's matrix function are we dealing with?
    void (*matrix_vector)(complex<double>*, complex<double>*, void*);
    void* matrix_extra_data; 
};

// MG preconditioner!
void mg_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data);

#endif // MG_REAL
