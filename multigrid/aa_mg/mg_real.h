
#include "mg.h"

#ifndef MG_REAL
#define MG_REAL

// General multigrid projector function!
void coarse_square_op(double* lhs, double* rhs, void* extra_data); 

// Useful mg functions.
struct mg_operator_struct_real;

// Normalize the projectors based on the coarse blocksize. 
void block_normalize(mg_operator_struct_real* mgstruct);

// Prolong from the coarse lattice to the fine lattice. 
void prolong(double* x_fine, double* x_coarse, mg_operator_struct_real* mgstruct);

// Restrict a fine vector to a coarse vector using the info in mgstruct.
void restrict(double* x_coarse, double* x_fine, mg_operator_struct_real* mgstruct);

// Multigrid operator struct.
struct mg_operator_struct_real
{
    int x_fine; // Fine x dimension.
    int y_fine; // Fine y dimension. 
    int blocksize_x; // How much to block in x direction.
    int blocksize_y; // How much to block in y direction. 
    int n_vector; // Number of vectors. 
    double** projectors; // Holds the projectors. First index n_vector, second size.
    void (*matrix_vector)(double*, double*, void*);
    void* matrix_extra_data; 
};

// Preconditioning struct
struct mg_precond_struct_real
{
    // How many MinRes pre-smooth steps?
    int n_pre_smooth;
    
    // How many MinRes post-smooth steps?
    int n_post_smooth;
    
    // What inner solver should we use?
    inner_solver in_solve_type; // MINRES, CG, or GCR
    int n_step; // Max steps for inner solver?
    double rel_res; // Rel_res for inner solver?
    
    // What's the mg_info?
    mg_operator_struct_real* mgstruct; 
    
    // What's matrix function are we dealing with?
    void (*matrix_vector)(double*, double*, void*);
    void* matrix_extra_data; 
};

// MG preconditioner!
void mg_preconditioner(double* lhs, double* rhs, int size, void* extra_data);

#endif // MG_REAL
