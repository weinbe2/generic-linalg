
#include <iostream>
#include <iomanip> // to set output precision.
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

using namespace std;

#include "generic_inverters.h"
#include "generic_inverters_precond.h"
#include "generic_vector.h"
#include "mg.h"
#include "mg_complex.h"

// General multigrid projector function!
void coarse_square_laplace(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{   
    // Grab the mg_precond_struct.
    mg_operator_struct_complex* mgstruct = (mg_operator_struct_complex*)extra_data; 
    
    // lhs and rhs are of size coarse_size. mgstruct.matrix_vector expects
    // fine_size. 
    int x_fine = mgstruct->x_fine;
    int y_fine = mgstruct->y_fine;
    int fine_size = x_fine*y_fine;
    int x_coarse = x_fine/mgstruct->blocksize_x[0]; // how many coarse sites are there in the x dir?
    int y_coarse = y_fine/mgstruct->blocksize_y[0]; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 
    int coarse_length = coarse_size*mgstruct->n_vector;
    
    // Okay... how the hell are we going to do this. 
    complex<double>* Px; // Holds prolonged current solution.
    complex<double>* APx; // Holds A times prolonged current solution.
    
    Px = new complex<double>[fine_size];
    APx = new complex<double>[fine_size];
    
    zero<double>(Px, fine_size); zero<double>(APx, fine_size); 
    
    // Prolong. 
    prolong(Px, rhs, mgstruct);
    
    // Apply the original matrix.
    (*mgstruct->matrix_vector)(APx, Px, mgstruct->matrix_extra_data);
    
    // Restrict. 
    zero<double>(lhs, coarse_length);
    restrict(lhs, APx, mgstruct); 
    
    delete[] Px;
    delete[] APx; 
    
}

// General multigrid projector function!
// Apply the current coarse operator.
void coarse_square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Grab the mg_precond_struct.
    mg_operator_struct_complex* mgstruct = (mg_operator_struct_complex*)extra_data; 
    
    // We need to prolong all the way up to the top.
    int curr_level = mgstruct->curr_level+1;
    
    // Prepare space for everything.
    complex<double>** Px; // holds prolonged current solutions.
    Px = new complex<double>*[curr_level+1];
    
    // Copy over.
    Px[curr_level] = new complex<double>[mgstruct->curr_coarse_size];
    copy<double>(Px[curr_level], rhs, mgstruct->curr_coarse_size);
    
    for (int i = curr_level-1; i >=0; i--)
    {
        Px[i] = new complex<double>[mgstruct->curr_fine_size];
        zero<double>(Px[i], mgstruct->curr_fine_size);
        
        // Prolong.
        prolong(Px[i], Px[i+1], mgstruct);
        
        // Level up! 
        if (i != 0)
        {
            level_up(mgstruct); 
        }
    }
    
    // Apply A.
    complex<double>** APx; // holds restricted current solutions.
    APx = new complex<double>*[curr_level+1];
    APx[0] = new complex<double>[mgstruct->curr_fine_size];
    
    zero<double>(APx[0], mgstruct->curr_fine_size);
    (*mgstruct->matrix_vector)(APx[0], Px[0], mgstruct->matrix_extra_data);
    
    // Bring it down.
    for (int i = 1; i <= curr_level; i++)
    {
        APx[i] = new complex<double>[mgstruct->curr_coarse_size];
        zero<double>(APx[i], mgstruct->curr_coarse_size);
        
        // Restrict.
        restrict(APx[i], APx[i-1], mgstruct);
        
        // Level down!
        if (i != curr_level)
        {
            level_down(mgstruct);
        }
    }
    
    copy<double>(lhs, APx[curr_level], mgstruct->curr_coarse_size); 
    for (int i = 0; i < curr_level+1; i++)
    {
        delete[] APx[i];
        delete[] Px[i];
    }
    
    
}

// General multigrid projector function!
// Apply the current fine operator.
void fine_square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Grab the mg_precond_struct.
    mg_operator_struct_complex* mgstruct = (mg_operator_struct_complex*)extra_data; 
    
    if (mgstruct->curr_level == 0)
    {
        zero<double>(lhs, mgstruct->curr_fine_size);
        (*mgstruct->matrix_vector)(lhs, rhs, mgstruct->matrix_extra_data);
    }
    else
    {
        level_up(mgstruct);
        coarse_square_staggered(lhs, rhs, extra_data);
        level_down(mgstruct);
    }
    
    
}


// Properly normalize the P vectors.
void block_normalize(mg_operator_struct_complex* mgstruct)
{
    int i;
    
    // Grab the current null vectors.
    complex<double>** null_vectors = mgstruct->null_vectors[mgstruct->curr_level];
    
    // Hold current sites.
    int curr_x_coarse, curr_y_coarse, curr_dof_coarse; 
    
    // Loop over every coarse site and generalized color index.
    for (i = 0; i < mgstruct->curr_coarse_size; i++)
    {
        int tmp = i;
        // Decompose into color, x, y.
        curr_dof_coarse = tmp % mgstruct->curr_dof_coarse;
        tmp = (tmp - curr_dof_coarse)/mgstruct->curr_dof_coarse;
        curr_x_coarse = tmp % mgstruct->curr_x_coarse;
        tmp = (tmp - curr_x_coarse)/mgstruct->curr_x_coarse;
        curr_y_coarse = tmp;
        
        // Convert to 'unit' corner of fine lattice.
        int curr_dof_fine = 0;
        int curr_x_fine = curr_x_coarse*mgstruct->blocksize_x[mgstruct->curr_level];
        int curr_y_fine = curr_y_coarse*mgstruct->blocksize_y[mgstruct->curr_level];
        
        double norm = 0.0;
        
        // Loop over all relevant fine sites.
        for (int y = curr_y_fine; y < curr_y_fine+mgstruct->blocksize_y[mgstruct->curr_level]; y++)
        {
            for (int x = curr_x_fine; x < curr_x_fine+mgstruct->blocksize_x[mgstruct->curr_level]; x++)
            {
                for (int dof = curr_dof_fine; dof < mgstruct->curr_dof_fine; dof++)
                {
                    // Construct fine site.
                    int fine_site = y*mgstruct->curr_x_fine*mgstruct->curr_dof_fine + x*mgstruct->curr_dof_fine + dof;
                    
                    // Build norm.
                    norm += real(conj(null_vectors[curr_dof_coarse][fine_site])*null_vectors[curr_dof_coarse][fine_site]);
                }
            }
        }
        
        norm = sqrt(norm);
        
        // Loop over all relevant fine sites.
        for (int y = curr_y_fine; y < curr_y_fine+mgstruct->blocksize_y[mgstruct->curr_level]; y++)
        {
            for (int x = curr_x_fine; x < curr_x_fine+mgstruct->blocksize_x[mgstruct->curr_level]; x++)
            {
                for (int dof = curr_dof_fine; dof < mgstruct->curr_dof_fine; dof++)
                {
                    // Construct fine site.
                    int fine_site = y*mgstruct->curr_x_fine*mgstruct->curr_dof_fine + x*mgstruct->curr_dof_fine + dof;
                    
                    // Normalize.
                    null_vectors[curr_dof_coarse][fine_site] /= norm;
                }
            }
        }
    }
        
    
    
}

// Properly orthonormalize the P vectors.
void block_orthonormalize(mg_operator_struct_complex* mgstruct)
{
    int  m, i;

    int x_coarse = mgstruct->curr_x_coarse; // how many coarse sites are there in the x dir?
    int y_coarse = mgstruct->curr_y_coarse; // how many coarse sites are there in the y dir?
    
    // Grab the current null vectors.
    complex<double>** null_vectors = mgstruct->null_vectors[mgstruct->curr_level];
    
    // Hold current sites.
    int curr_x_coarse, curr_y_coarse, curr_dof_coarse; 
    
    double norm;
    complex<double> dot_prod;
    
    // Loop over every coarse site separately from each color index.
    for (i = 0; i < x_coarse*y_coarse; i++)
    {
        // Get the current coarse sites.
        int tmp = i;
        curr_x_coarse = tmp % mgstruct->curr_x_coarse;
        tmp = (tmp - curr_x_coarse)/mgstruct->curr_x_coarse;
        curr_y_coarse = tmp;
        
        // Convert to 'unit' corner of fine lattice.
        int curr_dof_fine = 0;
        int curr_x_fine = curr_x_coarse*mgstruct->blocksize_x[mgstruct->curr_level];
        int curr_y_fine = curr_y_coarse*mgstruct->blocksize_y[mgstruct->curr_level];
        
        for (curr_dof_coarse = 1; curr_dof_coarse < mgstruct->curr_dof_coarse; curr_dof_coarse++)
        {
            norm = 0.0;
            // First, normalize the previous vector.
            
            for (int y = curr_y_fine; y < curr_y_fine+mgstruct->blocksize_y[mgstruct->curr_level]; y++)
            {
                for (int x = curr_x_fine; x < curr_x_fine+mgstruct->blocksize_x[mgstruct->curr_level]; x++)
                {
                    for (int dof = curr_dof_fine; dof < mgstruct->curr_dof_fine; dof++)
                    {
                        // Construct fine site.
                        int fine_site = y*mgstruct->curr_x_fine*mgstruct->curr_dof_fine + x*mgstruct->curr_dof_fine + dof;

                        // Build norm.
                        norm += real(conj(null_vectors[curr_dof_coarse-1][fine_site])*null_vectors[curr_dof_coarse-1][fine_site]);
                    }
                }
            }
            
            norm = sqrt(norm);
            
            // Loop over all relevant fine sites.
            for (int y = curr_y_fine; y < curr_y_fine+mgstruct->blocksize_y[mgstruct->curr_level]; y++)
            {
                for (int x = curr_x_fine; x < curr_x_fine+mgstruct->blocksize_x[mgstruct->curr_level]; x++)
                {
                    for (int dof = curr_dof_fine; dof < mgstruct->curr_dof_fine; dof++)
                    {
                        // Construct fine site.
                        int fine_site = y*mgstruct->curr_x_fine*mgstruct->curr_dof_fine + x*mgstruct->curr_dof_fine + dof;

                        // Normalize.
                        null_vectors[curr_dof_coarse-1][fine_site] /= norm;
                    }
                }
            }
            
            // Orthogonalize current vector against all previous vectors.
            for (m = 0; m < curr_dof_coarse; m++)
            {
                dot_prod = 0.0;
                
                // Compute inner product...
                for (int y = curr_y_fine; y < curr_y_fine+mgstruct->blocksize_y[mgstruct->curr_level]; y++)
                {
                    for (int x = curr_x_fine; x < curr_x_fine+mgstruct->blocksize_x[mgstruct->curr_level]; x++)
                    {
                        for (int dof = curr_dof_fine; dof < mgstruct->curr_dof_fine; dof++)
                        {
                            // Construct fine site.
                            int fine_site = y*mgstruct->curr_x_fine*mgstruct->curr_dof_fine + x*mgstruct->curr_dof_fine + dof;

                            // Dot product.
                            dot_prod = conj(null_vectors[m][fine_site])*null_vectors[curr_dof_coarse][fine_site];
                        }
                    }
                }
                
                // Project off.
                for (int y = curr_y_fine; y < curr_y_fine+mgstruct->blocksize_y[mgstruct->curr_level]; y++)
                {
                    for (int x = curr_x_fine; x < curr_x_fine+mgstruct->blocksize_x[mgstruct->curr_level]; x++)
                    {
                        for (int dof = curr_dof_fine; dof < mgstruct->curr_dof_fine; dof++)
                        {
                            // Construct fine site.
                            int fine_site = y*mgstruct->curr_x_fine*mgstruct->curr_dof_fine + x*mgstruct->curr_dof_fine + dof;

                            // Dot product.
                            null_vectors[curr_dof_coarse][fine_site] -= dot_prod*null_vectors[m][fine_site];
                        }
                    }
                }
            }
        }
           
    }
    
    block_normalize(mgstruct); 
}

// Prolong a coarse vector to a fine vector using the info in mgstruct.
void prolong(complex<double>* vec_fine, complex<double>* vec_coarse, mg_operator_struct_complex* mgstruct)
{
    int i;
    
    // Grab the current null vectors.
    complex<double>** null_vectors = mgstruct->null_vectors[mgstruct->curr_level];
    
    zero<double>(vec_fine, mgstruct->curr_fine_size); 
    
    // Hold current sites.
    int curr_x_coarse, curr_y_coarse, curr_dof_coarse; 
    
    // Loop over every coarse site and generalized color index.
    for (i = 0; i < mgstruct->curr_coarse_size; i++)
    {
        int tmp = i;
        // Decompose into color, x, y.
        curr_dof_coarse = tmp % mgstruct->curr_dof_coarse;
        tmp = (tmp - curr_dof_coarse)/mgstruct->curr_dof_coarse;
        curr_x_coarse = tmp % mgstruct->curr_x_coarse;
        tmp = (tmp - curr_x_coarse)/mgstruct->curr_x_coarse;
        curr_y_coarse = tmp;
        
        // Convert to 'unit' corner of fine lattice.
        int curr_dof_fine = 0;
        int curr_x_fine = curr_x_coarse*mgstruct->blocksize_x[mgstruct->curr_level];
        int curr_y_fine = curr_y_coarse*mgstruct->blocksize_y[mgstruct->curr_level];
        
        // Loop over all relevant fine sites.
        for (int y = curr_y_fine; y < curr_y_fine+mgstruct->blocksize_y[mgstruct->curr_level]; y++)
        {
            for (int x = curr_x_fine; x < curr_x_fine+mgstruct->blocksize_x[mgstruct->curr_level]; x++)
            {
                for (int dof = curr_dof_fine; dof < mgstruct->curr_dof_fine; dof++)
                {
                    // Construct fine site.
                    int fine_site = y*mgstruct->curr_x_fine*mgstruct->curr_dof_fine + x*mgstruct->curr_dof_fine + dof;
                    
                    // Update the fine with the coarse. 
                    vec_fine[fine_site] += null_vectors[curr_dof_coarse][fine_site]*vec_coarse[i];
                    
                }
            }
        }
    }
}

// Restrict a fine vector to a coarse vector using the info in mgstruct.
// Hermitian conjugate of prolong. 
void restrict(complex<double>* vec_coarse, complex<double>* vec_fine, mg_operator_struct_complex* mgstruct)
{
    int i;
    
    // Grab the current null vectors.
    complex<double>** null_vectors = mgstruct->null_vectors[mgstruct->curr_level];
    
    zero<double>(vec_coarse, mgstruct->curr_coarse_size);
    
    // Hold current sites.
    int curr_x_coarse, curr_y_coarse, curr_dof_coarse; 
    
    // Loop over every coarse site and generalized color index.
    for (i = 0; i < mgstruct->curr_coarse_size; i++)
    {
        int tmp = i;
        // Decompose into color, x, y.
        curr_dof_coarse = tmp % mgstruct->curr_dof_coarse;
        tmp = (tmp - curr_dof_coarse)/mgstruct->curr_dof_coarse;
        curr_x_coarse = tmp % mgstruct->curr_x_coarse;
        tmp = (tmp - curr_x_coarse)/mgstruct->curr_x_coarse;
        curr_y_coarse = tmp;
        
        // Convert to 'unit' corner of fine lattice.
        int curr_dof_fine = 0;
        int curr_x_fine = curr_x_coarse*mgstruct->blocksize_x[mgstruct->curr_level];
        int curr_y_fine = curr_y_coarse*mgstruct->blocksize_y[mgstruct->curr_level];
        
        // Loop over all relevant fine sites.
        for (int y = curr_y_fine; y < curr_y_fine+mgstruct->blocksize_y[mgstruct->curr_level]; y++)
        {
            for (int x = curr_x_fine; x < curr_x_fine+mgstruct->blocksize_x[mgstruct->curr_level]; x++)
            {
                for (int dof = curr_dof_fine; dof < mgstruct->curr_dof_fine; dof++)
                {
                    // Construct fine site.
                    int fine_site = y*mgstruct->curr_x_fine*mgstruct->curr_dof_fine + x*mgstruct->curr_dof_fine + dof;
                    
                    // Update the coarse with the fine.
                    
                    vec_coarse[i] += conj(null_vectors[curr_dof_coarse][fine_site])*vec_fine[fine_site];
                    
                }
            }
        }
    }
}

// Push the operator struct down a level, updating curr_* variables.
void level_down(mg_operator_struct_complex* mgstruct)
{
    if (mgstruct->curr_level < mgstruct->n_refine-1) // Can't go lower than the number of refinements!
    {
        // Update the number of d.o.f.
        if (mgstruct->curr_level == 0)
        {
            mgstruct->curr_dof_fine *= mgstruct->n_vector/mgstruct->Nc; // Nc for square laplace only. 
        }
        
        // Update the current fine.
        mgstruct->curr_x_fine /= mgstruct->blocksize_x[mgstruct->curr_level];
        mgstruct->curr_y_fine /= mgstruct->blocksize_y[mgstruct->curr_level];
        
        mgstruct->curr_fine_size = mgstruct->curr_y_fine*mgstruct->curr_x_fine*mgstruct->curr_dof_fine;
        
        // Update the current coarse.
        mgstruct->curr_x_coarse /= mgstruct->blocksize_x[mgstruct->curr_level+1];
        mgstruct->curr_y_coarse /= mgstruct->blocksize_y[mgstruct->curr_level+1];
        
        mgstruct->curr_coarse_size = mgstruct->curr_y_coarse*mgstruct->curr_x_coarse*mgstruct->curr_dof_coarse;
        
        // Update the level.
        mgstruct->curr_level++;
    }
}

// Pull the operator struct up a level, updating curr_* variables.
void level_up(mg_operator_struct_complex* mgstruct)
{
    if (mgstruct->curr_level > 0) // Can't go lower than the number of refinements!
    {
        // Update the level.
        mgstruct->curr_level--;
        
        // Update the number of d.o.f.
        if (mgstruct->curr_level == 0)
        {
            mgstruct->curr_dof_fine /= mgstruct->n_vector/mgstruct->Nc;  // Nc for square laplace only. 
        }
        
        // Update the current fine.
        mgstruct->curr_x_fine *= mgstruct->blocksize_x[mgstruct->curr_level];
        mgstruct->curr_y_fine *= mgstruct->blocksize_y[mgstruct->curr_level];
        
        mgstruct->curr_fine_size = mgstruct->curr_y_fine*mgstruct->curr_x_fine*mgstruct->curr_dof_fine;
        
        // Update the current coarse.
        mgstruct->curr_x_coarse *= mgstruct->blocksize_x[mgstruct->curr_level+1];
        mgstruct->curr_y_coarse *= mgstruct->blocksize_y[mgstruct->curr_level+1];
        
        mgstruct->curr_coarse_size = mgstruct->curr_y_coarse*mgstruct->curr_x_coarse*mgstruct->curr_dof_coarse;
    }
}


// MG preconditioner!! (Man, I'm excited!
void mg_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data, inversion_verbose_struct* verb)
{
    // Iterator.
    int i;
    
    cout << "[MG]: Entered mg_preconditioner.\n";
    mg_precond_struct_complex* mgprecond = (mg_precond_struct_complex*)extra_data; 
    
    // Standard defines.
    int x_fine = mgprecond->mgstruct->curr_x_fine;
    int y_fine = mgprecond->mgstruct->curr_y_fine;
    int fine_size = mgprecond->mgstruct->curr_fine_size;
    int x_coarse = mgprecond->mgstruct->curr_x_coarse; // how many coarse sites are there in the x dir?
    int y_coarse = mgprecond->mgstruct->curr_y_coarse; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 
    int coarse_length = coarse_size*mgprecond->mgstruct->n_vector; 
    
    // Store inversion info.
    inversion_info invif;
    
    // GET EXCITED! First off, let's do some pre-smoothing. 
    // 1. z1 = smoothed rhs. 
    complex<double>* z1 = new complex<double>[fine_size];
    copy<double>(z1, rhs, fine_size); 
    complex<double>* r1 = new complex<double>[fine_size];
    zero<double>(r1, fine_size);
    if (mgprecond->n_pre_smooth > 0 && mgprecond->in_smooth_type != NONE)
    {
        switch (mgprecond->in_smooth_type)
        {
            case NONE: // can't reach here, anyway.
                break;
            case CG:
                invif = minv_vector_cg(z1, rhs, fine_size, mgprecond->n_pre_smooth, 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case MR:
                invif = minv_vector_mr(z1, rhs, fine_size, mgprecond->n_pre_smooth, 1e-20, mgprecond->omega_smooth, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case GCR: 
                invif = minv_vector_gcr(z1, rhs, fine_size, mgprecond->n_pre_smooth, 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break; 
            case BICGSTAB:
                invif = minv_vector_bicgstab(z1, rhs, fine_size, mgprecond->n_pre_smooth, 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
        }
        printf("[L%d Presmooth]: Iterations %d Res %.8e Err N Algorithm %s\n", mgprecond->mgstruct->curr_level+1, invif.iter, sqrt(invif.resSq), invif.name.c_str()); fflush(stdout);
        
        // Compute r1 = r - A z1
        (*mgprecond->fine_matrix_vector)(r1, z1, mgprecond->matrix_extra_data); // Temporarily store Az1 in r1. 
        for (int i = 0; i < fine_size; i++)
        {
            r1[i] = rhs[i] - r1[i];
        }
        
    }
    else
    {
        copy<double>(z1, rhs, fine_size);
        // Compute r1 = r - A z1
        (*mgprecond->fine_matrix_vector)(r1, z1, mgprecond->matrix_extra_data); // Temporarily store Az1 in r1. 
        for (int i = 0; i < fine_size; i++)
        {
            r1[i] = rhs[i] - r1[i];
        }
        
    }
    
    // Next, solve z2 = P(P^\dag A P)^(-1) r1, r2 = r1 - A z2
    complex<double>* z2 = new complex<double>[fine_size]; zero<double>(z2, fine_size);
    complex<double>* r2 = new complex<double>[fine_size]; zero<double>(r2, fine_size); 
    if (mgprecond->in_solve_type != NONE)
    {
        
        // z2 = P ( P^\dag A P )^(-1) P^\dag r1 
        
        // Restrict z_smooth. This forms P^\dag r1. 
        complex<double>* rhs_coarse = new complex<double>[coarse_length];
        zero<double>(rhs_coarse, coarse_length); 
        restrict(rhs_coarse, r1, mgprecond->mgstruct);
        
        // Apply (P^\dag A P)^(-1). This forms  (P^\dag A P)^(-1) P^\dag r1. 
        // This gets modified if we're going to do multi-level MG, see commented out material below.
        complex<double>* lhs_coarse = new complex<double>[coarse_length];
        zero<double>(lhs_coarse, coarse_length);
        if (mgprecond->mgstruct->curr_level+1 == mgprecond->mgstruct->n_refine) // We're already on the coarsest level.
        {
            switch (mgprecond->in_solve_type)
            {
                case NONE: // The code can't reach here, anyway.
                    break;
                case MR:
                    invif = minv_vector_mr(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res, mgprecond->coarse_matrix_vector, mgprecond->matrix_extra_data, verb);
                    break;
                case CG:
                    invif = minv_vector_cg(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res, mgprecond->coarse_matrix_vector, mgprecond->matrix_extra_data, verb);
                    break;
                case GCR:
                    invif = minv_vector_gcr_restart(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res, mgprecond->n_restart, mgprecond->coarse_matrix_vector, mgprecond->matrix_extra_data, verb);
                    break;
                case BICGSTAB:
                    invif = minv_vector_bicgstab(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res, mgprecond->coarse_matrix_vector, mgprecond->matrix_extra_data, verb);
                    break;
            }
            
            printf("[L%d]: Iterations %d RelRes %.8e Err N Algorithm %s\n", mgprecond->mgstruct->curr_level+2, invif.iter, sqrt(invif.resSq)/sqrt(norm2sq<double>(rhs_coarse, coarse_length)), invif.name.c_str());
        }
        else // Apply the fine operator preconditioned with the coarse op.
        {
            printf("About to enter coarser solve.\n"); fflush(stdout);
            level_down(mgprecond->mgstruct);
            mg_preconditioner(lhs_coarse, rhs_coarse, coarse_length, extra_data, verb);
            level_up(mgprecond->mgstruct);
            printf("Exited coarser solve.\n"); fflush(stdout);
        }
            
        
        // Prolong lhs_coarse. This forms P ( ( P^\dag A P )^(-1) P^\dag r1
        prolong(z2, lhs_coarse, mgprecond->mgstruct); 
        
        // Add prolonged lhs to original rhs. This forms r + P ( (P^\dag A P)^(-1) - 1 ) P^\dag r.
        // While we're at it, add the solution from the first part (if we presmoothed). 
        for (i = 0; i < fine_size; i++)
        {
            lhs[i] = z1[i] + z2[i]; 
        }
        
        // Compute an updated residual from z2. 
        (*mgprecond->fine_matrix_vector)(r2, z2, mgprecond->matrix_extra_data); // Temporarily store Az2 in r2. 
        for (int i = 0; i < fine_size; i++)
        {
            r2[i] = r1[i] - r2[i];
        }
        
        // And we're done!
        
        // Clean up!
        delete[] rhs_coarse;
        delete[] lhs_coarse;
        
    }
    else // no inner solver, set lhs to z1, copy r1 into r2 (since z2 is trivially 0).
    {
        for (i = 0; i < fine_size; i++)
        {
            lhs[i] = z1[i];
            r2[i] = r1[i];
        }
        
    }
    
    complex<double>* z3 = new complex<double>[fine_size]; zero<double>(z3, fine_size); 
    // Almost done! Do some post-smoothing. on z3 = A^(-1) r2
    // 7. Post smooth.
    if (mgprecond->n_post_smooth > 0 && mgprecond->in_smooth_type != NONE)
    {
        switch (mgprecond->in_smooth_type)
        {
            case NONE: // Can't reach here, anyway.
                break;
            case CG:
                invif = minv_vector_cg(z3, r2, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case MR:
                invif = minv_vector_mr(z3, r2, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->omega_smooth, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case GCR:
                invif = minv_vector_gcr(z3, r2, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case BICGSTAB:
                invif = minv_vector_bicgstab(z3, r2, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
        }
        printf("[L%d Postsmooth]: Iterations %d Res %.8e Err N Algorithm %s\n", mgprecond->mgstruct->curr_level+1, invif.iter, sqrt(invif.resSq), invif.name.c_str());
    }
    // else z3 = 0 is fine.
    
    // Update lhs with z3. 
    for (i = 0; i < fine_size; i++)
    {
        lhs[i] += z3[i]; 
    }
    
    // Clean up!
    delete[] z1;
    delete[] z2;
    delete[] z3; 
    delete[] r1;
    delete[] r2; 
    
    cout << "[MG]: Exited mg_preconditioner.\n";
    
}

