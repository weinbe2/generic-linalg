
#include <iostream>
#include <iomanip> // to set output precision.
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

using namespace std;


#include "generic_vector.h"

#include "verbosity.h"

#include "generic_bicgstab.h"
#include "generic_bicgstab_l.h"
#include "generic_cg.h"
#include "generic_cr.h"
#include "generic_gcr.h"
#include "generic_minres.h"
#include "generic_cg_flex_precond.h"
#include "generic_gcr_var_precond.h"
#include "generic_bicgstab_precond.h"


#include "mg.h"
#include "mg_complex.h"
#include "lattice.h"

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
    
    // Check if we have a stencil!
    if (mgstruct->stencils[mgstruct->curr_level+1] != 0 && mgstruct->stencils[mgstruct->curr_level+1]->generated == true)
    {
        // Just apply the stencil! Life is so good!
        apply_stencil_2d(lhs, rhs, (void*)mgstruct->stencils[mgstruct->curr_level+1]);
    }
    else // We don't have a stencil, do it the old fashioned way!
    {
        // prolong, apply fine stencil, restrict.
        complex<double>* Prhs = new complex<double>[mgstruct->latt[mgstruct->curr_level]->get_lattice_size()];
        complex<double>* APrhs = new complex<double>[mgstruct->latt[mgstruct->curr_level]->get_lattice_size()];
        
        zero<double>(Prhs, mgstruct->latt[mgstruct->curr_level]->get_lattice_size());
        zero<double>(APrhs, mgstruct->latt[mgstruct->curr_level]->get_lattice_size());
        
        prolong(Prhs, rhs, mgstruct);
        
        fine_square_staggered(APrhs, Prhs, extra_data);
        
        restrict(lhs, APrhs, mgstruct);
        
        delete[] Prhs;
        delete[] APrhs;
    }
    
    /*
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
    }*/
    
    
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
                            dot_prod += conj(null_vectors[m][fine_site])*null_vectors[curr_dof_coarse][fine_site];
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
        // Update the level.
        mgstruct->curr_level++;
        
        // Update curr dof.
        mgstruct->curr_dof_fine = mgstruct->latt[mgstruct->curr_level]->get_nc(); // Top level has only one d.o.f. per site. 
        mgstruct->curr_x_fine = mgstruct->latt[mgstruct->curr_level]->get_lattice_dimension(0);
        mgstruct->curr_y_fine = mgstruct->latt[mgstruct->curr_level]->get_lattice_dimension(1);
        mgstruct->curr_fine_size = mgstruct->latt[mgstruct->curr_level]->get_lattice_size();

        mgstruct->curr_dof_coarse = mgstruct->latt[mgstruct->curr_level+1]->get_nc();
        mgstruct->curr_x_coarse = mgstruct->latt[mgstruct->curr_level+1]->get_lattice_dimension(0);
        mgstruct->curr_y_coarse = mgstruct->latt[mgstruct->curr_level+1]->get_lattice_dimension(1);
        mgstruct->curr_coarse_size = mgstruct->latt[mgstruct->curr_level+1]->get_lattice_size();
        
    }
}

// Pull the operator struct up a level, updating curr_* variables.
void level_up(mg_operator_struct_complex* mgstruct)
{
    if (mgstruct->curr_level > 0) // Can't go lower than the number of refinements!
    {
        // Update the level.
        mgstruct->curr_level--;
        
        // Update curr dof.
        mgstruct->curr_dof_fine = mgstruct->latt[mgstruct->curr_level]->get_nc(); // Top level has only one d.o.f. per site. 
        mgstruct->curr_x_fine = mgstruct->latt[mgstruct->curr_level]->get_lattice_dimension(0);
        mgstruct->curr_y_fine = mgstruct->latt[mgstruct->curr_level]->get_lattice_dimension(1);
        mgstruct->curr_fine_size = mgstruct->latt[mgstruct->curr_level]->get_lattice_size();

        mgstruct->curr_dof_coarse = mgstruct->latt[mgstruct->curr_level+1]->get_nc();
        mgstruct->curr_x_coarse = mgstruct->latt[mgstruct->curr_level+1]->get_lattice_dimension(0);
        mgstruct->curr_y_coarse = mgstruct->latt[mgstruct->curr_level+1]->get_lattice_dimension(1);
        mgstruct->curr_coarse_size = mgstruct->latt[mgstruct->curr_level+1]->get_lattice_size();
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
    if (mgprecond->n_pre_smooth[mgprecond->mgstruct->curr_level] > 0 && mgprecond->in_smooth_type != NONE)
    {
        switch (mgprecond->in_smooth_type)
        {
            case NONE: // can't reach here, anyway.
                break;
            case CG:
                invif = minv_vector_cg(z1, rhs, fine_size, mgprecond->n_pre_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case MINRES:
                invif = minv_vector_minres(z1, rhs, fine_size, mgprecond->n_pre_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->omega_smooth, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case GCR: 
                invif = minv_vector_gcr(z1, rhs, fine_size, mgprecond->n_pre_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break; 
            case BICGSTAB:
                invif = minv_vector_bicgstab(z1, rhs, fine_size, mgprecond->n_pre_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case CR:
                invif = minv_vector_cr(z1, rhs, fine_size, mgprecond->n_pre_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break; 
            case BICGSTAB_L:
                invif = minv_vector_bicgstab_l(z1, rhs, fine_size, mgprecond->n_pre_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->n_pre_smooth[mgprecond->mgstruct->curr_level], mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
        }
        printf("[L%d Presmooth]: Iterations %d Res %.8e Err N Algorithm %s\n", mgprecond->mgstruct->curr_level+1, invif.iter, sqrt(invif.resSq), invif.name.c_str()); fflush(stdout);
        mgprecond->mgstruct->dslash_count->presmooth[mgprecond->mgstruct->curr_level] += invif.ops_count; 
        
        // Compute r1 = r - A z1
        (*mgprecond->fine_matrix_vector)(r1, z1, mgprecond->matrix_extra_data); // Temporarily store Az1 in r1. 
        mgprecond->mgstruct->dslash_count->residual[mgprecond->mgstruct->curr_level]++;
        
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
        mgprecond->mgstruct->dslash_count->residual[mgprecond->mgstruct->curr_level]++;
        for (int i = 0; i < fine_size; i++)
        {
            r1[i] = rhs[i] - r1[i];
        }
        
    }
    
    // Next, solve z2 = P(P^\dag A P)^(-1) r1, r2 = r1 - A z2
    complex<double>* z2 = new complex<double>[fine_size]; zero<double>(z2, fine_size);
    complex<double>* r2 = new complex<double>[fine_size]; zero<double>(r2, fine_size); 
    //if (mgprecond->in_solve_type != NONE)
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
                    invif.resSq = 1;
                    invif.ops_count = 0; 
                    invif.iter = 0;
                    invif.success = true;
                    invif.name = "None";
                    break;
                case MINRES:
                    invif = minv_vector_minres(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res[mgprecond->mgstruct->curr_level], mgprecond->coarse_matrix_vector, mgprecond->matrix_extra_data, verb);
                    break;
                case CG:
                    invif = minv_vector_cg(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res[mgprecond->mgstruct->curr_level], mgprecond->coarse_matrix_vector, mgprecond->matrix_extra_data, verb);
                    break;
                case GCR:
                    invif = minv_vector_gcr_restart(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res[mgprecond->mgstruct->curr_level], mgprecond->n_restart, mgprecond->coarse_matrix_vector, mgprecond->matrix_extra_data, verb);
                    break;
                case BICGSTAB:
                    invif = minv_vector_bicgstab(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res[mgprecond->mgstruct->curr_level], mgprecond->coarse_matrix_vector, mgprecond->matrix_extra_data, verb);
                    break;
                case CR:
                    invif = minv_vector_cr_restart(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res[mgprecond->mgstruct->curr_level], mgprecond->n_restart, mgprecond->coarse_matrix_vector, mgprecond->matrix_extra_data, verb);
                    break;
            }
            
            printf("[L%d]: Iterations %d RelRes %.8e Err N Algorithm %s\n", mgprecond->mgstruct->curr_level+2, invif.iter, sqrt(invif.resSq)/sqrt(norm2sq<double>(rhs_coarse, coarse_length)), invif.name.c_str());
            mgprecond->mgstruct->dslash_count->krylov[mgprecond->mgstruct->curr_level+1] += invif.ops_count; 
        
        }
        else // Apply the fine operator preconditioned with the coarse op.
        {
            printf("About to enter coarser solve.\n"); fflush(stdout);
            level_down(mgprecond->mgstruct);
            
            switch (mgprecond->mlevel_type)
            {
                case MLEVEL_SMOOTH:
                    mg_preconditioner(lhs_coarse, rhs_coarse, coarse_length, extra_data, verb);
                    break;
                case MLEVEL_RECURSIVE:
                    switch (mgprecond->in_solve_type)
                    {
                        case NONE: // The code can't reach here, anyway. 
                            invif.resSq = 1;
                            invif.ops_count = 0;
                            invif.iter = 0;
                            invif.success = true;
                            invif.name = "None";
                            break;
                        case MINRES: // Really just pretend this is MLEVEL_SMOOTH.
                            mg_preconditioner(lhs_coarse, rhs_coarse, coarse_length, extra_data, verb); // Need some way to count these...
                            invif.resSq = 1;
                            invif.ops_count = 0; // This gets handled lower down.
                            invif.iter = 0;
                            invif.success = true;
                            invif.name = "None";
                            break;
                        case CG:
                            invif = minv_vector_cg_flex_precond_restart(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res[mgprecond->mgstruct->curr_level-1], mgprecond->n_restart, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data, mg_preconditioner, (void*)mgprecond, verb); 
                            break;
                        case GCR:
                        case CR: // Since I don't have flexible CR...
                            invif = minv_vector_gcr_var_precond_restart(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res[mgprecond->mgstruct->curr_level-1], mgprecond->n_restart, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data, mg_preconditioner, (void*)mgprecond, verb); 
                            break;
                        case BICGSTAB:
                            invif = minv_vector_bicgstab_precond(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res[mgprecond->mgstruct->curr_level-1], mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data, mg_preconditioner, (void*)mgprecond, verb);  // Should probably add a flag for restarting, but meh, it's bicg, it's already flexible. 
                            break;
                            
                    }
                    printf("[L%d]: Iterations %d RelRes %.8e Err N Algorithm %s\n", mgprecond->mgstruct->curr_level+1, invif.iter, sqrt(invif.resSq)/sqrt(norm2sq<double>(rhs_coarse, coarse_length)), invif.name.c_str());
                    mgprecond->mgstruct->dslash_count->krylov[mgprecond->mgstruct->curr_level] += invif.ops_count; 
                    break;
            }
            
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
        mgprecond->mgstruct->dslash_count->residual[mgprecond->mgstruct->curr_level]++; 
        
        for (int i = 0; i < fine_size; i++)
        {
            r2[i] = r1[i] - r2[i];
        }
        
        // And we're done!
        
        // Clean up!
        delete[] rhs_coarse;
        delete[] lhs_coarse;
        
    }
    /*else // no inner solver, set lhs to z1, copy r1 into r2 (since z2 is trivially 0).
    {
        for (i = 0; i < fine_size; i++)
        {
            lhs[i] = z1[i];
            r2[i] = r1[i];
        }
        
    }*/
    
    complex<double>* z3 = new complex<double>[fine_size]; zero<double>(z3, fine_size); 
    // Almost done! Do some post-smoothing. on z3 = A^(-1) r2
    // 7. Post smooth.
    if (mgprecond->n_post_smooth[mgprecond->mgstruct->curr_level] > 0 && mgprecond->in_smooth_type != NONE)
    {
        switch (mgprecond->in_smooth_type)
        {
            case NONE: // Can't reach here, anyway.
                break;
            case CG:
                invif = minv_vector_cg(z3, r2, fine_size, mgprecond->n_post_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case MINRES:
                invif = minv_vector_minres(z3, r2, fine_size, mgprecond->n_post_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->omega_smooth, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case GCR:
                invif = minv_vector_gcr(z3, r2, fine_size, mgprecond->n_post_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case BICGSTAB:
                invif = minv_vector_bicgstab(z3, r2, fine_size, mgprecond->n_post_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case CR:
                invif = minv_vector_gcr(z3, r2, fine_size, mgprecond->n_post_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
            case BICGSTAB_L:
                invif = minv_vector_bicgstab_l(z3, r2, fine_size, mgprecond->n_post_smooth[mgprecond->mgstruct->curr_level], 1e-20, mgprecond->n_post_smooth[mgprecond->mgstruct->curr_level], mgprecond->fine_matrix_vector, mgprecond->matrix_extra_data); 
                break;
        }
        printf("[L%d Postsmooth]: Iterations %d Res %.8e Err N Algorithm %s\n", mgprecond->mgstruct->curr_level+1, invif.iter, sqrt(invif.resSq), invif.name.c_str());
        mgprecond->mgstruct->dslash_count->postsmooth[mgprecond->mgstruct->curr_level] += invif.ops_count; 
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

// Generate a coarse stencil from a fine stencil. This takes advantage of the prolong and restrict functions
// explicitly, and also depends on the "sdir" variable the stencil object includes.
void generate_coarse_from_fine_stencil(stencil_2d* stenc_coarse, stencil_2d* stenc_fine, mg_operator_struct_complex* mgstruct)
{   
    if (stenc_coarse->generated || stenc_fine->stencil_size > 2 || stenc_coarse->stencil_size > 2)
    {
        return;
    }
    
    Lattice* latt_fine = stenc_fine->lat;
    Lattice* latt_coarse = stenc_coarse->lat; 
    
    int i,x,y;
    int coord[2];
    int color, c;
    
    int stencil_distance = stenc_coarse->stencil_size; 
    
    int fine_size = latt_fine->get_lattice_size(); // Only required to allocate temporary vectors. 
    int coarse_size = latt_coarse->get_lattice_size();
    int coarse_nc = latt_coarse->get_nc();
    
    // Coarse vectors.
    complex<double>* tmp_rhs = new complex<double>[coarse_size];
    complex<double>* tmp_lhs = new complex<double>[coarse_size];
    
    // Fine vectors.
    complex<double>* tmp_Prhs = new complex<double>[fine_size];
    complex<double>* tmp_APrhs = new complex<double>[fine_size];
    
    
    // Save the state of the fine stencil.
    stencil_dir saved_dir = stenc_fine->sdir; 
    
    // We build the coarse stencil one piece at a time: the fine clover, the fine hopping terms, the fine two-link terms.
    
    // Piece 1: Fine Clover.
    // For the fine clover, we learn about some of the coarse clover. This takes one application per coarse color.
    stenc_fine->sdir = DIR_0;
    for (color = 0; color < coarse_nc; color++)
    {
        zero<double>(tmp_rhs, coarse_size);
        
        // Set a 1 at each coarse site.
        for (i = 0; i < latt_coarse->get_volume(); i++)
        {
            tmp_rhs[i*coarse_nc+color] = 1.0;
        }
        
        // Apply P^\dag A_0 P
        prolong(tmp_Prhs, tmp_rhs, mgstruct);
        apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
        restrict(tmp_lhs, tmp_APrhs, mgstruct);
        
        // Loop over each coarse site, update the clover.
        for (i = 0; i < latt_coarse->get_lattice_size(); i++)
        {
            latt_coarse->index_to_coord(i, (int*)coord, c);
            
            stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index(coord, c)];
        }
    }
    
    // Piece 2: The Hopping Term.
    // For the fine hopping term, we learn about more of the coarse clover, and some of the coarse hopping term.
    // In the case of a distance = 1 fine stencil, we learn about the entire coarse hopping term,
    // and also finish learning about the coarse clover.
    // We learn about the hopping term in two hits per direction: one for even sites, one for odd sites.
    for (color = 0; color < coarse_nc; color++)
    {   
        // First step: Learn about the +x part of the stencil.
        // One for the even x sites, one for the odd x sites.
        {
            stenc_fine->sdir = DIR_XP1;


            // Even sites first.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 0; i < latt_coarse->get_volume(); i+=2)
            {
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_+x P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[0] % 2 == 0) // even, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index(coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // odd, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index(coord,c)+0*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index(coord, c)];
                }
            }

            // Odd sites second.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 1; i < latt_coarse->get_volume(); i+=2)
            {
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_+x P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[0] % 2 == 1) // odd, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // even, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index(coord,c)+0*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index(coord, c)];
                }
            }
        }
        // Done with +x.
        //cout << "Generated +x stencil.\n" << flush; 
        
        // Second step: Learn about the -x part of the stencil.
        // One for the even x sites, one for the odd x sites.
        {
            stenc_fine->sdir = DIR_XM1;


            // Even sites first.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 0; i < latt_coarse->get_volume(); i+=2)
            {
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_-x P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[0] % 2 == 0) // even, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // odd, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+2*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }

            // Odd sites second.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 1; i < latt_coarse->get_volume(); i+=2)
            {
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_-x P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[0] % 2 == 1) // odd, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // even, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+2*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }
        }
        // Done with -x.
        //cout << "Generated -x stencil.\n" << flush; 
        
        // Third step: Learn about the +y part of the stencil.
        // One for the even y sites, one for the odd y sites.
        {
            stenc_fine->sdir = DIR_YP1;


            // Even sites first.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 0; i < latt_coarse->get_volume()-latt_coarse->get_lattice_dimension(0); i++)
            {
                if ((i/latt_coarse->get_lattice_dimension(0)) % 2 == 1) // if it has an odd y...
                {
                    i += latt_coarse->get_lattice_dimension(0);
                }
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_+y P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[1] % 2 == 0) // even, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // odd, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+1*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }

            // Odd sites second.
            zero<double>(tmp_rhs, coarse_size);
            for (i = latt_coarse->get_lattice_dimension(0); i < latt_coarse->get_volume(); i++)
            {
                if ((i/latt_coarse->get_lattice_dimension(0)) % 2 == 0) // if it has an even y...
                {
                    i += latt_coarse->get_lattice_dimension(0);
                }
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_+y P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[1] % 2 == 1) // odd, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // even, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+1*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }
        }
        // Done with +y.
        //cout << "Generated +y stencil.\n" << flush; 
        
        
        // Second step: Learn about the -y part of the stencil.
        // One for the even y sites, one for the odd y sites.
        {
            stenc_fine->sdir = DIR_YM1;


            // Even sites first.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 0; i < latt_coarse->get_volume()-latt_coarse->get_lattice_dimension(0); i++)
            {
                if ((i/latt_coarse->get_lattice_dimension(0)) % 2 == 1) // if it has an odd y...
                {
                    i += latt_coarse->get_lattice_dimension(0);
                }
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_-y P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[1] % 2 == 0) // even, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // odd, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+3*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }

            // Odd sites second.
            zero<double>(tmp_rhs, coarse_size);
            for (i = latt_coarse->get_lattice_dimension(0); i < latt_coarse->get_volume(); i++)
            {
                if ((i/latt_coarse->get_lattice_dimension(0)) % 2 == 0) // if it has an even y...
                {
                    i += latt_coarse->get_lattice_dimension(0);
                }
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_-y P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[1] % 2 == 1) // odd, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // even, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+3*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }
        }
        // Done with -y.
        //cout << "Generated -y stencil.\n" << flush;
        
        
        // Need to add other pieces...
        
        //cout << "Finished color " << color << "\n" << flush; 
    }
         
    delete[] tmp_lhs;
    delete[] tmp_rhs; 
    
    delete[] tmp_Prhs;
    delete[] tmp_APrhs; 
    
    stenc_coarse->generated = true; 
    
    // Restore the state of the fine stencil.
    stenc_fine->sdir = saved_dir; 
    
}




