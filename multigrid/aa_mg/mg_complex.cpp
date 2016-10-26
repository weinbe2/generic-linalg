
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
#include "operators.h"

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

// General multigrid projector function!
// Apply the current coarse operator dagger.
void coarse_square_staggered_dagger(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Grab the mg_precond_struct.
    mg_operator_struct_complex* mgstruct = (mg_operator_struct_complex*)extra_data; 
    
    // Check if we have a stencil!
    if (mgstruct->dagger_stencils != 0 && mgstruct->dagger_stencils[mgstruct->curr_level+1] != 0 && mgstruct->dagger_stencils[mgstruct->curr_level+1]->generated == true)
    {
        // Just apply the stencil! Life is so good!
   
        apply_stencil_2d(lhs, rhs, (void*)mgstruct->dagger_stencils[mgstruct->curr_level+1]);
    }
    else // We don't have a stencil, do it the old fashioned way!
    {   
        // prolong, apply fine stencil, restrict.
        complex<double>* Prhs = new complex<double>[mgstruct->latt[mgstruct->curr_level]->get_lattice_size()];
        complex<double>* APrhs = new complex<double>[mgstruct->latt[mgstruct->curr_level]->get_lattice_size()];
        
        zero<double>(Prhs, mgstruct->latt[mgstruct->curr_level]->get_lattice_size());
        zero<double>(APrhs, mgstruct->latt[mgstruct->curr_level]->get_lattice_size());
        
        prolong(Prhs, rhs, mgstruct);
        
        fine_square_staggered_dagger(APrhs, Prhs, extra_data);
        
        restrict(lhs, APrhs, mgstruct);
        
        delete[] Prhs;
        delete[] APrhs;
    }
}

// General multigrid projector function!
// Apply the current fine operator dagger.
void fine_square_staggered_dagger(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Grab the mg_precond_struct.
    mg_operator_struct_complex* mgstruct = (mg_operator_struct_complex*)extra_data; 
    
    if (mgstruct->curr_level == 0)
    {
        zero<double>(lhs, mgstruct->curr_fine_size);
        (*mgstruct->matrix_vector_dagger)(lhs, rhs, mgstruct->matrix_extra_data);
    }
    else
    {
        level_up(mgstruct);
        coarse_square_staggered_dagger(lhs, rhs, extra_data);
        level_down(mgstruct);
    }
}

// General multigrid projector function!
// Apply the current coarse normal op.
void coarse_square_staggered_normal(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Grab the mg_precond_struct.
    mg_operator_struct_complex* mgstruct = (mg_operator_struct_complex*)extra_data; 
    
    complex<double> *tmp = new complex<double>[mgstruct->latt[mgstruct->curr_level+1]->get_lattice_size()];
    zero<double>(tmp, mgstruct->latt[mgstruct->curr_level+1]->get_lattice_size());
    zero<double>(lhs, mgstruct->latt[mgstruct->curr_level+1]->get_lattice_size());
    
    // Apply D
    coarse_square_staggered(tmp, rhs, extra_data);

    // Apply D dagger.
    coarse_square_staggered_dagger(lhs, tmp, extra_data);
    
    // Clean up.
    delete[] tmp; 
}

// General multigrid projector function!
// Apply the current fine normal op.
void fine_square_staggered_normal(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Grab the mg_precond_struct.
    mg_operator_struct_complex* mgstruct = (mg_operator_struct_complex*)extra_data; 
    
    complex<double> *tmp = new complex<double>[mgstruct->latt[mgstruct->curr_level]->get_lattice_size()];
    zero<double>(tmp, mgstruct->latt[mgstruct->curr_level]->get_lattice_size());
    zero<double>(lhs, mgstruct->latt[mgstruct->curr_level]->get_lattice_size());
    
    // Apply D
    fine_square_staggered(tmp, rhs, extra_data);
    
    // Apply D dagger.
    fine_square_staggered_dagger(lhs, tmp, extra_data);
    
    // Clean up.
    delete[] tmp; 
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
    
    // Create a tmp for smoothing as needed.
    complex<double>* tmp_smooth = new complex<double>[fine_size];
    
    // What operator are we smoothing with?
    void (*smooth_op)(complex<double>*,complex<double>*,void*);
    if (mgprecond->normal_eqn_smooth) // Smooth on D^dag D x = D^dag b, CGNR. 
    {
        smooth_op = mgprecond->fine_matrix_vector_normal;
    }
    else // Smooth on D x = b
    {
        smooth_op = mgprecond->fine_matrix_vector;
    }
    
    // GET EXCITED! First off, let's do some pre-smoothing. 
    // 1. z1 = smoothed rhs. 
    // z1^0 = 0, rhs is... well, rhs. 
    complex<double>* z1 = new complex<double>[fine_size];
    zero<double>(z1, fine_size);
    complex<double>* r1 = new complex<double>[fine_size];
    zero<double>(r1, fine_size);
    if (mgprecond->n_pre_smooth[mgprecond->mgstruct->curr_level] > 0 && mgprecond->in_smooth_type != MINV_INVALID)
    {
        // Are we smoothing via the normal equation?
        if (mgprecond->normal_eqn_smooth && !mgprecond->normal_eqn_mg) 
        {
            // We smooth via b = A^\dag rhs.
            zero<double>(tmp_smooth, fine_size);
            mgprecond->fine_matrix_vector_dagger(tmp_smooth, rhs, mgprecond->matrix_extra_data);
        }
        else
        {
            // smooth via b = rhs.
            copy<double>(tmp_smooth, rhs, fine_size);
        }
            
        // Populate an inverter struct.
        minv_inverter_params pre_solve;
        pre_solve.tol = 1e-20; 
        pre_solve.max_iters = mgprecond->n_pre_smooth[mgprecond->mgstruct->curr_level];
        pre_solve.restart = false;
        pre_solve.restart_freq = -1;
        pre_solve.minres_omega = 1.0; // should expose.
        pre_solve.bicgstabl_l = mgprecond->n_pre_smooth[mgprecond->mgstruct->curr_level];
        
        invif = minv_unpreconditioned(z1, tmp_smooth, fine_size, mgprecond->in_smooth_type, pre_solve, smooth_op, mgprecond->matrix_extra_data);

        printf("[L%d Presmooth]: Iterations %d Res %.8e Err N Algorithm %s\n", mgprecond->mgstruct->curr_level+1, invif.iter, sqrt(invif.resSq), invif.name.c_str()); fflush(stdout);
        mgprecond->mgstruct->dslash_count->presmooth[mgprecond->mgstruct->curr_level] += ((mgprecond->normal_eqn_smooth || mgprecond->normal_eqn_mg) ? 2 : 1)*invif.ops_count; 
        
        // Compute r1 = r - A z1.
        (*mgprecond->fine_matrix_vector)(r1, z1, mgprecond->matrix_extra_data); // Temporarily store Az1 in r1. 
        mgprecond->mgstruct->dslash_count->residual[mgprecond->mgstruct->curr_level]++;
        
        for (int i = 0; i < fine_size; i++)
        {
            r1[i] = rhs[i] - r1[i];
        }
        
    }
    else
    {
        // z1 = 0, r1 = rhs.
        copy<double>(r1, rhs, fine_size);
        zero<double>(z1, fine_size);
    }
    
    // Next, solve z2 = P(P^\dag A P)^(-1) r1. Current lhs is z1 + z2. 
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
                case BICGSTAB_L: // not really, BiCGstab(l) doesn't have a faithful precond form.
                    invif = minv_vector_bicgstab(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res[mgprecond->mgstruct->curr_level], mgprecond->coarse_matrix_vector, mgprecond->matrix_extra_data, verb);
                    break;
                case CR:
                    invif = minv_vector_cr_restart(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_max, mgprecond->rel_res[mgprecond->mgstruct->curr_level], mgprecond->n_restart, mgprecond->coarse_matrix_vector, mgprecond->matrix_extra_data, verb);
                    break;
            }
            
            printf("[L%d]: Iterations %d RelRes %.8e Err N Algorithm %s\n", mgprecond->mgstruct->curr_level+2, invif.iter, sqrt(invif.resSq)/sqrt(norm2sq<double>(rhs_coarse, coarse_length)), invif.name.c_str());
            mgprecond->mgstruct->dslash_count->krylov[mgprecond->mgstruct->curr_level+1] += (mgprecond->normal_eqn_mg ? 2 : 1)*invif.ops_count; 
        
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
                        case BICGSTAB_L: // kind of a lie. BiCGstab(l) doesn't have a faithful precond form.
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
        
        // And we're done!
        
        // Clean up!
        delete[] rhs_coarse;
        delete[] lhs_coarse;
        
    }
    else // if we're in a V cycle and not at the coarsest level, recurse down.
    {
        if (mgprecond->mgstruct->curr_level+1 != mgprecond->mgstruct->n_refine) // We're not on the coarsest level.
        {
            // Restrict z_smooth. This forms P^\dag r1. 
            complex<double>* rhs_coarse = new complex<double>[coarse_length];
            zero<double>(rhs_coarse, coarse_length); 
            restrict(rhs_coarse, r1, mgprecond->mgstruct);

            // Apply (P^\dag A P)^(-1). This forms  (P^\dag A P)^(-1) P^\dag r1. 
            // This gets modified if we're going to do multi-level MG, see commented out material below.
            complex<double>* lhs_coarse = new complex<double>[coarse_length];
            zero<double>(lhs_coarse, coarse_length);
            
            // Level down.
            printf("[L%d]: About to enter coarser solve.\n", mgprecond->mgstruct->curr_level+1); fflush(stdout);
            level_down(mgprecond->mgstruct);
            
            // Recurse.
            mg_preconditioner(lhs_coarse, rhs_coarse, coarse_length, extra_data, verb);
            
            // Level up.
            level_up(mgprecond->mgstruct);
            printf("[L%d]: Exited coarse solve.\n", mgprecond->mgstruct->curr_level+1); fflush(stdout);
            
            // Prolong.
            prolong(z2, lhs_coarse, mgprecond->mgstruct); 
            
            // While we're at it, add the solution from the first part (if we presmoothed). 
            for (i = 0; i < fine_size; i++)
            {
                lhs[i] = z1[i] + z2[i]; 
            }
            
            // Clean up!
            delete[] rhs_coarse;
            delete[] lhs_coarse; 
        }
        else
        {
            // no inner solver, set lhs to z1, set z2 to 0.
            copy<double>(lhs, z1, fine_size);
            zero<double>(z2, fine_size);   
        }
    }
    
    complex<double>* z3 = new complex<double>[fine_size]; zero<double>(z3, fine_size); 
    // Almost done! Do some post-smoothing. on z3 = A^(-1) r2
    // 7. Post smooth.
    if (mgprecond->n_post_smooth[mgprecond->mgstruct->curr_level] > 0 && mgprecond->in_smooth_type != MINV_INVALID)
    {
        // Form the residual equation.
        // Compute r2 = rhs - A lhs.
        (*mgprecond->fine_matrix_vector)(r2, lhs, mgprecond->matrix_extra_data); // Temporarily store A*lhs in r2.
        mgprecond->mgstruct->dslash_count->residual[mgprecond->mgstruct->curr_level]++;
        
        for (int i = 0; i < fine_size; i++)
        {
            r2[i] = rhs[i] - r2[i];
        }
        
        // Are we smoothing via the normal equation?
        if (mgprecond->normal_eqn_smooth && !mgprecond->normal_eqn_mg) 
        {
            // We smooth via b = A^\dag r2.
            zero<double>(tmp_smooth, fine_size);
            mgprecond->fine_matrix_vector_dagger(tmp_smooth, r2, mgprecond->matrix_extra_data); // CGNR
        }
        else
        {
            // smooth via b = rhs.
            copy<double>(tmp_smooth, r2, fine_size);
        }
        
        // Populate an inverter struct.
        minv_inverter_params post_solve;
        post_solve.tol = 1e-20; 
        post_solve.max_iters = mgprecond->n_post_smooth[mgprecond->mgstruct->curr_level];
        post_solve.restart = false;
        post_solve.restart_freq = -1;
        post_solve.minres_omega = 1.0; // should expose.
        post_solve.bicgstabl_l = mgprecond->n_post_smooth[mgprecond->mgstruct->curr_level];
        
        invif = minv_unpreconditioned(z3, tmp_smooth, fine_size, mgprecond->in_smooth_type, post_solve, smooth_op, mgprecond->matrix_extra_data);
        
        printf("[L%d Postsmooth]: Iterations %d Res %.8e Err N Algorithm %s\n", mgprecond->mgstruct->curr_level+1, invif.iter, sqrt(invif.resSq), invif.name.c_str());
        mgprecond->mgstruct->dslash_count->postsmooth[mgprecond->mgstruct->curr_level] += ((mgprecond->normal_eqn_smooth || mgprecond->normal_eqn_mg) ? 2 : 1)*invif.ops_count; 
    }
    
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
    delete[] tmp_smooth; 
    
    cout << "[MG]: Exited mg_preconditioner.\n";
    
}

// Generate a coarse stencil from a fine stencil. This takes advantage of the prolong and restrict functions
// explicitly, and also depends on the "sdir" variable the stencil object includes.
// ignore_shifts = true -> set shifts to zero before building stencil. Otherwise, build shifts directly into new stencil.
void generate_coarse_from_fine_stencil(stencil_2d* stenc_coarse, stencil_2d* stenc_fine, mg_operator_struct_complex* mgstruct, bool ignore_shifts)
{   
    if (stenc_coarse->generated || stenc_fine->stencil_size > 2 || stenc_coarse->stencil_size > 2)
    {
        return;
    }
    
    Lattice* latt_fine = stenc_fine->lat;
    Lattice* latt_coarse = stenc_coarse->lat; 
    
    int i;
    int coord[2];
    int color, c;
    
    //int stencil_distance = stenc_coarse->stencil_size; 
    
    int fine_size = latt_fine->get_lattice_size(); // Only required to allocate temporary vectors. 
    int coarse_size = latt_coarse->get_lattice_size();
    int coarse_nc = latt_coarse->get_nc();
    
    // Coarse vectors.
    complex<double>* tmp_rhs = new complex<double>[coarse_size];
    complex<double>* tmp_lhs = new complex<double>[coarse_size];
    
    // Fine vectors.
    complex<double>* tmp_Prhs = new complex<double>[fine_size];
    complex<double>* tmp_APrhs = new complex<double>[fine_size];
    
    // Back up shifts.
    complex<double> shift = stenc_fine->shift;
    complex<double> eo_shift = stenc_fine->eo_shift;
    complex<double> dof_shift = stenc_fine->dof_shift;
    
    if (ignore_shifts)
    {
        stenc_fine->shift = 0.0;
        stenc_fine->eo_shift = 0.0;
        stenc_fine->dof_shift = 0.0;
    }
    
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
    
    // Put shifts back in
    if (ignore_shifts)
    {
        stenc_fine->shift = shift;
        stenc_fine->eo_shift = eo_shift;
        stenc_fine->dof_shift = dof_shift;
    }
         
    delete[] tmp_lhs;
    delete[] tmp_rhs; 
    
    delete[] tmp_Prhs;
    delete[] tmp_APrhs; 
    
    stenc_coarse->generated = true; 
    
    // Restore the state of the fine stencil.
    stenc_fine->sdir = saved_dir; 
    
}

void generate_coarse_from_fine_stencil(stencil_2d* stenc_coarse, stencil_2d* stenc_fine, mg_operator_struct_complex* mgstruct)
{
    generate_coarse_from_fine_stencil(stenc_coarse, stenc_fine, mgstruct, false);
}


