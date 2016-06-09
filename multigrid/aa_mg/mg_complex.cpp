
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
void coarse_square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data)
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
    
    // Apply gamma_5.
    /*for (i = 0; i < fine_size; i++)
    {
        int x = fine_size % mgstruct->x_fine;
        int y = fine_size / mgstruct->x_fine; 
        
        if ((x+y)%2 == 1)
        {
            APx[i] = -APx[i];
        }
    }*/
    
    // Restrict. 
    zero<double>(lhs, coarse_length);
    restrict(lhs, APx, mgstruct); 
    
    delete[] Px;
    delete[] APx; 
    
}


// Properly normalize the P vectors.
void block_normalize(mg_operator_struct_complex* mgstruct)
{
    int n, i;
    int x_fine = mgstruct->curr_x_fine;
    int y_fine = mgstruct->curr_y_fine;
    int fine_size = x_fine*y_fine;
    int x_coarse = mgstruct->curr_x_coarse; // how many coarse sites are there in the x dir?
    int y_coarse = mgstruct->curr_y_coarse; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 
    
    // Grab the current null vectors.
    complex<double>** null_vectors = mgstruct->null_vectors[mgstruct->curr_level];
    
    // Hold current sites.
    int curr_x, curr_y, curr_x_coarse, curr_y_coarse, curr_coarse; 
    
    // Build up norms in this array...
    double* norms = new double[coarse_size];
    
    // Loop over every vector.
    for (n = 0; n < mgstruct->n_vector; n++)
    {
        zero<double>(norms, coarse_size); 
        // Loop over the fine size.
        
        for (i = 0; i < fine_size; i++)
        {
            // What's the current coarse site? First, find the fine site.
            curr_x = i % x_fine;
            curr_y = i / x_fine; 
            
            // Now, find the coarse site. 
            curr_x_coarse = curr_x / mgstruct->blocksize_x[0];
            curr_y_coarse = curr_y / mgstruct->blocksize_y[0]; 
            curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 
            
            // Update the norm!
            norms[curr_coarse] += real(conj(null_vectors[n][i])*null_vectors[n][i]);
        }
        
        // Sqrt all of the norms.
        for (i = 0; i < coarse_size; i++)
        {
            norms[i] = sqrt(norms[i]);
        }
        
        // Normalize the null_vectors.
        for (i = 0; i < fine_size; i++)
        {
            // What's the current coarse site? First, find the fine site.
            curr_x = i % x_fine;
            curr_y = i / x_fine; 
            
            // Now, find the coarse site. 
            curr_x_coarse = curr_x / mgstruct->blocksize_x[0];
            curr_y_coarse = curr_y / mgstruct->blocksize_y[0]; 
            curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 
            
            // Update the norm!
            null_vectors[n][i] /= norms[curr_coarse];
        }
    }
    
    delete[] norms; 
}

// Properly orthonormalize the P vectors.
void block_orthonormalize(mg_operator_struct_complex* mgstruct)
{
    int n, m, i;
    int x_fine = mgstruct->curr_x_fine;
    int y_fine = mgstruct->curr_y_fine;
    int fine_size = x_fine*y_fine;
    int x_coarse = mgstruct->curr_x_coarse; // how many coarse sites are there in the x dir?
    int y_coarse = mgstruct->curr_y_coarse; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 
    
    // Grab the current null vectors.
    complex<double>** null_vectors = mgstruct->null_vectors[mgstruct->curr_level];
    
    // Hold current sites.
    int curr_x, curr_y, curr_x_coarse, curr_y_coarse, curr_coarse; 
    
    // Perform a modified gram schmidt. 
    
    if (mgstruct->n_vector > 1)
    {
        // Build up norms in this array...
        double* norms = new double[coarse_size];
        complex<double>* dot_prod = new complex<double>[coarse_size];

        
        // Loop over every vector.
        for (n = 1; n < mgstruct->n_vector; n++)
        {
            //cout << "\nn " << n << "\n";
            zero<double>(norms, coarse_size);

            // Get the new norm.
            for (i = 0; i < fine_size; i++)
            {
                // What's the current coarse site? First, find the fine site.
                curr_x = i % x_fine;
                curr_y = i / x_fine; 

                // Now, find the coarse site. 
                curr_x_coarse = curr_x / mgstruct->blocksize_x[0];
                curr_y_coarse = curr_y / mgstruct->blocksize_y[0]; 
                curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 

                // Update the norm of the previous vector,
                // Compute new dot products. 
                norms[curr_coarse] += real(conj(null_vectors[n-1][i])*null_vectors[n-1][i]);
            }
            
            // Sqrt all of the norms.
            /*cout << "norm ";
            for (i = 0; i < coarse_size; i++)
            {
                norms[i] = sqrt(norms[i]);
                cout << norms[i] << " ";
            }
            cout << "\n";*/

            // Normalize the null_vectors.
            for (i = 0; i < fine_size; i++)
            {
                // What's the current coarse site? First, find the fine site.
                curr_x = i % x_fine;
                curr_y = i / x_fine; 

                // Now, find the coarse site. 
                curr_x_coarse = curr_x / mgstruct->blocksize_x[0];
                curr_y_coarse = curr_y / mgstruct->blocksize_y[0]; 
                curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 

                // Update the norm!
                null_vectors[n-1][i] /= norms[curr_coarse];
            }
            
            for (m = 0; m < n; m++)
            {
                // Zero out the dot product. 
                
                zero<double>(dot_prod, coarse_size);
                // Pull off one piece at a time.
                for (i = 0; i < fine_size; i++)
                {
                    // What's the current coarse site? First, find the fine site.
                    curr_x = i % x_fine;
                    curr_y = i / x_fine; 

                    // Now, find the coarse site. 
                    curr_x_coarse = curr_x / mgstruct->blocksize_x[0];
                    curr_y_coarse = curr_y / mgstruct->blocksize_y[0]; 
                    curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 

                    dot_prod[curr_coarse] += conj(null_vectors[m][i])*null_vectors[n][i];
                }

            
                /*cout << "dot " << m << " ";
                for (i = 0; i < coarse_size; i++)
                {
                    cout << dot_prod[i] << " ";
                }
                cout << "\n";*/
                

                // Normalize the null_vectors.
                for (i = 0; i < fine_size; i++)
                {
                    // What's the current coarse site? First, find the fine site.
                    curr_x = i % x_fine;
                    curr_y = i / x_fine; 

                    // Now, find the coarse site. 
                    curr_x_coarse = curr_x / mgstruct->blocksize_x[0];
                    curr_y_coarse = curr_y / mgstruct->blocksize_y[0]; 
                    curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 

                    // Project off previous vector.
                    null_vectors[n][i] -= dot_prod[curr_coarse]*null_vectors[m][i];
                }
            }
        }

        delete[] norms; 
        delete[] dot_prod;
    }
    
    block_normalize(mgstruct); 
}

// Prolong a coarse vector to a fine vector using the info in mgstruct.
void prolong(complex<double>* vec_fine, complex<double>* vec_coarse, mg_operator_struct_complex* mgstruct)
{
    int n, i;
    int x_fine = mgstruct->curr_x_fine;
    int y_fine = mgstruct->curr_y_fine;
    int fine_size = x_fine*y_fine;
    int x_coarse = mgstruct->curr_x_coarse; // how many coarse sites are there in the x dir?
    //int y_coarse = mgstruct->curr_y_coarse; // how many coarse sites are there in the y dir?
    
    // Grab the current null vectors.
    complex<double>** null_vectors = mgstruct->null_vectors[mgstruct->curr_level];
    
    // Hold current sites.
    int curr_x, curr_y, curr_x_coarse, curr_y_coarse, curr_coarse; 
    
    zero<double>(vec_fine, fine_size); 
    
    // Loop over every vector.
    for (n = 0; n < mgstruct->n_vector; n++)
    {
        // Loop over the fine size.
        
        for (i = 0; i < fine_size; i++)
        {
            // What's the current coarse site? First, find the fine site.
            curr_x = i % x_fine;
            curr_y = i / x_fine; 
            
            // Now, find the coarse site. 
            curr_x_coarse = curr_x / mgstruct->blocksize_x[0];
            curr_y_coarse = curr_y / mgstruct->blocksize_y[0]; 
            curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 
            
            // Update the fine with the coarse. 
            vec_fine[i] += null_vectors[n][i]*vec_coarse[curr_coarse*mgstruct->n_vector+n];
        }
    }
}

// Restrict a fine vector to a coarse vector using the info in mgstruct.
// Hermitian conjugate of prolong. 
void restrict(complex<double>* vec_coarse, complex<double>* vec_fine, mg_operator_struct_complex* mgstruct)
{
    int n, i;
    int x_fine = mgstruct->curr_x_fine;
    int y_fine = mgstruct->curr_y_fine;
    int fine_size = x_fine*y_fine;
    int x_coarse = mgstruct->curr_x_coarse; // how many coarse sites are there in the x dir?
    int y_coarse = mgstruct->curr_y_coarse; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 
    
    // Grab the current null vectors.
    complex<double>** null_vectors = mgstruct->null_vectors[mgstruct->curr_level];
    
    // Hold current sites.
    int curr_x, curr_y, curr_x_coarse, curr_y_coarse, curr_coarse; 
    
    zero<double>(vec_coarse, mgstruct->n_vector*coarse_size); 
    
    // Loop over every vector.
    for (n = 0; n < mgstruct->n_vector; n++)
    {
        // Loop over the fine size.
        
        for (i = 0; i < fine_size; i++)
        {
            // What's the current coarse site? First, find the fine site.
            curr_x = i % x_fine;
            curr_y = i / x_fine; 
            
            // Now, find the coarse site. 
            curr_x_coarse = curr_x / mgstruct->blocksize_x[0];
            curr_y_coarse = curr_y / mgstruct->blocksize_y[0]; 
            curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 
            
            // Update the fine with the coarse. 
            vec_coarse[curr_coarse*mgstruct->n_vector+n] += conj(null_vectors[n][i])*vec_fine[i];
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
            mgstruct->curr_dof_fine *= mgstruct->n_vector;
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
            mgstruct->curr_dof_fine /= mgstruct->n_vector; 
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
void mg_preconditioner(complex<double>* lhs, complex<double>* rhs, int size, void* extra_data)
{
    cout << "[MG]: Entered mg_preconditioner.\n";
    mg_precond_struct_complex* mgprecond = (mg_precond_struct_complex*)extra_data; 
    
    // Standard defines.
    int x_fine = mgprecond->mgstruct->curr_x_fine;
    int y_fine = mgprecond->mgstruct->curr_y_fine;
    int fine_size = x_fine*y_fine;
    int x_coarse = mgprecond->mgstruct->curr_x_coarse; // how many coarse sites are there in the x dir?
    int y_coarse = mgprecond->mgstruct->curr_y_coarse; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 
    int coarse_length = coarse_size*mgprecond->mgstruct->n_vector; 
    
    // Store inversion info.
    inversion_info invif;
    
    // GET EXCITED! First off, let's do some pre-smoothing. 
    // 1. z_presmooth = smoothed rhs. 
    complex<double>* z_presmooth = new complex<double>[fine_size];
    zero<double>(z_presmooth, fine_size);
    if (mgprecond->n_pre_smooth > 0 && mgprecond->in_smooth_type != NONE)
    {
        switch (mgprecond->in_smooth_type)
        {
            case NONE: // can't reach here, anyway.
                break;
            case CG:
                invif = minv_vector_minres(z_presmooth, rhs, fine_size, mgprecond->n_pre_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
                break;
            case MINRES:
                invif = minv_vector_minres(z_presmooth, rhs, fine_size, mgprecond->n_pre_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
                break;
            case GCR:
                invif = minv_vector_gcr(z_presmooth, rhs, fine_size, mgprecond->n_pre_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
                break; 
        }
        printf("[L1 Presmooth]: Iterations %d Res %.8e Err N Algorithm %s\n", invif.iter, sqrt(invif.resSq), invif.name.c_str());
    }
    else
    {
        copy<double>(z_presmooth, rhs, fine_size);
    }
    
    if (mgprecond->in_solve_type != NONE)
    {
        // Compute updated residual. 
        // 2. r_pre = rhs - A z_presmooth
        complex<double>* Az_presmooth = new complex<double>[fine_size];
        zero<double>(Az_presmooth, fine_size);
        (*mgprecond->mgstruct->matrix_vector)(Az_presmooth, z_presmooth, mgprecond->mgstruct->matrix_extra_data);

        complex<double>* r_presmooth = new complex<double>[fine_size];
        for (int i = 0; i < fine_size; i++)
        {
            r_presmooth[i] = rhs[i] - Az_presmooth[i];
        }

        // Restrict r_presmooth. 
        // 3. rhs_coarse = restrict(r_presmooth)
        complex<double>* rhs_coarse = new complex<double>[coarse_length];
        zero<double>(rhs_coarse, coarse_length); 
        restrict(rhs_coarse, r_presmooth, mgprecond->mgstruct);

        // 4. Perform coarse solve.
        complex<double>* lhs_coarse = new complex<double>[coarse_length];
        zero<double>(lhs_coarse, coarse_length);
        switch (mgprecond->in_solve_type)
        {
            case NONE: // The code can't reach here, anyway.
                break;
            case MINRES:
                invif = minv_vector_minres(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_step, mgprecond->rel_res, mgprecond->matrix_vector, mgprecond->matrix_extra_data);
                break;
            case CG:
                invif = minv_vector_cg(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_step, mgprecond->rel_res, mgprecond->matrix_vector, mgprecond->matrix_extra_data);
                break;
            case GCR:
                invif = minv_vector_gcr(lhs_coarse, rhs_coarse, coarse_length, mgprecond->n_step, mgprecond->rel_res, mgprecond->matrix_vector, mgprecond->matrix_extra_data);
                break;
        }
        printf("[L2]: Iterations %d RelRes %.8e Err N Algorithm %s\n", invif.iter, sqrt(invif.resSq)/sqrt(norm2sq<double>(rhs_coarse, coarse_length)), invif.name.c_str());

        // Project the lhs to the fine vector.
        // 5. lhs_postsmooth = prolong(lhs_coarse)
        complex<double>* lhs_postsmooth = new complex<double>[fine_size];
        zero<double>(lhs_postsmooth, fine_size);
        prolong(lhs_postsmooth, lhs_coarse, mgprecond->mgstruct); 


        // Update the solution. 
        // 6. lhs = initial smooth (z_presmooth) + coarse solve (lhs_postsmooth)
        for (int i = 0; i < fine_size; i++)
        {
            lhs[i] = z_presmooth[i] + lhs_postsmooth[i];
        }
        
        delete[] Az_presmooth;
        delete[] r_presmooth; 
        delete[] rhs_coarse;
        delete[] lhs_coarse; 
        delete[] lhs_postsmooth;
    }
    else // no inner solver
    {
        copy<double>(lhs, z_presmooth, fine_size); 
    }
    
    // Almost done! Do some post-smoothing.
    // 7. Post smooth.
    if (mgprecond->n_post_smooth > 0 && mgprecond->in_smooth_type != NONE)
    {
        switch (mgprecond->in_smooth_type)
        {
            case NONE: // Can't reach here, anyway.
                break;
            case CG:
                invif = minv_vector_cg(lhs, rhs, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
                break;
            case MINRES:
                invif = minv_vector_minres(lhs, rhs, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
                break;
            case GCR:
                invif = minv_vector_gcr(lhs, rhs, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
                break;
        }
        //invif = minv_vector_minres(lhs, lhs_postsmooth, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
        printf("[L1 Postsmooth]: Iterations %d Res %.8e Err N Algorithm %s\n", invif.iter, sqrt(invif.resSq), invif.name.c_str());
    }
    
    // Clean up!
    delete[] z_presmooth;
    
    cout << "[MG]: Exited mg_preconditioner.\n";
    
}

