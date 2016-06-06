
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
#include "mg_real.h"

// General multigrid projector function!
void coarse_square_op(double* lhs, double* rhs, void* extra_data)
{
    // Iterators.
    int i, j, k, n; 
    int tmp; 
    
    // Grab the mg_precond_struct.
    mg_operator_struct_real* mgstruct = (mg_operator_struct_real*)extra_data; 
    
    // lhs and rhs are of size coarse_size. mgstruct.matrix_vector expects
    // fine_size. 
    int x_fine = mgstruct->x_fine;
    int y_fine = mgstruct->y_fine;
    int fine_size = x_fine*y_fine;
    int x_coarse = x_fine/mgstruct->blocksize_x; // how many coarse sites are there in the x dir?
    int y_coarse = y_fine/mgstruct->blocksize_y; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 
    int coarse_length = coarse_size*mgstruct->n_vector;
    
    // Okay... how the hell are we going to do this. 
    double* Px; // Holds prolonged current solution.
    double* APx; // Holds A times prolonged current solution.
    
    Px = new double[fine_size];
    APx = new double[fine_size];
    
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

// Properly normalize the P vectors.
void block_normalize(mg_operator_struct_real* mgstruct)
{
    int n, i;
    int x_fine = mgstruct->x_fine;
    int y_fine = mgstruct->y_fine;
    int fine_size = x_fine*y_fine;
    int x_coarse = x_fine/mgstruct->blocksize_x; // how many coarse sites are there in the x dir?
    int y_coarse = y_fine/mgstruct->blocksize_y; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 
    
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
            curr_x_coarse = curr_x / mgstruct->blocksize_x;
            curr_y_coarse = curr_y / mgstruct->blocksize_y; 
            curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 
            
            // Update the norm!
            norms[curr_coarse] += mgstruct->projectors[n][i]*mgstruct->projectors[n][i];
        }
        
        // Sqrt all of the norms.
        for (i = 0; i < coarse_size; i++)
        {
            norms[i] = sqrt(norms[i]);
        }
        
        // Normalize the projectors.
        for (i = 0; i < fine_size; i++)
        {
            // What's the current coarse site? First, find the fine site.
            curr_x = i % x_fine;
            curr_y = i / x_fine; 
            
            // Now, find the coarse site. 
            curr_x_coarse = curr_x / mgstruct->blocksize_x;
            curr_y_coarse = curr_y / mgstruct->blocksize_y; 
            curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 
            
            // Update the norm!
            mgstruct->projectors[n][i] /= norms[curr_coarse];
        }
    }
    
    delete[] norms; 
}

// Prolong a coarse vector to a fine vector using the info in mgstruct.
void prolong(double* vec_fine, double* vec_coarse, mg_operator_struct_real* mgstruct)
{
    int n, i;
    int x_fine = mgstruct->x_fine;
    int y_fine = mgstruct->y_fine;
    int fine_size = x_fine*y_fine;
    int x_coarse = x_fine/mgstruct->blocksize_x; // how many coarse sites are there in the x dir?
    int y_coarse = y_fine/mgstruct->blocksize_y; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 
    
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
            curr_x_coarse = curr_x / mgstruct->blocksize_x;
            curr_y_coarse = curr_y / mgstruct->blocksize_y; 
            curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 
            
            // Update the fine with the coarse. 
            vec_fine[i] += mgstruct->projectors[n][i]*vec_coarse[curr_coarse*mgstruct->n_vector+n];
        }
    }
}

// Restrict a fine vector to a coarse vector using the info in mgstruct.
void restrict(double* vec_coarse, double* vec_fine, mg_operator_struct_real* mgstruct)
{
    int n, i;
    int x_fine = mgstruct->x_fine;
    int y_fine = mgstruct->y_fine;
    int fine_size = x_fine*y_fine;
    int x_coarse = x_fine/mgstruct->blocksize_x; // how many coarse sites are there in the x dir?
    int y_coarse = y_fine/mgstruct->blocksize_y; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 
    
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
            curr_x_coarse = curr_x / mgstruct->blocksize_x;
            curr_y_coarse = curr_y / mgstruct->blocksize_y; 
            curr_coarse = curr_y_coarse*x_coarse + curr_x_coarse; 
            
            // Update the fine with the coarse. 
            vec_coarse[curr_coarse*mgstruct->n_vector+n] += mgstruct->projectors[n][i]*vec_fine[i];
        }
    }
}


// MG preconditioner!! (Man, I'm excited!
void mg_preconditioner(double* lhs, double* rhs, int size, void* extra_data)
{
    cout << "Entered mg_preconditioner.\n";
    mg_precond_struct_real* mgprecond = (mg_precond_struct_real*)extra_data; 
    
    // Standard defines.
    int x_fine = mgprecond->mgstruct->x_fine;
    int y_fine = mgprecond->mgstruct->y_fine;
    int fine_size = x_fine*y_fine;
    int x_coarse = x_fine/mgprecond->mgstruct->blocksize_x; // how many coarse sites are there in the x dir?
    int y_coarse = y_fine/mgprecond->mgstruct->blocksize_y; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 
    int coarse_length = coarse_size*mgprecond->mgstruct->n_vector; 
    
    // Store inversion info.
    inversion_info invif;
    
    // GET EXCITED! First off, let's do some pre-smoothing. 
    // 1. z_presmooth = smoothed rhs. 
    double* z_presmooth = new double[fine_size];
    zero<double>(z_presmooth, fine_size);
    if (mgprecond->n_pre_smooth > 0)
    {
        invif = minv_vector_minres(z_presmooth, rhs, fine_size, mgprecond->n_pre_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
        printf("Presmooth: Algorithm %s took %d iterations to reach a residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq));
    }
    else
    {
        copy<double>(z_presmooth, rhs, fine_size);
    }
    
    // Compute updated residual. 
    // 2. r_pre = rhs - A z_presmooth
    double* Az_presmooth = new double[fine_size];
    zero<double>(Az_presmooth, fine_size);
    (*mgprecond->mgstruct->matrix_vector)(Az_presmooth, z_presmooth, mgprecond->mgstruct->matrix_extra_data);
    
    double* r_presmooth = new double[fine_size];
    for (int i = 0; i < fine_size; i++)
    {
        r_presmooth[i] = rhs[i] - Az_presmooth[i];
    }
    
    // Restrict r_presmooth. 
    // 3. rhs_coarse = restrict(r_presmooth)
    double* rhs_coarse = new double[coarse_length];
    zero<double>(rhs_coarse, coarse_length); 
    restrict(rhs_coarse, r_presmooth, mgprecond->mgstruct);
    
    // 4. Perform coarse solve.
    double* lhs_coarse = new double[coarse_length];
    zero<double>(lhs_coarse, coarse_length);
    switch (mgprecond->in_solve_type)
    {
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
    printf("Coarse solve: Algorithm %s took %d iterations to reach a relative residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq)/sqrt(norm2sq<double>(rhs_coarse, coarse_length)));
    
    // Project the lhs to the fine vector.
    // 5. lhs_postsmooth = prolong(lhs_coarse)
    double* lhs_postsmooth = new double[fine_size];
    zero<double>(lhs_postsmooth, fine_size);
    prolong(lhs_postsmooth, lhs_coarse, mgprecond->mgstruct); 
    
    // Update the solution. 
    // 6. lhs = initial smooth (z_presmooth) + coarse solve (lhs_postsmooth)
    for (int i = 0; i < fine_size; i++)
    {
        lhs[i] = z_presmooth[i] + lhs_postsmooth[i];
    }
    
    // Almost done! Do some post-smoothing.
    // 7. Post smooth.
    if (mgprecond->n_post_smooth > 0)
    {
        invif = minv_vector_minres(lhs, rhs, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
        //invif = minv_vector_minres(lhs, lhs_postsmooth, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
        printf("Postsmooth: Algorithm %s took %d iterations to reach a residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq));
    }
    
    /*
    // Allocate the coarse rhs, lhs. 
    double* rhs_coarse = new double[coarse_length];
    double* lhs_coarse = new double[coarse_length];
    zero<double>(rhs_coarse, coarse_length);
    zero<double>(lhs_coarse, coarse_length);
    
    // We need to restrict A*(rhs_smooth) - rhs.
    double* tmp1 = new double[fine_size];
    double* tmp2 = new double[fine_size];
    // tmp1 = A*rhs_smooth
    (*mgprecond->mgstruct->matrix_vector)(tmp1, rhs_presmooth, mgprecond->mgstruct->matrix_extra_data);
    for (int i = 0; i < fine_size; i++)
    {
        tmp2[i] = rhs_presmooth[i];
        rhs_presmooth[i] = rhs[i] - tmp1[i];
        
    }
    // Now, rhs_presmooth = (A r_smooth - rhs). 
    // tmp2 = r_smooth. 
    
    // Restrict the rhs to the coarse vector. 
    restrict(rhs_coarse, rhs_presmooth, mgprecond->mgstruct);
    
    // Solve the coarse system.
    switch (mgprecond->in_solve_type)
    {
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
    
    printf("Coarse solve: Algorithm %s took %d iterations to reach a relative residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq)/sqrt(norm2sq<double>(rhs_coarse, coarse_length)));
    
    // Project the lhs to the fine vector.
    double* lhs_postsmooth = new double[fine_size];
    zero<double>(lhs_postsmooth, fine_size);
    prolong(lhs_postsmooth, lhs_coarse, mgprecond->mgstruct); 
    
    // We want to post-smooth, using
    // tmp2 - lhs_smooth as the initial guess,
    // rhs as the right hand side. 
    for (int i = 0; i < fine_size; i++)
    {
        lhs[i] = tmp2[i] + lhs_postsmooth[i];
    }
    
    // Almost done! Do some post-smoothing.
    if (mgprecond->n_post_smooth > 0)
    {
        invif = minv_vector_minres(lhs, rhs, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
        //invif = minv_vector_minres(lhs, lhs_postsmooth, fine_size, mgprecond->n_post_smooth, 1e-20, mgprecond->mgstruct->matrix_vector, mgprecond->mgstruct->matrix_extra_data); 
        printf("Postsmooth: Algorithm %s took %d iterations to reach a residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq));
    }
    else
    {
        // We're good.
        //copy<double>(lhs, lhs_postsmooth, fine_size);
    }*/
    
    
    // Clean up!
    delete[] z_presmooth;
    delete[] Az_presmooth;
    delete[] r_presmooth; 
    delete[] rhs_coarse;
    delete[] lhs_coarse; 
    delete[] lhs_postsmooth; 
    
    cout << "Exited mg_preconditioner.\n";
    
}

