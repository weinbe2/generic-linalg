// Sun Sep 18 11:44:35 EDT 2016
// Evan S Weinberg
// C++ file for the Block CG inverter.

// To do:
// 1. Template to support float, double. 

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_traits.h"
#include "generic_vector.h"

#include "generic_block_cg_mrhs.h"
#include "generic_gelim.h"

using namespace std;

// Solves multiple rhs using BlockCG: http://www.row1.ca/s/pdfs/courses/BlockCG.pdf
// "The Block Conjugate Gradient for Multiple Right Hand Sides in a Direct Current Resistivity Inversion"
inversion_info minv_vector_block_cg_mrhs(double **phi, double **phi0, int n_rhs, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{

  // Initialize vectors.
  double **r, **p, **Ap;
  double **lambda, **psi;
  double **rsq, **rsqNew, **pAp; 
  double *bsqrt; 
  double truersq; // Temporarily hold rsq values. 
  double tmp; 
  int tmp_int; 
  double* tmp_ptr; 
  
  int* mapping; // holds the mapping between vectors and the original vector ordering.
                // we do this because we implement the unequal convergence variant. (6.3)
  int n_rem = n_rhs; // counts how many systems aren't converged. The first 
                     // (n_rem) vectors in any array still need more iterations.
                     // n_rem stands for "number remaining"
  
  int k,i,j,l; // vector element iterators
  int n,m; // vector iterators
  inversion_info invif(n_rhs); 

  // Allocate memory.
  
  // Vectors. Initialize.
  r = new double*[n_rhs];
  p = new double*[n_rhs];
  Ap = new double*[n_rhs];
  for (n = 0; n < n_rhs; n++)
  {
    r[n] = new double[size]; zero<double>(r[n], size);
    p[n] = new double[size]; zero<double>(p[n], size);
    Ap[n] = new double[size]; zero<double>(Ap[n], size);
  }
  
  // Generalized alpha, beta, r. Initialize.
  pAp = new double*[n_rhs];
  lambda = new double*[n_rhs];
  psi = new double*[n_rhs];
  rsq = new double*[n_rhs];
  rsqNew = new double*[n_rhs];
  
  for (n = 0; n < n_rhs; n++)
  {
    pAp[n] = new double[n_rhs]; zero<double>(pAp[n], n_rhs);
    lambda[n] = new double[n_rhs]; zero<double>(lambda[n], n_rhs);
    psi[n] = new double[n_rhs]; zero<double>(psi[n], n_rhs);
    rsq[n] = new double[n_rhs]; zero<double>(rsq[n], n_rhs);
    rsqNew[n] = new double[n_rhs]; zero<double>(rsqNew[n], n_rhs);
  }
  
  // Prepare space for original norms.
  bsqrt = new double[n_rhs];

  // Initialize the mapping.
  mapping = new int[n_rhs];
  for (n = 0; n < n_rhs; n++)
  {
    mapping[n] = n; // All vectors are currently in order.
  }
  
  // End allocating memory.
  
  // Find norm of rhs.
  for (n = 0; n < n_rhs; n++)
  {
    bsqrt[n] = sqrt(norm2sq<double>(phi0[n], size));
  }
  
  // 1. Compute r = b - Ax
  for (n = 0; n < n_rhs; n++)
  {
    (*matrix_vector)(p[n], phi[n], extra_info); invif.ops_count++;
    for (i = 0; i < size; i++)
    {
      r[n][i] = phi0[n][i] - p[n][i];
    }
  }

  // Compute rsq matrix.
  for (n = 0; n < n_rhs; n++) // row index
  {
    rsq[n][n] = norm2sq<double>(r[n], size);
    rsqNew[n][n] = rsq[n][n]; // This is for convenience later.
    for (m = n+1; m < n_rhs; m++) // column index
    {
      rsq[n][m] = dot<double>(r[n], r[m], size); // 
      rsq[m][n] = rsq[n][m]; // needs a conj for complex case.
      
      rsqNew[n][m] = rsq[n][m];
      rsqNew[m][n] = rsq[m][n];
    }
  }
  
  for (n = 0; n < n_rem; n++)
  {
    if (sqrt(rsq[n][n])/bsqrt[n] < eps) // Check if this solution has converged.
    {
      // If so, permute it away. Need to update x, p, Ap, rsq, bnorm.
      n_rem--;
      
      // Pay attention to what swap we did.
      mapping[n] = n_rem;
      mapping[n_rem] = n;
      
      // Permute x.
      tmp_ptr = phi[n_rem];
      phi[n_rem] = phi[n];
      phi[n] = tmp_ptr;
      
      // Permute p.
      tmp_ptr = p[n_rem];
      p[n_rem] = p[n];
      p[n] = tmp_ptr;
      
      // Permute Ap.
      tmp_ptr = Ap[n_rem];
      Ap[n_rem] = Ap[n];
      Ap[n] = tmp_ptr;
      
      // Permute rsq, rows first.
      tmp_ptr = rsq[n_rem];
      rsq[n_rem] = rsq[n];
      rsq[n] = tmp_ptr;
      
      // Now columns.
      for (m = 0; m <= n_rem; m++)
      {
        tmp = rsq[m][n_rem];
        rsq[m][n_rem] = rsq[m][n];
        rsq[m][n] = tmp;
      }
      
      // Permute rsqNew, rows first.
      tmp_ptr = rsqNew[n_rem];
      rsqNew[n_rem] = rsqNew[n];
      rsqNew[n] = tmp_ptr;
      
      // Now columns.
      for (m = 0; m <= n_rem; m++)
      {
        tmp = rsqNew[m][n_rem];
        rsqNew[m][n_rem] = rsqNew[m][n];
        rsqNew[m][n] = tmp;
      }
      
      // Permute bnorm.
      tmp = bsqrt[n_rem];
      bsqrt[n_rem] = bsqrt[n];
      bsqrt[n] = tmp; 
      
      // We swapped with the end, so we need to recheck the end.
      n--; 
    }
  }
  
  // In the off chance n_rem already equals 0, we're done!
  if (n_rem > 0)
  {
    
    // 2. p_0 = -r_0.
    for (n = 0; n < n_rem; n++)
    {
      for (i = 0; i < size; i++)
      {
        p[n][i] = -r[n][i];
      }
    }

    // Compute Ap.
    for (n = 0; n < n_rem; n++)
    {
      zero<double>(Ap[n], size);
      (*matrix_vector)(Ap[n], p[n], extra_info); invif.ops_count++;
    }


    // Compute pAp matrix.
    for (n = 0; n < n_rem; n++) // row index
    {
      pAp[n][n] = dot<double>(p[n], Ap[n], size);
      for (m = n+1; m < n_rem; m++) // column index
      {
        pAp[n][m] = dot<double>(p[n], Ap[m], size); // 
        pAp[m][n] = pAp[n][m]; // needs a conj for complex case.
      }
    }

    /*
    // We need pAp inverse (and never pAp again), so do the in-place inverse.
    cout << "pAp\n";
    cout << pAp[0][0] << " " << pAp[1][0] << "\n";
    cout << pAp[0][1] << " " << pAp[1][1] << "\n\n";
    gaussian_elimination_matrix_inverse(pAp, pAp, n_rem);
    cout << pAp[0][0] << " " << pAp[1][0] << "\n";
    cout << pAp[0][1] << " " << pAp[1][1] << "\n\n";
    */
    
    // iterate till convergence
    for(k = 0; k< max_iter; k++) {
      

      // Lambda = (pAp)^(-1) rsq. (Not anymore) pAp contains the inverse already.
      
      // Technically I need to transpose rsq, but it's already symmetric.
      gaussian_elimination_multi_rhs(lambda, rsq, pAp, n_rhs, n_rhs);
      
      // Transpose Lambda.
      for (i = 0; i < n_rhs-1; i++)
      {
        for (j = i+1; j < n_rhs; j++)
        {
          tmp = lambda[i][j];
          lambda[i][j] = lambda[j][i];
          lambda[j][i] = tmp;
        }
      }
      
      /*for (n = 0; n < n_rem; n++)
      {
        for (m = 0; m < n_rem; m++)
        {
          lambda[m][n] = 0.0;
          for (l = 0; l < n_rem; l++)
          {
            lambda[m][n] += pAp[m][l]*rsq[l][n];
          }
        }
      }*/

      // Update remaining X = X + P*Lambda, R = R + AP * Lambda.
      // Keep Lambda in cache as much as possible.
      for (i = 0; i < size; i++) 
      {
        for (n = 0; n < n_rem; n++)
        {
          for (m = 0; m < n_rem; m++)
          {
            phi[n][i] = phi[n][i] + p[m][i]*lambda[m][n];
            r[n][i] = r[n][i] + Ap[m][i]*lambda[m][n];
          }
        }
      }
      
      // Update residuals matrix.
      for (n = 0; n < n_rhs; n++) // row index
      {
        rsqNew[n][n] = norm2sq<double>(r[n], size);
        for (m = n+1; m < n_rhs; m++) // column index
        {
          rsqNew[n][m] = dot<double>(r[n], r[m], size); // 
          rsqNew[m][n] = rsqNew[n][m]; // needs a conj for complex case.
        }
      }
      
      // Need to create a new verbosity routine.
      //print_verbosity_resid(verb, "CG", k+1, invif.ops_count, sqrt(rsqNew)/bsqrt); 
      if (verb != 0)
      {
        if (verb->verbosity == VERB_DETAIL)
        {
          std::cout << verb->verb_prefix << "BlockCG " << " Iter " << k+1 << " Ops " << invif.ops_count << " RelRes ";
          for (n = 0; n < n_rhs; n++)
          {
            std::cout << sqrt(rsqNew[mapping[n]][mapping[n]])/bsqrt[mapping[n]] << " ";
          }
          
          std::cout << "\n";
        }
      }
      
      // Check residuals. Permute solutions away if we can. 
      for (n = 0; n < n_rem; n++)
      {
        if (sqrt(rsqNew[n][n])/bsqrt[n] < eps) // Check if this solution has converged.
        {
          // If so, permute it away. Need to update x, p, Ap, rsq, rsqNew. 
          // Don't need to worry about pAp, Lambda. 
          n_rem--;

          // Pay attention to what swap we did.
          mapping[n] = n_rem;
          mapping[n_rem] = n;

          // Permute x.
          tmp_ptr = phi[n_rem];
          phi[n_rem] = phi[n];
          phi[n] = tmp_ptr;

          // Permute p.
          tmp_ptr = p[n_rem];
          p[n_rem] = p[n];
          p[n] = tmp_ptr;

          // Permute Ap.
          tmp_ptr = Ap[n_rem];
          Ap[n_rem] = Ap[n];
          Ap[n] = tmp_ptr;

          // Permute rsq, rows first.
          tmp_ptr = rsq[n_rem];
          rsq[n_rem] = rsq[n];
          rsq[n] = tmp_ptr;

          // Now columns.
          for (m = 0; m <= n_rem; m++)
          {
            tmp = rsq[m][n_rem];
            rsq[m][n_rem] = rsq[m][n];
            rsq[m][n] = tmp;
          }
          
          // Permute rsqNew, rows first.
          tmp_ptr = rsqNew[n_rem];
          rsqNew[n_rem] = rsqNew[n];
          rsqNew[n] = tmp_ptr;

          // Now columns.
          for (m = 0; m <= n_rem; m++)
          {
            tmp = rsqNew[m][n_rem];
            rsqNew[m][n_rem] = rsqNew[m][n];
            rsqNew[m][n] = tmp;
          }

          
          // Permute bnorm.
          tmp = bsqrt[n_rem];
          bsqrt[n_rem] = bsqrt[n];
          bsqrt[n] = tmp; 

          // We swapped with the end, so we need to recheck the end.
          n--; 
        }
      }

      cout << "n_rem = " << n_rem << "\n\n"; 
      // Exit if all residuals are small enough.
      if (n_rem == 0 || k == max_iter-1)
      {
        break;
      }
      
      // Compute rsq^(-1) in place.
      /*cout << "rsq\n";
      cout << rsq[0][0] << " " << rsq[1][0] << "\n";
      cout << rsq[0][1] << " " << rsq[1][1] << "\n\n";
      gaussian_elimination_matrix_inverse(rsq, rsq, n_rem);
      cout << rsq[0][0] << " " << rsq[1][0] << "\n";
      cout << rsq[0][1] << " " << rsq[1][1] << "\n\n";*/
      
      // Psi = (rsq)^(-1) rsqNew. 
      
      
      // Technically I need to transpose rsq, but it's already symmetric.
      gaussian_elimination_multi_rhs(psi, rsqNew, rsq, n_rhs, n_rhs);
      
      // Transpose Lambda.
      for (i = 0; i < n_rhs-1; i++)
      {
        for (j = i+1; j < n_rhs; j++)
        {
          tmp = psi[i][j];
          psi[i][j] = psi[j][i];
          psi[j][i] = tmp;
        }
      }
      
      /*
      for (n = 0; n < n_rem; n++)
      {
        for (m = 0; m < n_rem; m++)
        {
          psi[n][m] = 0.0;
          for (l = 0; l < n_rem; l++)
          {
            psi[n][m] += rsq[n][l]*rsqNew[l][m];
          }
        }
      }*/
      
      
      // We're done with rsq. Copy rsqNew in.
      for (n = 0; n < n_rem; n++)
      {
        for (m = 0; m < n_rem; m++)
        {
          rsq[n][m] = rsqNew[n][m];
        }
      }
      
      // Update remaining P = -R + P*Psi
      // Keep Lambda in cache as much as possible.
      for (i = 0; i < size; i++) 
      {
        for (n = 0; n < n_rem; n++)
        {
          for (m = 0; m < n_rem; m++)
          {
            p[n][i] = -r[n][i] + p[m][i]*psi[m][n];
          }
        }
      }

      // Compute Ap.
      for (n = 0; n < n_rem; n++)
      {
        zero<double>(Ap[n], size);
        (*matrix_vector)(Ap[n], p[n], extra_info); invif.ops_count++;
      }
      
      // Compute updated pAp matrix.
      for (n = 0; n < n_rem; n++) // row index
      {
        pAp[n][n] = dot<double>(p[n], Ap[n], size);
        for (m = n+1; m < n_rem; m++) // column index
        {
          pAp[n][m] = dot<double>(p[n], Ap[m], size); // 
          pAp[m][n] = pAp[n][m]; // needs a conj for complex case.
        }
      }
      
      
      // We need pAp inverse (and never pAp again), so do the in-place inverse.
      /*cout << "pAp\n";
      cout << pAp[0][0] << " " << pAp[1][0] << "\n";
      cout << pAp[0][1] << " " << pAp[1][1] << "\n\n";
      gaussian_elimination_matrix_inverse(pAp, pAp, n_rem);
      cout << pAp[0][0] << " " << pAp[1][0] << "\n";
      cout << pAp[0][1] << " " << pAp[1][1] << "\n\n";*/
    } 
  } // end check n_rem > 0;
    
  if(k == max_iter-1) {
    invif.success = false;
  }
  else
  {
     invif.success = true;
  }
  k++;
  
  // Undo the permutation damage I did.
  // Only need to permute solutions phi.
  for (n = 0; n < n_rhs; n++)
  {
    // Find the true n'th vector.
    if (mapping[n] != n)
    {
      for (m = n+1; m < n_rhs; m++)
      {
        if (mapping[m] == n) // Match, swap!
        {
          tmp_ptr = phi[m];
          phi[m] = phi[n];
          phi[n] = tmp_ptr;
          
          mapping[m] = mapping[n];
          mapping[n] = n;
          
          break;
        }
      } 
    }
  }
  
  // Calculate explicit rsqs.
  for (n = 0; n < n_rhs; n++)
  {
    (*matrix_vector)(Ap[n],phi[n],extra_info); invif.ops_count++;
    invif.resSqmrhs[n] = diffnorm2sq<double>(Ap[n], phi0[n], size);
  }
  
  // Need to update verbosity type things.
  //print_verbosity_summary(verb, "CG", invif.success, k, invif.ops_count, sqrt(truersq)/bsqrt);
  
  if (verb != 0)
  {
    if (verb->verbosity == VERB_SUMMARY || verb->verbosity == VERB_RESTART_DETAIL || verb->verbosity == VERB_DETAIL)
    {
      std::cout << verb->verb_prefix << "BlockCG " << " Success " << (invif.success ? "Y" : "N") << " Iter " << k+1 << " Ops " << invif.ops_count << " RelRes ";
      for (n = 0; n < n_rhs; n++)
      {
        std::cout << sqrt(invif.resSqmrhs[n])/bsqrt[mapping[n]] << " ";
      }

      std::cout << "\n";
    }
  }
  
  // Free all the things!
  
  // Vectors.
  for (n = 0; n < n_rhs; n++)
  {
    delete[] r[n];
    delete[] p[n];
    delete[] Ap[n];
  }
  delete[] r;
  delete[] p;
  delete[] Ap;
  
  // Generalized alpha, beta, r.
  for (n = 0; n < n_rhs; n++)
  {
    delete[] pAp[n];
    delete[] lambda[n];
    delete[] psi[n];
    delete[] rsq[n];
    delete[] rsqNew[n];
  }
  delete[] pAp;
  delete[] lambda;
  delete[] psi;
  delete[] rsq;
  delete[] rsqNew;
  
  // Space for original norms.
  delete[] bsqrt;

  // Delete the mapping.
  delete[] mapping;
  
  // End freeing all the things!
  
  invif.iter = k;
  invif.name = "BlockCG";
  return invif; // Convergence 
} 

/*
// Performs CG(restart_freq) with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_cg_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  int ops_count; 
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  stringstream ss;
  ss << "CG(" << restart_freq << ")";
  
  iter = 0; ops_count = 0; 
  do
  {
    invif = minv_vector_cg(phi, phi0, size, restart_freq, res, matrix_vector, extra_info, &verb_rest);
    iter += invif.iter;
    ops_count += invif.ops_count; 
    
    print_verbosity_restart(verb, ss.str(), iter, ops_count, sqrt(invif.resSq)/bsqrt);
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq)/bsqrt > res);
  
  invif.iter = iter; invif.ops_count = ops_count; 
  
  print_verbosity_summary(verb, ss.str(), invif.success, iter, invif.ops_count, sqrt(invif.resSq)/bsqrt);
  
  invif.name = ss.str();
  // invif.resSq is good.
  if (sqrt(invif.resSq)/bsqrt > res)
  {
    invif.success = false;
  }
  else
  {
    invif.success = true;
  }
  
  return invif;
}


inversion_info minv_vector_cg(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{
// CG solutions to Mphi = b 
//  see http://en.wikipedia.org/wiki/Conjugate_gradient_method

  
  // Initialize vectors.
  complex<double> *r, *p, *Ap;
  complex<double> alpha, beta, denom;
  double rsq, rsqNew, bsqrt, truersq;
  int k,i;
  inversion_info invif;

  // Allocate memory.
  r = new complex<double>[size];
  p = new complex<double>[size];
  Ap = new complex<double>[size];

  // Initialize values.
  rsq = 0.0; rsqNew = 0.0; bsqrt = 0.0; truersq = 0.0; k=0;

  // Zero vectors;
  zero<double>(r, size); 
  zero<double>(p, size); zero<double>(Ap, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // 1. Compute r = b - Ax
  (*matrix_vector)(p, phi, extra_info); invif.ops_count++;
  for (i = 0; i < size; i++)
  {
    r[i] = phi0[i] - p[i];
  }
  
  // 2. p_0 = r_0.
  copy<double>(p, r, size);
  
  // Compute Ap.
  zero<double>(Ap, size);
  (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;
  
  // Compute rsq.
  rsq = norm2sq<double>(r, size);

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    alpha = rsq/dot<double>(p, Ap, size);

    for (i = 0; i < size; i++)
    {
      phi[i] = phi[i] + alpha*p[i];
      r[i] = r[i] - alpha*Ap[i];
    }
    
    // Exit if new residual is small enough
    rsqNew = norm2sq<double>(r, size);
      
    print_verbosity_resid(verb, "CG", k+1, invif.ops_count, sqrt(rsqNew)/bsqrt);

    if (sqrt(rsqNew) < eps*bsqrt || k == max_iter - 1) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
  
    // Update vec using new residual
    beta = rsqNew / rsq;
    rsq = rsqNew; 
    
    for (i = 0; i < size; i++) {
      p[i] = r[i] + beta * p[i];
    }
    
    // Compute the new Ap.
    (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;
  } 
    
  if(k == max_iter-1) {
    invif.success = false;
  }
  else
  {
     invif.success = true;
  }
	
  k++; 
  
  (*matrix_vector)(Ap,phi,extra_info); invif.ops_count++;
  truersq = diffnorm2sq<double>(Ap, phi0, size);
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;

  
  print_verbosity_summary(verb, "CG", invif.success, k, invif.ops_count, sqrt(truersq)/bsqrt);
  
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "CG";
  return invif; // Convergence 
} 


// Performs CG(restart_freq) with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_cg_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  int iter; // counts total number of iterations.
  int ops_count; 
  inversion_info invif;
  double bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  inversion_verbose_struct verb_rest;
  shuffle_verbosity_restart(&verb_rest, verb);
  
  stringstream ss;
  ss << "CG(" << restart_freq << ")";
  
  iter = 0; ops_count = 0; 
  do
  {
    invif = minv_vector_cg(phi, phi0, size, restart_freq, res, matrix_vector, extra_info, &verb_rest);
    iter += invif.iter;
    ops_count += invif.ops_count; 
    
    print_verbosity_restart(verb, ss.str(), iter, ops_count, sqrt(invif.resSq)/bsqrt);
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq)/bsqrt > res);
  
  invif.iter = iter; invif.ops_count = ops_count; 
  
  print_verbosity_summary(verb, ss.str(), invif.success, iter, invif.ops_count, sqrt(invif.resSq)/bsqrt);
  
  invif.name = ss.str();
  // invif.resSq is good.
  if (sqrt(invif.resSq)/bsqrt > res)
  {
    invif.success = false;
  }
  else
  {
    invif.success = true;
  }
  
  return invif;
}
*/
