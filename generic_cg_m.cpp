// Fri Apr 29 11:54:22 EDT 2016
// Evan S Weinberg
// C++ file for CG-M inverter.

// To do:
// 1. Template to support float, double. 

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_vector.h"

#include "generic_cg_m.h"

using namespace std;

// Solves lhs = A^(-1) rhs using multishift CG as defined in http://arxiv.org/pdf/hep-lat/9612014.pdf
// Assumes there are n_shift values in "shifts", and that they are sorted.
// resid_freq_check is how often to check the residual of other solutions. This lets us stop iterating on converged systems. 
inversion_info minv_vector_cg_m(double **phi, double *phi0, int n_shift, int size, int resid_freq_check, int max_iter, double eps, double* shifts, void (*matrix_vector)(double*,double*,void*), void* extra_info, inversion_verbose_struct* verb)
{
  
  // Initialize vectors.
  double *r, *p, *Ap;
  double **p_s;
  double alpha, beta, beta_prev, rsq, rsqNew, bsqrt, truersq, tmp; 
  double *alpha_s, *beta_s, *zeta_s, *zeta_s_prev;
  int k,i,n;
  int n_shift_rem = n_shift; // number of systems to still iterate on. 
  int* mapping; // holds the mapping between vectors and the original vector ordering.
                // this is because some vectors may converge before others. 
  double* tmp_ptr; // temporary pointer for swaps.
  double tmp_dbl; // temporary double for swaps. 
  
  // Prepare an inversion_info for multiple residuals.
  inversion_info invif(n_shift); 

  // Allocate memory.
  alpha_s = new double[n_shift];
  beta_s = new double[n_shift];
  zeta_s = new double[n_shift];
  zeta_s_prev = new double[n_shift];
  
  p_s = new double*[n_shift];
  for (n = 0; n < n_shift; n++)
  {
    p_s[n] = new double[size];
  }
  
  r = new double[size];
  p = new double[size];
  Ap = new double[size];
  
  // Initialize mapping.
  mapping = new int[n_shift];
  for (n = 0; n < n_shift; n++)
  {
    mapping[n] = n; // All vectors are currently in order.
  }

  // Initialize values.
  rsq = 0.0; rsqNew = 0.0; bsqrt = 0.0; truersq = 0.0; k=0;
  for (n = 0; n < n_shift; n++)
  {
    // beta_0, zeta_0, zeta_-1
    beta_s[n] = zeta_s[n] = zeta_s_prev[n] = 1.0;
    // alpha_0. 
    alpha_s[n] = 0.0;
  }
  beta = 1.0; alpha = 0.0;

  // Zero vectors;
  zero<double>(r, size); 
  zero<double>(p, size); zero<double>(Ap, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // There can't be an initial guess... though it is sort of possible, in reference to:
  // http://arxiv.org/pdf/0810.1081v1.pdf
  
  // 1. x_sigma = 0, r = p_sigma = b.
  for (n = 0; n < n_shift; n++)
  {
    copy<double>(p_s[n], phi0, size);
    zero<double>(phi[n], size);
  }
  copy<double>(p, phi0, size);
  copy<double>(r, phi0, size);
  
  // Compute Ap.
  zero<double>(Ap, size);
  (*matrix_vector)(Ap, p, extra_info); invif.ops_count++;
  
  // Compute rsq.
  rsq = norm2sq<double>(r, size);

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // 2. beta_i = - rsq / pAp. Which is a weird switch from the normal notation, but whatever.
    beta_prev = beta; 
    beta = -rsq/dot<double>(p, Ap, size);
    //cout << "beta = " << beta << "\n";
    
    for (n = 0; n < n_shift_rem; n++)
    {
      // 3. Calculate beta_i^sigma, zeta_i+1^sigma according to 2.42 to 2.44.
      // zeta_{i+1}^sigma = complicated...
      tmp = zeta_s[n]; // Save zeta_i to pop into the prev zeta.
      zeta_s[n] = (zeta_s[n]*zeta_s_prev[n]*beta_prev)/(beta*alpha*(zeta_s_prev[n]-zeta_s[n]) + zeta_s_prev[n]*beta_prev*(1-shifts[n]*beta));
      zeta_s_prev[n] = tmp; 
      
      //cout << "zeta_n = " << zeta_s[n] << ", zeta_{n-1} = " << zeta_s_prev[n];
      
      // beta_i^sigma = beta_i zeta_{n+1}^sigma / zeta_n^sigma
      beta_s[n] = beta*zeta_s[n]/zeta_s_prev[n];
      
      // 4. x_s = x_s - beta_s p_s
      for (i = 0; i < size; i++)
      {
        phi[n][i] = phi[n][i] - beta_s[n]*p_s[n][i];
      }
      
      //cout << ", beta_n = " << beta_s[n] << "\n"; 
    }

    // 5. r = r + beta Ap
    for (i = 0; i < size; i++)
    {
      r[i] = r[i] + beta*Ap[i];
    }
    
    // Exit if new residual is small enough
    rsqNew = norm2sq<double>(r, size);
    
    print_verbosity_resid(verb, "CG-M", k+1, invif.ops_count, sqrt(rsqNew)/bsqrt); 
    
    // The residual of the shifted systems is zeta_s[n]*sqrt(rsqNew). Stop iterating on converged systems.
    if (k % resid_freq_check == 0)
    {
      for (n = 0; n < n_shift_rem; n++)
      {
        if (zeta_s[n]*sqrt(rsqNew) < eps*bsqrt) // if the residual of vector 'n' is sufficiently small...
        {
          // Permute it out.
          n_shift_rem--;
          
          cout << "Vector " << n << " has converged.\n" << flush;
          
          if (n_shift_rem != n) // Reorder in the case of out-of-order convergence. 
          {
            // Update mapping.
            mapping[n] = n_shift_rem;
            mapping[n_shift_rem] = n;

            // Permute phi, p_s, alpha_s, beta_s, zeta_s, zeta_s_prev.
            tmp_ptr = phi[n_shift_rem];
            phi[n_shift_rem] = phi[n];
            phi[n] = tmp_ptr;
            
            tmp_ptr = p_s[n_shift_rem];
            p_s[n_shift_rem] = p_s[n];
            p_s[n] = tmp_ptr;
            
            tmp_dbl = alpha_s[n_shift_rem];
            alpha_s[n_shift_rem] = alpha_s[n];
            alpha_s[n] = tmp_dbl;
            
            tmp_dbl = beta_s[n_shift_rem];
            beta_s[n_shift_rem] = beta_s[n];
            beta_s[n] = tmp_dbl;
            
            tmp_dbl = zeta_s[n_shift_rem];
            zeta_s[n_shift_rem] = zeta_s[n];
            zeta_s[n] = tmp_dbl;
            
            tmp_dbl = zeta_s_prev[n_shift_rem];
            zeta_s_prev[n_shift_rem] = zeta_s_prev[n];
            zeta_s_prev[n] = tmp_dbl;
            
            // We swapped with the end, so we need to recheck the end.
            n--;
          }
        }
      }
    }

    if (sqrt(rsqNew) < eps*bsqrt || n_shift_rem == 0 || k == max_iter-1) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    
  
    // 6. alpha = rsqNew / rsq.
    alpha = rsqNew / rsq;
    rsq = rsqNew; 
    
    //cout << "alpha = " << alpha << "\n";  
    
    for (n = 0; n < n_shift_rem; n++)
    {
      // 7. alpha_s = alpha * zeta_s * beta_s / (zeta_s_prev * beta)
      alpha_s[n] = alpha*zeta_s[n]*beta_s[n]/(zeta_s_prev[n] * beta);
      //cout << "alpha_n = " << alpha_s[n] << "\n";
      
      // 8. p_s = zeta_s_prev r + alpha_s p_s
      for (i = 0; i < size; i++)
      {
        p_s[n][i] = zeta_s[n]*r[i] + alpha_s[n]*p_s[n][i]; // Note, there's a typo in the paper where they have zeta_s_prev.
      }
    }
    
    // Compute the new Ap.
    for (i = 0; i < size; i++)
    {
      p[i] = r[i] + alpha*p[i];
    }
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
  
  // Undo the permutation damage.
  // Only need to permute phi.
  for (n = 0; n < n_shift; n++)
  {
    // Find the true n'th vector.
    if (mapping[n] != n)
    {
      for (int m = n+1; m < n_shift; m++)
      {
        if (mapping[m] == n) // Match, swap.
        {
          tmp_ptr = phi[m];
          phi[m] = phi[n];
          phi[n] = tmp_ptr;
          
          mapping[m] = mapping[n];
          mapping[n] = n;
          
          n--;
          break;
        }
      }
    }
  }
  
  // Calculate explicit rsqs.
  for (n = 0; n < n_shift; n++)
  {
    zero<double>(Ap, size);
    (*matrix_vector)(Ap, phi[n], extra_info); invif.ops_count++;
    for (i = 0; i < size; i++)
    {
      Ap[i] = Ap[i] + (shifts[n]*phi[n][i]);
    }
    invif.resSqmrhs[n] = diffnorm2sq<double>(Ap, phi0, size);
  }
  
  
  // Free all the things!
  delete[] r;
  delete[] p;
  delete[] Ap;
  
  for (i = 0; i < n_shift; i++)
  {
    delete[] p_s[i];
  }
  delete[] p_s;
  
  delete[] alpha_s;
  delete[] beta_s;
  delete[] zeta_s;
  delete[] zeta_s_prev;
  
  delete[] mapping; 

  // Need to update verbosity type things.
  //print_verbosity_summary(verb, "CG", invif.success, k, invif.ops_count, sqrt(truersq)/bsqrt);
  
  if (verb != 0)
  {
    if (verb->verbosity == VERB_SUMMARY || verb->verbosity == VERB_RESTART_DETAIL || verb->verbosity == VERB_DETAIL)
    {
      std::cout << verb->verb_prefix << "CG-M " << " Success " << (invif.success ? "Y" : "N") << " Iter " << k+1 << " Ops " << invif.ops_count << " RelRes ";
      for (n = 0; n < n_shift; n++)
      {
        std::cout << sqrt(invif.resSqmrhs[n])/bsqrt << " ";
      }

      std::cout << "\n";
    }
  }
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "CG";
  return invif; // Convergence 
} 


/*inversion_info minv_vector_cg(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
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
} */

