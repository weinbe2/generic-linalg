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
//#include <random>

#include "generic_vector.h"

#include "generic_bicgstab_m.h"

using namespace std;

// Solves lhs = A^(-1) rhs using multishift BiCGstab as defined in http://arxiv.org/pdf/hep-lat/9612014.pdf
// Assumes there are n_shift values in "shifts", and that they are sorted.
// resid_freq_check is how often to check the residual of other solutions. This lets us stop iterating on converged systems. 
// REMARK: Since we use 'phi' for the soln vector, we use 'psi' as a variable instead of 'phi' in the context
//         of the source paper. 
inversion_info minv_vector_bicgstab_m(double **phi, double *phi0, int n_shift, int size, int resid_freq_check, int max_iter, double eps, double* shifts, void (*matrix_vector)(double*,double*,void*), void* extra_info, bool worst_first, inversion_verbose_struct* verb)
{
  
  // Initialize vectors.
  double *r, *r_prev, *s, *As, *w, *Aw, *wdag;
  double **s_s;
  double alpha, beta, beta_prev, psi, delta, delta_prev, chi, rsqNew, bsqrt, truersq, tmp; 
  double *alpha_s, *beta_s, *zeta_s, *zeta_s_prev, *chi_s, *rho_s, *rho_s_prev;
  int k,i,n;
  int n_shift_rem = n_shift; // number of systems to still iterate on. 
  int* mapping; // holds the mapping between vectors and the original vector ordering.
                // this is because some vectors may converge before others. 
  double* tmp_ptr; // temporary pointer for swaps.
  double tmp_dbl; // temporary double for swaps. 
  int tmp_int; // temporary int for swaps. 
  bool breakdown = false; // check for breakdown
  
  //std::mt19937 generator (1338u); // RNG, 1337u is the seed. 
  
  // Prepare an inversion_info for multiple residuals.
  inversion_info invif(n_shift); 

  // Allocate memory.
  alpha_s = new double[n_shift];
  beta_s = new double[n_shift];
  zeta_s = new double[n_shift];
  zeta_s_prev = new double[n_shift];
  chi_s = new double[n_shift];
  rho_s = new double[n_shift];
  rho_s_prev = new double[n_shift];
  
  s_s = new double*[n_shift];
  for (n = 0; n < n_shift; n++)
  {
    s_s[n] = new double[size];
  }
  
  r = new double[size];
  r_prev = new double[size];
  s = new double[size];
  As = new double[size];
  w = new double[size];
  Aw = new double[size];
  wdag = new double[size];
  
  // Initialize mapping.
  mapping = new int[n_shift];
  for (n = 0; n < n_shift; n++)
  {
    mapping[n] = n; // All vectors are currently in order.
  }

  // Initialize values.
  rsqNew = 0.0; bsqrt = 0.0; truersq = 0.0; k=0;
  for (n = 0; n < n_shift; n++)
  {
    // beta_0, zeta_0, zeta_-1
    beta_s[n] = zeta_s[n] = zeta_s_prev[n] = rho_s[n] = rho_s_prev[n] = 1.0;
    // alpha_0. 
    alpha_s[n] = chi_s[n] = 0.0;
  }
  beta = 1.0; alpha = 0.0;

  // Zero vectors;
  zero<double>(r, size); 
  zero<double>(s, size); zero<double>(As, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // There can't be an initial guess... though it is sort of possible, in reference to:
  // http://arxiv.org/pdf/0810.1081v1.pdf
  
  // 1. x_sigma = 0, r = s_sigma = b.
  for (n = 0; n < n_shift; n++)
  {
    copy<double>(s_s[n], phi0, size);
    zero<double>(phi[n], size);
  }
  copy<double>(s, phi0, size);
  copy<double>(r, phi0, size);
  copy<double>(r_prev, r, size); 
  
  // Compute As.
  zero<double>(As, size);
  (*matrix_vector)(As, s, extra_info); invif.ops_count++;
  
  // Assign wdag.
  //gaussian<double>(wdag, size, generator); 
  copy<double>(wdag, r, size);
  
  // 2. delta_0 = <wdag, r>, psi_0 = <wdag, As>/delta_0
  delta = delta_prev = dot<double>(wdag, r, size);
  psi = dot<double>(wdag, As, size)/delta;

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // 3. beta_i = - 1.0/psi. Which is a weird switch from the normal notation, but whatever.
    beta_prev = beta; 
    beta = -1.0/psi;
    
    // 5. w = r + beta As
    for (i = 0; i < size; i++)
    {
      w[i] = r[i] + beta*As[i];
    }
    
    // 6. Compute Aw. chi = <Aw, w>/|Aw|^2.
    zero<double>(Aw, size);
    (*matrix_vector)(Aw, w, extra_info); invif.ops_count++;
    chi = dot<double>(Aw, w, size)/norm2sq<double>(Aw, size);
    
    // 8. Update r. 
    copy<double>(r_prev, r, size); 
    for (i = 0; i < size; i++)
    {
      r[i] = w[i] - chi*Aw[i];
    }
    
    
    // Update beta_s, zeta_s, chi_s, rho_s.
    for (n = 0; n < n_shift_rem; n++)
    {
      // 4. Calculate beta_i^sigma, zeta_i+1^sigma according to 2.42 to 2.44.
      // zeta_{i+1}^sigma = complicated...
      tmp = zeta_s[n]; // Save zeta_i to pop into the prev zeta.
      zeta_s[n] = (zeta_s[n]*zeta_s_prev[n]*beta_prev)/(beta*alpha*(zeta_s_prev[n]-zeta_s[n]) + zeta_s_prev[n]*beta_prev*(1.0-shifts[n]*beta));
      zeta_s_prev[n] = tmp; 
      

      // beta_i^sigma = beta_i zeta_{n+1}^sigma / zeta_n^sigma
      beta_s[n] = beta*zeta_s[n]/zeta_s_prev[n];
      
      // 7. Calculate chi_i^sigma, rho_{i+1}^sigma according to 2.23 to 2.24.
      chi_s[n] = chi/(1.0+chi*shifts[n]);
      tmp = rho_s[n];
      rho_s[n] = rho_s[n]/(1.0+chi*shifts[n]);
      rho_s_prev[n] = tmp; 
      
      // 9. x_s = x_s - beta_s s_s + chi_s rho_s^prev zeta_s w (Maybe typos in paper?)
      for (i = 0; i < size; i++)
      {
        phi[n][i] = phi[n][i] - beta_s[n]*s_s[n][i] + chi_s[n]*rho_s_prev[n]*zeta_s[n]*w[i]; // ?
        //phi[n][i] = phi[n][i] - beta_s[n]*s_s[n][i] + chi_s[n]*rho_s[n]*zeta_s[n]*w[i]; // ?
      }
      
      //cout << ", beta_n = " << beta_s[n] << "\n"; 
    }
    
    // Exit if new residual is small enough
    rsqNew = norm2sq<double>(r, size);
    
    print_verbosity_resid(verb, "BICGSTAB-M", k+1, invif.ops_count, sqrt(rsqNew)/bsqrt); 
    
    
    // The residual of the shifted systems is zeta_s[n]*sqrt(rsqNew). Stop iterating on converged systems.
    if (k % resid_freq_check == 0)
    {
      for (n = 0; n < n_shift_rem; n++)
      {
        if (abs(zeta_s[n]*rho_s[n])*sqrt(rsqNew) < eps*bsqrt) // if the residual of vector 'n' is sufficiently small...
        {
          // Permute it out.
          n_shift_rem--;
          
          //cout << "Vector " << mapping[n] << " has converged.\n" << flush;
          
          if (n_shift_rem != n) // Reorder in the case of out-of-order convergence. 
          {
            // Update mapping.
            tmp_int = mapping[n_shift_rem];
            mapping[n_shift_rem] = mapping[n];
            mapping[n] = tmp_int;

            // Permute phi, s_s, alpha_s, beta_s, zeta_s, zeta_s_prev, 
            //                   chi_s, rho_s, rho_s_prev, shifts. 
            tmp_ptr = phi[n_shift_rem];
            phi[n_shift_rem] = phi[n];
            phi[n] = tmp_ptr;
            
            tmp_ptr = s_s[n_shift_rem];
            s_s[n_shift_rem] = s_s[n];
            s_s[n] = tmp_ptr;
            
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
            
            tmp_dbl = chi_s[n_shift_rem];
            chi_s[n_shift_rem] = chi_s[n];
            chi_s[n] = tmp_dbl;
            
            tmp_dbl = rho_s[n_shift_rem];
            rho_s[n_shift_rem] = rho_s[n];
            rho_s[n] = tmp_dbl;
            
            tmp_dbl = rho_s_prev[n_shift_rem];
            rho_s_prev[n_shift_rem] = rho_s_prev[n];
            rho_s_prev[n] = tmp_dbl;
            
            tmp_dbl = shifts[n_shift_rem];
            shifts[n_shift_rem] = shifts[n];
            shifts[n] = tmp_dbl;
            
            // We swapped with the end, so we need to recheck the end.
            n--;
          }
        }
      }
    }

    if (/*sqrt(rsqNew) < eps*bsqrt || */(worst_first && (abs(zeta_s[0]*rho_s[0])*sqrt(rsqNew) < eps*bsqrt)) || rsqNew != rsqNew || n_shift_rem == 0 || k == max_iter-1 ) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    
    // 10. delta = <wdag, r>
    delta_prev = delta;
    delta = dot<double>(wdag, r, size);
    
    // 11. alpha = -(beta_i*delta_i)/(delta_i_prev*chi)
    alpha = - beta*delta/(delta_prev*chi);

    //cout << "alpha = " << alpha << "\n";  
    
    // 12. alpha_s = alpha * zeta_s * beta_s / (zeta_s_prev * beta)
    for (n = 0; n < n_shift_rem; n++)
    {
      alpha_s[n] = alpha*zeta_s[n]*beta_s[n]/(zeta_s_prev[n] * beta);
      //cout << "alpha_n = " << alpha_s[n] << "\n";
    }
    
    // 13. s = r + alpha*(s - chi As)
    for (i = 0; i < size; i++)
    {
      s[i] = r[i] + alpha*(s[i] - chi*As[i]);
    }
    
    // 14. s_sigma = ... well, what's in the paper.
    for (n = 0; n < n_shift_rem; n++)
    {
      if (abs(shifts[n]) != 0.0) 
      {
        for (i = 0; i < size; i++)
        {
          s_s[n][i] = zeta_s[n]*rho_s[n]*r[i] + alpha_s[n]*(s_s[n][i] - chi_s[n]/beta_s[n]*(zeta_s[n]*rho_s_prev[n]*w[i] - zeta_s_prev[n]*rho_s_prev[n]*r_prev[i]));
        }
      }
      else
      {
        copy<double>(s_s[n], s, size);
      }
    }
    
    // 15. Compute the new As. 
    zero<double>(As, size);
    (*matrix_vector)(As, s, extra_info); invif.ops_count++;
    
    // 16. psi = <wdag, As>/delta
    psi = dot<double>(wdag, As, size)/delta;
    
    // Check for breakdown.
    if (psi == 0)
    {
      breakdown = true; 
      break;
    }
      
  } 
    
  if(k == max_iter-1 || breakdown) {
    invif.success = false;
  }
  else
  {
     invif.success = true;
  }
  k++;
  
  // Undo the permutation damage.
  // Only need to permute phi, shifts. 
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
          
          tmp_dbl = shifts[m];
          shifts[m] = shifts[n];
          shifts[n] = tmp_dbl;
          
          mapping[m] = mapping[n];
          mapping[n] = n;
          
          n--;
          break;
        }
      }
    }
  }
  
  // Calculate explicit rsqs.
  double* relres = new double[n_shift];
  for (n = 0; n < n_shift; n++)
  {
    zero<double>(As, size);
    (*matrix_vector)(As, phi[n], extra_info); invif.ops_count++;
    for (i = 0; i < size; i++)
    {
      As[i] = As[i] + (shifts[n]*phi[n][i]);
    }
    invif.resSqmrhs[n] = diffnorm2sq<double>(As, phi0, size);
    relres[n] = sqrt(invif.resSqmrhs[n])/bsqrt;
  }
  
  
  // Free all the things!
  delete[] r;
  delete[] r_prev; 
  delete[] s;
  delete[] As;
  delete[] w;
  delete[] Aw;
  delete[] wdag;
  
  for (i = 0; i < n_shift; i++)
  {
    delete[] s_s[i];
  }
  delete[] s_s;
  
  delete[] alpha_s;
  delete[] beta_s;
  delete[] zeta_s;
  delete[] zeta_s_prev;
  delete[] chi_s;
  delete[] rho_s;
  delete[] rho_s_prev;
  
  delete[] mapping; 

  print_verbosity_summary_multi(verb, "BICGSTAB-M", invif.success, k, invif.ops_count, relres, n_shift);
  delete[] relres; 
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "BICGSTAB-M";
  return invif; // Convergence 
} 

// Solves lhs = A^(-1) rhs using multishift BiCGstab as defined in http://arxiv.org/pdf/hep-lat/9612014.pdf
// Assumes there are n_shift values in "shifts", and that they are sorted.
// resid_freq_check is how often to check the residual of other solutions. This lets us stop iterating on converged systems. 
// REMARK: Since we use 'phi' for the soln vector, we use 'psi' as a variable instead of 'phi' in the context
//         of the source paper. 
inversion_info minv_vector_bicgstab_m(complex<double> **phi, complex<double> *phi0, int n_shift, int size, int resid_freq_check, int max_iter, double eps, double* shifts, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, bool worst_first, inversion_verbose_struct* verb)
{
  
  // Initialize vectors.
  complex<double> *r, *r_prev, *s, *As, *w, *Aw, *wdag;
  complex<double> **s_s;
  complex<double> alpha, beta, beta_prev, psi, delta, delta_prev, chi, tmp;
  double rsqNew, bsqrt, truersq;
  complex<double> *alpha_s, *beta_s, *zeta_s, *zeta_s_prev, *chi_s, *rho_s, *rho_s_prev;
  int k,i,n;
  int n_shift_rem = n_shift; // number of systems to still iterate on. 
  int* mapping; // holds the mapping between vectors and the original vector ordering.
                // this is because some vectors may converge before others. 
  complex<double>* tmp_ptr; // temporary pointer for swaps.
  double tmp_dbl_real; // temporary
  complex<double> tmp_dbl; // temporary double for swaps. 
  int tmp_int; // temporary int for swaps. 
  bool breakdown = false; // check for breakdown
  
  //std::mt19937 generator (1338u); // RNG, 1337u is the seed. 
  
  // Prepare an inversion_info for multiple residuals.
  inversion_info invif(n_shift); 

  // Allocate memory.
  alpha_s = new complex<double>[n_shift];
  beta_s = new complex<double>[n_shift];
  zeta_s = new complex<double>[n_shift];
  zeta_s_prev = new complex<double>[n_shift];
  chi_s = new complex<double>[n_shift];
  rho_s = new complex<double>[n_shift];
  rho_s_prev = new complex<double>[n_shift];
  
  s_s = new complex<double>*[n_shift];
  for (n = 0; n < n_shift; n++)
  {
    s_s[n] = new complex<double>[size];
  }
  
  r = new complex<double>[size];
  r_prev = new complex<double>[size];
  s = new complex<double>[size];
  As = new complex<double>[size];
  w = new complex<double>[size];
  Aw = new complex<double>[size];
  wdag = new complex<double>[size];
  
  // Initialize mapping.
  mapping = new int[n_shift];
  for (n = 0; n < n_shift; n++)
  {
    mapping[n] = n; // All vectors are currently in order.
  }

  // Initialize values.
  rsqNew = 0.0; bsqrt = 0.0; truersq = 0.0; k=0;
  for (n = 0; n < n_shift; n++)
  {
    // beta_0, zeta_0, zeta_-1
    beta_s[n] = zeta_s[n] = zeta_s_prev[n] = rho_s[n] = rho_s_prev[n] = 1.0;
    // alpha_0. 
    alpha_s[n] = chi_s[n] = 0.0;
  }
  beta = 1.0; alpha = 0.0;

  // Zero vectors;
  zero<double>(r, size); 
  zero<double>(s, size); zero<double>(As, size);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2sq<double>(phi0, size));
  
  // There can't be an initial guess... though it is sort of possible, in reference to:
  // http://arxiv.org/pdf/0810.1081v1.pdf
  
  // 1. x_sigma = 0, r = s_sigma = b.
  for (n = 0; n < n_shift; n++)
  {
    copy<double>(s_s[n], phi0, size);
    zero<double>(phi[n], size);
  }
  copy<double>(s, phi0, size);
  copy<double>(r, phi0, size);
  copy<double>(r_prev, r, size); 
  
  // Compute As.
  zero<double>(As, size);
  (*matrix_vector)(As, s, extra_info); invif.ops_count++;
  
  // Assign wdag.
  //gaussian<double>(wdag, size, generator); 
  copy<double>(wdag, r, size);
  conj<double>(wdag, size); 
  
  // 2. delta_0 = <wdag, r>, psi_0 = <wdag, As>/delta_0
  delta = delta_prev = dot<double>(wdag, r, size);
  psi = dot<double>(wdag, As, size)/delta;

  // iterate till convergence
  for(k = 0; k< max_iter; k++) {
    
    // 3. beta_i = - 1.0/psi. Which is a weird switch from the normal notation, but whatever.
    beta_prev = beta; 
    beta = -1.0/psi;
    
    // 5. w = r + beta As
    for (i = 0; i < size; i++)
    {
      w[i] = r[i] + beta*As[i];
    }
    
    // 6. Compute Aw. chi = <Aw, w>/|Aw|^2.
    zero<double>(Aw, size);
    (*matrix_vector)(Aw, w, extra_info); invif.ops_count++;
    chi = dot<double>(Aw, w, size)/norm2sq<double>(Aw, size);
    
    // 8. Update r. 
    copy<double>(r_prev, r, size); 
    for (i = 0; i < size; i++)
    {
      r[i] = w[i] - chi*Aw[i];
    }
    
    
    // Update beta_s, zeta_s, chi_s, rho_s.
    for (n = 0; n < n_shift_rem; n++)
    {
      // 4. Calculate beta_i^sigma, zeta_i+1^sigma according to 2.42 to 2.44.
      // zeta_{i+1}^sigma = complicated...
      tmp = zeta_s[n]; // Save zeta_i to pop into the prev zeta.
      zeta_s[n] = (zeta_s[n]*zeta_s_prev[n]*beta_prev)/(beta*alpha*(zeta_s_prev[n]-zeta_s[n]) + zeta_s_prev[n]*beta_prev*(1.0-shifts[n]*beta));
      zeta_s_prev[n] = tmp; 
      

      // beta_i^sigma = beta_i zeta_{n+1}^sigma / zeta_n^sigma
      beta_s[n] = beta*zeta_s[n]/zeta_s_prev[n];
      
      // 7. Calculate chi_i^sigma, rho_{i+1}^sigma according to 2.23 to 2.24.
      chi_s[n] = chi/(1.0+chi*shifts[n]);
      tmp = rho_s[n];
      rho_s[n] = rho_s[n]/(1.0+chi*shifts[n]);
      rho_s_prev[n] = tmp; 
      
      // 9. x_s = x_s - beta_s s_s + chi_s rho_s^prev zeta_s w (Maybe typos in paper?)
      for (i = 0; i < size; i++)
      {
        phi[n][i] = phi[n][i] - beta_s[n]*s_s[n][i] + chi_s[n]*rho_s_prev[n]*zeta_s[n]*w[i]; // ?
        //phi[n][i] = phi[n][i] - beta_s[n]*s_s[n][i] + chi_s[n]*rho_s[n]*zeta_s[n]*w[i]; // ?
      }
      
      //cout << ", beta_n = " << beta_s[n] << "\n"; 
    }
    
    // Exit if new residual is small enough
    rsqNew = norm2sq<double>(r, size);
    
    print_verbosity_resid(verb, "BICGSTAB-M", k+1, invif.ops_count, sqrt(rsqNew)/bsqrt); 
    
    
    // The residual of the shifted systems is zeta_s[n]*sqrt(rsqNew). Stop iterating on converged systems.
    if (k % resid_freq_check == 0)
    {
      for (n = 0; n < n_shift_rem; n++)
      {
        if (abs(zeta_s[n]*rho_s[n])*sqrt(rsqNew) < eps*bsqrt) // if the residual of vector 'n' is sufficiently small...
        {
          // Permute it out.
          n_shift_rem--;
          
          //cout << "Vector " << mapping[n] << " has converged.\n" << flush;
          
          if (n_shift_rem != n) // Reorder in the case of out-of-order convergence. 
          {
            // Update mapping.
            tmp_int = mapping[n_shift_rem];
            mapping[n_shift_rem] = mapping[n];
            mapping[n] = tmp_int;

            // Permute phi, s_s, alpha_s, beta_s, zeta_s, zeta_s_prev, 
            //                   chi_s, rho_s, rho_s_prev, shifts. 
            tmp_ptr = phi[n_shift_rem];
            phi[n_shift_rem] = phi[n];
            phi[n] = tmp_ptr;
            
            tmp_ptr = s_s[n_shift_rem];
            s_s[n_shift_rem] = s_s[n];
            s_s[n] = tmp_ptr;
            
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
            
            tmp_dbl = chi_s[n_shift_rem];
            chi_s[n_shift_rem] = chi_s[n];
            chi_s[n] = tmp_dbl;
            
            tmp_dbl = rho_s[n_shift_rem];
            rho_s[n_shift_rem] = rho_s[n];
            rho_s[n] = tmp_dbl;
            
            tmp_dbl = rho_s_prev[n_shift_rem];
            rho_s_prev[n_shift_rem] = rho_s_prev[n];
            rho_s_prev[n] = tmp_dbl;
            
            tmp_dbl_real = shifts[n_shift_rem];
            shifts[n_shift_rem] = shifts[n];
            shifts[n] = tmp_dbl_real;
            
            // We swapped with the end, so we need to recheck the end.
            n--;
          }
        }
      }
    }

    if (/*sqrt(rsqNew) < eps*bsqrt || */(worst_first && (abs(zeta_s[0]*rho_s[0])*sqrt(rsqNew) < eps*bsqrt)) || rsqNew != rsqNew || n_shift_rem == 0 || k == max_iter-1 ) {
      //        printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    
    // 10. delta = <wdag, r>
    delta_prev = delta;
    delta = dot<double>(wdag, r, size);
    
    // 11. alpha = -(beta_i*delta_i)/(delta_i_prev*chi)
    alpha = - beta*delta/(delta_prev*chi);

    //cout << "alpha = " << alpha << "\n";  
    
    // 12. alpha_s = alpha * zeta_s * beta_s / (zeta_s_prev * beta)
    for (n = 0; n < n_shift_rem; n++)
    {
      alpha_s[n] = alpha*zeta_s[n]*beta_s[n]/(zeta_s_prev[n] * beta);
      //cout << "alpha_n = " << alpha_s[n] << "\n";
    }
    
    // 13. s = r + alpha*(s - chi As)
    for (i = 0; i < size; i++)
    {
      s[i] = r[i] + alpha*(s[i] - chi*As[i]);
    }
    
    // 14. s_sigma = ... well, what's in the paper.
    for (n = 0; n < n_shift_rem; n++)
    {
      if (abs(shifts[n]) != 0.0) 
      {
        for (i = 0; i < size; i++)
        {
          s_s[n][i] = zeta_s[n]*rho_s[n]*r[i] + alpha_s[n]*(s_s[n][i] - chi_s[n]/beta_s[n]*(zeta_s[n]*rho_s_prev[n]*w[i] - zeta_s_prev[n]*rho_s_prev[n]*r_prev[i]));
        }
      }
      else
      {
        copy<double>(s_s[n], s, size);
      }
    }
    
    // 15. Compute the new As. 
    zero<double>(As, size);
    (*matrix_vector)(As, s, extra_info); invif.ops_count++;
    
    // 16. psi = <wdag, As>/delta
    psi = dot<double>(wdag, As, size)/delta;
    
    // Check for breakdown.
    if (abs(psi) == 0)
    {
      breakdown = true; 
      break;
    }
      
  } 
    
  if(k == max_iter-1 || breakdown) {
    invif.success = false;
  }
  else
  {
     invif.success = true;
  }
  k++;
  
  // Undo the permutation damage.
  // Only need to permute phi, shifts. 
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
          
          tmp_dbl_real = shifts[m];
          shifts[m] = shifts[n];
          shifts[n] = tmp_dbl_real;
          
          mapping[m] = mapping[n];
          mapping[n] = n;
          
          n--;
          break;
        }
      }
    }
  }
  
  // Calculate explicit rsqs.
  double* relres = new double[n_shift];
  for (n = 0; n < n_shift; n++)
  {
    zero<double>(As, size);
    (*matrix_vector)(As, phi[n], extra_info); invif.ops_count++;
    for (i = 0; i < size; i++)
    {
      As[i] = As[i] + (shifts[n]*phi[n][i]);
    }
    invif.resSqmrhs[n] = diffnorm2sq<double>(As, phi0, size);
    relres[n] = sqrt(invif.resSqmrhs[n])/bsqrt;
  }
  
  
  // Free all the things!
  delete[] r;
  delete[] r_prev; 
  delete[] s;
  delete[] As;
  delete[] w;
  delete[] Aw;
  delete[] wdag;
  
  for (i = 0; i < n_shift; i++)
  {
    delete[] s_s[i];
  }
  delete[] s_s;
  
  delete[] alpha_s;
  delete[] beta_s;
  delete[] zeta_s;
  delete[] zeta_s_prev;
  delete[] chi_s;
  delete[] rho_s;
  delete[] rho_s_prev;
  
  delete[] mapping; 

  print_verbosity_summary_multi(verb, "BICGSTAB-M", invif.success, k, invif.ops_count, relres, n_shift);
  delete[] relres; 
  
  invif.resSq = truersq;
  invif.iter = k;
  invif.name = "BICGSTAB-M";
  return invif; // Convergence 
} 
