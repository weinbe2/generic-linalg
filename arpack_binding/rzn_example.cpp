
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <string.h>

using std::complex;

#include "arpack_interface.h"
#include "generic_vector.h"
#include "generic_bicgstab.h"

// Various maxima. This is used to allocate more than
// enough memory.

// Maximum dimension of the matrix.
#define MAXN 1024

// Maximum number of eigenvalues requested.
#define MAXNEV 256

// Maximum number of basis vectors for Arnoldi iterations.
// Remember: MAXNEV+2 < MAXNCV < MAXN
#define MAXNCV 512

// Unique to this problem: define size of 1D laplace chain.
#define NDIM 200
#define MASS 0.1
#define PI 3.14159265358979323846
typedef struct {
   double ndim;
   double mass;
} laplace_1d_t;
   

// This function does the multiply!
void laplace_1d(complex<double>* lhs, complex<double>* rhs, void* extra_data);


int main(int argc, char** argv)
{
   // Declare variables.
   int i, j; // Iterators.
   int n, nev, ncv; // matrix dimension, eigenvals, internal vals.
   int maxitr; // maximum number of iterations.
   char which[3]; // Sets which part of spectrum to get.
   double tol; // Set tolerance. 0 = machine.
   complex<double> sigma; // Set to zero for now.
   laplace_1d_t sys_info; // Set the mass term and dim of matrix.
   arpack_solve_t info_solve; // Hold info on the solve.
   complex<double> **evec, *eval; // Will hold eigenvalues, eigenvectors.
   complex<double> *resid; // Hold residuals of eigenvalues.
   complex<double> *ax, *r; // Temporary space to hold vectors to compute
                   // residuals.
   int maxitr_cg; // Max number of iterations for cg.
   double tol_cg; // Tolerance for cg.
                   
   arpack_dcn_t* ar_strc; // Holds arpack info.
   
   // Set dimension of matrix, number of eigenvalues, etc.
   n = NDIM;
   nev = 20; // nev < n
   ncv = 40; // nev+2 < ncv < n.
   sys_info.mass = MASS;
   sys_info.ndim = NDIM;
	
   
   // Configure calculating eigenvectors.
   maxitr = 4000; // 4000 iterations!
   strcpy(which, "SM"); // This could be SM (smallest magnitude).
                 // This could also be LM (largest magnitude)
                 // SR (smallest real), SI (smallest imaginary),
                 // and similar for largest.
   tol = 1e-7; // 0 = machine precision.
   sigma = MASS*0.99; // Zero for now.
	
   // Allocate space.
   evec = new complex<double>*[nev]; // eigenvals
	for (i = 0; i < nev; i++)
	{
		evec[i] = new complex<double>[n];
	}
   eval = new complex<double>[nev]; // eigenvecs
   resid = new complex<double>[nev]; // residuals
   ax = new complex<double>[n]; // temp vector.
   r = new complex<double>[n]; // temp vector.
   
   printf("Allocated.\n"); fflush(stdout);
   
   // Initialize arpack. Largely just allocates internal memory.
   ar_strc = arpack_dcn_init(MAXN, MAXNEV, MAXNCV);
   
   printf("Init arpack.\n"); fflush(stdout);
   
   // Get some eigenvalues and vectors!
   info_solve = arpack_dcn_getev(ar_strc, eval, evec, n, nev, ncv,
                              maxitr, which, tol, sigma,
                              &laplace_1d, (void*)(&sys_info));
   
   printf("Got evals/evecs.\n"); fflush(stdout);         
   
   if (info_solve.is_error == 0 && info_solve.nconv > 0)
   {
      // Print info about the eigensolve.
      printf("Number of converged eigenvalues: %d\n", info_solve.nconv);
      printf("Number of iteration steps:       %d\n", info_solve.niter);
      printf("Number of matrix multiplies:     %d\n", info_solve.nops);
      
      // Get the residual norm:
      // || A*x - lambda*x ||
      // for the nconv accurately computed eigenvals/vecs.
      // Put the residual here:
      
      // Loop over all good eigenvalues.
      for (i=0;i<info_solve.nconv;i++)
      {
         // Get pointer to start of i'th eigenvector.
         complex<double> *evec_local = evec[i];
         
         // If you want to print it out, uncomment this.
         //for (int j=0;j<NDIM;j++)
         //{
         //   printf("%d\t%.8e\n", i, evec_local[j]);
         //}
         
         // Get A*x, where x is an eigenvector.
         laplace_1d(ax, evec_local, (void*)(&sys_info));
         
         resid[i] = 0;
         

         // Get the norm squared.
         for (j=0;j<NDIM;j++)
         {
			 resid[i] += real(conj(ax[j]-eval[i]*evec_local[j])*(ax[j]-eval[i]*evec_local[j]));
         }
         resid[i] = sqrt(resid[i])/(abs(eval[i])*abs(eval[i]));
      }
      
      printf("Got residuals.\n"); fflush(stdout);
   
      // Print eigenvalues, residual.
      // Also, print what mode we get for 1D Laplace
      for (i=0;i<info_solve.nconv;i++)
      {
         // We compute the mode. This can overflow on some
         // ops, though, so we're careful of that.
         double tmp,tmp2; 
         tmp = (real(eval[i])-MASS)/4; 
         tmp2 = (tmp > 0 ? sqrt(tmp) : 0); // make sure tmp > 0
         tmp = (tmp2 < 1 ? (tmp2 > -1 ? asin(tmp2) : asin(-1.0)) : asin(1.0))*NDIM/PI; // make sure -1 < tmp2 < 1
         // tmp = NDIM/PI*asin(sqrt((d[i]-MASS)/4));
         
         printf("%d\t%.8e\t%.8e\t%.8e\t%.8e\n", i+1, real(eval[i]), imag(eval[i]), real(resid[i]), tmp);
      }
   
   }
   
   // Now redo it with shift invert!
   strcpy(which, "LM"); // This could be SM (smallest magnitude).
                 // This could also be LM (largest magnitude)
                 // SR (smallest real), SI (smallest imaginary),
                 // and similar for largest.
   maxitr_cg = 4000;
   tol_cg = 1e-9;
   
   // Get some eigenvalues and vectors using the shift-invert op. 
   info_solve = arpack_dcn_getev_sinv(ar_strc, eval, evec, n, nev, ncv,
                              maxitr, which, tol, sigma,
                              &laplace_1d, &minv_vector_bicgstab, maxitr_cg, tol_cg, (void*)(&sys_info));
   
   printf("Got evals/evecs.\n"); fflush(stdout);         
   
   if (info_solve.is_error == 0 && info_solve.nconv > 0)
   {
      // Print info about the eigensolve.
      printf("Number of converged eigenvalues: %d\n", info_solve.nconv);
      printf("Number of iteration steps:       %d\n", info_solve.niter);
      printf("Number of matrix multiplies:     %d\n", info_solve.nops);
      
      // Get the residual norm:
      // || A*x - lambda*x ||
      // for the nconv accurately computed eigenvals/vecs.
      // Put the residual here:
      
      // Loop over all good eigenvalues.
      for (i=0;i<info_solve.nconv;i++)
      {
         // Get pointer to start of i'th eigenvector.
         complex<double> *evec_local = evec[i];
         
         // If you want to print it out, uncomment this.
         //for (int j=0;j<NDIM;j++)
         //{
         //   printf("%d\t%.8e\n", i, evec_local[j]);
         //}
         
         // Get A*x, where x is an eigenvector.
         laplace_1d((complex<double>*)ax, evec_local, (void*)(&sys_info));
         
         resid[i] = 0;
         

         // Get the norm squared.
         for (j=0;j<NDIM;j++)
         {
			 resid[i] += real(conj(ax[j]-eval[i]*evec_local[j])*(ax[j]-eval[i]*evec_local[j]));
         }
         resid[i] = sqrt(resid[i])/(abs(eval[i])*abs(eval[i]));
      }
      
      printf("Got residuals.\n"); fflush(stdout);
   
      // Print eigenvalues, residual.
      // Also, print what mode we get for 1D Laplace
      for (i=0;i<info_solve.nconv;i++)
      {
         // We compute the mode. This can overflow on some
         // ops, though, so we're careful of that.
         double tmp,tmp2; 
         tmp = (real(eval[i])-MASS)/4; 
         tmp2 = (tmp > 0 ? sqrt(tmp) : 0); // make sure tmp > 0
         tmp = (tmp2 < 1 ? (tmp2 > -1 ? asin(tmp2) : asin(-1.0)) : asin(1.0))*NDIM/PI; // make sure -1 < tmp2 < 1
         // tmp = NDIM/PI*asin(sqrt((d[i]-MASS)/4));
         
         printf("%d\t%.8e\t%.8e\t%.8e\t%.8e\n", i+1, real(eval[i]), imag(eval[i]), real(resid[i]), tmp);
      }
   
   }
   
   // Clean up!
	for (i = 0; i < nev; i++)
	{
		delete[] evec[i];
	}
   delete[] evec;
	delete[] eval;
	delete[] resid;
	delete[] ax;
	delete[] r;
   arpack_dcn_free(&ar_strc);
   
}

// matrix multiplication function!
// Return lhs = A*rhs;

void laplace_1d(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
   // Apply the 1D Laplace eqn with a mass!
   
   // Declare variables.
   int ndim, i;
   double mass;
   laplace_1d_t* system_info;
   
   // Get input data.
   system_info = (laplace_1d_t*)extra_data;
   mass = system_info->mass;
   ndim = system_info->ndim;
         
   for (i=0;i<ndim;i++)
   {
      lhs[i] = (2+mass)*rhs[i]-rhs[(i-1+ndim)%ndim]-rhs[(i+1)%ndim];
   }
}
