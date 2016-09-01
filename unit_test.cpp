// Fri Apr 29 12:03:02 EDT 2016
// Evan S Weinberg 2016

// This is a set of testing routines that make sure the inverters (CG, BiCGStab, GMRES)
// always work! Need to update with imaginary value tests. 
    
// Need to put in a test for power iterations. Prob requires returning eigenvector. 

#include <iostream>
#include <iomanip> // to set output precision.
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_inverters.h"
#include "generic_inverters_precond.h"
#include "generic_eigenvalues.h"
#include "generic_vector.h"
#include "generic_traits.h"
#include "verbosity.h"

using namespace std; 

// For now, define the length in a direction.
#define N 128

// Define pi.
#define PI 3.141592653589793

// Define mass.
#define MASS 0.1*0.1

// Square laplacian function.
void square_laplacian(double* lhs, double* rhs, void* extra_data);

// Zero out vectors, set a point source.
void initialize_test(double* lattice, double* lhs, double* rhs, double* check, int size);

// Check solution.
double check_test(double* lhs, double* rhs, double* check, int size, void (*matrix_vector)(double*,double*,void*), void* extra_info);


int main(int argc, char** argv)
{  
   // Declare some variables.
   int i, j;
   double *lattice; // At some point, I'll have a generic (template) lattice class.
   double *lhs, *rhs, *check; // For some Kinetic terms.
   double eig = 0.0; // To test power iterations.
   double explicit_resid = 0.0;
   double bnorm = 0.0;
   inversion_info invif;
   eigenvalue_info eigif;
    
   // Set output precision.
   cout << setiosflags(ios::scientific) << setprecision(6);
    
   // Create a verbosity struct.
   inversion_verbose_struct verb;
   verb.verbosity = VERB_DETAIL;
   verb.verb_prefix = "";
   verb.precond_verbosity = VERB_PASS_THROUGH;
   verb.precond_verb_prefix = "Prec ";
 
   // Initialize the lattice. Indexing: index = y*N + x.
   lattice = new double[N*N];
   lhs = new double[N*N];
   rhs = new double[N*N];   
   check = new double[N*N]; 
   
   printf("Begin Check CG.\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_cg(lhs, rhs, N*N, 4000, 1e-6, square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check CG.\n");
   printf("\n\n\n");
   
   printf("Begin Check CR.\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_cr(lhs, rhs, N*N, 4000, 1e-6, square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check CR.\n");
   printf("\n\n\n");
   
   printf("Begin Check restarted CR(8).\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_cr_restart(lhs, rhs, N*N, 4000, 1e-6, 8, square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check restarted CR(8).\n");
   printf("\n\n\n");
    
   printf("Begin Check BiCGStab.\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_bicgstab(lhs, rhs, N*N, 4000, 1e-6, square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check BiCGStab.\n");
   printf("\n\n\n");
	
   printf("Begin Check restarted BiCGStab(8).\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_bicgstab_restart(lhs, rhs, N*N, 4000, 1e-6, 8, square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check restarted BiCGStab(8).\n");
   printf("\n\n\n");
	
	return 0;
    
   printf("Begin Check GCR.\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_gcr(lhs, rhs, N*N, 4000, 1e-6, square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check GCR.\n");
   printf("\n\n\n");
    
   printf("Begin Check restarted GCR(8).\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_gcr_restart(lhs, rhs, N*N, 4000, 1e-6, 8, square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check restarted GCR(8).\n");
   printf("\n\n\n");
   
   printf("Begin Check unrestarted GMRES.\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_gmres(lhs, rhs, N*N, 4000, 1e-6, square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check unrestarted GMRES.\n");
   printf("\n\n\n");
   
   
   printf("Begin Check restarted GMRES(8).\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_gmres_restart(lhs, rhs, N*N, 4000, 1e-6, 8,  square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check restarted GMRES(8).\n");
   printf("\n\n\n");
    
   printf("Begin Check SOR with omega = 0.1.\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_sor(lhs, rhs, N*N, 10000, 1e-6, 0.01, square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check SOR with omega = 0.1.\n");
   printf("\n\n\n");
    
   printf("Begin Check MinRes.\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_minres(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check MinRes.\n");
   printf("\n\n\n");
    
   printf("Begin Check MinRes (Relaxation param 0.8).\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   invif = minv_vector_minres(lhs, rhs, N*N, 10000, 1e-6, 0.8, square_laplacian, NULL, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check MinRes (Relaxation param 0.8).\n");
   printf("\n\n\n");
    
   printf("Begin Check Preconditioned CG (3 iter MinRes).\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   // Prepare MR preconditioner.
   minres_precond_struct_real mps; 
   mps.n_step = 3;
   mps.rel_res = 1e-15; // Make n_step the dominant factor. 
   mps.matrix_vector = square_laplacian; 
   mps.matrix_extra_data = NULL;
   // End Prepare MR preconditioner.
   invif = minv_vector_cg_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, minres_preconditioner, (void*)&mps, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check Preconditioned CG (3 iter MinRes).\n");
   printf("\n\n\n");
    
   printf("Begin Check Flexibly Preconditioned CG (1e-1 rel resid MinRes).\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   // Prepare MR preconditioner. 
   mps.n_step = 10000; // make rel_res the dominant factor.
   mps.rel_res = 1e-1; 
   mps.matrix_vector = square_laplacian; 
   mps.matrix_extra_data = NULL;
   // End Prepare MR preconditioner.
   invif = minv_vector_cg_flex_precond(lhs, rhs, N*N, 4000, 1e-6, square_laplacian, NULL, minres_preconditioner, (void*)&mps, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check Flexibly Preconditioned CG (1e-1 rel resid MinRes).\n");
   printf("\n\n\n");
    
   printf("Begin Check Restarted Flexibly Preconditioned CG(12) (0.8 rel resid MinRes).\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   // Prepare MR preconditioner. 
   mps.n_step = 10000; // make rel_res the dominant factor.
   mps.rel_res = 0.8; 
   mps.matrix_vector = square_laplacian; 
   mps.matrix_extra_data = NULL;
   // End Prepare MR preconditioner.
   invif = minv_vector_cg_flex_precond_restart(lhs, rhs, N*N, 4000, 1e-6, 12, square_laplacian, NULL, minres_preconditioner, (void*)&mps, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check Restarted Flexibly Preconditioned CG(12) (0.8 rel resid MinRes).\n");
   printf("\n\n\n");
    
   printf("Begin Check Variably Preconditioned GCR (1e-1 rel resid MinRes).\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   // Prepare MR preconditioner.
   mps.n_step = 10000; // Make rel_res the dominant factor. 
   mps.rel_res = 1e-1; // Make n_step the dominant factor. 
   mps.matrix_vector = square_laplacian; 
   mps.matrix_extra_data = NULL;
   // End Prepare MinRes preconditioner.
   invif = minv_vector_gcr_var_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, minres_preconditioner, (void*)&mps, &verb);
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check Variably Preconditioned GCR (1e-1 rel resid MinRes).\n");
   printf("\n\n\n");
    
   printf("Begin Check Variably Preconditioned GCR (8 iter GCR).\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   //Prepare  GCR preconditioner. 
   gcr_precond_struct_real gps; 
   gps.n_step = 8; 
   gps.rel_res = 1e-20; // Make n_step the dominant factor. 
   gps.matrix_vector = square_laplacian; 
   gps.matrix_extra_data = NULL;
   invif = minv_vector_gcr_var_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, gcr_preconditioner, (void*)&gps, &verb); /**/
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check Variably Preconditioned GCR (8 iter GCR).\n");
   printf("\n\n\n");
    
   printf("Begin Check Restarted Variably Preconditioned GCR(12) (8 iter GCR).\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   //Prepare  GCR preconditioner. 
   gps.n_step = 8; 
   gps.rel_res = 1e-20; // Make n_step the dominant factor. 
   gps.matrix_vector = square_laplacian; 
   gps.matrix_extra_data = NULL;
   invif = minv_vector_gcr_var_precond_restart(lhs, rhs, N*N, 10000, 1e-6, 12, square_laplacian, NULL, gcr_preconditioner, (void*)&gps, &verb); /**/
   if (invif.success == true)
   {
      printf("GOOD Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   else
   {
      printf("FAIL Iter: %d Resid: %.15e.\n", invif.iter, sqrt(invif.resSq));
   }
   explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check Restarted Variably Preconditioned GCR(12) (8 iter GCR).\n");
   printf("\n\n\n");
    
   printf("Begin Check Power Iteration.\n");
   initialize_test(lattice, lhs, rhs, check, N*N);
   eigif = eig_vector_poweriter(&eig, rhs, N*N, 10000, 1e-7, square_laplacian, NULL);
   if (eigif.success == true)
   {
      printf("GOOD Iter: %d RelResid: %.15e Eval: %.15e.\n", eigif.iter, eigif.relative_diff, eig);
   }
   else
   {
      printf("FAIL Iter: %d RelResid: %.15e Eval: %.15e.\n", eigif.iter, eigif.relative_diff, eig);
   }
   //explicit_resid = check_test(lhs, rhs, check, N*N, square_laplacian, NULL); 
   //printf("Explicit Resid: %.15e.\n", explicit_resid);
   printf("End Check Power Iteration.\n");
   printf("\n\n\n");
   
   // Free the lattice.
   delete[] lattice;
   delete[] lhs;
   delete[] rhs;
   delete[] check;
}

// Square lattice.
// Kinetic term for a 2D laplacian w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" allows us to generalize these functions later.
// It would become an internal structure in C++ code.
void square_laplacian(double* lhs, double* rhs, void* extra_data)
{
   // Declare variables.
   int i;
   int x,y;

   // For a 2D square lattice, the stencil is:
   //     |  0 -1  0 |
   //     | -1 +4 -1 |
   //     |  0 -1  0 |
   //
   // e2 = yhat
   // ^
   // | 
   // |-> e1 = xhat

   // Apply the stencil.
   for (i = 0; i < N*N; i++)
   {
      lhs[i] = 0.0;
      x = i%N; // integer mod.
      y = i/N; // integer divide.
      
      // + e1.
      lhs[i] = lhs[i]-rhs[y*N+((x+1)%N)];
     
      // - e1.
      lhs[i] = lhs[i]-rhs[y*N+((x+N-1)%N)]; // The extra +N is because of the % sign convention.
      
      // + e2.
      lhs[i] = lhs[i]-rhs[((y+1)%N)*N+x];
    
      // - e2.
      lhs[i] = lhs[i]-rhs[((y+N-1)%N)*N+x];

      // 0
      // Added mass term here.
      lhs[i] = lhs[i]+(4+MASS)*rhs[i];
   }
       
}


// Zero out vectors, set a point source.
void initialize_test(double* lattice, double* lhs, double* rhs, double* check, int size)
{
   int i;
   int half_size = (int)(sqrt(size)/2+0.5);
   
   zero<double>(lattice,size);
   zero<double>(lhs,size);
   zero<double>(rhs,size);
   zero<double>(check,size);

   // Set a point on the rhs.
   rhs[half_size+half_size*half_size*2] = 1.0;
   //rhs[N/2+(N/2)*N] = 1.0;
   
   // Set a point on the lhs.
   lhs[half_size+half_size*half_size*2] = 1.0;
}

// Check solution.
double check_test(double* lhs, double* rhs, double* check, int size, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
   int i;
   double explicit_resid = 0.0;
   
   // Check and make sure we get the right answer.
   matrix_vector(check, lhs, extra_info);
   
   for (i = 0; i < size; i++)
   {
      explicit_resid += (rhs[i] - check[i])*(rhs[i] - check[i]);
   }
   explicit_resid = sqrt(explicit_resid);
   
   return explicit_resid; 
}


