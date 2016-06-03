// This is a sample code that constructs general kinetic terms!

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_inverters.h"
#include "generic_inverters_precond.h"
#include "generic_vector.h"


using namespace std; 

// For now, define the length in a direction.
#define N 128

// Define pi.
#define PI 3.141592653589793

// Define mass.
#define MASS 0.1*0.1

// Square laplacian function.
void square_laplacian(double* lhs, double* rhs, void* extra_data);



int main(int argc, char** argv)
{  
   // Declare some variables.
   int i, j;
   double *lattice; // At some point, I'll have a generic (template) lattice class.
   double *lhs, *rhs, *check; // For some Kinetic terms.
   double explicit_resid = 0.0;
   double bnorm = 0.0;
   inversion_info invif;
 
   // Initialize the lattice. Indexing: index = y*N + x.
   lattice = new double[N*N];
   lhs = new double[N*N];
   rhs = new double[N*N];   
   check = new double[N*N];   
   // Zero it out.
   for (i=0; i < N*N; i++)
   {
      lattice[i] = lhs[i] = rhs[i] = check[i] = 0;
   }

   // Set a point on the rhs.
   rhs[N/2+(N/2)*N] = 1.0;
   
   // Get norm for rhs.
   for (i=0;i<N*N;i++)
   {
     bnorm = bnorm + rhs[i]*rhs[i];
   }
   bnorm = sqrt(bnorm);
   
   // Set a point on the lhs.
   lhs[N/2+(N/2)*N+1] = 1.0;
   
   printf("Solving A [lhs] = [rhs] for lhs, using a point source.\n");
   
   // lhs = A^(-1) rhs
   // Arguments:
   // 1: lhs
   // 2: rhs
   // 3: size of vector
   // 4: maximum iterations
   // 5: residual
   // 5a for gmres_restart: how often to restart.
   // 5a for sor: overrelaxation parameter. Set to omega*largest eigenvalue < 1. 
   // 6: function pointer
   // 7: "extra data": can set this to not-null to pass in gauge fields, etc.
   
   //invif = minv_vector_cg(lhs, rhs, N*N, 4000, 1e-6, square_laplacian, NULL);
   //invif = minv_vector_bicgstab(lhs, rhs, N*N, 4000, 1e-8, square_laplacian, NULL);
   //invif = minv_vector_gmres_norestart(lhs, rhs, N*N, 4000, 1e-8, square_laplacian, NULL);
   //invif = minv_vector_gmres_restart(lhs, rhs, N*N, 4000, 1e-8, 8, square_laplacian, NULL);
   //invif = minv_vector_gmres_restart(lhs, rhs, N*N, 4000, 1e-8, 20, square_laplacian, NULL);
   //invif = minv_vector_sor(lhs, rhs, N*N, 10000, 1e-6, 0.1, square_laplacian, NULL);
    
   // Indentity preconditioner.
   //invif = minv_vector_cg_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, identity_preconditioner, NULL); 
    
   // Minres preconditioner.
   /*minres_precond_struct_real mps; 
   mps.n_step = 10;
   mps.rel_res = 1e-15; // Make n_step the dominant factor. 
   mps.matrix_vector = square_laplacian; 
   mps.matrix_extra_data = NULL;
   invif = minv_vector_cg_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, minres_preconditioner, (void*)&mps); /**/
    
   // Indentity preconditioner. 
   //invif = minv_vector_gcr_var_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, identity_preconditioner, NULL); 
   
   // Minres preconditioner. 
   minres_precond_struct_real mps; 
   mps.n_step = 10000; // Make rel_res the dominant factor. 
   mps.rel_res = 1e-1; // Make n_step the dominant factor. 
   mps.matrix_vector = square_laplacian; 
   mps.matrix_extra_data = NULL;
   invif = minv_vector_gcr_var_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, minres_preconditioner, (void*)&mps); /**/
    
   if (invif.success == true)
   {
     printf("Algorithm %s took %d iterations to reach a residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq));
   }
   else // failed, maybe.
   {
     printf("Potential error! Algorithm %s took %d iterations to reach a residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq));
     printf("This may be because the max iterations was reached.\n");
   }
   

   
   printf("Computing [check] = A [lhs] as a confirmation.\n");

   // Check and make sure we get the right answer.
   square_laplacian(check, lhs, NULL);
   
   for (i = 0; i < N*N; i++)
   {
      explicit_resid += (rhs[i] - check[i])*(rhs[i] - check[i]);
   }
   explicit_resid = sqrt(explicit_resid);
   
   printf("[check] should equal [rhs]. The residual is %15.20e.\n", explicit_resid);

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


