// This is a sample code that constructs general kinetic terms!

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_inverters.h"
#include "generic_inverters_precond.h"
#include "generic_vector.h"

using namespace std; 

// For now, define the length in a direction.
#define N 16

// Define pi.
#define PI 3.141592653589793

// Define mass.
#define MASS 0.1*0.1

// Square laplacian w/out u1 function.
void square_laplacian(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square laplacian w/ u1. 
void square_laplacian_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Load complex gauge field from file. 
// Based on Rich Brower's u1 gauge routines. 
// Reads in a U1 phase lattice from file, returns complex fields. 
// Rich's code has 'y' as the fast direction. Need to transpose!
void read_lattice_u1(complex<double>* gauge_field, string input_file)
{
   double phase_tmp;   
   fstream in_file;
   
   in_file.open(input_file,ios::in); 
   for(int x =0;x< N;x++)
   {
      for(int y =0;y< N;y++)
      {
         for(int mu=0; mu<2; mu++)
         {
            in_file >> phase_tmp;
            gauge_field[y*2*N+x*2+mu] = polar(1.0,phase_tmp);
            //cout << polar(1.0, phase_tmp) << "\n";
         }
      }
   }
   in_file.close(); 
   
   return;
}

int main(int argc, char** argv)
{  
   // Declare some variables.
   int i, j;
   complex<double> *lattice, *lhs, *rhs, *check; // For some Kinetic terms.
   double explicit_resid = 0.0;
   double bnorm = 0.0;
   inversion_info invif;
 
   // Initialize the gauge field. Indexing: index = y*2*N + x*2 + mu.
   lattice = new complex<double>[N*N*2]; 
   // Load lattice
   read_lattice_u1(lattice, "./cfg/phase16b16.dat"); 
    
   // Unit gauge.
   //for (i=0; i < N*N; i++) { lattice[2*i] = lattice[2*i+1] = 1.0; }
    
   // Initialize the fermions. Indexing: index = y*N + x.
   lhs = new complex<double>[N*N];
   rhs = new complex<double>[N*N];   
   check = new complex<double>[N*N];   
   // Zero it out.
   for (i=0; i < N*N; i++)
   {
      lhs[i] = rhs[i] = check[i] = 0;
   }

   // Set a point on the rhs.
   rhs[N/2+(N/2)*N] = 1.0;
   
   // Get norm for rhs.
   bnorm = sqrt(norm2sq<double>(rhs, N*N));
   
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
   // 6: function pointer
   // 7: "extra data": can set this to not-null to pass in gauge fields, etc.
   
   //invif = minv_vector_cg(lhs, rhs, N*N, 4000, 1e-6, square_laplacian_u1, (void*)lattice); 
   //invif = minv_vector_bicgstab(lhs, rhs, N*N, 4000, 1e-6, square_laplacian_u1, (void*)lattice); 
   //invif = minv_vector_gcr(lhs, rhs, N*N, 4000, 1e-6, square_laplacian_u1, (void*)lattice); 
   //invif = minv_vector_gmres_norestart(lhs, rhs, N*N, 4000, 1e-6, square_laplacian_u1, (void*)lattice);
   //invif = minv_vector_minres(lhs, rhs, N*N, 4000, 1e-6, square_laplacian_u1, (void*)lattice); 
    
    
   // Identity preconditioned.
   //invif = minv_vector_cg_precond(lhs, rhs, N*N, 1000, 1e-6, square_laplacian_u1, (void*)lattice, identity_preconditioner, (void*)&mps);
    
   // MinRes preconditioned.
   
   minres_precond_struct_complex mps; 
   mps.n_step = 3;
   mps.rel_res = 1e-15; // Fix to 3 iters. 
   mps.matrix_vector = square_laplacian_u1; 
   mps.matrix_extra_data = (void*)lattice;
   invif = minv_vector_cg_precond(lhs, rhs, N*N, 1000, 1e-6, square_laplacian_u1, (void*)lattice, minres_preconditioner, (void*)&mps); /**/
    
   
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
   //square_laplacian(check, lhs, NULL);
   square_laplacian_u1(check, lhs, (void*)lattice); 
   
   for (i = 0; i < N*N; i++)
   {
      explicit_resid += real(conj(rhs[i] - check[i])*(rhs[i] - check[i]));
   }
   explicit_resid = sqrt(explicit_resid)/bnorm;
   
   printf("[check] should equal [rhs]. The relative residual is %15.20e.\n", explicit_resid);

   // Free the lattice.
   delete[] lattice;
   delete[] lhs;
   delete[] rhs;
   delete[] check;
}

// Square lattice.
// Kinetic term for a 2D laplacian w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" doesn't include anything.
void square_laplacian(complex<double>* lhs, complex<double>* rhs, void* extra_data)
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

// Square lattice.
// Kinetic term for a 2D laplacian w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" is a cast to a complex gauge_field[N*N*2], 
//    loaded by the function read_lattice_u1. 
void square_laplacian_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
   // Declare variables.
   int i;
   int x,y;
    
   complex<double>* lattice = (complex<double>*)extra_data; 

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
      lhs[i] = lhs[i]-conj(lattice[y*N*2+x*2])*rhs[y*N+((x+1)%N)];
     
      // - e1.
      lhs[i] = lhs[i]-lattice[y*N*2+((x+N-1)%N)*2]*rhs[y*N+((x+N-1)%N)]; // The extra +N is because of the % sign convention.
      
      // + e2.
      lhs[i] = lhs[i]-conj(lattice[y*N*2+x*2+1])*rhs[((y+1)%N)*N+x];
    
      // - e2.
      lhs[i] = lhs[i]-lattice[((y+N-1)%N)*N*2+x*2+1]*rhs[((y+N-1)%N)*N+x];

      // 0
      // Added mass term here.
      lhs[i] = lhs[i]+(4+MASS)*rhs[i];
   }
       
}
