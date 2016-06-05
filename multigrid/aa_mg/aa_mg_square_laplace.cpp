// This is a sample code that constructs general kinetic terms!

#include <iostream>
#include <iomanip> // to set output precision.
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_inverters.h"
#include "generic_inverters_precond.h"
#include "generic_vector.h"

// Do restrict/prolong test?
#define PDAGP_TEST

// Do two vector restrict/prolong?
#define PDAGP_2TEST

using namespace std; 

// For now, define the length in a direction.
#define N 8

// Define pi.
#define PI 3.141592653589793

// Define mass.
#define MASS 0.1*0.1

// What's the X blocksize?
#define X_BLOCKSIZE 2
// What's the Y blocksize?
#define Y_BLOCKSIZE 2

// Square laplacian function.
void square_laplacian(double* lhs, double* rhs, void* extra_data);

// General multigrid projector function!
void coarse_square_laplacian(double* lhs, double* rhs, void* extra_data); 

// Useful mg functions.
struct mg_operator_struct_real;

// Normalize the projectors based on the coarse blocksize. 
void block_normalize(mg_operator_struct_real* mgstruct);

// Prolong from the coarse lattice to the fine lattice. 
void prolong(double* x_fine, double* x_coarse, mg_operator_struct_real* mgstruct);

// Restrict a fine vector to a coarse vector using the info in mgstruct.
void restrict(double* x_coarse, double* x_fine, mg_operator_struct_real* mgstruct);

// Multigrid operator struct.
struct mg_operator_struct_real
{
    int blocksize_x; // How much to block in x direction.
    int blocksize_y; // How much to block in y direction. 
    int n_vector; // Number of vectors. 
    double** projectors; // Holds the projectors. First index n_vector, second size.
    void (*matrix_vector)(double*, double*, void*);
    void* matrix_extra_data; 
};

enum inner_solver
{
    MINRES = 0,
    CG = 1,
    GCR = 2
};

// Preconditioning struct
struct mg_precond_struct_real;

// MG preconditioner!
void mg_preconditioner(double* lhs, double* rhs, int size, void* extra_data);

int main(int argc, char** argv)
{  
    // Declare some variables.
    cout << setiosflags(ios::fixed) << setprecision(6);
    int i, j, x, y;
    double *lattice; // At some point, I'll have a generic (template) lattice class.
    double *lhs, *rhs, *check; // For some Kinetic terms.
    double explicit_resid = 0.0;
    double bnorm = 0.0;
    inversion_info invif;
    
    // Describe the fine lattice. 
    int x_fine = N;
    int y_fine = N;
    int fine_size = x_fine*y_fine;
    
    // Build an mg_struct.
    mg_operator_struct_real mgstruct;
    mgstruct.blocksize_x = X_BLOCKSIZE;
    mgstruct.blocksize_y = Y_BLOCKSIZE;
    mgstruct.n_vector = 1;
    mgstruct.matrix_vector = square_laplacian;
    mgstruct.matrix_extra_data = NULL; 
    
    // Describe the coarse lattice. 
    int x_coarse = x_fine/mgstruct.blocksize_x; // how many coarse sites are there in the x dir?
    int y_coarse = y_fine/mgstruct.blocksize_y; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 

    // Initialize the lattice. Indexing: index = y*N + x.
    lattice = new double[fine_size];
    lhs = new double[fine_size];
    rhs = new double[fine_size];   
    check = new double[fine_size];   
    // Zero it out.
    for (i=0; i < fine_size; i++)
    {
      lattice[i] = lhs[i] = rhs[i] = check[i] = 0;
    }

    // Set a point on the rhs.
    rhs[x_fine/2+(y_fine/2)*x_fine] = 1.0;

    // Get norm for rhs.
    for (i=0;i<fine_size;i++)
    {
     bnorm = bnorm + rhs[i]*rhs[i];
    }
    bnorm = sqrt(bnorm);

    // Set a point on the lhs.
    lhs[x_fine/2+(y_fine/2)*x_fine+1] = 1.0;
    
    // Create a projector.
    mgstruct.projectors = new double*[mgstruct.n_vector];
    mgstruct.projectors[0] = new double[fine_size];

    cout << "Creating " << mgstruct.n_vector << " projector(s).\n";
    // Make a constant projector.
    for (i = 0; i < N*N; i++)
    {
        mgstruct.projectors[0][i] = 1;
    }
    
#ifdef PDAGP_TEST
    {
        // Begin PdagP test.

        // Print the fine rhs.
        cout << "Check fine point src:\n"; 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << rhs[x+y*x_fine] << " ";
            }
            cout << "\n";
        }

        cout << "Check projector:\n"; 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << mgstruct.projectors[0][x+y*x_fine] << " ";
            }
            cout << "\n";
        }

        // Test block normalizing the projector.
        cout << "Block normalize the projector(s). X blocks: " << x_coarse << " Y blocks: " << y_coarse << "\n";
        block_normalize(&mgstruct); 

        cout << "Check block normalized vector:\n"; 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << mgstruct.projectors[0][x+y*x_fine] << " ";
            }
            cout << "\n";
        }


        // Restrict the original source. 
        double* rhs_coarse = new double[coarse_size*mgstruct.n_vector];
        cout << "Restricting the source.\n";
        restrict(rhs_coarse, rhs, &mgstruct);

        // Check coarse source.
        cout << "Check coarse src:\n"; 
        for (int y = 0; y < y_coarse; y++)
        {
            for (int x = 0; x < x_coarse; x++)
            {
                cout << rhs_coarse[x+y*x_coarse] << " ";
            }
            cout << "\n";
        }

        // Prolong the restricted source. 
        double *rhs_PdagP = new double[fine_size];
        cout << "Prolonging the restricted source.\n";
        prolong(rhs_PdagP, rhs_coarse, &mgstruct);

        cout << "Check re-prolonged source vector:\n"; 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << rhs_PdagP[x+y*x_fine] << " ";
            }
            cout << "\n";
        }


        // Try applying the coarse Laplace operator!
        double* rhs_A_coarse = new double[coarse_size*mgstruct.n_vector]; 
        zero<double>(rhs_A_coarse, coarse_size*mgstruct.n_vector); 

        cout << "Applying the coarse Laplace operator to coarse source.\n";
        coarse_square_laplacian(rhs_A_coarse, rhs_coarse, (void*)&mgstruct);

        // Check A coarse source.
        cout << "Check A times coarse src:\n"; 
        for (int y = 0; y < y_coarse; y++)
        {
            for (int x = 0; x < x_coarse; x++)
            {
                cout << rhs_A_coarse[x+y*x_coarse] << " ";
            }
            cout << "\n";
        }

        // Prolong rhs_A_coarse.
        cout << "Prolong A times coarse source.\n";
        double* rhs_PAP_fine = new double[fine_size];
        zero<double>(rhs_PAP_fine, fine_size);
        prolong(rhs_PAP_fine, rhs_A_coarse, &mgstruct);

        // Check PAP. 
        cout << "Check PAP on source.\n";
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << rhs_PAP_fine[x+y*x_fine] << " ";
            }
            cout << "\n";
        }
        
        delete[] rhs_coarse;
        delete[] rhs_A_coarse;
        delete[] rhs_PdagP;
        delete[] rhs_PAP_fine;
    }
#endif
    
    
#ifdef PDAGP_2TEST
    
    // Test adding a second projector.
    {
        mgstruct.n_vector = 2; 
        double* tmp_store = mgstruct.projectors[0];
        delete[] mgstruct.projectors;
        mgstruct.projectors = new double*[mgstruct.n_vector];
        mgstruct.projectors[0] = tmp_store; 
        mgstruct.projectors[1] = new double[fine_size];

        // Add an even/odd vector. 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                if ((x+y)%2 == 0)
                    mgstruct.projectors[1][x+y*x_fine] = 1;
                else
                    mgstruct.projectors[1][x+y*x_fine] = -1;
            }
        }

        cout << "Check projector:\n"; 
        for (int n = 0; n < mgstruct.n_vector; n++)
        {
            cout << "Vector " << n << "\n";
            for (int y = 0; y < y_fine; y++)
            {
                for (int x = 0; x < x_fine; x++)
                {
                    cout << mgstruct.projectors[n][x+y*x_fine] << " ";
                }
                cout << "\n";
            }
        }

        // Test block normalizing the projector.
        cout << "Block normalize the projector(s). X blocks: " << x_coarse << " Y blocks: " << y_coarse << "\n";
        block_normalize(&mgstruct); 

        cout << "Check block normalized vector:\n"; 
        for (int n = 0; n < mgstruct.n_vector; n++)
        {
            cout << "Vector " << n << "\n";
            for (int y = 0; y < y_fine; y++)
            {
                for (int x = 0; x < x_fine; x++)
                {
                    cout << mgstruct.projectors[n][x+y*x_fine] << " ";
                }
                cout << "\n";
            }
        }

        // Restrict the original source. 
        double* rhs_coarse = new double[coarse_size*mgstruct.n_vector];
        cout << "Restricting the source.\n";
        restrict(rhs_coarse, rhs, &mgstruct);

        // Check coarse source.
        cout << "Check coarse src:\n"; 
        for (int y = 0; y < y_coarse; y++)
        {
            for (int x = 0; x < x_coarse; x++)
            {
                cout << "(";
                for (int n = 0; n < mgstruct.n_vector; n++)
                {
                    cout << rhs_coarse[(x+y*x_coarse)*mgstruct.n_vector+n] << ",";
                }
                cout << ") ";
            }
            cout << "\n";
        }

        // Prolong the restricted source. 
        double* rhs_PdagP = new double[fine_size];
        cout << "Prolonging the restricted source.\n";
        prolong(rhs_PdagP, rhs_coarse, &mgstruct);

        cout << "Check re-prolonged source vector:\n"; 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << rhs_PdagP[x+y*x_fine] << " ";
            }
            cout << "\n";
        }
        
        delete[] rhs_coarse;
        delete[] rhs_PdagP;
    }
    
#endif
    
    for (i = 0; i < mgstruct.n_vector; i++)
    {
        delete[] mgstruct.projectors[i];
    }
    delete[] mgstruct.projectors; 


    
    // Let's actually test a multigrid solve!
    
    
    
    return 0; 
        

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
    //invif = minv_vector_gcr(lhs, rhs, N*N, 4000, 1e-8, square_laplacian, NULL);
    //invif = minv_vector_gcr_restart(lhs, rhs, N*N, 4000, 1e-8, 20, square_laplacian, NULL);

    // Indentity preconditioner.
    //invif = minv_vector_cg_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, identity_preconditioner, NULL); 

    // Minres preconditioner.
    /*minres_precond_struct_real mps; 
    mps.n_step = 10;
    mps.rel_res = 1e-15; // Make n_step the dominant factor. 
    mps.matrix_vector = square_laplacian; 
    mps.matrix_extra_data = NULL;
    invif = minv_vector_cg_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, minres_preconditioner, (void*)&mps); /**/

    // Identity preconditioner to flex cg.
    //invif = minv_vector_cg_flex_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, identity_preconditioner, NULL);

    // Minres preconditioner to flex cg. 
    /*minres_precond_struct_real mps; 
    mps.n_step = 10000; // Make rel_res the dominant factor. 
    mps.rel_res = 1e-1; 
    mps.matrix_vector = square_laplacian; 
    mps.matrix_extra_data = NULL;
    invif = minv_vector_cg_flex_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, minres_preconditioner, (void*)&mps); /**/

    // Restarted CG with GCR preconditioner. 
    /*gcr_precond_struct_real gps; 
    gps.n_step = 10000; // Make rel_res the dominant factor.
    gps.rel_res = 0.8; // Make n_step the dominant factor. 
    gps.matrix_vector = square_laplacian; 
    gps.matrix_extra_data = NULL;
    invif = minv_vector_cg_flex_precond_restart(lhs, rhs, N*N, 10000, 1e-6, 12, square_laplacian, NULL, gcr_preconditioner, (void*)&gps); /**/


    // Indentity preconditioner. 
    //invif = minv_vector_gcr_var_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, identity_preconditioner, NULL); 

    // Minres preconditioner. 
    /*minres_precond_struct_real mps; 
    mps.n_step = 10000; // Make rel_res the dominant factor. 
    mps.rel_res = 1e-1; 
    mps.matrix_vector = square_laplacian; 
    mps.matrix_extra_data = NULL;
    invif = minv_vector_gcr_var_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, minres_preconditioner, (void*)&mps); /**/

    // GCR preconditioner. 
    /*gcr_precond_struct_real gps; 
    gps.n_step = 8; 
    gps.rel_res = 1e-20; // Make n_step the dominant factor. 
    gps.matrix_vector = square_laplacian; 
    gps.matrix_extra_data = NULL;
    invif = minv_vector_gcr_var_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, gcr_preconditioner, (void*)&gps); /**/

    // Restarted GCR with preconditioner. 
    /*gcr_precond_struct_real gps; 
    gps.n_step = 8; 
    gps.rel_res = 1e-20; // Make n_step the dominant factor. 
    gps.matrix_vector = square_laplacian; 
    gps.matrix_extra_data = NULL;
    invif = minv_vector_gcr_var_precond_restart(lhs, rhs, N*N, 10000, 1e-6, 12, square_laplacian, NULL, gcr_preconditioner, (void*)&gps); /**/

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

// Properly normalize the P vectors.
void block_normalize(mg_operator_struct_real* mgstruct)
{
    int n, i;
    int x_fine = N;
    int y_fine = N;
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
    int x_fine = N;
    int y_fine = N;
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
    int x_fine = N;
    int y_fine = N;
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

// General multigrid projector function!
void coarse_square_laplacian(double* lhs, double* rhs, void* extra_data)
{
    // Iterators.
    int i, j, k, n; 
    int tmp; 
    
    // Grab the mg_precond_struct.
    mg_operator_struct_real mgstruct = *(mg_operator_struct_real*)extra_data; 
    
    // lhs and rhs are of size coarse_size. mgstruct.matrix_vector expects
    // fine_size. 
    int fine_size = N*N;
    int coarse_size = fine_size*mgstruct.n_vector/(mgstruct.blocksize_x*mgstruct.blocksize_y); 
    
    // Okay... how the hell are we going to do this. 
    double* Px; // Holds prolonged current solution.
    double* APx; // Holds A times prolonged current solution.
    
    Px = new double[fine_size];
    APx = new double[fine_size];
    
    zero<double>(Px, fine_size); zero<double>(APx, fine_size); 
    
    // Prolong. 
    prolong(Px, rhs, &mgstruct);
    
    // Apply the original matrix.
    (*mgstruct.matrix_vector)(APx, Px, mgstruct.matrix_extra_data);
    
    // Restrict. 
    zero<double>(lhs, coarse_size);
    restrict(lhs, APx, &mgstruct); 
    
    delete[] Px;
    delete[] APx; 
    
}

struct mg_precond_struct_real
{
    // How many MinRes pre-smooth steps?
    int n_pre_smooth;
    
    // How many MinRes post-smooth steps?
    int n_post_smooth;
    
    // What inner solver should we use?
    inner_solver in_solve_type; // MINRES, CG, or GCR
    int n_step; // Max steps for inner solver?
    double rel_res; // Rel_res for inner solver?
    
    // What's the mg_info?
    mg_operator_struct_real* mgstruct; 
    
    // What's matrix function are we dealing with?
    void (*matrix_vector)(double*, double*, void*);
    void* matrix_extra_data; 
};

// MG preconditioner!! (Man, I'm excited!
void mg_preconditioner(double* lhs, double* rhs, int size, void* extra_data)
{
    mg_precond_struct_real* mgprecond = (mg_precond_struct_real*)extra_data; 
    
    // GET EXCITED! First off, let's do some pre-smoothing. 
    
}


