// This is a sample code that constructs general kinetic terms!

#include <iostream>
#include <iomanip> // to set output precision.
#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <random>

#include "generic_inverters.h"
#include "generic_inverters_precond.h"
#include "generic_vector.h"
#include "mg.h"
#include "mg_real.h"

// Do restrict/prolong test?
//#define PDAGP_TEST

// Do two vector restrict/prolong?
//#define PDAGP_2TEST

using namespace std; 

// For now, define the length in a direction.
#define N 8

// Define pi.
#define PI 3.141592653589793

// Define mass.
#define MASS 0.01*0.01

// What's the X blocksize?
#define X_BLOCKSIZE 2
// What's the Y blocksize?
#define Y_BLOCKSIZE 2

// Print null vectors?
//#define PRINT_NULL_VECTOR

// Do null vector generation? Currently uses BiCGStab
#define GEN_NULL_VECTOR

// How many BiCGStab iterations do we use?
#define GEN_NULL_VECTOR_STEP 10000
#define GEN_NULL_VECTOR_REL_RESID 1e-8

// 1 for just const vector, 2 for const + even/odd vector. Specifies
// the number to generate if GEN_NULL_VECTOR is defined
#define VECTOR_COUNT 2

// Square laplacian function.
void square_laplacian(double* lhs, double* rhs, void* extra_data);


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
    mgstruct.x_fine = N;
    mgstruct.y_fine = N; 
    mgstruct.blocksize_x = X_BLOCKSIZE;
    mgstruct.blocksize_y = Y_BLOCKSIZE;
    mgstruct.n_vector = VECTOR_COUNT;
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
    for (i = 0; i < mgstruct.n_vector; i++)
    {
        mgstruct.projectors[i] = new double[fine_size];
    }

    cout << "Creating " << mgstruct.n_vector << " projector(s).\n";
#ifndef GEN_NULL_VECTOR
    // Make a constant projector.
    for (i = 0; i < N*N; i++)
    {
        mgstruct.projectors[0][i] = 1;
        /*if (mgstruct.n_vector > 1)
        {
            x = i % N;
            y = i / N;
            mgstruct.projectors[1][i] = ((x+y)%2 == 0) ? 1 : -1;
        }*/
        // For block orthonormalize test.
        if (mgstruct.n_vector > 1)
        {
            x = i % N;
            y = i / N;
            mgstruct.projectors[1][i] = (x%2 == 0) ? 1 : -0.5;
        }
    }
    
#ifdef PRINT_NULL_VECTOR
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
#endif // PRINT_NULL_VECTOR

#else // generate null vector
    // Create an RNG. 
    std::mt19937 generator (1337u); // 1337u is the seed.
    
    // We generate null vectors by solving Ax = 0, with a
    // gaussian initial guess.
    // For sanity with the residual, we really solve Ax = -Ax_0,
    // where x has a zero initial guess, x_0 is a random vector.
    double* rand_guess = new double[fine_size];
    double* Arand_guess = new double[fine_size];
    
    for (i = 0; i < mgstruct.n_vector; i++)
    {
        // Use mgstruct.projectors[i] has a temporary holder.
        gaussian<double>(rand_guess, fine_size, generator);
        
        zero<double>(Arand_guess, fine_size); 
        square_laplacian(Arand_guess, rand_guess, NULL);
        for (j = 0; j < fine_size; j++)
        {
           Arand_guess[j] = -Arand_guess[j]; 
        }
        zero<double>(mgstruct.projectors[i], fine_size);
        
        minv_vector_bicgstab(mgstruct.projectors[i], Arand_guess, fine_size, GEN_NULL_VECTOR_STEP, GEN_NULL_VECTOR_REL_RESID, square_laplacian, NULL); 
        for (j = 0; j < fine_size; j++)
        {
            mgstruct.projectors[i][j] += rand_guess[j];
        }
        
        normalize(mgstruct.projectors[i], fine_size); 
    }
    
    delete[] rand_guess; 
    delete[] Arand_guess; 
    
    // Normalize projectors.
    for (i = 0; i < mgstruct.n_vector; i++)
    {
        normalize(mgstruct.projectors[i], fine_size);
    }
    
#ifdef PRINT_NULL_VECTOR
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
#endif // PRINT_NULL_VECTOR
    
#endif // generate null vector. 
    cout << "Performing block orthonormalize of null vectors...\n";
    block_orthonormalize(&mgstruct); 
    
    #ifdef PRINT_NULL_VECTOR
    cout << "\nCheck projector:\n"; 
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
#endif // PRINT_NULL_VECTOR
    
    
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
        coarse_square_laplace(rhs_A_coarse, rhs_coarse, (void*)&mgstruct);

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
    
    cout << "Solving coarse system only.\n";
    double* rhs_coarse = new double[coarse_size*mgstruct.n_vector];
    zero<double>(rhs_coarse, coarse_size*mgstruct.n_vector);
    restrict(rhs_coarse, rhs, &mgstruct);
    
    double* lhs_coarse = new double[coarse_size*mgstruct.n_vector];
    zero<double>(lhs_coarse, coarse_size*mgstruct.n_vector);
    
    invif = minv_vector_cg(lhs_coarse, rhs_coarse, coarse_size*mgstruct.n_vector, 10000, 1e-6, coarse_square_laplace, (void*)&mgstruct);
    
    
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
    double* A_lhs_coarse = new double[coarse_size*mgstruct.n_vector];
    zero<double>(A_lhs_coarse, coarse_size*mgstruct.n_vector);
    
    coarse_square_laplace(A_lhs_coarse, lhs_coarse, (void*)&mgstruct);

    for (i = 0; i < coarse_size*mgstruct.n_vector; i++)
    {
      explicit_resid += (A_lhs_coarse[i] - rhs_coarse[i])*(A_lhs_coarse[i] - rhs_coarse[i]);
    }
    explicit_resid = sqrt(explicit_resid);

    printf("[check] should equal [rhs]. The residual is %15.20e.\n", explicit_resid);
    
    double* pro_lhs_coarse = new double[N*N];
    zero<double>(pro_lhs_coarse, N*N);
    prolong(pro_lhs_coarse, lhs_coarse, &mgstruct); 
    double* pro_rhs_coarse = new double[N*N];
    zero<double>(pro_rhs_coarse, N*N);
    square_laplacian(pro_rhs_coarse, pro_lhs_coarse, NULL);
    
    // Try a direct solve.
    cout << "\nSolve fine system.\n";
    
    invif = minv_vector_cg(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL);
    
    if (invif.success == true)
    {
     printf("Algorithm %s took %d iterations to reach a relative residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq)/bnorm);
    }
    else // failed, maybe.
    {
     printf("Potential error! Algorithm %s took %d iterations to reach a relative residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq)/bnorm);
     printf("This may be because the max iterations was reached.\n");
    }
    
    // Compare PAP solution to real solution. 
    cout << "\nCompare solutions.\n";
    double comparison = 0;
    double resid_comparison = 0;
    for (i = 0; i < N*N; i++)
    {
        comparison += (pro_lhs_coarse[i]-lhs[i])*(pro_lhs_coarse[i]-lhs[i]);
        resid_comparison += (pro_rhs_coarse[i]-rhs[i])*(pro_rhs_coarse[i]-rhs[i]);
    }
    comparison = sqrt(explicit_resid);
    printf("The solutions deviate by %15.20e.\n", comparison);
    printf("The projected residual has a rel res of %15.20e.\n", sqrt(resid_comparison)/bnorm);
    
    delete[] rhs_coarse; 
    delete[] lhs_coarse;
    delete[] A_lhs_coarse; 
    delete[] pro_lhs_coarse; 
    delete[] pro_rhs_coarse; 
    

    // Let's actually test a multigrid solve!
    cout << "\nTest MG solve.\n";
    
    // Block normalize the null vectors.
    block_normalize(&mgstruct); 
    
    // Set up the MG preconditioner. 
    mg_precond_struct_real mgprecond;
    
    mgprecond.n_pre_smooth = 6; // 3 MinRes smoother steps before coarsening.
    mgprecond.n_post_smooth = 6; // 3 MinRes smoother steps after refining.
    mgprecond.in_solve_type = CG; // What inner solver? MINRES, CG, or GCR.
    mgprecond.n_step = 10000; // max number of steps to use for inner solver.
    mgprecond.rel_res = 1e-1; // Maximum relative residual for inner solver.
    mgprecond.mgstruct = &mgstruct; // Contains null vectors, fine operator. (Since we don't construct the fine op.)
    mgprecond.matrix_vector = coarse_square_laplace; // Function which applies the coarse operator. 
    mgprecond.matrix_extra_data = (void*)&mgstruct; // What extra_data the coarse operator expects. 
    
    // Well, maybe this will work?
    zero<double>(lhs, N*N);
    invif = minv_vector_cg_flex_precond(lhs, rhs, N*N, 10000, 1e-6, square_laplacian, NULL, mg_preconditioner, (void*)&mgprecond); /**/
    
    if (invif.success == true)
    {
     printf("Algorithm %s took %d iterations to reach a relative residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq)/bnorm);
    }
    else // failed, maybe.
    {
     printf("Potential error! Algorithm %s took %d iterations to reach a relative residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq)/bnorm);
     printf("This may be because the max iterations was reached.\n");
    }


    printf("Computing [check] = A [lhs] as a confirmation.\n");

    // Check and make sure we get the right answer.
    square_laplacian(check, lhs, NULL);

    explicit_resid = 0.0;
    for (i = 0; i < N*N; i++)
    {
      explicit_resid += (rhs[i] - check[i])*(rhs[i] - check[i]);
    }
    explicit_resid = sqrt(explicit_resid)/bnorm;

    printf("[check] should equal [rhs]. The relative residual is %15.20e.\n", explicit_resid);

    // Free the lattice.
    delete[] lattice;
    delete[] lhs;
    delete[] rhs;
    delete[] check;
    
    // Clean up!
    for (i = 0; i < mgstruct.n_vector; i++)
    {
        delete[] mgstruct.projectors[i];
    }
    delete[] mgstruct.projectors; 
    
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


