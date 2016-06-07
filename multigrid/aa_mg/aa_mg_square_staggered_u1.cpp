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
#include "mg_complex.h"
#include "u1_utils.h"

// Do restrict/prolong test?
//#define PDAGP_TEST

// Do two vector restrict/prolong?
//#define PDAGP_2TEST

// Try solving just the coarse solver. 
//#define COARSE_ONLY

using namespace std; 

// For now, define the length in a direction.
#define N 8

// Define pi.
#define PI 3.141592653589793

// Define mass.
#define MASS 0.1

// standard deviation of angle for random gauge field.
//#define BETA 100.0

// What's the X blocksize?
#define X_BLOCKSIZE 4
// What's the Y blocksize?
#define Y_BLOCKSIZE 4

// Print null vectors?
#define PRINT_NULL_VECTOR

// Do null vector generation? Currently uses BiCGStab
#define GEN_NULL_VECTOR

// How many BiCGStab iterations do we use?
#define GEN_NULL_VECTOR_STEP 10
#define GEN_NULL_VECTOR_REL_RESID 1e-8

// 1 for just const vector, 2 for const + even/odd vector. Specifies
// the number to generate if GEN_NULL_VECTOR is defined
#define VECTOR_COUNT 4

// Are we testing a random gauge rotation?
//#define TEST_RANDOM_GAUGE

// Are we testing a random field?
//#define TEST_RANDOM_FIELD

// Square staggered 2d operator w/out u1 function.
void square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square staggered 2d operator w/ u1 function.
void square_staggered_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);  

int main(int argc, char** argv)
{  
    // Declare some variables.
    cout << setiosflags(ios::fixed) << setprecision(6);
    int i, j, x, y;
    complex<double> *lattice; // Holds the gauge field.
    complex<double> *lhs, *rhs, *check; // For some Kinetic terms.
    double explicit_resid = 0.0;
    double bnorm = 0.0;
    inversion_info invif;
    
    // Create an RNG. 
    std::mt19937 generator (1337u); // 1337u is the seed. 
    
    // Describe the fine lattice. 
    int x_fine = N;
    int y_fine = N;
    int fine_size = x_fine*y_fine;
    
    // Inverter information.
    double outer_precision = 1e-8; 
    int outer_restart = 16; 
    double inner_precision = 1e-1;
    
    cout << "[VOL]: X " << x_fine << " Y " << y_fine << " Volume " << fine_size << "\n";
    
    // Initialize the lattice. Indexing: index = y*N + x.
    lattice = new complex<double>[2*fine_size];
    lhs = new complex<double>[fine_size];
    rhs = new complex<double>[fine_size];   
    check = new complex<double>[fine_size];   
    // Zero it out.
    zero<double>(lattice, 2*fine_size);
    zero<double>(rhs, fine_size);
    zero<double>(lhs, fine_size);
    zero<double>(check, fine_size);
    
    // Create a free lattice.
    cout << "[GAUGE]: Creating a gauge field.\n";
    unit_gauge_u1(lattice, x_fine, y_fine);
    
#ifdef TEST_RANDOM_GAUGE
    // Generate and perform a random gauge transformation.
    complex<double>* gauge_trans = new complex<double>[fine_size];
    rand_trans_u1(gauge_trans, x_fine, y_fine, generator);
    apply_gauge_trans_u1(lattice, gauge_trans, x_fine, y_fine);
    cout << "[GAUGE]: Performed a random gauge rotation.\n";
#endif
#ifdef TEST_RANDOM_FIELD
    gauss_gauge_u1(lattice, x_fine, y_fine, generator, BETA);
    cout << "[GAUGE]: Created a U(1) gauge field with angle standard deviation " << 1.0/sqrt(BETA) << "\n";
#endif
    
    cout << "[GAUGE]: The average plaquette is " << get_plaquette_u1(lattice, x_fine, y_fine) << ".\n";
    
    // Build an mg_struct.
    mg_operator_struct_complex mgstruct;
    mgstruct.x_fine = N;
    mgstruct.y_fine = N; 
    mgstruct.blocksize_x = X_BLOCKSIZE;
    mgstruct.blocksize_y = Y_BLOCKSIZE;
    mgstruct.n_vector = VECTOR_COUNT;
    mgstruct.matrix_vector = square_staggered_u1;
    mgstruct.matrix_extra_data = (void*)lattice; 
    
    cout << "[MG]: X_Block " << X_BLOCKSIZE << " Y_Block " << Y_BLOCKSIZE << " NullVectors " << VECTOR_COUNT << "\n";
    
    // Describe the coarse lattice. 
    int x_coarse = x_fine/mgstruct.blocksize_x; // how many coarse sites are there in the x dir?
    int y_coarse = y_fine/mgstruct.blocksize_y; // how many coarse sites are there in the y dir?
    int coarse_size = x_coarse*y_coarse; 


    // Set a point on the rhs.
    rhs[x_fine/2+(y_fine/2)*x_fine] = 1.0;

    // Get norm for rhs.
    bnorm = sqrt(norm2sq<double>(rhs, fine_size));

    // Set a point on the lhs.
    //lhs[x_fine/2+(y_fine/2)*x_fine+1] = 1.0;
    
    // Create a projector.
    mgstruct.projectors = new complex<double>*[mgstruct.n_vector];
    for (i = 0; i < mgstruct.n_vector; i++)
    {
        mgstruct.projectors[i] = new complex<double>[fine_size];
        zero<double>(mgstruct.projectors[i], fine_size);
    }
    

    cout << "[MG]: Creating " << mgstruct.n_vector << " projector(s).\n";
#ifndef GEN_NULL_VECTOR
    // Make a constant projector.
    if (mgstruct.n_vector == 1)
    {
        cout << "[MG]: Null vector 1 is a constant.\n";
        for (i = 0; i < fine_size; i++)
        {
            mgstruct.projectors[0][i] = 1;
#ifdef TEST_RANDOM_GAUGE
            mgstruct.projectors[0][i] *= conj(gauge_trans[i]);
#endif
        }
    }
    else if (mgstruct.n_vector == 2) // constant, even/odd phase. 
    {
        cout << "[MG]: Null vector 1 is a constant.\n";
        cout << "[MG]: Null vector 2 is an even/odd phase.\n";
        for (i = 0; i < fine_size; i++)
        {
            mgstruct.projectors[0][i] = 1;
            x = i % N;
            y = i / N;
            mgstruct.projectors[1][i] = ((x+y)%2 == 0) ? complex<double>(0.0,1.0) : complex<double>(0.0,-1.0);
#ifdef TEST_RANDOM_GAUGE
            mgstruct.projectors[0][i] *= conj(gauge_trans[i]);
            mgstruct.projectors[1][i] *= conj(gauge_trans[i]);
#endif
        }
    }
    else if (mgstruct.n_vector == 4) // 4 corners of hypercube.
    {
        cout << "[MG]: Null vector 1 is a constant on unit corner (0,0).\n";
        cout << "[MG]: Null vector 2 is a constant on unit corner (1,0).\n";
        cout << "[MG]: Null vector 3 is a constant on unit corner (0,1).\n";
        cout << "[MG]: Null vector 4 is a constant on unit corner (1,1).\n";
        for (i = 0; i < fine_size; i++)
        {
            x = i % N;
            y = i / N;
            mgstruct.projectors[2*(y%2)+(x%2)][i] = 1.0;
            
#ifdef TEST_RANDOM_GAUGE
            mgstruct.projectors[2*(y%2)+(x%2)][i] *= conj(gauge_trans[i]);
#endif
        }
    }
    else // invalid.
    {
        cout << "Unless you are generating null vectors, you can only use 1, 2, or 4 null vectors.\n";
        return 0;
    }

    
#ifdef TEST_RANDOM_GAUGE
    delete[] gauge_trans;
#endif
    
    
#ifdef PRINT_NULL_VECTOR
    cout << "[MG]: Check projector:\n"; 
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
    
    // We generate null vectors by solving Ax = 0, with a
    // gaussian initial guess.
    // For sanity with the residual, we really solve Ax = -Ax_0,
    // where x has a zero initial guess, x_0 is a random vector.
    complex<double>* rand_guess = new complex<double>[fine_size];
    complex<double>* Arand_guess = new complex<double>[fine_size];
    
    for (i = 0; i < mgstruct.n_vector; i++)
    {
        gaussian<double>(rand_guess, fine_size, generator);
        
        zero<double>(Arand_guess, fine_size); 
        square_staggered_u1(Arand_guess, rand_guess, (void*)lattice);
        for (j = 0; j < fine_size; j++)
        {
           Arand_guess[j] = -Arand_guess[j]; 
        }
        zero<double>(mgstruct.projectors[i], fine_size);
        
        minv_vector_gcr(mgstruct.projectors[i], Arand_guess, fine_size, GEN_NULL_VECTOR_STEP, GEN_NULL_VECTOR_REL_RESID, square_staggered_u1, (void*)lattice); 
        
        for (j = 0; j < fine_size; j++)
        {
            mgstruct.projectors[i][j] += rand_guess[j];
        }
        
        //minv_vector_gcr(mgstruct.projectors[i], rand_guess, fine_size, GEN_NULL_VECTOR_STEP, GEN_NULL_VECTOR_REL_RESID, square_staggered_u1, (void*)lattice); 
        
        normalize(mgstruct.projectors[i], fine_size); 
    }
    
    // This causes a segfault related to the RNG when
    // the vector is initialized.
    delete[] rand_guess; 
    delete[] Arand_guess; 
    
    // Normalize projectors.
    /*for (i = 0; i < mgstruct.n_vector; i++)
    {
        normalize(mgstruct.projectors[i], fine_size);
    }*/
    
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
    cout << "[MG]: Performing block orthonormalize of null vectors...\n";
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
        complex<double>* rhs_coarse = new complex<double>[coarse_size*mgstruct.n_vector];
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
        complex<double> *rhs_PdagP = new complex<double>[fine_size];
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
        complex<double>* rhs_A_coarse = new complex<double>[coarse_size*mgstruct.n_vector]; 
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
        complex<double>* rhs_PAP_fine = new complex<double>[fine_size];
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
        complex<double>* tmp_store = mgstruct.projectors[0];
        delete[] mgstruct.projectors;
        mgstruct.projectors = new complex<double>*[mgstruct.n_vector];
        mgstruct.projectors[0] = tmp_store; 
        mgstruct.projectors[1] = new complex<double>[fine_size];

        // Add an even/odd vector. 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                if ((x+y)%2 == 0)
                    mgstruct.projectors[1][x+y*x_fine] = complex<double>(0.0,1.0);
                else
                    mgstruct.projectors[1][x+y*x_fine] = complex<double>(0.0,-1.0);
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
        complex<double>* rhs_coarse = new complex<double>[coarse_size*mgstruct.n_vector];
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
        complex<double>* rhs_PdagP = new complex<double>[fine_size];
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
    
#ifdef COARSE_ONLY
    cout << "Solving coarse system only.\n";
    complex<double>* rhs_coarse = new complex<double>[coarse_size*mgstruct.n_vector];
    zero<double>(rhs_coarse, coarse_size*mgstruct.n_vector);
    restrict(rhs_coarse, rhs, &mgstruct);
    
    complex<double>* lhs_coarse = new complex<double>[coarse_size*mgstruct.n_vector];
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
    complex<double>* A_lhs_coarse = new complex<double>[coarse_size*mgstruct.n_vector];
    zero<double>(A_lhs_coarse, coarse_size*mgstruct.n_vector);
    
    coarse_square_laplace(A_lhs_coarse, lhs_coarse, (void*)&mgstruct);

    for (i = 0; i < coarse_size*mgstruct.n_vector; i++)
    {
      explicit_resid += real(conj(A_lhs_coarse[i] - rhs_coarse[i])*(A_lhs_coarse[i] - rhs_coarse[i]));
    }
    explicit_resid = sqrt(explicit_resid);

    printf("[check] should equal [rhs]. The residual is %15.20e.\n", explicit_resid);
    
    complex<double>* pro_lhs_coarse = new complex<double>[N*N];
    zero<double>(pro_lhs_coarse, N*N);
    prolong(pro_lhs_coarse, lhs_coarse, &mgstruct); 
    complex<double>* pro_rhs_coarse = new complex<double>[N*N];
    zero<double>(pro_rhs_coarse, N*N);
    square_laplacian(pro_rhs_coarse, pro_lhs_coarse, NULL);
    
#endif // COARSE_ONLY
    
    // Try a direct solve.
    cout << "\n[ORIG]: Solve fine system.\n";
    
    invif = minv_vector_gcr_restart(lhs, rhs, N*N, 100000, outer_precision, outer_restart, square_staggered_u1, (void*)lattice);
    
    if (invif.success == true)
    {
     printf("[ORIG]: Iterations %d RelRes %.8e Err N Algorithm %s\n", invif.iter, sqrt(invif.resSq)/bnorm, invif.name.c_str());
    }
    else // failed, maybe.
    {
     printf("[ORIG]: Iterations %d RelRes %.8e Err Y Algorithm %s\n", invif.iter, sqrt(invif.resSq)/bnorm, invif.name.c_str());
     printf("[ORIG]: This may be because the max iterations was reached.\n");
    }
    
    printf("[ORIG]: Computing [check] = A [lhs] as a confirmation.\n");

    // Check and make sure we get the right answer.
    zero<double>(check, fine_size);
    square_staggered_u1(check, lhs, (void*)lattice);

    explicit_resid = 0.0;
    for (i = 0; i < N*N; i++)
    {
      explicit_resid += real(conj(rhs[i] - check[i])*(rhs[i] - check[i]));
    }
    explicit_resid = sqrt(explicit_resid)/bnorm;

    printf("[ORIG]: [check] should equal [rhs]. The relative residual is %15.20e.\n", explicit_resid);
    
#ifdef COARSE_ONLY
    // Compare PAP solution to real solution. 
    cout << "\nCompare solutions.\n";
    double comparison = 0;
    double resid_comparison = 0;
    for (i = 0; i < N*N; i++)
    {
        comparison += real(conj(pro_lhs_coarse[i]-lhs[i])*(pro_lhs_coarse[i]-lhs[i]));
        resid_comparison += real(conj(pro_rhs_coarse[i]-rhs[i])*(pro_rhs_coarse[i]-rhs[i]));
    }
    comparison = sqrt(explicit_resid);
    printf("The solutions deviate by %15.20e.\n", comparison);
    printf("The projected residual has a rel res of %15.20e.\n", sqrt(resid_comparison)/bnorm);
    
    delete[] rhs_coarse; 
    delete[] lhs_coarse;
    delete[] A_lhs_coarse; 
    delete[] pro_lhs_coarse; 
    delete[] pro_rhs_coarse; 
#endif // COARSE_ONLY 

    // Let's actually test a multigrid solve!
    cout << "\n[MG]: Test MG solve.\n";
    
    // Block normalize the null vectors.
    block_normalize(&mgstruct); 
    
    // Set up the MG preconditioner. 
    mg_precond_struct_complex mgprecond;
    
    mgprecond.n_pre_smooth = 6; // 6 MinRes smoother steps before coarsening.
    mgprecond.n_post_smooth = 6; // 6 MinRes smoother steps after refining.
    mgprecond.in_solve_type = GCR; // What inner solver? MINRES, CG, or GCR.
    mgprecond.n_step = 10000; // max number of steps to use for inner solver.
    mgprecond.rel_res = inner_precision; // Maximum relative residual for inner solver.
    mgprecond.mgstruct = &mgstruct; // Contains null vectors, fine operator. (Since we don't construct the fine op.)
    mgprecond.matrix_vector = coarse_square_staggered; // Function which applies the coarse operator. 
    mgprecond.matrix_extra_data = (void*)&mgstruct; // What extra_data the coarse operator expects. 
    
    // Well, maybe this will work?
    zero<double>(lhs, N*N);
    invif = minv_vector_gcr_var_precond_restart(lhs, rhs, N*N, 10000, outer_precision, outer_restart, square_staggered_u1, (void*)lattice, mg_preconditioner, (void*)&mgprecond); /**/
    
    if (invif.success == true)
    {
     printf("[L1]: Iterations %d RelRes %.8e Err N Algorithm %s\n", invif.iter, sqrt(invif.resSq)/bnorm, invif.name.c_str());
    }
    else // failed, maybe.
    {
     printf("[L1]: Iterations %d RelRes %.8e Err Y Algorithm %s\n", invif.iter, sqrt(invif.resSq)/bnorm, invif.name.c_str());
     printf("[L1]: This may be because the max iterations was reached.\n");
    }


    printf("[MG]: Computing [check] = A [lhs] as a confirmation.\n");

    // Check and make sure we get the right answer.
    square_staggered_u1(check, lhs, (void*)lattice);

    explicit_resid = 0.0;
    for (i = 0; i < N*N; i++)
    {
      explicit_resid += real(conj(rhs[i] - check[i])*(rhs[i] - check[i]));
    }
    explicit_resid = sqrt(explicit_resid)/bnorm;

    printf("[MG]: [check] should equal [rhs]. The relative residual is %15.20e.\n", explicit_resid);

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
}



// Square lattice.
// Kinetic term for a 2D staggered w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" doesn't include anything.
// Apsi[x][y] = m psi[x,y] - U[y][x,x+1] 
void square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
   // Declare variables.
   int i;
   int x,y;
   double eta1;

   // For a 2D square lattice, the stencil is:
   //   1 |  0 -eta1  0 |
   //   - | +1    0  -1 |  , where eta1 = (-1)^x
   //   2 |  0 +eta1  0 |
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
      eta1 = 1 - 2*(x%2);
      
      // + e1.
      lhs[i] = lhs[i]-rhs[y*N+((x+1)%N)];
     
      // - e1.
      lhs[i] = lhs[i]+ rhs[y*N+((x+N-1)%N)]; // The extra +N is because of the % sign convention.
      
      // + e2.
      lhs[i] = lhs[i]- eta1*rhs[((y+1)%N)*N+x];
    
      // - e2.
      lhs[i] = lhs[i]+ eta1*rhs[((y+N-1)%N)*N+x];

      // Normalization.
      lhs[i] = 0.5*lhs[i];

      // 0
      // Added mass term here.
      lhs[i] = lhs[i]+ MASS*rhs[i];
   }
       
}

// Square lattice.
// Kinetic term for a 2D staggered w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" is a cast to a complex gauge_field[N*N*2], 
//    loaded by the function read_lattice_u1. 
// Apsi[x][y] = m psi[x,y] - U[y][x,x+1] 
void square_staggered_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
   // Declare variables.
   int i;
   int x,y;
   double eta1; 
    
   complex<double>* lattice = (complex<double>*)extra_data; 

   // For a 2D square lattice, the stencil is:
   //   1 |  0 -eta1  0 |
   //   - | +1    0  -1 |  , where eta1 = (-1)^x
   //   2 |  0 +eta1  0 |
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
      eta1 = 1 - 2*(x%2);
      
      // + e1.
      lhs[i] = lhs[i]-conj(lattice[y*N*2+x*2])*rhs[y*N+((x+1)%N)];
     
      // - e1.
      lhs[i] = lhs[i]+ lattice[y*N*2+((x+N-1)%N)*2]*rhs[y*N+((x+N-1)%N)]; // The extra +N is because of the % sign convention.
      
      // + e2.
      lhs[i] = lhs[i]- eta1*conj(lattice[y*N*2+x*2+1])*rhs[((y+1)%N)*N+x];
    
      // - e2.
      lhs[i] = lhs[i]+ eta1*lattice[((y+N-1)%N)*N*2+x*2+1]*rhs[((y+N-1)%N)*N+x];

      // Normalization.
      lhs[i] = 0.5*lhs[i];

      // 0
      // Added mass term here.
      lhs[i] = lhs[i]+ MASS*rhs[i];
   }
       
}

